#include <mayo.h>
#include <mem.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

// ---------------- GF(16) arithmetic ----------------

static unsigned char mul_f(unsigned char a, unsigned char b)
{
    unsigned char p;
    p = (a & 1) * b;
    p ^= (a & 2) * b;
    p ^= (a & 4) * b;
    p ^= (a & 8) * b;

    // reduce mod x^4 + x + 1
    unsigned char top_p = p & 0xf0;
    unsigned char out = (p ^ (top_p >> 4) ^ (top_p >> 3)) & 0x0f;
    return out;
}

// ---------------- Fault Simulation ----------------
static void decode_to_nibbles(const unsigned char *in, int len, unsigned char *out, int n)
{
    int j = 0;
    for (int i = 0; i < len && j < n; ++i)
    {
        out[j++] = in[i] & 0xF;
        if (j < n)
            out[j++] = in[i] >> 4;
    }
}

// static void encode_from_nibbles(unsigned char *sig, const unsigned char *nibbles, int n)
// {
//     memset(sig, 0, (n + 1) / 2);
//     for (int i = 0; i < n; i++)
//     {
//         if (i & 1)
//             sig[i >> 1] |= (nibbles[i] & 0xF) << 4;
//         else
//             sig[i >> 1] |= (nibbles[i] & 0xF);
//     }
// }

unsigned char compute_Ox0(const sk_t *esk, const unsigned char *x, int o)
{
    unsigned char val = 0;
    for (int j = 0; j < o; j++)
    {
        unsigned char O0j = esk->O[j] & 0xF; // row 0 → O[0*o + j]
        val ^= mul_f(O0j, x[j]);
    }
    return val & 0xF;
}

static int example_fault_sim(const mayo_params_t *p)
{
    printf("Running MAYO fault simulation using %s\n", PARAM_name(p));

    unsigned char msg[32] = {0xAB};

    /* ---- Move ALL declarations here ---- */
    unsigned char *pk = NULL;
    unsigned char *sk = NULL;
    unsigned char *sig = NULL;
    unsigned char *nibbles = NULL;
    sk_t *esk = NULL;
    uint64_t *epk = NULL;

    int o = PARAM_o(p);
    int v = PARAM_v(p);
    int k = PARAM_k(p);
    int n = PARAM_n(p);

    size_t bufsize = PARAM_sig_bytes(p) + sizeof(msg);
    size_t siglen;

    unsigned char *s0 = NULL;
    unsigned char *s1 = NULL;
    unsigned char *x0 = NULL;
    unsigned char *x1 = NULL;
    unsigned char *delta_s = NULL;
    unsigned char *dx = NULL;
    unsigned char *Ox_dx = NULL;
    int ok = 1;

    /* ---- Allocate memory ---- */
    pk = calloc(1, PARAM_cpk_bytes(p));
    sk = calloc(1, PARAM_csk_bytes(p));
    esk = calloc(1, sizeof(sk_t));
    epk = calloc(1, sizeof(pk_t));
    sig = malloc(bufsize);
    x0 = calloc(o, sizeof(unsigned char));
    x1 = calloc(o, sizeof(unsigned char));
    delta_s = calloc(v, sizeof(unsigned char));
    dx = calloc(o, sizeof(unsigned char));
    Ox_dx = calloc(v, sizeof(unsigned char));

    siglen = bufsize;

    /* ---- Keygen ---- */
    mayo_keypair(p, pk, sk);
    mayo_expand_sk(p, sk, esk);
    mayo_expand_pk(p, pk, epk);

    /* ---- Sign (fault occurs internally) ---- */
    if (mayo_sign(p, sig, &siglen, msg, sizeof(msg), sk) != MAYO_OK)
    {
        printf("mayo_sign failed\n");
        goto cleanup;
    }

    // printf("Signature (%zu bytes):\n", siglen);
    // for (size_t i = 0; i < siglen; i++)
    // {
    //     printf("%02x", sig[i]);
    //     // Optional: Add a newline every 32 bytes for readability
    //     if ((i + 1) % 32 == 0)
    //         printf("\n");
    // }
    // printf("\n\n");

    nibbles = calloc(n * k, 1);
    decode_to_nibbles(sig, (int)siglen, nibbles, n * k);

    s0 = nibbles;
    s1 = nibbles + n;

    for (int i = 0; i < o; i++)
    {
        x0[i] = s0[v + i] & 0xF;
        x1[i] = s1[v + i] & 0xF;
    }

    for (int i = 0; i < v; i++)
    {
        delta_s[i] = (s0[i] ^ s1[i]) & 0xF;
    }

    printf("delta_s = s0 (xor) s1 = O · (x0 (xor) x1)\n");

    for (int i = 0; i < o; i++)
        dx[i] = x0[i] ^ x1[i];

    memset(Ox_dx, 0, v);

    for (int row = 0; row < v; row++)
    {
        for (int j = 0; j < o; j++)
        {
            unsigned char Oij = esk->O[row * o + j] & 0xF;
            Ox_dx[row] ^= mul_f(Oij, dx[j]);
        }
    }

    ok = 1;
    for (int i = 0; i < v; i++)
    {
        if (Ox_dx[i] != delta_s[i])
        {
            ok = 0;
            break;
        }
    }

    printf("delta_s:\n");
    for (int i = 0; i < v; i++)
    {
        printf("%02x ", delta_s[i] & 0xF);
        if ((i + 1) % 32 == 0)
            printf("\n");
    }
    printf("\n");

    printf("Ox_dx:\n");
    for (int i = 0; i < v; i++)
    {
        printf("%02x ", Ox_dx[i] & 0xF);
        if ((i + 1) % 32 == 0)
            printf("\n");
    }
    printf("\n");

    if (ok)
        printf(" \n [SUCCESS] delta_s matches O · (x0 xor x1)\n");
    else
        printf("\n [FAIL] delta_s mismatch\n");

cleanup:
    if (sig)
        free(sig);
    if (nibbles)
        free(nibbles);
    if (pk)
        free(pk);
    if (sk)
        mayo_secure_free(sk, PARAM_csk_bytes(p));
    if (esk)
        free(esk);
    if (epk)
        free(epk);
    if (x0)
        free(x0);
    if (x1)
        free(x1);
    if (delta_s)
        free(delta_s);
    if (dx)
        free(dx);
    if (Ox_dx)
        free(Ox_dx);

    return 0;
}

int main(void)
{
    example_fault_sim(NULL);
}
