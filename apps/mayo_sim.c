// SPDX-License-Identifier: Apache-2.0
//
// Fault simulation harness for MAYO

// CHanges in mayo : skip computation for v[0][0] and
// the zeroth col of L to get entire VL[0][0] as zero

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

static void encode_from_nibbles(unsigned char *sig, const unsigned char *nibbles, int n)
{
    memset(sig, 0, (n + 1) / 2);
    for (int i = 0; i < n; i++)
    {
        if (i & 1)
            sig[i >> 1] |= (nibbles[i] & 0xF) << 4;
        else
            sig[i >> 1] |= (nibbles[i] & 0xF);
    }
}

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
    // int v = PARAM_v(p);
    // int k = PARAM_k(p);
    int n = PARAM_n(p);

    unsigned char xvec[o]; // VLA; safe because no 'goto' jumps over it
    size_t bufsize = PARAM_sig_bytes(p) + sizeof(msg);
    size_t siglen;

    /* ---- Allocate memory ---- */
    pk = calloc(1, PARAM_cpk_bytes(p));
    sk = calloc(1, PARAM_csk_bytes(p));
    esk = calloc(1, sizeof(sk_t));
    epk = calloc(1, sizeof(pk_t));
    sig = malloc(bufsize);

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

    printf("Signature (%zu bytes):\n", siglen);
    for (size_t i = 0; i < siglen; i++)
    {
        printf("%02x", sig[i]);
        // Optional: Add a newline every 32 bytes for readability
        if ((i + 1) % 32 == 0)
            printf("\n");
    }
    printf("\n\n");

    nibbles = calloc(n, 1);
    decode_to_nibbles(sig, (int)siglen, nibbles, n);

    printf("Guessing missing (O x)[0] using verify oracle\n");

    unsigned char recovered = 0xFF;
    unsigned char orig = nibbles[0];

    for (unsigned char i = 0; i < 16; i++)
    {
        nibbles[0] = orig ^ i;
        encode_from_nibbles(sig, nibbles, n);

        int res = mayo_verify(p, msg, sizeof(msg), sig, pk);
        if (res == MAYO_OK)
        {
            recovered = i;
            // printf("[Recovered Ox[0] = 0x%x (verify passed)\n", i);
            break;
        }
    }

    nibbles[0] = orig; // restore

    if (recovered == 0xFF)
    {
        printf("[FAIL] No guess validated\n");
        goto cleanup;
    }

    printf("[RECOVERED] (O x)[0] = 0x%x\n", recovered);
    /* ---- Verify ---- */

    // printf("[FAULT] Signature VALID\n");

    // /* ---- Decode signature ---- */
    // nibbles = calloc(n, 1);
    // decode_to_nibbles(sig, (int)siglen, nibbles, n);

    // printf("[SIM] decoded s[0] = 0x%x\n", nibbles[0]);

    /* ---- Extract oil ---- */
    int oil_offset = n - o;
    for (int i = 0; i < o; i++)
        xvec[i] = nibbles[oil_offset + i] & 0xF;

    /* ---- Compute true (O x)[0] ---- */
    unsigned char trueOx0 = compute_Ox0(esk, xvec, o);
    printf("True (O x)[0] = 0x%x\n", trueOx0);

    if (trueOx0 == recovered)
        printf("MATCH: recovered (O x)[0] matches\n");
    else
        printf("FAIL: mismatch\n");
    printf("\n");
    // if (res == MAYO_OK)
    //     return 1;
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

    return 0;
}
