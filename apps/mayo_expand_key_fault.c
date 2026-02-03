// Simulation: compute L, compute L' with one faulted m_vec_mul_add, recover oil scalar O[j]
// target = L = P1P1t_times_O
// fault in multable tab : zero the lower byte

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <mayo.h>
#include <arithmetic.h>
#include <simple_arithmetic.h>
#include <randombytes.h>
#include <mem.h>


static void decode_to_nibbles(const unsigned char *in, size_t len,
                              unsigned char *out, size_t n)
{
    size_t j = 0;
    for (size_t i = 0; i < len && j < n; i++)
    {
        out[j++] = in[i] & 0xF;
        if (j < n)
            out[j++] = in[i] >> 4;
    }
}

// static unsigned char mul_f(unsigned char a, unsigned char b)
// {
//     unsigned char p;
//     p = (a & 1) * b;
//     p ^= (a & 2) * b;
//     p ^= (a & 4) * b;
//     p ^= (a & 8) * b;

//     // reduce mod x^4 + x + 1
//     unsigned char top_p = p & 0xf0;
//     unsigned char out = (p ^ (top_p >> 4) ^ (top_p >> 3)) & 0x0f;
//     return out;
// }

static inline uint32_t mul_table_fault(uint8_t a)
{
    uint32_t x = (uint32_t)a * 0x08040201;
    uint32_t hm = 0xf0f0f0f0;
    uint32_t hh = x & hm;
    uint32_t tab = x ^ (hh >> 4) ^ (hh >> 3);

    // FAULT
    tab &= ~0x000000FF; // zero low byte → remove γ0*(a·1)
    return tab;
}

// static unsigned char inverse_f(unsigned char a)
// {
//     unsigned char a2 = mul_f(a, a);
//     unsigned char a4 = mul_f(a2, a2);
//     unsigned char a8 = mul_f(a4, a4);
//     unsigned char a6 = mul_f(a2, a4);
//     unsigned char a14 = mul_f(a8, a6);

//     return a14;
// }

static inline void m_vec_mul_add_fault(int m_vec_limbs,
                                       const uint64_t *in,
                                       unsigned char a,
                                       uint64_t *acc)
{
    (void)m_vec_limbs;
    uint32_t tab = mul_table_fault(a);
    uint64_t ask = 0x1111111111111111ULL;

    for (int i = 0; i < M_VEC_LIMBS_MAX; i++)
    {
        acc[i] ^= (in[i] & ask) * (tab & 0xff) ^ ((in[i] >> 1) & ask) * ((tab >> 8) & 0xf) ^ ((in[i] >> 2) & ask) * ((tab >> 16) & 0xf) ^ ((in[i] >> 3) & ask) * ((tab >> 24) & 0xf);
    }
}

static void P1P1t_times_O_fault(const mayo_params_t *p, const uint64_t *P1,
                                const unsigned char *O, uint64_t *acc,
                                int bs_target)
{
    const int o = PARAM_o(p);
    const int v = PARAM_v(p);
    const int m_vec_limbs = PARAM_m_vec_limbs(p);

    int bs_mat_entries_used = 0;

    for (int r = 0; r < v; r++)
    {
        for (int c = r; c < v; c++)
        {

            if (c == r)
            {
                bs_mat_entries_used++;
                continue;
            }

            const uint64_t *Pblk =
                P1 + m_vec_limbs * bs_mat_entries_used;

            for (int k = 0; k < o; k++)
            {
                uint64_t *dst_rc = acc + m_vec_limbs * (r * o + k);
                uint64_t *dst_cr = acc + m_vec_limbs * (c * o + k);

                if (bs_mat_entries_used == bs_target)
                {
                    // FAULTY multiplication
                    m_vec_mul_add_fault(m_vec_limbs, Pblk, O[c * o + k], dst_rc);
                    m_vec_mul_add_fault(m_vec_limbs, Pblk, O[r * o + k], dst_cr);
                }
                else
                {
                    // normal
                    m_vec_mul_add(m_vec_limbs, Pblk, O[c * o + k], dst_rc);
                    m_vec_mul_add(m_vec_limbs, Pblk, O[r * o + k], dst_cr);
                }
            }
            bs_mat_entries_used++;
        }
    }
}

int mayo_sign_fault(const mayo_params_t *p, unsigned char *sig,
                    size_t *siglen, const unsigned char *m, size_t mlen,
                    const unsigned char *csk,
                    int bs_target)
{

    const unsigned char *seed_sk;
    sk_t sk __attribute__((aligned(32)));
    int ret = mayo_expand_sk(p, csk, &sk);
    if (ret != MAYO_OK)
        return ret;

    seed_sk = csk;

    uint64_t *P1 = sk.p;
    uint64_t *L = P1 + PARAM_P1_limbs(p);
    uint64_t Ltmp[K_MAX * O_MAX * M_VEC_LIMBS_MAX];

    unsigned char tenc[M_BYTES_MAX], t[M_MAX];                               // no secret data
    unsigned char y[M_MAX];                                                  // secret data
    unsigned char salt[SALT_BYTES_MAX];                                      // not secret data
    unsigned char V[K_MAX * V_BYTES_MAX + R_BYTES_MAX], Vdec[V_MAX * K_MAX]; // secret data
    unsigned char A[((M_MAX + 7) / 8 * 8) * (K_MAX * O_MAX + 1)] = {0};      // secret data
    unsigned char x[K_MAX * N_MAX];                                          // not secret data
    unsigned char rvec[K_MAX * O_MAX + 1] = {0};                             // secret data
    unsigned char s[K_MAX * N_MAX];
    unsigned char tmp[DIGEST_BYTES_MAX + SALT_BYTES_MAX + SK_SEED_BYTES_MAX + 1];
    unsigned char *vi;
    unsigned char Ox[V_MAX];

    memset(A, 0, sizeof(A));
    memset(rvec, 0, sizeof(rvec));

    const int param_m = PARAM_m(p);
    const int param_n = PARAM_n(p);
    const int param_o = PARAM_o(p);
    const int param_k = PARAM_k(p);
    const int param_v = PARAM_v(p);
    const int param_m_bytes = PARAM_m_bytes(p);
    const int param_v_bytes = PARAM_v_bytes(p);
    const int param_r_bytes = PARAM_r_bytes(p);
    const int param_sig_bytes = PARAM_sig_bytes(p);
    const int param_A_cols = PARAM_A_cols(p);
    const int param_digest_bytes = PARAM_digest_bytes(p);
    const int param_sk_seed_bytes = PARAM_sk_seed_bytes(p);
    const int param_salt_bytes = PARAM_salt_bytes(p);

    shake256(tmp, param_digest_bytes, m, mlen);
    memcpy(tmp + param_digest_bytes, salt, param_salt_bytes);

    randombytes(tmp + param_digest_bytes, param_salt_bytes);
    memcpy(tmp + param_digest_bytes + param_salt_bytes, seed_sk,
           param_sk_seed_bytes);

    shake256(salt, param_salt_bytes, tmp,
             param_digest_bytes + param_salt_bytes + param_sk_seed_bytes);

    memcpy(tmp + param_digest_bytes, salt, param_salt_bytes);

    shake256(tenc, param_m_bytes, tmp, param_digest_bytes + param_salt_bytes);

    decode(tenc, t, param_m);

    unsigned char *ctrbyte = tmp + param_digest_bytes + param_salt_bytes + param_sk_seed_bytes;
    for (int ctr = 0; ctr < 256; ++ctr)
    {
        *ctrbyte = (unsigned char)ctr;
        shake256(V, param_k * param_v_bytes + param_r_bytes, tmp,
                 param_digest_bytes + param_salt_bytes + param_sk_seed_bytes + 1);

        // decode the v_i vectors
        for (int i = 0; i <= param_k - 1; ++i)
        {
            decode(V + i * param_v_bytes, Vdec + i * param_v, param_v);
        }

        // fault in L

        memset(Ltmp, 0, sizeof(Ltmp));

        memset(L, 0, PARAM_P2_limbs(p) * sizeof(uint64_t));
        P1P1t_times_O_fault(p, P1, sk.O, L, bs_target);

        // continue with normal flow

        compute_M_and_VPV(p, Vdec, L, P1, Ltmp, (uint64_t *)A);

        compute_rhs(p, (uint64_t *)A, t, y);
        compute_A(p, Ltmp, A);

        for (int i = 0; i < param_m; i++)
        {
            A[(1 + i) * (param_k * param_o + 1) - 1] = 0;
        }
        decode(V + param_k * param_v_bytes, rvec,
               param_k *
                   param_o);
        if (sample_solution(p, A, y, rvec, x, param_k, param_o, param_m, param_A_cols))
        {
            goto success;
        }
        else
        {
            memset(Ltmp, 0, sizeof(Ltmp));
            memset(A, 0, sizeof(A));
        }
    }

    goto err;

success:
    for (int i = 0; i <= param_k - 1; ++i)
    {
        vi = Vdec + i * (param_n - param_o);
        mat_mul(sk.O, x + i * param_o, Ox, param_o, param_n - param_o, 1);
        mat_add(vi, Ox, s + i * param_n, param_n - param_o, 1);
        memcpy(s + i * param_n + (param_n - param_o), x + i * param_o, param_o);
    }
    encode(s, sig, param_n * param_k);

    memcpy(sig + param_sig_bytes - param_salt_bytes, salt, param_salt_bytes);
    *siglen = param_sig_bytes;
    return MAYO_OK;

err:
    mayo_secure_clear(V, sizeof(V));
    mayo_secure_clear(Vdec, sizeof(Vdec));
    mayo_secure_clear(A, sizeof(A));
    mayo_secure_clear(rvec, sizeof(rvec));
    mayo_secure_clear(sk.O, sizeof(sk.O));
    mayo_secure_clear(&sk, sizeof(sk_t));
    mayo_secure_clear(Ox, sizeof(Ox));
    mayo_secure_clear(tmp, sizeof(tmp));
    mayo_secure_clear(Ltmp, sizeof(Ltmp));
    return MAYO_ERR;
}

int main(void)
{
    const mayo_params_t *p = &MAYO_1;

    const int o = PARAM_o(p);
    const int v = PARAM_v(p);
    const int k = PARAM_k(p);
    const int n = PARAM_n(p);

    printf("PARAMS: n=%d m=%d v=%d o=%d k=%d\n",
           n, PARAM_m(p), v, o, k);

    //------------------------------------------------------------
    // 1. Generate keypair
    //------------------------------------------------------------
    unsigned char pk[PARAM_cpk_bytes(p)];
    unsigned char sk[PARAM_csk_bytes(p)];
    mayo_keypair(p, pk, sk);

    //------------------------------------------------------------
    // Expand SK so we know true O
    //------------------------------------------------------------
    sk_t E;
    mayo_expand_sk(p, sk, &E);

    //------------------------------------------------------------
    // Choose an oil nibble to attack
    //------------------------------------------------------------
    int oil_index = 0; // attack O[0]
    unsigned char trueOil = E.O[oil_index] & 0xF;
    printf("True O[%d] = 0x%x\n", oil_index, trueOil);

    //------------------------------------------------------------
    // 2. Produce baseline signature
    //------------------------------------------------------------
    unsigned char msg[32] = {0x42};
    unsigned char sig[512], sig_fault[512];
    size_t slen = sizeof(sig), slen2 = sizeof(sig_fault);

    mayo_sign(p, sig, &slen, msg, sizeof(msg), sk);

    //------------------------------------------------------------
    // 3. Produce faulted signature (fault on bs_entry = 0)
    //------------------------------------------------------------
    int bs_target = 0;
    mayo_sign_fault(p, sig_fault, &slen2, msg, sizeof(msg), sk, bs_target);

    //------------------------------------------------------------
    // 4. Compute Δsig in nibble space
    //------------------------------------------------------------
    size_t nibcount = k * n;
    unsigned char *n0 = calloc(nibcount, 1);
    unsigned char *nF = calloc(nibcount, 1);
    decode_to_nibbles(sig, slen, n0, nibcount);
    decode_to_nibbles(sig_fault, slen2, nF, nibcount);

    unsigned char *Delta = calloc(nibcount, 1);
    for (size_t i = 0; i < nibcount; i++)
        Delta[i] = n0[i] ^ nF[i];

    //------------------------------------------------------------
    // 5. Compute Γ_j propagated to signature:
    // Make a fake key with O[j]=1, all other O=0.
    //------------------------------------------------------------
    unsigned char O_unit[v * o];
    memset(O_unit, 0, sizeof(O_unit));
    O_unit[oil_index] = 1;

    // Build SK_unit = SK but with O replaced by O_unit
    sk_t E_unit = E;
    memcpy(E_unit.O, O_unit, v * o);

    // Re-sign with modified O
    unsigned char sig_unit[512];
    size_t slenU = sizeof(sig_unit);
    mayo_sign(p, sig_unit, &slenU, msg, sizeof(msg),
              (unsigned char *)(&E_unit));

    unsigned char *nU = calloc(nibcount, 1);
    decode_to_nibbles(sig_unit, slenU, nU, nibcount);

    //------------------------------------------------------------
    // 6. Recover oil nibble: try a in 0..15 such that:
    //   nU * a == Delta   (component-wise GF(16) multiply)
    //------------------------------------------------------------
    int recovered = -1;
    for (size_t i = 0; i < nibcount; i++)
    {
        unsigned char g = nU[i] & 0xF; // fn(j)[i]
        if (g == 0)
            continue;

        unsigned char d = Delta[i] & 0xF; // delta_sig[i]

        // Solve O_j = g^{-1} * d
        unsigned char inv = inverse_f(g);
        unsigned char a = mul_f(inv, d);

        recovered = a;
        break;
    }

    printf("Recovered O[%d] = 0x%x\n", oil_index, recovered);

    free(n0);
    free(nF);
    free(Delta);
    free(nU);

    return 0;
}
