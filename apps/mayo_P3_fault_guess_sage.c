#include <inttypes.h>
#include <mayo.h>

#include <aes_ctr.h>
#include <arithmetic.h>
#include <fips202.h>
#include <mem.h>
#include <randombytes.h>
#include <simple_arithmetic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#define MAX_UNK 7000
#define MAX_EQ  6000
#define PK_PRF  AES_128_CTR


unsigned char seed_sk_fixed[24] = {
    0x9a, 0xf8, 0xcb, 0x5c, 0xb5, 0xb7, 0xb1, 0x2c, 0x7a, 0x15, 0xac, 0x4e,
    0x0a, 0x26, 0xae, 0x99, 0x6c, 0x6c, 0xf8, 0x9a, 0x95, 0xb5, 0xfb, 0x64};

// GF(16) arithmetic




static unsigned char extract_m_element(const uint64_t *vec, int ell,
                                       int m_vec_limbs) {
  (void)m_vec_limbs;
  int limb = ell / 16;
  int pos = ell % 16;
  uint64_t w = vec[limb];
  int base = pos * 4;
  unsigned char val = 0;
  if (w & (1ULL << (base + 0))) val |= 1;
  if (w & (1ULL << (base + 1))) val |= 2;
  if (w & (1ULL << (base + 2))) val |= 4;
  if (w & (1ULL << (base + 3))) val |= 8;
  return val;
}

static void reconstruct_full_P3(const mayo_params_t *p, const uint64_t *epk,
                                uint64_t *P3_full) {
  int o = PARAM_o(p);
  int m_vec_limbs = PARAM_m_vec_limbs(p);
  const uint64_t *P3_upper = epk + PARAM_P1_limbs(p) + PARAM_P2_limbs(p);
  int idx = 0;

  for (int r = 0; r < o; r++) {
    for (int c = r; c < o; c++) {
      const uint64_t *src = P3_upper + m_vec_limbs * idx;
      for (int l = 0; l < m_vec_limbs; l++) {
        P3_full[m_vec_limbs * (r * o + c) + l] = src[l];
        if (r != c)
          P3_full[m_vec_limbs * (c * o + r) + l] = src[l];
      }
      idx++;
    }
  }
}



static inline void P1_times_O_correct(const mayo_params_t *p, const uint64_t *P1,
                                      const unsigned char *O, uint64_t *acc) {
#ifndef ENABLE_PARAMS_DYNAMIC
  (void)p;
#endif
  mul_add_m_upper_triangular_mat_x_mat(PARAM_m_vec_limbs(p), P1, O, acc,
                                       PARAM_v(p), PARAM_v(p), PARAM_o(p), 1);
}

static inline void compute_P3_correct(const mayo_params_t *p, const uint64_t *P1,
                                      uint64_t *P2, const unsigned char *O,
                                      uint64_t *P3) {
  const int m_vec_limbs = PARAM_m_vec_limbs(p);
  const int param_v = PARAM_v(p);
  const int param_o = PARAM_o(p);

  P1_times_O_correct(p, P1, O, P2);
  mul_add_mat_trans_x_m_mat(m_vec_limbs, O, P2, P3, param_v, param_o, param_o);
}

static void dump_hex(const char *label, const unsigned char *buf, size_t len) {
  printf("%s (%zu bytes):\n", label, len);
  for (size_t i = 0; i < len; i++) {
    printf("%02x", buf[i]);
    if ((i + 1) % 32 == 0) printf("\n");
    else if ((i + 1) % 2 == 0) printf(" ");
  }
  if (len % 32) printf("\n");
}

/*
 * Print the polynomial system induced by the row-0-only P1*O fault.
 *
 * Variables:
 *
 *      x_k_j = O[k][j]
 *
 * for:
 *
 *      0 <= k < v
 *      0 <= j < o
 *
 * Output format is SageMath-compatible.
 *
 * Each equation is printed as:
 *
 *      polynomial == 0
 *
 * Coefficients are printed as F(n), where the Sage script defines
 * F(n) as the conversion from a MAYO GF(16) nibble to GF(16).
 */
static void dump_linear_equations(
    const mayo_params_t *p,
    const uint64_t *epk)
{
    int v   = PARAM_v(p);
    int o   = PARAM_o(p);
    int m   = PARAM_m(p);
    int mvl = PARAM_m_vec_limbs(p);

    /*
     * Expanded public key layout:
     *
     *     P1 || P2 || P3
     *
     * Since P1_times_O() was skipped completely,
     * P2 was never modified:
     *
     *     P2_fault = P2_old
     *
     * and:
     *
     *     P3_fault = O^T * P2_old
     */
    const uint64_t *P2 =
        epk + PARAM_P1_limbs(p);

    /*
     * Reconstruct the stored upper-triangular P3
     * into a full o x o matrix.
     */
    uint64_t *P3_full =
        calloc((size_t)o * o * mvl, sizeof(uint64_t));

    if (!P3_full) {
        perror("calloc P3_full");
        return;
    }

    reconstruct_full_P3(p, epk, P3_full);

    FILE *fp = fopen("../mayo_equations_linear.txt", "w");

    if (!fp) {
        perror("mayo_equations_linear.txt");
        free(P3_full);
        return;
    }

    fprintf(fp,
            "# MAYO linear equations from fault:\n");
    fprintf(fp,
            "# P1_times_O() completely skipped\n");
    fprintf(fp,
            "# Faulty relation: P3 = O^T * P2\n\n");

    fprintf(fp, "# v = %d\n", v);
    fprintf(fp, "# o = %d\n", o);
    fprintf(fp, "# m = %d\n", m);

    fprintf(fp,
            "# unknowns = v*o = %d\n",
            v * o);

    fprintf(fp,
            "# equations = m*o*(o+1)/2 = %d\n",
            m * (o * (o + 1) / 2));

    fprintf(fp,
            "# variable x_k_j means O[k][j]\n\n");

    int eq = 0;

    /*
     * For every GF(16) coefficient/component ell.
     */
    for (int ell = 0; ell < m; ell++) {

        /*
         * P3 is stored as an upper-triangular matrix,
         * so only generate equations for i <= j.
         */
        for (int i = 0; i < o; i++) {

            for (int j = i; j < o; j++) {

                /*
                 * Public RHS:
                 *
                 *     y = P3[i][j][ell]
                 */
                unsigned char y =
                    extract_m_element(
                        P3_full +
                        (size_t)mvl * (i * o + j),
                        ell,
                        mvl);

                fprintf(fp,
                        "# eq=%d ell=%d i=%d j=%d\n",
                        eq, ell, i, j);

                int first = 1;

#define PRINT_TERM(...)                    \
    do {                                   \
        if (!first)                        \
            fprintf(fp, " + ");            \
        fprintf(fp, __VA_ARGS__);          \
        first = 0;                         \
    } while (0)

                /*
                 * ==================================================
                 * LINEAR EQUATION FROM:
                 *
                 *       P3 = O^T * P2
                 *
                 * ==================================================
                 *
                 * For off-diagonal i < j:
                 *
                 * P3[i,j] =
                 *
                 *   sum_k O[k,i] * P2[k,j]
                 *
                 * + sum_k O[k,j] * P2[k,i]
                 *
                 *
                 * For diagonal i == j:
                 *
                 * P3[i,i] =
                 *
                 *   sum_k O[k,i] * P2[k,i]
                 *
                 * We use only one direction on the diagonal,
                 * consistent with MAYO's upper-triangular
                 * representation.
                 */

                for (int k = 0; k < v; k++) {

                    /*
                     * First direction:
                     *
                     *     O[k,i] * P2[k,j]
                     *
                     * Unknown:
                     *
                     *     x_k_i = O[k][i]
                     */
                    unsigned char coeff =
                        extract_m_element(
                            P2 +
                            (size_t)mvl * (k * o + j),
                            ell,
                            mvl);

                    if (coeff != 0) {

                        PRINT_TERM(
                            "F(%u)*x_%d_%d",
                            coeff,
                            k,
                            i);
                    }

                    /*
                     * For off-diagonal entries only:
                     *
                     *     O[k,j] * P2[k,i]
                     *
                     * Unknown:
                     *
                     *     x_k_j = O[k][j]
                     */
                    if (i != j) {

                        coeff =
                            extract_m_element(
                                P2 +
                                (size_t)mvl * (k * o + i),
                                ell,
                                mvl);

                        if (coeff != 0) {

                            PRINT_TERM(
                                "F(%u)*x_%d_%d",
                                coeff,
                                k,
                                j);
                        }
                    }
                }

                /*
                 * Move P3[i,j] to the left-hand side:
                 *
                 *     linear_expression = y
                 *
                 * Over GF(16), characteristic = 2, so:
                 *
                 *     -y = +y
                 *
                 * Therefore:
                 *
                 *     linear_expression + y = 0
                 */
                if (y != 0) {

                    PRINT_TERM(
                        "F(%u)",
                        y);
                }

                /*
                 * Handle an all-zero equation.
                 */
                if (first)
                    fprintf(fp, "F(0)");

                fprintf(fp, " == 0\n\n");

#undef PRINT_TERM

                eq++;
            }
        }
    }

    fclose(fp);
    free(P3_full);

    printf("\n");
    printf("Linear system written to mayo_equations_linear.txt\n");
    printf("Fault model : P1_times_O() completely skipped\n");
    printf("Relation    : P3 = O^T * P2\n");
    printf("Unknowns    : %d\n", v * o);
    printf("Equations   : %d\n", eq);
}
// =====================================================================
// Top-level: run genuinely faulty keygen, then recover.
// =====================================================================

static void example_fault_row0_only(const mayo_params_t *p) {
  printf("Fault sim: bs_mat_rows forced to 1 (only P2 row 0 = P1O + P2_old)\n");

  unsigned char *pk = calloc(1, PARAM_cpk_bytes(p));
  unsigned char *sk = seed_sk_fixed;
  sk_t *esk = calloc(1, sizeof(sk_t));
  uint64_t *epk = calloc(1, sizeof(pk_t));

  // Linked against the faulted compute_P3 above, so this keygen
  // genuinely produces a faulty pk (only P2 row 0 updated correctly).
  mayo_keypair(p, pk, sk);

  mayo_expand_sk(p, sk, esk);
  mayo_expand_pk(p, pk, epk); // re-derives P1/P2_old fresh from the seed --
                              // these are correct/public regardless of the
                              // fault, since only the stored P3 bytes in pk
                              // reflect the faulted compute_P3 call.

  dump_hex("Faulty pk", pk, PARAM_cpk_bytes(p));

  // Optional: show how far the faulty P3 diverges from the correct one,
  // purely informational (uses a private, correct-computation copy of O).
  {
    int  o = PARAM_o(p), mvl = PARAM_m_vec_limbs(p);
    int param_P1_limbs = PARAM_P1_limbs(p);

    uint64_t *P2_copy = calloc(PARAM_P2_limbs(p), sizeof(uint64_t));
    memcpy(P2_copy, epk + param_P1_limbs, PARAM_P2_limbs(p) * sizeof(uint64_t));
    uint64_t P3_correct[O_MAX * O_MAX * M_VEC_LIMBS_MAX] = {0};
    compute_P3_correct(p, epk, P2_copy, esk->O, P3_correct);

    uint64_t *P3_faulty_full = calloc((size_t)o * o * mvl, sizeof(uint64_t));
    reconstruct_full_P3(p, epk, P3_faulty_full);

    int diffs = 0;
    for (int i = 0; i < o * o * mvl; i++)
      if (P3_faulty_full[i] != P3_correct[i * 0 + (i % (o*o*mvl))]) {} // placeholder guard
    // simpler direct compare against P3_correct laid out (o,o,mvl) same as reconstruct_full_P3
    diffs = 0;
    for (int r = 0; r < o; r++)
      for (int c = 0; c < o; c++)
        for (int l = 0; l < mvl; l++)
          if (P3_faulty_full[mvl*(r*o+c)+l] != P3_correct[mvl*(r*o+c)+l])
            diffs++;
    printf("Faulty P3 vs correct P3: %d/%d limb entries differ (expected: nonzero)\n",
           diffs, o*o*mvl);

    free(P2_copy);
    free(P3_faulty_full);
  }

  // The actual attack: guess row 0, linearize, solve.
  dump_linear_equations(p, epk);

  free(pk);
  free(esk);
  free(epk);
}

int main(void) {
  example_fault_row0_only(NULL);
  return 0;
}