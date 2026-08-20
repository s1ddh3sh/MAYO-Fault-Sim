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
#define MAX_EQ 6000
#define PK_PRF AES_128_CTR

unsigned char seed_sk_fixed[24] = {
    0x9a, 0xf8, 0xcb, 0x5c, 0xb5, 0xb7, 0xb1, 0x2c, 0x7a, 0x15, 0xac, 0x4e,
    0x0a, 0x26, 0xae, 0x99, 0x6c, 0x6c, 0xf8, 0x9a, 0x95, 0xb5, 0xfb, 0x64};

static unsigned char gf_add(unsigned char a, unsigned char b) {
  return (a ^ b) & 0xF;
}

static unsigned char gf_mul(unsigned char a, unsigned char b) {
  unsigned char p;
  p = (a & 1) * b;
  p ^= (a & 2) * b;
  p ^= (a & 4) * b;
  p ^= (a & 8) * b;

  unsigned char top_p = p & 0xf0;
  unsigned char out = (p ^ (top_p >> 4) ^ (top_p >> 3)) & 0x0f;
  return out;
}
static unsigned char gf_inv(unsigned char a) {
  unsigned char a2 = gf_mul(a, a);
  unsigned char a4 = gf_mul(a2, a2);
  unsigned char a8 = gf_mul(a4, a4);
  unsigned char a6 = gf_mul(a2, a4);
  unsigned char a14 = gf_mul(a8, a6);
  return a14;
}
/*
 * ============================================================
 * Hex parser
 * ============================================================
 */

static int hex2bytes(const char *hex, unsigned char *out, int expected_len) {
  size_t hex_len = strlen(hex);

  if (hex_len != (size_t)(2 * expected_len))
    return -1;

  for (int i = 0; i < expected_len; i++) {

    unsigned int x;

    if (sscanf(hex + 2 * i, "%2x", &x) != 1)
      return -1;

    out[i] = (unsigned char)x;
  }

  return expected_len;
}

static unsigned char extract_m_element(const uint64_t *vec, int ell) {
  int limb = ell / 16;
  int pos = ell % 16;
  uint64_t w = vec[limb];
  int base = pos * 4;
  unsigned char val = 0;
  if (w & (1ULL << (base + 0)))
    val |= 1;
  if (w & (1ULL << (base + 1)))
    val |= 2;
  if (w & (1ULL << (base + 2)))
    val |= 4;
  if (w & (1ULL << (base + 3)))
    val |= 8;
  return val;
}

/*
 * ============================================================
 * Reconstruct full P3 from the compact upper triangular P3
 * ============================================================
 */

static void reconstruct_full_P3(const mayo_params_t *p, const uint64_t *epk,
                                uint64_t *P3_full) {
  int o = PARAM_o(p);

  int m_vec_limbs = PARAM_m_vec_limbs(p);

  const uint64_t *P3_upper = epk + PARAM_P1_limbs(p) + PARAM_P2_limbs(p);

  int idx = 0;

  for (int i = 0; i < o; i++) {

    for (int j = i; j < o; j++) {

      for (int l = 0; l < m_vec_limbs; l++) {

        uint64_t x = P3_upper[idx * m_vec_limbs + l];

        /*
         * (i,j)
         */
        P3_full[(i * o + j) * m_vec_limbs + l] = x;

        /*
         * (j,i)
         */
        if (i != j) {

          P3_full[(j * o + i) * m_vec_limbs + l] = x;
        }
      }

      idx++;
    }
  }
}

/*
 * ============================================================
 * Build the linear system for Attack 2
 * ============================================================
 *
 * After fixing O[0,:], the equation
 *
 * P3(r,j)
 *
 * becomes
 *
 *   sum_k O[k,r] P2[k,j]
 *
 *   +
 *
 *   O[0,r] *
 *       sum_t P1[0,t] O[t,j]
 *
 * which is linear in O[t,j], t >= 1.
 *
 * Unknown layout:
 *
 *   x[(k-1)*o + j] = O[k,j]
 *
 * for
 *
 *   k = 1,...,v-1.
 *
 */

static int build_partial_system(const mayo_params_t *p, const uint64_t *P1,
                                const uint64_t *P2, const uint64_t *P3_full,
                                const unsigned char *row0, unsigned char *A,
                                unsigned char *b) {
  int v = PARAM_v(p);

  int o = PARAM_o(p);

  int m = PARAM_m(p);

  int m_vec_limbs = PARAM_m_vec_limbs(p);

  int unknowns = (v - 1) * o;

  int equations = m * (o * (o + 1) / 2);

  memset(A, 0, (size_t)equations * unknowns);

  memset(b, 0, equations);

  int eq = 0;

  for (int ell = 0; ell < m; ell++) {

    for (int i = 0; i < o; i++) {

      for (int j = i; j < o; j++) {

        /*
         * RHS = observed faulty P3(i,j)
         */
        unsigned char rhs =
            extract_m_element(P3_full + (i * o + j) * m_vec_limbs, ell);

        unsigned char known = 0;

        /*
         * ==================================================
         * First contribution:
         *
         * O[0,i] * P2[0,j]
         *
         * + O[0,i] *
         *       sum_t P1[0,t] O[t,j]
         * ==================================================
         */

        unsigned char P2_0j = extract_m_element(P2 + j * m_vec_limbs, ell);

        /*
         * O[0,i] * P2[0,j]
         */
        known = gf_add(known, gf_mul(row0[i], P2_0j));

        /*
         * t = 0:
         *
         * O[0,i] *
         * P1[0,0] *
         * O[0,j]
         *
         * This is completely known because row 0
         * has been fixed.
         */
        unsigned char P1_00 = extract_m_element(P1, ell);

        known = gf_add(known, gf_mul(P1_00, gf_mul(row0[i], row0[j])));

        /*
         * t >= 1:
         *
         * O[0,i] *
         * P1[0,t] *
         * O[t,j]
         *
         * coefficient of O[t,j] is known.
         */
        for (int t = 1; t < v; t++) {

          unsigned char P1_0t = extract_m_element(P1 + t * m_vec_limbs, ell);

          unsigned char coeff = gf_mul(row0[i], P1_0t);

          int idx = (t - 1) * o + j;

          A[eq * unknowns + idx] = gf_add(A[eq * unknowns + idx], coeff);
        }

        /*
         * ==================================================
         * Second contribution for off-diagonal entries:
         *
         * O[0,j] * P2[0,i]
         *
         * + O[0,j] *
         *       sum_t P1[0,t] O[t,i]
         *
         * This is needed because the public P3 representation
         * is symmetric and only its upper triangle is stored.
         * ==================================================
         */

        if (i != j) {

          unsigned char P2_0i = extract_m_element(P2 + i * m_vec_limbs, ell);

          known = gf_add(known, gf_mul(row0[j], P2_0i));

          /*
           * t = 0
           */
          known = gf_add(known, gf_mul(P1_00, gf_mul(row0[j], row0[i])));

          /*
           * t >= 1
           */
          for (int t = 1; t < v; t++) {

            unsigned char P1_0t = extract_m_element(P1 + t * m_vec_limbs, ell);

            unsigned char coeff = gf_mul(row0[j], P1_0t);

            int idx = (t - 1) * o + i;

            A[eq * unknowns + idx] = gf_add(A[eq * unknowns + idx], coeff);
          }
        }

        /*
         * ==================================================
         * Linear P2 contributions for k >= 1:
         *
         * O[k,i] * P2[k,j]
         *
         * ==================================================
         */

        for (int k = 1; k < v; k++) {

          unsigned char P2_kj =
              extract_m_element(P2 + (k * o + j) * m_vec_limbs, ell);

          int idx = (k - 1) * o + i;

          A[eq * unknowns + idx] = gf_add(A[eq * unknowns + idx], P2_kj);

          /*
           * Symmetric contribution:
           *
           * O[k,j] * P2[k,i]
           */
          if (i != j) {

            unsigned char P2_ki =
                extract_m_element(P2 + (k * o + i) * m_vec_limbs, ell);

            int idx2 = (k - 1) * o + j;

            A[eq * unknowns + idx2] = gf_add(A[eq * unknowns + idx2], P2_ki);
          }
        }

        /*
         * Move known terms to RHS.
         *
         * In GF(16):
         *
         *      a - b = a + b
         */
        b[eq] = gf_add(rhs, known);

        eq++;
      }
    }
  }

  return eq;
}

/*
 * ============================================================
 * Gaussian elimination over GF(16)
 * ============================================================
 *
 * This follows the same approach as the existing complete
 * attack:
 *
 *      solve A x = b
 *
 * and then explicitly check A*x == b.
 *
 * Returns the rank.
 */

static int solve_linear_system(unsigned char *A, unsigned char *b,
                               unsigned char *x, int rows, int cols,
                               unsigned char *M) {
  /*
   * Augmented matrix:
   *
   *      [ A | b ]
   */
  for (int r = 0; r < rows; r++) {

    for (int c = 0; c < cols; c++) {

      M[r * (cols + 1) + c] = A[r * cols + c];
    }

    M[r * (cols + 1) + cols] = b[r];
  }

  int rank = 0;

  for (int col = 0; col < cols && rank < rows; col++) {

    /*
     * Find pivot.
     */
    int pivot = -1;

    for (int r = rank; r < rows; r++) {

      if (M[r * (cols + 1) + col] != 0) {

        pivot = r;
        break;
      }
    }

    if (pivot == -1)
      continue;

    /*
     * Swap pivot row into position.
     */
    if (pivot != rank) {

      for (int c = 0; c <= cols; c++) {

        unsigned char tmp = M[pivot * (cols + 1) + c];

        M[pivot * (cols + 1) + c] = M[rank * (cols + 1) + c];

        M[rank * (cols + 1) + c] = tmp;
      }
    }

    /*
     * Normalize pivot to 1.
     */
    unsigned char pivot_val = M[rank * (cols + 1) + col];

    unsigned char inv = gf_inv(pivot_val);

    for (int c = col; c <= cols; c++) {

      M[rank * (cols + 1) + c] = gf_mul(M[rank * (cols + 1) + c], inv);
    }

    /*
     * Eliminate this column from every other row.
     */
    for (int r = 0; r < rows; r++) {

      if (r == rank)
        continue;

      unsigned char factor = M[r * (cols + 1) + col];

      if (factor == 0)
        continue;

      for (int c = col; c <= cols; c++) {

        M[r * (cols + 1) + c] = gf_add(
            M[r * (cols + 1) + c], gf_mul(factor, M[rank * (cols + 1) + c]));
      }
    }

    rank++;
  }

  /*
   * Check consistency:
   *
   *      0 ... 0 | nonzero
   *
   */
  for (int r = rank; r < rows; r++) {

    int zero_row = 1;

    for (int c = 0; c < cols; c++) {

      if (M[r * (cols + 1) + c] != 0) {

        zero_row = 0;
        break;
      }
    }

    if (zero_row && M[r * (cols + 1) + cols] != 0) {

      /*
       * Inconsistent.
       */
      return -1;
    }
  }

  /*
   * Extract solution.
   *
   * The system should have full column rank for the
   * intended attack.
   */
  memset(x, 0, cols);

  for (int r = 0; r < rank; r++) {

    int pivot_col = -1;

    for (int c = 0; c < cols; c++) {

      if (M[r * (cols + 1) + c] == 1) {

        pivot_col = c;
        break;
      }
    }

    if (pivot_col >= 0) {

      x[pivot_col] = M[r * (cols + 1) + cols];
    }
  }

  return rank;
}

/*
 * ============================================================
 * Check A*x = b
 *
 * Used exactly as an experimental validation check.
 * ============================================================
 */

static int check_solution(const unsigned char *A, const unsigned char *b,
                          const unsigned char *x, int rows, int cols) {
  for (int r = 0; r < rows; r++) {

    unsigned char acc = 0;

    for (int c = 0; c < cols; c++) {

      unsigned char a = A[r * cols + c];

      if (a != 0) {

        acc = gf_add(acc, gf_mul(a, x[c]));
      }
    }

    if (acc != b[r])
      return 0;
  }

  return 1;
}

/*
 * ============================================================
 * Recover O for Attack 2
 *
 * IMPORTANT:
 *
 * The first row is taken from esk->O only for validation.
 *
 * No enumeration is performed.
 *
 * ============================================================
 */

static int recover_O_partial(const mayo_params_t *p, const uint64_t *epk,
                             const sk_t *esk, unsigned char *O_rec) {
  int v = PARAM_v(p);

  int o = PARAM_o(p);

  int m = PARAM_m(p);

  int m_vec_limbs = PARAM_m_vec_limbs(p);

  int unknowns = (v - 1) * o;

  int equations = m * (o * (o + 1) / 2);

  //   printf("\n");
  //   printf("============================================\n");
  //   printf(" Attack 2: Partial-Fault Oil Recovery\n");
  //   printf("============================================\n");

  //   printf("v          = %d\n", v);
  //   printf("o          = %d\n", o);
  //   printf("m          = %d\n", m);
  //   printf("Variables  = %d\n", unknowns);
  //   printf("Equations  = %d\n", equations);

  /*
   * ---------------------------------------------------------
   * P1, P2 and faulty P3.
   * ---------------------------------------------------------
   */

  const uint64_t *P1 = epk;

  const uint64_t *P2 = epk + PARAM_P1_limbs(p);

  uint64_t *P3_full = calloc((size_t)o * o * m_vec_limbs, sizeof(uint64_t));

  if (!P3_full)
    return 0;

  reconstruct_full_P3(p, epk, P3_full);

  /*
   * ---------------------------------------------------------
   * Allocate system.
   * ---------------------------------------------------------
   */

  unsigned char *A = calloc((size_t)equations * unknowns, 1);

  unsigned char *b = calloc(equations, 1);

  unsigned char *x = calloc(unknowns, 1);

  unsigned char *M = calloc((size_t)equations * (unknowns + 1), 1);

  unsigned char *row0 = calloc(o, 1);

  if (!A || !b || !x || !M || !row0) {

    free(P3_full);
    free(A);
    free(b);
    free(x);
    free(M);
    free(row0);

    return 0;
  }

  /*
   * ---------------------------------------------------------
   * IMPORTANT:
   *
   * For validation we use the true O[0,:].
   *
   * This avoids 16^o exhaustive enumeration.
   *
   * The attack itself does NOT have access to this row.
   * ---------------------------------------------------------
   */

  for (int i = 0; i < o; i++) {

    row0[i] = esk->O[i] & 0xF;
  }

  //   printf("\n");
  //   printf("Row-0 guess used for validation:\n");

  //   for (int i = 0; i < o; i++) {

  //     printf("%x ", row0[i]);
  //   }

  //   printf("\n");

  /*
   * ---------------------------------------------------------
   * Construct linearized system.
   * ---------------------------------------------------------
   */

  build_partial_system(p, P1, P2, P3_full, row0, A, b);

  /*
   * ---------------------------------------------------------
   * Solve.
   * ---------------------------------------------------------
   */

  int rank = solve_linear_system(A, b, x, equations, unknowns, M);

  if (rank < 0) {

    // printf("\nRESULT: FAIL -- system is inconsistent\n");

    free(P3_full);
    free(A);
    free(b);
    free(x);
    free(M);
    free(row0);

    return 0;
  }

  //   printf("\n");
  //   printf("Rank = %d\n", rank);

  /*
   * Check Ax=b explicitly.
   */
  int solution_ok = check_solution(A, b, x, equations, unknowns);

  if (!solution_ok) {

    // printf("RESULT: FAIL -- Ax != b\n");

    free(P3_full);
    free(A);
    free(b);
    free(x);
    free(M);
    free(row0);

    return 0;
  }

  //   printf("Linear system: OK\n");

  /*
   * ---------------------------------------------------------
   * Reconstruct complete O.
   *
   * O[0,:] = row0
   *
   * O[k,i] = x[(k-1)*o+i]
   *
   * ---------------------------------------------------------
   */

  memset(O_rec, 0, (size_t)v * o);

  for (int i = 0; i < o; i++) {

    O_rec[i] = row0[i];
  }

  for (int k = 1; k < v; k++) {

    for (int i = 0; i < o; i++) {

      O_rec[k * o + i] = x[(k - 1) * o + i] & 0xF;
    }
  }

  /*
   * ---------------------------------------------------------
   * Compare recovered O with ground truth.
   *
   * Validation only.
   * ---------------------------------------------------------
   */

  int mismatches = 0;

  for (int k = 0; k < v; k++) {

    for (int i = 0; i < o; i++) {

      if ((O_rec[k * o + i] & 0xF) != (esk->O[k * o + i] & 0xF)) {

        mismatches++;
      }
    }
  }

  //   printf("Oil matrix comparison: %s\n", mismatches == 0 ? "PASS" : "FAIL");

  //   printf("Oil mismatches: %d / %d\n", mismatches, v * o);

  /*
   * Print recovered O.
   */
  //   printf("\nRecovered O:\n");

  //   for (int k = 0; k < v; k++) {

  //     printf("O[%d] = [ ", k);

  //     for (int i = 0; i < o; i++) {

  //       printf("%x ", O_rec[k * o + i] & 0xF);
  //     }

  //     printf("]\n");
  //   }
  //   printf("Actual O:\n");

  //   for (int k = 0; k < v; k++) {

  //     printf("O[%d] = [ ", k);

  //     for (int i = 0; i < o; i++) {

  //       printf("%x ", esk->O[k * o + i] & 0xF);
  //     }

  //     printf("]\n");
  //   }
  free(P3_full);
  free(A);
  free(b);
  free(x);
  free(M);
  free(row0);

  return mismatches == 0;
}

/*
 * ============================================================
 * Correct P3 using recovered O
 * ============================================================
 */

static void compute_P3_correct(const mayo_params_t *p, const uint64_t *P1,
                               const uint64_t *P2_input, const unsigned char *O,
                               uint64_t *P3) {
  int m_vec_limbs = PARAM_m_vec_limbs(p);

  int v = PARAM_v(p);

  int o = PARAM_o(p);

  /*
   * Make writable copy of P2.
   */
  uint64_t *P2 = calloc(PARAM_P2_limbs(p), sizeof(uint64_t));

  memcpy(P2, P2_input, PARAM_P2_limbs(p) * sizeof(uint64_t));

  /*
   * Correct MAYO operation:
   *
   *      P2 <- P2 + P1 O
   */
  P1_times_O(p, P1, O, P2);

  /*
   *      P3 = O^T P2
   */
  mul_add_mat_trans_x_m_mat(m_vec_limbs, O, P2, P3, v, o, o);

  free(P2);
}

/*
 * ============================================================
 * Pack bitsliced vectors into compact P3 representation.
 *
 * Same layout as compact MAYO public key.
 * ============================================================
 */

static void pack_m_vecs(const uint64_t *in, unsigned char *out, int vecs,
                        int m) {
  int m_vec_limbs = (m + 15) / 16;

  const unsigned char *_in = (const unsigned char *)in;

  for (int i = 0; i < vecs; i++) {

    memmove(out + i * m / 2, _in + i * m_vec_limbs * sizeof(uint64_t), m / 2);
  }
}

/*
 * ============================================================
 * Rebuild corrected public key from recovered O
 * ============================================================
 */

static unsigned char *rebuild_corrected_pk(const mayo_params_t *p,
                                           const uint64_t *epk,
                                           const unsigned char *faulty_pk,
                                           const unsigned char *O_rec) {
  //   int v = PARAM_v(p);

  int o = PARAM_o(p);

  int m = PARAM_m(p);

  int m_vec_limbs = PARAM_m_vec_limbs(p);

  int pk_seed_bytes = PARAM_pk_seed_bytes(p);

  int P3_limbs = PARAM_P3_limbs(p);

  const uint64_t *P1 = epk;

  const uint64_t *P2 = epk + PARAM_P1_limbs(p);

  /*
   * Full P3 matrix.
   */
  uint64_t *P3 = calloc((size_t)o * o * m_vec_limbs, sizeof(uint64_t));

  /*
   * Recompute correct P3.
   */
  compute_P3_correct(p, P1, P2, O_rec, P3);

  /*
   * Extract upper triangular P3.
   */
  uint64_t *P3_upper = calloc(P3_limbs, sizeof(uint64_t));

  m_upper(p, P3, P3_upper, o);

  /*
   * Construct compact public key.
   */
  unsigned char *corrected_pk = calloc(PARAM_cpk_bytes(p), 1);

  /*
   * Public key seed remains unchanged.
   */
  memcpy(corrected_pk, faulty_pk, pk_seed_bytes);

  /*
   * Encode P3.
   */
  pack_m_vecs(P3_upper, corrected_pk + pk_seed_bytes, P3_limbs / m_vec_limbs,
              m);

  free(P3);
  free(P3_upper);

  return corrected_pk;
}

/*
 * ============================================================
 * Main
 * ============================================================
 */

int main(int argc, char *argv[]) {
#ifdef ENABLE_PARAMS_DYNAMIC
  const mayo_params_t *p = &MAYO_1;
#else
  const mayo_params_t *p = NULL;
#endif

  if (argc != 2) {

    fprintf(stderr, "Usage: %s <faulty_pk_hex>\n", argv[0]);

    return 1;
  }

  int v = PARAM_v(p);

  int o = PARAM_o(p);

  int pk_bytes = PARAM_cpk_bytes(p);

  //   printf("============================================\n"
  //          " MAYO Attack 2 -- Partial Fault Solver\n"
  //          "============================================\n");

  //   printf("v = %d\n", v);

  //   printf("o = %d\n", o);

  /*
   * ---------------------------------------------------------
   * 1. Read faulty compact public key.
   * ---------------------------------------------------------
   */

  unsigned char *faulty_pk = calloc(pk_bytes, 1);

  if (!faulty_pk) {

    // fprintf(stderr, "Allocation failure\n");

    return 1;
  }

  if (hex2bytes(argv[1], faulty_pk, pk_bytes) < 0) {

    // fprintf(stderr, "Invalid faulty public key\n");

    free(faulty_pk);

    return 1;
  }

  //   printf("Faulty public key received: %d bytes\n", pk_bytes);

  /*
   * ---------------------------------------------------------
   * 2. Expand keys.
   *
   * esk is used ONLY as ground truth for validation.
   * ---------------------------------------------------------
   */

  unsigned char *sk = seed_sk_fixed;

  sk_t *esk = calloc(1, sizeof(sk_t));

  uint64_t *epk = calloc(1, sizeof(pk_t));

  if (!esk || !epk) {

    // fprintf(stderr, "Allocation failure\n");

    free(faulty_pk);
    free(esk);
    free(epk);

    return 1;
  }

  /*
   * Expand secret key only for experimental validation.
   */
  mayo_expand_sk(p, sk, esk);

  /*
   * Expand the supplied faulty public key.
   */
  mayo_expand_pk(p, faulty_pk, epk);

  /*
   * ---------------------------------------------------------
   * 3. Recover O using Attack 2.
   *
   * NO 16^o enumeration.
   *
   * Correct row 0 is supplied from esk only to validate
   * the linearized system.
   * ---------------------------------------------------------
   */

  unsigned char *O_rec = calloc((size_t)v * o, 1);

  if (!O_rec) {

    free(faulty_pk);
    free(esk);
    free(epk);

    return 1;
  }

  int recovery_ok = recover_O_partial(p, epk, esk, O_rec);

  if (!recovery_ok) {

    // printf("\n"
    //        "============================================\n"
    //        " ATTACK 2: OIL RECOVERY FAILED\n"
    //        "============================================\n");

    free(faulty_pk);
    free(esk);
    free(epk);
    free(O_rec);
    return 1;
  }

  //   printf("\n"
  //          "============================================\n"
  //          " ATTACK 2: OIL RECOVERY SUCCESS\n"
  //          "============================================\n");

  /*
   * ---------------------------------------------------------
   * 4. Correct the public key.
   * ---------------------------------------------------------
   */

  //   printf("\n"
  //          "Reconstructing correct P3 and public key...\n");

  unsigned char *corrected_pk = rebuild_corrected_pk(p, epk, faulty_pk, O_rec);

  if (!corrected_pk) {

    // fprintf(stderr, "Public-key correction failed\n");

    free(faulty_pk);
    free(esk);
    free(epk);
    free(O_rec);
    return 1;
  }

  /*
   * ---------------------------------------------------------
   * 5. Generate genuine public key for validation.
   *
   * This uses the fixed seed only for experimental
   * comparison.
   * ---------------------------------------------------------
   */

  unsigned char *correct_pk = calloc(pk_bytes, 1);

  if (!correct_pk) {

    free(faulty_pk);
    free(esk);
    free(epk);
    free(O_rec);
    free(corrected_pk);
    return 1;
  }

  mayo_keypair(p, correct_pk, sk);

  /*
   * ---------------------------------------------------------
   * 6. Compare corrected PK with genuine PK.
   * ---------------------------------------------------------
   */

  int pk_match = memcmp(corrected_pk, correct_pk, pk_bytes) == 0;

  //   printf("\n"
  //          "============================================\n"
  //          " Public-Key Correction\n"
  //          "============================================\n");

  //   printf("Corrected PK == genuine PK: %s\n", pk_match ? "PASS" : "FAIL");

  //   if (!pk_match) {

  //     int first_mismatch = -1;

  //     for (int i = 0; i < pk_bytes; i++) {

  //       if (corrected_pk[i] != correct_pk[i]) {

  //         first_mismatch = i;
  //         break;
  //       }
  //     }

  // if (first_mismatch >= 0) {

  //   printf("First mismatch at byte %d\n", first_mismatch);

  //   printf("Corrected = %02x\n", corrected_pk[first_mismatch]);

  //   printf("Genuine   = %02x\n", correct_pk[first_mismatch]);
  // }
  //   }

  int success = recovery_ok && pk_match;

  //   printf("\n"
  //          "============================================\n"
  //          " ATTACK 2 RESULT\n"
  //          "============================================\n");

  //   printf("Oil recovery       : %s\n", recovery_ok ? "PASS" : "FAIL");

  //   printf("Public-key correct : %s\n", pk_match ? "PASS" : "FAIL");

  //   printf("Overall             : %s\n", success ? "SUCCESS" : "FAIL");

  free(faulty_pk);

  free(esk);

  free(epk);

  free(O_rec);

  free(corrected_pk);

  free(correct_pk);
  printf("%d", success);
  return 0;
}   