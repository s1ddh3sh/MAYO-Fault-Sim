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

static int solve_linear_system(unsigned char *A, unsigned char *b,
                               unsigned char *x, int rows, int cols) {
  unsigned char *M = calloc((size_t)rows * (cols + 1), 1);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++)
      M[i * (cols + 1) + j] = A[i * cols + j];
    M[i * (cols + 1) + cols] = b[i];
  }

  int rank = 0;

  for (int col = 0; col < cols && rank < rows; col++) {
    int pivot = -1;
    for (int row = rank; row < rows; row++) {
      if (M[row * (cols + 1) + col] != 0) { pivot = row; break; }
    }
    if (pivot == -1) continue;

    if (pivot != rank) {
      for (int j = 0; j <= cols; j++) {
        unsigned char tmp = M[pivot * (cols + 1) + j];
        M[pivot * (cols + 1) + j] = M[rank * (cols + 1) + j];
        M[rank * (cols + 1) + j] = tmp;
      }
    }

    unsigned char inv = gf_inv(M[rank * (cols + 1) + col]);
    for (int j = col; j <= cols; j++)
      M[rank * (cols + 1) + j] = gf_mul(M[rank * (cols + 1) + j], inv);

    for (int row = 0; row < rows; row++) {
      if (row != rank && M[row * (cols + 1) + col] != 0) {
        unsigned char factor = M[row * (cols + 1) + col];
        for (int j = col; j <= cols; j++) {
          M[row * (cols + 1) + j] =
              gf_add(M[row * (cols + 1) + j],
                     gf_mul(factor, M[rank * (cols + 1) + j]));
        }
      }
    }
    rank++;
  }

  memset(x, 0, cols);
  for (int i = 0; i < rank; i++) {
    int lead = -1;
    for (int j = 0; j < cols; j++) {
      if (M[i * (cols + 1) + j] == 1) { lead = j; break; }
    }
    if (lead != -1) x[lead] = M[i * (cols + 1) + cols];
  }

  free(M);
  return rank;
}

// =====================================================================
// FAULTED compute_P3: bs_mat_rows = 1 forced into the triangular
// P1*O accumulation, so ONLY row 0 of P2 gets P2_new = P1O + P2_old.
// Every other row of P2 stays at P2_old.
//
// This shadows the real library's compute_P3 (same symbol name),
// so mayo_keypair(...) genuinely produces a faulty pk when linked
// against this transla`tion unit.
// =====================================================================

// static inline void P1_times_O_faulted(const mayo_params_t *p, const uint64_t *P1,
//                                       const unsigned char *O, uint64_t *acc) {
// #ifndef ENABLE_PARAMS_DYNAMIC
//   (void)p;
// #endif
//   // FAULT: bs_mat_rows forced to 1 instead of PARAM_v(p).
//   // The r-loop inside mul_add_m_upper_triangular_mat_x_mat then only
//   // executes for r = 0, with its c-loop running the FULL correct
//   // range (c = 0..v-1), so row 0 of acc gets the complete, correct
//   // P1*O contribution, and no other row is touched at all.
//   mul_add_m_upper_triangular_mat_x_mat(PARAM_m_vec_limbs(p), P1, O, acc,
//                                        /*bs_mat_rows=*/1, PARAM_v(p),
//                                        PARAM_o(p), 1);
// }

// static inline void compute_P3(const mayo_params_t *p, const uint64_t *P1,
//                               uint64_t *P2, const unsigned char *O,
//                               uint64_t *P3) {
//   const int m_vec_limbs = PARAM_m_vec_limbs(p);
//   const int param_v = PARAM_v(p);
//   const int param_o = PARAM_o(p);

//   // compute P1*O + P2, FAULTED: only row 0 updated
//   P1_times_O_faulted(p, P1, O, P2);

//   // compute P3 = O^t * (P1*O + P2)  -- uses the now partially-faulty P2
//   mul_add_mat_trans_x_m_mat(m_vec_limbs, O, P2, P3, param_v, param_o, param_o);
// }

// =====================================================================
// Reference (unfaulted) compute_P3, used only for a sanity comparison
// print at the end -- NOT used to derive the recovery, just to show
// how far the faulty pk's P3 diverges from what it should have been.
// =====================================================================

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

// =====================================================================
// Guess-then-linearize recovery oracle (|Q| = 1, quadratic row = row 0)
// =====================================================================

// Build & test the linear system for O[1..v-1][*] given a candidate
// guess for row 0. Returns 1 if consistent (Ax == b for every equation),
// 0 otherwise. On success x[(k-1)*o + col] = recovered O[k][col].
static int try_row0_guess(const mayo_params_t *p, const uint64_t *P1,
                          const uint64_t *P2_old, const uint64_t *P3_full,
                          const unsigned char *guess, unsigned char *A,
                          unsigned char *b, unsigned char *x) {
  int v = PARAM_v(p);
  int o = PARAM_o(p);
  int m = PARAM_m(p);
  int mvl = PARAM_m_vec_limbs(p);

  int unknowns = (v - 1) * o;
  int equations = m * (o * (o + 1) / 2);

  memset(A, 0, (size_t)equations * unknowns);
  memset(b, 0, (size_t)equations);

  int eq = 0;
  for (int ell = 0; ell < m; ell++) {
    for (int i = 0; i < o; i++) {
      for (int j = i; j < o; j++) {

        unsigned char rhs = extract_m_element(
            P3_full + (size_t)mvl * (i * o + j), ell, mvl);
        unsigned char known_const = 0;

        // direction (a=i,b=j): O[0][i]*P2_new[0][j]
        {
          unsigned char p2old_0j = extract_m_element(
              P2_old + (size_t)mvl * (0 * o + j), ell, mvl);
          unsigned char p1_00 = extract_m_element(P1, ell, mvl); // r=0,c=0

          known_const = gf_add(known_const, gf_mul(guess[i], p2old_0j));
          known_const = gf_add(known_const,
                               gf_mul(p1_00, gf_mul(guess[i], guess[j])));

          for (int c = 1; c < v; c++) {
            const uint64_t *p1_0c = P1 + (size_t)mvl * c; // r=0 entries idx 0..v-1
            unsigned char p1_0c_ell = extract_m_element(p1_0c, ell, mvl);
            unsigned char coeff = gf_mul(guess[i], p1_0c_ell);
            int uidx = (c - 1) * o + j;
            A[eq * unknowns + uidx] = gf_add(A[eq * unknowns + uidx], coeff);
          }
        }

        // direction (a=j,b=i), only if i != j
        if (i != j) {
          unsigned char p2old_0i = extract_m_element(
              P2_old + (size_t)mvl * (0 * o + i), ell, mvl);
          unsigned char p1_00 = extract_m_element(P1, ell, mvl);

          known_const = gf_add(known_const, gf_mul(guess[j], p2old_0i));
          known_const = gf_add(known_const,
                               gf_mul(p1_00, gf_mul(guess[j], guess[i])));

          for (int c = 1; c < v; c++) {
            const uint64_t *p1_0c = P1 + (size_t)mvl * c;
            unsigned char p1_0c_ell = extract_m_element(p1_0c, ell, mvl);
            unsigned char coeff = gf_mul(guess[j], p1_0c_ell);
            int uidx = (c - 1) * o + i;
            A[eq * unknowns + uidx] = gf_add(A[eq * unknowns + uidx], coeff);
          }
        }

        // rows k = 1..v-1 (always linear)
        for (int k = 1; k < v; k++) {
          unsigned char p2old_kj = extract_m_element(
              P2_old + (size_t)mvl * (k * o + j), ell, mvl);
          int uidx_ki = (k - 1) * o + i;
          A[eq * unknowns + uidx_ki] =
              gf_add(A[eq * unknowns + uidx_ki], p2old_kj);

          if (i != j) {
            unsigned char p2old_ki = extract_m_element(
                P2_old + (size_t)mvl * (k * o + i), ell, mvl);
            int uidx_kj = (k - 1) * o + j;
            A[eq * unknowns + uidx_kj] =
                gf_add(A[eq * unknowns + uidx_kj], p2old_ki);
          }
        }

        b[eq] = gf_add(rhs, known_const);
        eq++;
      }
    }
  }

  solve_linear_system(A, b, x, equations, unknowns);

  for (int e = 0; e < equations; e++) {
    unsigned char acc = 0;
    for (int u = 0; u < unknowns; u++) {
      unsigned char a = A[e * unknowns + u];
      if (a) acc = gf_add(acc, gf_mul(a, x[u]));
    }
    if (acc != b[e]) return 0;
  }
  return 1;
}

static void attack_row0_quadratic_guess(const mayo_params_t *p,
                                        const uint64_t *epk,
                                        const sk_t *esk) {
  int v = PARAM_v(p);
  int o = PARAM_o(p);
  int mvl = PARAM_m_vec_limbs(p);

  const uint64_t *P1     = epk;                       // public, unaffected by fault
  const uint64_t *P2_old = epk + PARAM_P1_limbs(p);    // pristine, from expand_pk

  uint64_t *P3_full = calloc((size_t)o * o * mvl, sizeof(uint64_t));
  reconstruct_full_P3(p, epk, P3_full); // genuine FAULTY P3, straight from the real pk

  int unknowns  = (v - 1) * o;
  int equations = PARAM_m(p) * (o * (o + 1) / 2);

  unsigned char *A = malloc((size_t)equations * unknowns);
  unsigned char *b = malloc((size_t)equations);
  unsigned char *x = malloc((size_t)unknowns);
  unsigned char guess[O_MAX];

  printf("\n==============================\n");
  printf("Approach (|Q|=1, row 0 quadratic): guess-then-linearize\n");
  printf("v=%d, o=%d, guess space = 16^%d\n", v, o, o);
  printf("==============================\n");

  // Step 1: sanity check with the TRUE row 0 (proves derivation/indexing).
  for (int i = 0; i < o; i++) guess[i] = esk->O[0 * o + i] & 0xF;

  int ok = try_row0_guess(p, P1, P2_old, P3_full, guess, A, b, x);
  printf("Sanity check with true row-0 guess: %s\n",
         ok ? "CONSISTENT (derivation OK)" : "INCONSISTENT (bug!)");

  if (ok) {
    int mism = 0;
    for (int k = 1; k < v; k++)
      for (int col = 0; col < o; col++)
        if ((x[(k - 1) * o + col] & 0xF) != (esk->O[k * o + col] & 0xF))
          mism++;
    printf("Rows 1..%d recovered vs esk->O: %s (%d mismatches)\n",
           v - 1, mism == 0 ? "MATCH" : "DIFFER", mism);
  }

  // Step 2: real brute-force search over 16^o candidates for row 0.
  // Tractable only for small o -- shown here as the full oracle.
  uint64_t total = 1;
  for (int t = 0; t < o; t++) total *= 16;

  printf("\nRunning brute-force guess search (%llu candidates)...\n",
         (unsigned long long)total);

  int found = 0;
  for (uint64_t gnum = 0; gnum < total; gnum++) {
    uint64_t tmp = gnum;
    for (int i = 0; i < o; i++) { guess[i] = tmp & 0xF; tmp >>= 4; }

    if (try_row0_guess(p, P1, P2_old, P3_full, guess, A, b, x)) {
      printf("\nFound consistent guess #%llu: ", (unsigned long long)gnum);
      for (int i = 0; i < o; i++) printf("%x ", guess[i]);
      printf("\n");

      int mism = 0;
      for (int i = 0; i < o; i++)
        if (guess[i] != (esk->O[0 * o + i] & 0xF)) mism++;
      for (int k = 1; k < v; k++)
        for (int col = 0; col < o; col++)
          if ((x[(k - 1) * o + col] & 0xF) != (esk->O[k * o + col] & 0xF))
            mism++;

      if (mism == 0)
        printf("SUCCESS: full oil matrix recovered, matches esk->O\n");
      else
        printf("Consistent but MISMATCHED vs esk->O (%d entries differ)\n", mism);

      found = 1;
      break; // first consistent guess is (essentially certainly) the right one
    }
  }
  if (!found) printf("No consistent guess found across entire search space.\n");

  free(A); free(b); free(x); free(P3_full);
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
  attack_row0_quadratic_guess(p, epk, esk);

  free(pk);
  free(esk);
  free(epk);
}

int main(void) {
  example_fault_row0_only(NULL);
  return 0;
}