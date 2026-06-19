#include <mayo.h>
#include <mem.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#define MAX_UNK 7000
#define MAX_EQ 6000

// ---------------- GF(16) arithmetic ----------------

static unsigned char gf_add(unsigned char a, unsigned char b)
{
    return (a ^ b) & 0xF;
}

static unsigned char gf_mul(unsigned char a, unsigned char b)
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

static unsigned char gf_inv(unsigned char a)
{
    unsigned char a2 = gf_mul(a, a);
    unsigned char a4 = gf_mul(a2, a2);
    unsigned char a8 = gf_mul(a4, a4);
    unsigned char a6 = gf_mul(a2, a4);
    unsigned char a14 = gf_mul(a8, a6);

    return a14;
}

// Gaussian Elimination GF(16)

static unsigned char extract_m_element(const uint64_t *vec, int ell, int m_vec_limbs)
{
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

static void reconstruct_full_P3(const mayo_params_t *p, const uint64_t *epk, uint64_t *P3_full)
{
    int o = PARAM_o(p);
    int m_vec_limbs = PARAM_m_vec_limbs(p);

    const uint64_t *P3_upper =
        epk + PARAM_P1_limbs(p) + PARAM_P2_limbs(p);

    int idx = 0;

    for (int r = 0; r < o; r++)
    {
        for (int c = r; c < o; c++)
        {
            const uint64_t *src =
                P3_upper + m_vec_limbs * idx;

            for (int l = 0; l < m_vec_limbs; l++)
            {
                P3_full[m_vec_limbs * (r * o + c) + l] = src[l];
                if (r != c)
                    P3_full[m_vec_limbs * (c * o + r) + l] = src[l];
            }

            idx++;
        }
    }
}

static int solve_linear_system(
    unsigned char *A,
    unsigned char *b,
    unsigned char *x,
    int rows,
    int cols)
{
    unsigned char *M = calloc(rows * (cols + 1), 1);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            M[i * (cols + 1) + j] = A[i * cols + j];
        M[i * (cols + 1) + cols] = b[i];
    }

    int rank = 0;

    for (int col = 0; col < cols && rank < rows; col++)
    {
        int pivot = -1;
        for (int row = rank; row < rows; row++)
        {
            if (M[row * (cols + 1) + col] != 0)
            {
                pivot = row;
                break;
            }
        }

        if (pivot == -1)
            continue;

        if (pivot != rank)
        {
            for (int j = 0; j <= cols; j++)
            {
                unsigned char tmp = M[pivot * (cols + 1) + j];
                M[pivot * (cols + 1) + j] = M[rank * (cols + 1) + j];
                M[rank * (cols + 1) + j] = tmp;
            }
        }

        unsigned char inv = gf_inv(M[rank * (cols + 1) + col]);

        for (int j = col; j <= cols; j++)
            M[rank * (cols + 1) + j] = gf_mul(M[rank * (cols + 1) + j], inv);

        for (int row = 0; row < rows; row++)
        {
            if (row != rank && M[row * (cols + 1) + col] != 0)
            {
                unsigned char factor = M[row * (cols + 1) + col];

                for (int j = col; j <= cols; j++)
                {
                    M[row * (cols + 1) + j] = gf_add(M[row * (cols + 1) + j], gf_mul(factor, M[rank * (cols + 1) + j]));
                }
            }
        }

        rank++;
    }

    memset(x, 0, cols);

    for (int i = 0; i < rank; i++)
    {
        int lead = -1;
        for (int j = 0; j < cols; j++)
        {
            if (M[i * (cols + 1) + j] == 1)
            {
                lead = j;
                break;
            }
        }

        if (lead != -1)
            x[lead] = M[i * (cols + 1) + cols];
    }

    free(M);
    return rank;
}

static void verify_fault_equation(
    const mayo_params_t *p,
    const uint64_t *epk,
    const sk_t *esk)
{
    int v = PARAM_v(p);
    int o = PARAM_o(p);
    int m = PARAM_m(p);
    int m_vec_limbs = PARAM_m_vec_limbs(p);

    uint64_t *P2 =
        (uint64_t *)(epk + PARAM_P1_limbs(p));

    uint64_t *P3_full =
        calloc(o * o * m_vec_limbs, sizeof(uint64_t));

    reconstruct_full_P3(p, epk, P3_full);

    printf("\nVerifying P3 = O^T P2 + P2^T O\n");

    for (int ell = 0; ell < m; ell++)
    {
        for (int i = 0; i < o; i++)
        {
            for (int j = i; j < o; j++)
            {
                unsigned char rhs = 0;

                for (int k = 0; k < v; k++)
                {
                    unsigned char Oki =
                        esk->O[k * o + i] & 0xF;

                    unsigned char Okj =
                        esk->O[k * o + j] & 0xF;

                    const uint64_t *p2_kj =
                        P2 + m_vec_limbs * (k * o + j);

                    const uint64_t *p2_ki =
                        P2 + m_vec_limbs * (k * o + i);

                    unsigned char p2kj =
                        extract_m_element(
                            p2_kj,
                            ell,
                            m_vec_limbs);

                    unsigned char p2ki =
                        extract_m_element(
                            p2_ki,
                            ell,
                            m_vec_limbs);

                    if (i == j)
                    {
                        rhs ^= gf_mul(Oki, p2kj);
                    }
                    else
                    {
                        rhs ^= gf_mul(Oki, p2kj);
                        rhs ^= gf_mul(Okj, p2ki);
                    }
                }

                const uint64_t *p3_entry =
                    P3_full +
                    m_vec_limbs * (i * o + j);

                unsigned char lhs =
                    extract_m_element(
                        p3_entry,
                        ell,
                        m_vec_limbs);

                if (lhs != rhs)
                {
                    printf("Mismatch at poly %d, (%d,%d)\n",
                           ell, i, j);
                    printf("Stored P3 = %x, Computed RHS = %x\n",
                           lhs, rhs);

                    free(P3_full);
                    return;
                }
            }
        }
    }

    printf("Equation holds for all entries.\n");

    free(P3_full);
}

static inline void compute_P3(const mayo_params_t *p, const uint64_t *P1, uint64_t *P2, const unsigned char *O, uint64_t *P3)
{

    const int m_vec_limbs = PARAM_m_vec_limbs(p);
    const int param_v = PARAM_v(p);
    const int param_o = PARAM_o(p);
    // printf("Param_v : %d\n", param_v);
    // printf("Param_o : %d\n", param_o);
    // compute P1*O + P2
    P1_times_O(p, P1, O, P2);

    // compute P3 = O^t * (P1*O + P2)
    mul_add_mat_trans_x_m_mat(m_vec_limbs, O, P2, P3, param_v, param_o, param_o);
}
void m_upper(const mayo_params_t *p, const uint64_t *in, uint64_t *out, int size)
{
#ifndef ENABLE_PARAMS_DYNAMIC
    (void)p;
#endif
    // Look into AVX2'ing this
    const int m_vec_limbs = PARAM_m_vec_limbs(p);
    int m_vecs_stored = 0;
    for (int r = 0; r < size; r++)
    {
        for (int c = r; c < size; c++)
        {
            m_vec_copy(m_vec_limbs, in + m_vec_limbs * (r * size + c), out + m_vec_limbs * m_vecs_stored);
            if (r != c)
            {
                m_vec_add(m_vec_limbs, in + m_vec_limbs * (c * size + r), out + m_vec_limbs * m_vecs_stored);
            }
            m_vecs_stored++;
        }
    }
}

static void pack_m_vecs(const uint64_t *in, unsigned char *out, int vecs, int m){
    const int m_vec_limbs = (m + 15) / 16;
    unsigned char *_in = (unsigned char *) in;
    for (int i = 0; i < vecs; i++)
    {
        memmove(out + (i*m/2), _in + i*m_vec_limbs*sizeof(uint64_t), m/2);
    }
}

// ---------------- Rebuild public key from recovered O ----------------

// Re-derives P3 from the recovered oil secret and the (faulted) epk's P1/P2,
// then repacks cpk exactly like mayo_keypair_compact does after compute_P3,
// so the result can be diffed against the genuine public key.
static void rebuild_pk_from_recovered_O(
    const mayo_params_t *p,
    const uint64_t *epk,
    const unsigned char *recovered_x, // laid out as x[i * v + k], i in [0,o), k in [0,v)
    const unsigned char *real_pk,
    int real_pk_len)
{
    int v = PARAM_v(p);
    int o = PARAM_o(p);
    int m_vec_limbs = PARAM_m_vec_limbs(p);

    int param_pk_seed_bytes = PARAM_pk_seed_bytes(p);
    int param_P1_limbs = PARAM_P1_limbs(p);
    int param_P2_limbs = PARAM_P2_limbs(p);
    int param_P3_limbs = PARAM_P3_limbs(p);

    // Re-layout recovered oil matrix into O[k * o + i] form, matching esk->O / sk_t,
    // which is what compute_P3 / P1_times_O expect.
    unsigned char *O_rec = calloc((size_t)v * o, 1);
    for (int i = 0; i < o; i++)
        for (int k = 0; k < v; k++)
            O_rec[k * o + i] = recovered_x[i * v + k] & 0xF;

    // P1 is untouched by the fault; P2 will be overwritten in place by compute_P3
    // (it computes P1*O + P2 into P2), so work on a fresh copy taken from epk.
    const uint64_t *P1 = epk;
    uint64_t *P2_work = calloc(param_P2_limbs, sizeof(uint64_t));
    memcpy(P2_work, epk + param_P1_limbs, param_P2_limbs * sizeof(uint64_t));

    uint64_t *P3 = calloc((size_t)o * o * m_vec_limbs, sizeof(uint64_t));

    compute_P3(p, P1, P2_work, O_rec, P3);

    uint64_t *P3_upper = calloc(param_P3_limbs, sizeof(uint64_t));
    m_upper(p, P3, P3_upper, o);

    unsigned char *rebuilt_pk = calloc(1, real_pk_len);

    // seed_pk: the fault sim doesn't expose the original seed_pk directly here,
    // so just copy it from the real pk for this check (it's public anyway and
    // untouched by the fault / by O-recovery).
    memcpy(rebuilt_pk, real_pk, param_pk_seed_bytes);

    pack_m_vecs(P3_upper, rebuilt_pk + param_pk_seed_bytes,
                param_P3_limbs / m_vec_limbs, PARAM_m(p));

    printf("\n==============================\n");
    printf("Recomputing public key from recovered O\n");
    printf("==============================\n");

    int diff = memcmp(rebuilt_pk, real_pk, real_pk_len);

    if (diff == 0)
    {
        printf("SUCCESS: rebuilt public key matches the real public key.\n");
    }
    else
    {
        printf("FAIL: rebuilt public key differs from the real public key.\n");

        int first_mismatch = -1;
        for (int i = 0; i < real_pk_len; i++)
        {
            if (rebuilt_pk[i] != real_pk[i])
            {
                first_mismatch = i;
                break;
            }
        }
        printf("First differing byte at offset %d (real=0x%02x, rebuilt=0x%02x)\n",
               first_mismatch,
               real_pk[first_mismatch],
               rebuilt_pk[first_mismatch]);
    }

    free(O_rec);
    free(P2_work);
    free(P3);
    free(P3_upper);
    free(rebuilt_pk);
}

// ---------------- Fault Simulation ----------------

static void example_fault_P3_OtP2(const mayo_params_t *p)
{
    printf("Fault sim: P3 = O^T P2\n");

    int v = PARAM_v(p);
    int o = PARAM_o(p);
    int m = PARAM_m(p);
    int m_vec_limbs = PARAM_m_vec_limbs(p);

    unsigned char *pk = calloc(1, PARAM_cpk_bytes(p));
    unsigned char *sk = calloc(1, PARAM_csk_bytes(p));
    sk_t *esk = calloc(1, sizeof(sk_t));
    uint64_t *epk = calloc(1, sizeof(pk_t));

    mayo_keypair(p, pk, sk);
    mayo_expand_sk(p, sk, esk);
    mayo_expand_pk(p, pk, epk);

    // int o = PARAM_o(p);
    // int m_vec_limbs = PARAM_m_vec_limbs(p);

    uint64_t *P3_upper =
        epk + PARAM_P1_limbs(p) + PARAM_P2_limbs(p);

    printf("P3_upper = 0x%016" PRIx64 "\n", *P3_upper);

    // printf("\nRaw P3_upper limbs (first few entries):\n");

    // for (int i = 0; i < 3; i++) // first 3 upper entries
    // {
    //     printf("Entry %d:\n", i);
    //     for (int l = 0; l < m_vec_limbs; l++)
    //         printf("  limb %d = %016lx\n", l,
    //                P3_upper[i * m_vec_limbs + l]);
    // }

    verify_fault_equation(p, epk, esk);

    uint64_t *P2 = epk + PARAM_P1_limbs(p);
    uint64_t *P3_full = calloc(o * o * m_vec_limbs, sizeof(uint64_t));

    reconstruct_full_P3(p, epk, P3_full);
    printf("P3_full = 0x%016" PRIx64 "\n", *P3_full);

    // printf("\nReconstructed P3_full (0,0) limbs:\n");

    // for (int l = 0; l < m_vec_limbs; l++)
    // {
    //     printf("  limb %d = %016lx\n",
    //            l,
    //            P3_full[l]);
    // }

    int unknowns = v * o;
    int equations = m * (o * (o + 1) / 2);

    printf("SOlving %d equations in %d variables : \n", equations, unknowns);

    unsigned char *A = calloc(equations * unknowns, 1);
    unsigned char *b = calloc(equations, 1);
    unsigned char *x = calloc(unknowns, 1);

    int eq = 0;

    // printf("\nRecovering oil columns via P2^T x = b\n");

    for (int ell = 0; ell < m; ell++)
    {
        for (int i = 0; i < o; i++)
        {
            for (int j = i; j < o; j++)
            {
                const uint64_t *p3_entry =
                    P3_full +
                    m_vec_limbs * (i * o + j);

                b[eq] =
                    extract_m_element(
                        p3_entry,
                        ell,
                        m_vec_limbs);

                for (int k = 0; k < v; k++)
                {
                    int idx_i =
                        i * v + k;
                    int idx_j =
                        j * v + k;

                    const uint64_t *p2_kj =
                        P2 +
                        m_vec_limbs * (k * o + j);

                    const uint64_t *p2_ki =
                        P2 +
                        m_vec_limbs * (k * o + i);

                    A[eq * unknowns + idx_i] =
                        gf_add(
                            A[eq * unknowns + idx_i],
                            extract_m_element(
                                p2_kj,
                                ell,
                                m_vec_limbs));

                    if (i != j)
                        A[eq * unknowns + idx_j] =
                            gf_add(
                                A[eq * unknowns + idx_j],
                                extract_m_element(
                                    p2_ki,
                                    ell,
                                    m_vec_limbs));
                }

                eq++;
            }
        }
    }
    int rank = solve_linear_system(A, b, x, equations, unknowns);

    printf("System rank = %d\n", rank);

    // printf("Rank = %d\n", rank);

    printf("\n==============================\n");
    printf("Recovered Oil Matrix (column-wise)\n");
    printf("==============================\n");

    for (int i = 0; i < o; i++)
    {
        printf("\nColumn %d:\n", i);
        for (int k = 0; k < v; k++)
        {
            printf("%x ", x[i * v + k] & 0xF);
        }
        printf("\n");
    }

    printf("\n==============================\n");
    printf("Real Oil Matrix (column-wise)\n");
    printf("==============================\n");

    for (int i = 0; i < o; i++)
    {
        printf("\nColumn %d:\n", i);
        for (int k = 0; k < v; k++)
        {
            printf("%x ", esk->O[k * o + i] & 0xF);
        }
        printf("\n");
    }

    // Verification
    int ok = 1;

    for (int i = 0; i < o; i++)
    {
        for (int k = 0; k < v; k++)
        {
            unsigned char rec =
                x[i * v + k] & 0xF;

            unsigned char real =
                esk->O[k * o + i] & 0xF;

            if (rec != real)
            {
                ok = 0;
                break;
            }
        }
    }

    printf("\n==============================\n");
    if (ok)
        printf("SUCCESS: Oil fully recovered\n");
    else
        printf("FAIL: Oil mismatch\n");
    printf("==============================\n");

    // Now recompute P3 from the recovered O and rebuild cpk, to check that
    // a fresh keypair_compact-style packing reproduces the real public key.
    rebuild_pk_from_recovered_O(p, epk, x, pk, PARAM_cpk_bytes(p));

    free(pk);
    mayo_secure_free(sk, PARAM_csk_bytes(p));
    free(esk);
    free(epk);
    free(P3_full);
    free(A);
    free(b);
    free(x);
}

int main(void)
{
    example_fault_P3_OtP2(NULL);
}