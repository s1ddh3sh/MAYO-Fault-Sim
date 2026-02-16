#include <mayo.h>
#include <mem.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

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

static unsigned char extract_m_element(
    const uint64_t *vec,
    int ell,
    int m_vec_limbs)
{
    int limb = ell / 16;
    int pos = ell % 16;

    uint64_t word = vec[limb];

    unsigned char val = 0;

    for (int bit = 0; bit < 4; bit++)
    {
        uint64_t mask = 1ULL << (pos + 16 * bit);
        if (word & mask)
            val |= (1 << bit);
    }

    return val & 0xF;
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
            M[rank * (cols + 1) + j] =
                gf_mul(M[rank * (cols + 1) + j], inv);

        for (int row = 0; row < rows; row++)
        {
            if (row != rank &&
                M[row * (cols + 1) + col] != 0)
            {
                unsigned char factor =
                    M[row * (cols + 1) + col];

                for (int j = col; j <= cols; j++)
                {
                    M[row * (cols + 1) + j] =
                        gf_add(M[row * (cols + 1) + j],
                               gf_mul(factor,
                                      M[rank * (cols + 1) + j]));
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

// ---------------- Fault Simulation ----------------

static int example_fault_P3_OtP2(const mayo_params_t *p)
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

    uint64_t *P2 = epk + PARAM_P1_limbs(p);
    uint64_t *P3_fault = epk + PARAM_P1_limbs(p) + PARAM_P2_limbs(p);

    printf("0x%016" PRIx64 "\n", *P3_fault);

    // printf("\nRecovering oil columns via P2^T x = b\n");

    for (int col = 0; col < o; col++)
    {
        printf("\n--- Oil column %d ---\n", col);

        int rows = m * o;
        int cols = v;

        unsigned char *A = calloc(rows * cols, 1);
        unsigned char *b = calloc(rows, 1);
        unsigned char *x = calloc(cols, 1);

        int eq = 0;

        for (int ell = 0; ell < m; ell++)
        {
            for (int j = 0; j < o; j++)
            {
                const uint64_t *P3_entry =
                    P3_fault +
                    m_vec_limbs * (col * o + j);

                b[eq] =
                    extract_m_element(P3_entry,
                                      ell,
                                      m_vec_limbs);

                for (int r = 0; r < v; r++)
                {
                    const uint64_t *P2_entry =
                        P2 +
                        m_vec_limbs * (r * o + j);

                    A[eq * cols + r] =
                        extract_m_element(P2_entry,
                                          ell,
                                          m_vec_limbs);
                }

                eq++;
            }
        }

        int rank =
            solve_linear_system(A, b, x, rows, cols);

        printf("Rank = %d\n", rank);

        printf("Recovered:\n");
        for (int i = 0; i < v; i++)
            printf("%x ", x[i]);
        printf("\n\nReal:\n");
        for (int i = 0; i < v; i++)
            printf("%x ", esk->O[i * o + col] & 0xF);
        printf("\n");

        free(A);
        free(b);
        free(x);
    }

    free(pk);
    mayo_secure_free(sk, PARAM_csk_bytes(p));
    free(esk);
    free(epk);

    return 0;
}

int main(void)
{
    example_fault_P3_OtP2(NULL);
}
