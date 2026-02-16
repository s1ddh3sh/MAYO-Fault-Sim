#include <mayo.h>
#include <mem.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// ---------------- GF(16) arithmetic ----------------

static unsigned char gf_add(unsigned char a, unsigned char b)
{
    return (a ^ b) & 0xF;
}

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

static unsigned char gf_inv(unsigned char a)
{
    unsigned char a2 = mul_f(a, a);
    unsigned char a4 = mul_f(a2, a2);
    unsigned char a8 = mul_f(a4, a4);
    unsigned char a6 = mul_f(a2, a4);
    unsigned char a14 = mul_f(a8, a6);

    return a14;
}

// Gaussian Elimination GF(16)

static void solve_linear_system(
    unsigned char *A, // o x v  (row-major)
    unsigned char *b, // o
    unsigned char *x, // v (output)
    int o,
    int v)
{
    unsigned char M[o][v + 1];

    // Build augmented matrix
    for (int i = 0; i < o; i++)
    {
        for (int j = 0; j < v; j++)
            M[i][j] = A[i * v + j];
        M[i][v] = b[i];
    }

    int rank = 0;

    for (int col = 0; col < v && rank < o; col++)
    {
        int pivot = -1;
        for (int row = rank; row < o; row++)
        {
            if (M[row][col] != 0)
            {
                pivot = row;
                break;
            }
        }

        if (pivot == -1)
            continue;

        // swap
        if (pivot != rank)
        {
            for (int j = 0; j <= v; j++)
            {
                unsigned char tmp = M[pivot][j];
                M[pivot][j] = M[rank][j];
                M[rank][j] = tmp;
            }
        }

        // normalize
        unsigned char inv = gf_inv(M[rank][col]);
        for (int j = col; j <= v; j++)
            M[rank][j] = mul_f(M[rank][j], inv);

        // eliminate
        for (int row = 0; row < o; row++)
        {
            if (row != rank && M[row][col] != 0)
            {
                unsigned char factor = M[row][col];
                for (int j = col; j <= v; j++)
                {
                    M[row][j] =
                        gf_add(M[row][j],
                               mul_f(factor, M[rank][j]));
                }
            }
        }

        rank++;
    }

    // Set free variables = 0
    for (int j = 0; j < v; j++)
        x[j] = 0;

    // Back substitute
    for (int i = 0; i < rank; i++)
    {
        int lead = -1;
        for (int j = 0; j < v; j++)
            if (M[i][j] == 1)
            {
                lead = j;
                break;
            }

        if (lead != -1)
            x[lead] = M[i][v];
    }

    printf("System rank = %d, solution space dim ≈ %d\n",
           rank, v - rank);
}

// ---------------- Fault Simulation ----------------

static int example_fault_P3_OtP2(const mayo_params_t *p)
{
    printf("Fault sim: P3 = O^T P2\n");

    int v = PARAM_v(p);
    int o = PARAM_o(p);

    unsigned char *pk = calloc(1, PARAM_cpk_bytes(p));
    unsigned char *sk = calloc(1, PARAM_csk_bytes(p));
    sk_t *esk = calloc(1, sizeof(sk_t));
    uint64_t *epk = calloc(1, sizeof(pk_t));

    mayo_keypair(p, pk, sk);
    mayo_expand_sk(p, sk, esk);
    mayo_expand_pk(p, pk, epk);

    uint64_t *P2 = epk + PARAM_P1_limbs(p);
    uint64_t *P3_fault = epk + PARAM_P1_limbs(p) + PARAM_P2_limbs(p);

    printf("\nRecovering oil columns via P2^T x = b\n");

    for (int col = 0; col < o; col++)
    {
        printf("\n--- Oil column %d ---\n", col);

        unsigned char A[o * v];
        unsigned char b_vec[o];
        unsigned char x_sol[v];

        // Build A = P2^T
        for (int i = 0; i < o; i++)
        {
            for (int j = 0; j < v; j++)
            {
                // decode P2(j,i)
                A[i * v + j] =
                    ((unsigned char *)P2)[j * o + i] & 0xF;
            }
        }

        // Build RHS from faulty P3 row
        for (int i = 0; i < o; i++)
        {
            b_vec[i] =
                ((unsigned char *)P3_fault)[col * o + i] & 0xF;
        }

        solve_linear_system(A, b_vec, x_sol, o, v);

        printf("Recovered (one solution):\n");
        for (int i = 0; i < v; i++)
            printf("%x ", x_sol[i]);
        printf("\n");

        printf("Real O column:\n");
        for (int i = 0; i < v; i++)
            printf("%x ", esk->O[i * o + col] & 0xF);
        printf("\n");
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
