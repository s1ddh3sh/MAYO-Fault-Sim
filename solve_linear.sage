# ============================================================
# Solve linear MAYO fault equations:
#
#       P3 = O^T * P2
#
# Unknowns:
#
#       x_k_j = O[k][j]
#
# The input file mayo_equations.txt contains equations like:
#
# F(14)*x_0_0 + F(2)*x_1_0 + ... + F(4) == 0
#
# All arithmetic is over GF(16).
# ============================================================


# -------------------------------------------------------
# 1. Construct MAYO's GF(16)
#
# GF(16) = GF(2)[a] / (a^4 + a + 1)
# -------------------------------------------------------

B.<z> = PolynomialRing(GF(2))

K.<a> = GF(
    2^4,
    modulus=z^4 + z + 1
)


# -------------------------------------------------------
# 2. Convert MAYO nibble 0..15 to GF(16)
#
# MAYO nibble:
#
#   b0 + b1*a + b2*a^2 + b3*a^3
#
# -------------------------------------------------------

def F(n):

    n = int(n)

    return K(
        ((n >> 0) & 1)
        + ((n >> 1) & 1) * a
        + ((n >> 2) & 1) * a^2
        + ((n >> 3) & 1) * a^3
    )


# -------------------------------------------------------
# 3. MAYO parameters
# -------------------------------------------------------

v = 78
o = 8

nvars = v * o

print("v =", v)
print("o =", o)
print("Expected unknowns =", nvars)


# -------------------------------------------------------
# 4. Create variables
#
# x_k_j = O[k][j]
#
# Therefore:
#
# x_0_0 ... x_77_7
#
# Total = 78 * 8 = 624 variables
# -------------------------------------------------------

names = [

    "x_%d_%d" % (k, j)

    for k in range(v)

    for j in range(o)
]


R = PolynomialRing(
    K,
    names=names,
    order='degrevlex'
)


gens = R.gens()


# Map:
#
# X[(k,j)] = x_k_j

X = {}

idx = 0

for k in range(v):

    for j in range(o):

        X[(k,j)] = gens[idx]

        idx += 1


# -------------------------------------------------------
# 5. Environment for eval()
#
# Makes names such as:
#
# x_0_0
# x_45_3
#
# available when parsing the equation file.
# -------------------------------------------------------

env = {

    str(g): g

    for g in gens

}

env["F"] = F


# -------------------------------------------------------
# 6. Load equations
# -------------------------------------------------------

equations = []


with open("mayo_equations_linear.txt", "r") as fp:

    for line_no, line in enumerate(fp, start=1):

        line = line.strip()

        # Ignore empty lines
        if not line:
            continue

        # Ignore comments
        if line.startswith("#"):
            continue

        if "==" not in line:

            print(
                "WARNING: skipping malformed line",
                line_no
            )

            continue


        lhs, rhs = line.split("==", 1)

        lhs = lhs.strip()

        rhs = rhs.strip()


        # Convert:
        #
        # lhs == rhs
        #
        # into:
        #
        # lhs - rhs == 0
        #
        # Usually rhs is already 0.

        lhs_expr = eval(
            lhs,
            {"__builtins__": {}},
            env
        )

        rhs_expr = eval(
            rhs,
            {"__builtins__": {}},
            env
        )


        f = R(lhs_expr - rhs_expr)

        equations.append(f)


print("\nLoaded equations:", len(equations))
print("Variables:", R.ngens())


# -------------------------------------------------------
# 7. Verify all equations are linear
# -------------------------------------------------------

nonlinear = []

for i, f in enumerate(equations):

    if f.degree() > 1:

        nonlinear.append(
            (i, f.degree(), f)
        )


if nonlinear:

    print("\nERROR: Nonlinear equations detected!")

    for i, degree, f in nonlinear[:10]:

        print(
            "Equation",
            i,
            "degree =",
            degree
        )

        print(f)

    raise ValueError(
        "Input system is not completely linear"
    )


print("All equations are linear.")


# -------------------------------------------------------
# 8. Convert polynomial equations into:
#
#              A*x = b
#
#
# Each equation has form:
#
#   a0*x0 + a1*x1 + ... + c = 0
#
# Therefore:
#
#   a0*x0 + a1*x1 + ... = -c
#
#
# Over GF(16), characteristic 2:
#
#              -c = c
#
# but we use -c for generality.
# -------------------------------------------------------

A = Matrix(
    K,
    len(equations),
    nvars,
    lambda i, j:
        equations[i].monomial_coefficient(
            gens[j]
        )
)


b = vector(
    K,

    [
        -f.constant_coefficient()

        for f in equations
    ]
)


print("\nLinear system constructed.")

print(
    "A dimensions:",
    A.nrows(),
    "x",
    A.ncols()
)

print(
    "b dimension:",
    len(b)
)


# -------------------------------------------------------
# 9. Compute rank
# -------------------------------------------------------

print("\nComputing coefficient rank...")

rank_A = A.rank()

print(
    "rank(A) =",
    rank_A
)


# -------------------------------------------------------
# 10. Check consistency / SAT
#
# System:
#
#       A*x = b
#
# is SAT iff:
#
#       rank(A) = rank([A | b])
#
# -------------------------------------------------------

Aug = A.augment(
    b.column()
)

rank_Aug = Aug.rank()


print(
    "rank([A|b]) =",
    rank_Aug
)


print("\n========================================")

if rank_A != rank_Aug:

    print("RESULT: UNSAT")

    print(
        "The linear equations are inconsistent."
    )

    print(
        "rank(A)      =",
        rank_A
    )

    print(
        "rank([A|b])  =",
        rank_Aug
    )

    print("========================================")

    quit()


print("RESULT: SAT")

print(
    "The linear equations are consistent."
)

print("========================================")


# -------------------------------------------------------
# 11. Analyze solution space
# -------------------------------------------------------

kernel_dim = nvars - rank_A


print("\nSystem information:")

print(
    "Number of equations :",
    len(equations)
)

print(
    "Number of unknowns  :",
    nvars
)

print(
    "Rank                :",
    rank_A
)

print(
    "Kernel dimension    :",
    kernel_dim
)


if kernel_dim == 0:

    print(
        "\nThe system has a UNIQUE solution."
    )

else:

    print(
        "\nThe system has multiple solutions."
    )

    print(
        "Solution-space dimension:",
        kernel_dim
    )


# -------------------------------------------------------
# 12. Find one solution
#
# solve_right() returns x such that:
#
#              A*x = b
#
# If A is rectangular but consistent, Sage may require
# solve_right() on the matrix directly.
# -------------------------------------------------------

print("\nSolving linear system...")


try:

    solution = A.solve_right(b)

except ValueError as e:

    print(
        "Could not directly solve system:",
        e
    )

    quit()


print("Solution found.")


# -------------------------------------------------------
# 13. Print recovered O matrix
#
# Variable ordering is:
#
# x_0_0, x_0_1, ..., x_0_7,
# x_1_0, ...
#
# Therefore:
#
# index = k*o + j
# -------------------------------------------------------

print("\nRecovered O matrix:")
print()


for k in range(v):

    row = []

    for j in range(o):

        index = k * o + j

        row.append(
            solution[index]
        )

    print(
        "O[%2d] =" % k,
        row
    )


# -------------------------------------------------------
# 14. Verify solution
#
# Check:
#
#       A*solution == b
#
# -------------------------------------------------------

print("\nVerifying solution...")


if A * solution == b:

    print(
        "PASS: A * solution == b"
    )

else:

    print(
        "ERROR: solution does not satisfy system"
    )


# -------------------------------------------------------
# 15. Final SAT result
# -------------------------------------------------------

print("\n========================================")

print("SAT")

print(
    "rank            =",
    rank_A
)

print(
    "unknowns        =",
    nvars
)

print(
    "free variables  =",
    kernel_dim
)

if kernel_dim == 0:

    print(
        "O is uniquely determined by the fault."
    )

else:

    print(
        "O is NOT uniquely determined."
    )

print("========================================")