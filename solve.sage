# -------------------------------------------------------
# Construct MAYO's GF(16)
#
# GF(16) = GF(2)[a] / (a^4 + a + 1)
# -------------------------------------------------------

B.<z> = PolynomialRing(GF(2))

K.<a> = GF(2**4, modulus=z^4 + z + 1)


# Convert MAYO nibble 0..15 into the corresponding
# polynomial-basis GF(16) element.
def F(n):

    n = int(n)

    return (
        ((n >> 0) & 1)
        + ((n >> 1) & 1) * a
        + ((n >> 2) & 1) * a^2
        + ((n >> 3) & 1) * a^3
    )

v = 78
o = 8

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

X = {}

idx = 0

for k in range(v):
    for j in range(o):
        X[(k,j)] = gens[idx]
        idx += 1

# Make names such as x_0_0 available to eval()
env = {
    str(g): g
    for g in gens
}

env["F"] = F

equations = []

with open("mayo_equations.txt", "r") as fp:

    for line in fp:

        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        lhs, rhs = line.split("==")

        lhs = lhs.strip()

        f = eval(lhs, {"__builtins__": {}}, env)

        equations.append(R(f))


print("Loaded equations:", len(equations))
print("Variables:", R.ngens())

from collections import Counter

degrees = Counter(f.degree() for f in equations)

print("Degree distribution:")

for d in sorted(degrees):
    print("degree", d, ":", degrees[d])

    
linear_eqs = [
    f for f in equations
    if f.degree() <= 1
]

quadratic_eqs = [
    f for f in equations
    if f.degree() == 2
]

print("Linear equations:", len(linear_eqs))
print("Quadratic equations:", len(quadratic_eqs))

quad_mons = sorted(
    set(
        mon
        for f in equations
        for mon in f.monomials()
        if mon.degree() == 2
    ),
    key=str
)

CQ = Matrix(
    K,
    len(equations),
    len(quad_mons),
    lambda i,j:
        equations[i].monomial_coefficient(quad_mons[j])
)

print("Quadratic coefficient matrix:")
print("rows =", CQ.nrows())
print("cols =", CQ.ncols())
print("rank =", CQ.rank())

N = CQ.left_kernel()

print("Quadratic-cancelling combinations:", N.dimension())

derived_linear = []

for lam in N.basis():
    g = sum(
        lam[i] * equations[i]
        for i in range(len(equations))
    )

    if g != 0 and g.degree() <= 1:
        derived_linear.append(g)

print("Derived linear equations:",
      len(derived_linear))


# -------------------------------------------------------
# Analyze derived linear equations
# -------------------------------------------------------



nvars = R.ngens()

AL = Matrix(
    K,
    len(derived_linear),
    nvars,
    lambda i, j:
        derived_linear[i].monomial_coefficient(gens[j])
)

# f = linear_part + constant = 0
# In characteristic 2, -constant == constant,
# but write the general expression anyway.
bL = vector(
    K,
    [
        -derived_linear[i].constant_coefficient()
        for i in range(len(derived_linear))
    ]
)

print("\nDerived linear system:")
print("Equations:", AL.nrows())
print("Variables:", AL.ncols())
print("Coefficient rank:", AL.rank())

Aug = AL.augment(bL.column())

print("Augmented rank:", Aug.rank())

if AL.rank() != Aug.rank():
    print("ERROR: derived linear system is inconsistent")
else:
    print(
        "Independent linear constraints:",
        AL.rank()
    )
    print(
        "Linear solution-space dimension:",
        nvars - AL.rank()
    )


print("\nFirst-row variables in derived equations:")

for j in range(o):
    var = X[(0,j)]

    count = sum(
        1
        for g in derived_linear
        if g.monomial_coefficient(var) != 0
    )

    print(var, "appears in", count, "derived equations")



first_row_vars = set(
    X[(0,j)]
    for j in range(o)
)

with_first_row = []
without_first_row = []

for mon in quad_mons:

    support = set()

    for var in gens:
        if mon.degree(var) > 0:
            support.add(var)

    if support & first_row_vars:
        with_first_row.append(mon)
    else:
        without_first_row.append(mon)

print("\nQuadratic monomial structure:")
print("Total:", len(quad_mons))
print("Containing first-row variable:",
      len(with_first_row))
print("WITHOUT first-row variable:",
      len(without_first_row))


qq_mons = []
qy_mons = []
other_mons = []

first_row_vars = set(X[(0,j)] for j in range(o))

for mon in quad_mons:

    support = []

    for var in gens:
        d = mon.degree(var)
        for _ in range(d):
            support.append(var)

    if all(var in first_row_vars for var in support):
        qq_mons.append(mon)

    elif any(var in first_row_vars for var in support):
        qy_mons.append(mon)

    else:
        other_mons.append(mon)

print("\nQuadratic classification:")
print("q*q monomials:", len(qq_mons))
print("q*y monomials:", len(qy_mons))
print("y*y monomials:", len(other_mons))


print("\nq*q monomials:")
for mon in qq_mons:
    print(mon)