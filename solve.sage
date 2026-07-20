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