from sage.all import *
import re
import sys
import time


# ============================================================
# PARAMETERS
# ============================================================

v = 78
o = 8

NUM_GF16_VARS = v * o

INPUT_FILE = (
    sys.argv[1]
    if len(sys.argv) > 1
    else "sage/mayo_equations_quadratic.txt"
)

# How many of the parsed equations to actually use (2808 available).
# You only strictly need ~616 independent ones to pin down l(g); the
# rest become consistency checks on g. Using all of them gives more
# redundancy (helps if some subset is rank-deficient) but costs more
# time in the substitution step. Start smaller while debugging.
MAX_EQUATIONS = (
    int(sys.argv[2])
    if len(sys.argv) > 2
    else 10**9
)


# ============================================================
# GF(16)
# ============================================================

PR = PolynomialRing(GF(2), "z")
z = PR.gen()

F = GF(
    2**4,
    name="a",
    modulus=z**4 + z + 1
)

a = F.gen()


# ============================================================
# GF(16) POLYNOMIAL RING (same layout as original script)
# ============================================================

gf_names = [
    "x_%d_%d" % (k, j)
    for k in range(v)
    for j in range(o)
]

R = PolynomialRing(
    F,
    names=gf_names,
    order="degrevlex"
)

Rvars = list(R.gens())

gf_var_map = {}

for k in range(v):
    for j in range(o):
        gf_var_map[(k, j)] = Rvars[k * o + j]


print("========================================")
print("MAYO BILINEAR-STRUCTURE SOLVER")
print("========================================")
print("Input file       :", INPUT_FILE)
print("GF(16) variables :", NUM_GF16_VARS)


# ============================================================
# INT -> GF(16), PARSER  (unchanged from original script)
# ============================================================

def int_to_gf16(c):
    c = int(c) & 0xF
    value = F(0)
    for bit in range(4):
        if (c >> bit) & 1:
            value += a**bit
    return value


VAR_RE = re.compile(r"^x_(\d+)_(\d+)(?:\^(\d+))?$")
COEFF_RE = re.compile(r"^F\(\s*(\d+)\s*\)$")
EQ_PREFIX_RE = re.compile(r"^\s*Eq\s+\d+\s*:\s*", re.IGNORECASE)


def parse_variable_factor(text):
    m = VAR_RE.fullmatch(text.strip())
    if not m:
        raise ValueError("Invalid variable factor: %r" % text)
    k = int(m.group(1))
    j = int(m.group(2))
    exponent = int(m.group(3)) if m.group(3) else 1
    if not (0 <= k < v):
        raise ValueError("Invalid row index: %d" % k)
    if not (0 <= j < o):
        raise ValueError("Invalid column index: %d" % j)
    if exponent < 1:
        raise ValueError("Invalid exponent: %d" % exponent)
    return gf_var_map[(k, j)], exponent


def parse_term(term):
    term = term.strip()
    if not term:
        return R(0)
    factors = [x.strip() for x in term.split("*") if x.strip()]
    coeff = F(1)
    variable_factors = []
    for factor in factors:
        m = COEFF_RE.fullmatch(factor)
        if m:
            coeff *= int_to_gf16(int(m.group(1)))
            continue
        if factor.isdigit():
            coeff *= int_to_gf16(int(factor))
            continue
        variable, exponent = parse_variable_factor(factor)
        variable_factors.append((variable, exponent))
    result = R(coeff)
    for variable, exponent in variable_factors:
        result *= variable**exponent
    return result


def parse_side(side):
    side = side.strip()
    if not side:
        return R(0)
    result = R(0)
    for raw_term in side.split("+"):
        term = raw_term.strip()
        if not term:
            continue
        result += parse_term(term)
    return result


def parse_equation(line):
    line = line.strip()
    line = EQ_PREFIX_RE.sub("", line, count=1)
    if "==" in line:
        lhs_text, rhs_text = line.split("==", 1)
    elif "=" in line:
        lhs_text, rhs_text = line.split("=", 1)
    else:
        raise ValueError("Equation has no '=' or '==': %s" % line)
    lhs = parse_side(lhs_text)
    rhs = parse_side(rhs_text)
    return lhs + rhs   # char 2: lhs=rhs  ->  lhs+rhs=0


# ============================================================
# READ SYSTEM
# ============================================================

gf16_system = []

with open(INPUT_FILE, "r") as fp:
    for line_no, line in enumerate(fp, start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            f = parse_equation(line)
        except Exception as e:
            print("\nERROR parsing line", line_no)
            print(line)
            print("Reason:", e)
            raise
        gf16_system.append(f)
        if len(gf16_system) >= MAX_EQUATIONS:
            break

if not gf16_system:
    raise RuntimeError("No GF(16) equations parsed")

print("Parsed equations :", len(gf16_system))


# ============================================================
# STEP 1: SPLIT VARIABLES INTO g (row 0 of O) AND l (everything else)
#
#   g_vars = x_0_0 .. x_0_{o-1}     (8 unknowns)
#   l_vars = everything else        (616 unknowns)
#
# Each parsed equation is assumed to have the shape:
#
#   eq_i(g, l) = Q_i(g)  +  B_i(g) . l   = 0
#
# where Q_i(g) has degree <= 2 in g only (no l), and each entry of
# B_i(g) (the coefficient multiplying a single l-variable) is affine
# (degree <= 1) in g. There is NO l*l term anywhere, by assumption.
# ============================================================

g_vars = Rvars[0:o]
l_vars = Rvars[o:]
num_l = len(l_vars)

l_index = {v_: i for i, v_ in enumerate(l_vars)}

S = PolynomialRing(F, o, "g")
g = S.gens()
FS = S.fraction_field()

print()
print("g-variables (guessed)   :", o)
print("l-variables (linear-out):", num_l)


# ============================================================
# STEP 2: EXTRACT B0, B1..B8 (numeric GF(16) matrices) AND Q (list
# of small polynomials in S, one per equation)
#
#   B(g)[i,j]  =  B0[i,j]  +  sum_a  g[a] * Ba[a][i,j]
# ============================================================

num_eq = len(gf16_system)

print()
print("Extracting bilinear structure from %d equations..." % num_eq)
t0 = time.time()

B0 = matrix(F, num_eq, num_l, sparse=True)
Ba = [matrix(F, num_eq, num_l, sparse=True) for _ in range(o)]
Q_list = [S(0)] * num_eq

for i, f in enumerate(gf16_system):

    for mon, coeff in zip(f.monomials(), f.coefficients()):

        exps = mon.exponents()[0]
        g_exp = exps[:o]
        l_exp = exps[o:]

        l_nz = [k for k, e in enumerate(l_exp) if e != 0]

        if len(l_nz) == 0:
            # pure function of g -> goes into Q_i(g)
            gm = S(1)
            for aidx in range(o):
                if g_exp[aidx]:
                    gm *= g[aidx] ** g_exp[aidx]
            Q_list[i] += S(coeff) * gm

        elif len(l_nz) == 1 and l_exp[l_nz[0]] == 1 and sum(g_exp) <= 1:
            # single l-variable, coefficient affine in g
            idx = l_nz[0]
            deg1 = [aidx for aidx in range(o) if g_exp[aidx] == 1]
            if not deg1:
                B0[i, idx] += coeff
            else:
                Ba[deg1[0]][i, idx] += coeff

        else:
            raise ValueError(
                "eq %d: term with exponents %s (coeff %s) does not "
                "match the assumed Q(g) + B(g).l structure -- an "
                "l*l or g-degree>1-times-l term was found."
                % (i, exps, coeff)
            )

print("Extraction done in %.2f s" % (time.time() - t0))


# ============================================================
# STEP 3: PICK A NUMERICALLY FULL-RANK BLOCK OF ROWS
#
# We evaluate B(g) at a random point g_test in F^8 and find a set of
# num_l linearly independent rows there. Generically (i.e. for all
# but a low-dimensional set of "bad" g), that same row set stays
# independent as a *symbolic* matrix in g -- this is what lets us
# invert only a 616x616 block instead of eliminating 2808 rows.
# ============================================================

import random


def random_point():
    return [F.random_element() for _ in range(o)]


def eval_B_numeric(gval):
    Bn = matrix(F, num_eq, num_l, sparse=True)
    Bn += B0
    for aidx in range(o):
        if gval[aidx] != F(0):
            Bn += gval[aidx] * Ba[aidx]
    return Bn


print()
print("Searching for a full-rank pivot block...")

pivot_rows = None

for attempt in range(10):

    gval = random_point()
    Bn = eval_B_numeric(gval)

    # pivot rows of Bn = pivot columns of Bn^T
    piv = Bn.transpose().pivots()

    print("  attempt %d: rank = %d / %d" % (attempt, len(piv), num_l))

    if len(piv) == num_l:
        pivot_rows = list(piv)
        break

if pivot_rows is None:
    raise RuntimeError(
        "Could not find %d independent rows out of %d equations -- "
        "either MAX_EQUATIONS is too small, or the l-part of the "
        "system is genuinely rank-deficient (underdetermined) and "
        "needs equations from another (ell,i) block combined in."
        % (num_l, num_eq)
    )

remaining_rows = [i for i in range(num_eq) if i not in set(pivot_rows)]

print("Pivot block size     :", len(pivot_rows))
print("Remaining (check) eqs:", len(remaining_rows))


# ============================================================
# STEP 4: BUILD THE 616x616 SYMBOLIC BLOCK AND SOLVE l(g)
# ============================================================

print()
print("Building symbolic pivot matrix and solving for l(g)...")
t0 = time.time()

Bsub_rows = []
rhs = []

for i in pivot_rows:
    row = []
    for j in range(num_l):
        c0 = B0[i, j]
        entry = S(c0)
        for aidx in range(o):
            ca = Ba[aidx][i, j]
            if ca != F(0):
                entry += S(ca) * g[aidx]
        row.append(FS(entry))
    Bsub_rows.append(row)
    rhs.append(FS(-Q_list[i]))

Bsub = matrix(FS, Bsub_rows)
rhs_vec = vector(FS, rhs)

l_solution = Bsub.solve_right(rhs_vec)

print("Solved l(g) in %.2f s" % (time.time() - t0))


# ============================================================
# STEP 5: SUBSTITUTE INTO REMAINING EQUATIONS -> CONSISTENCY SYSTEM
# PURELY IN g (8 VARIABLES)
# ============================================================

print()
print("Deriving consistency equations from %d leftover rows..."
      % len(remaining_rows))
t0 = time.time()

consistency_polys = []
seen = set()

for count, i in enumerate(remaining_rows):

    expr = FS(Q_list[i])

    for j in range(num_l):
        c0 = B0[i, j]
        ca_nonzero = [(aidx, Ba[aidx][i, j]) for aidx in range(o)
                      if Ba[aidx][i, j] != F(0)]
        if c0 == F(0) and not ca_nonzero:
            continue
        coeff_expr = S(c0)
        for aidx, ca in ca_nonzero:
            coeff_expr += S(ca) * g[aidx]
        expr += FS(coeff_expr) * l_solution[j]

    if expr == 0:
        continue

    poly = expr.numerator()   # clear denominator; assumes denom != 0
    poly = S(poly)

    key = poly.monomials(), tuple(poly.coefficients())
    if poly not in seen:
        seen.add(poly)
        consistency_polys.append(poly)

    if (count + 1) % 200 == 0:
        print("  processed %d / %d" % (count + 1, len(remaining_rows)))

print("Consistency extraction done in %.2f s" % (time.time() - t0))
print("Distinct nonzero consistency polynomials:", len(consistency_polys))


# ============================================================
# STEP 6: SOLVE THE SMALL 8-VARIABLE SYSTEM
#
# Add field equations g_a^16 - g_a = 0 to force solutions to lie in
# GF(16) (not an extension), then Groebner basis + variety().
# ============================================================

field_eqs = [g[aidx]**16 - g[aidx] for aidx in range(o)]

small_ideal = S.ideal(consistency_polys + field_eqs)

print()
print("========================================")
print("Groebner basis on 8-variable system")
print("========================================")

t0 = time.time()
Gsmall = small_ideal.groebner_basis()
print("Done in %.2f s, basis size %d" % (time.time() - t0, len(Gsmall)))

if S(1) in Gsmall:
    print("System is UNSAT (no solution for row 0 / g)")
    sys.exit(0)

print("Ideal does not contain 1 -- attempting to enumerate solutions")

try:
    sols = small_ideal.variety()
    print("Found %d solution(s) for g:" % len(sols))
    for sol in sols:
        gvals = [sol[g[aidx]] for aidx in range(o)]
        print("  g =", gvals)
except Exception as e:
    print("variety() failed (%s) -- falling back to printing the "
          "Groebner basis; solve the small system by hand / with an "
          "MQ solver on just 8 variables:" % e)
    for poly in Gsmall:
        print(" ", poly)
    sols = []


# ============================================================
# STEP 7: BACK-SUBSTITUTE g INTO l(g) TO RECOVER FULL SOLUTION
# ============================================================

for sol in sols:

    gvals = [sol[g[aidx]] for aidx in range(o)]

    full = {}
    for aidx in range(o):
        full[g_vars[aidx]] = gvals[aidx]

    ok = True
    for j in range(num_l):
        expr = l_solution[j]
        denom = expr.denominator()(*gvals)
        if denom == 0:
            print("  (degenerate point: denominator vanishes for this "
                  "g, skipping -- would need a perturbed pivot block)")
            ok = False
            break
        numer = expr.numerator()(*gvals)
        full[l_vars[j]] = numer / denom

    if not ok:
        continue

    print()
    print("Candidate full solution (g followed by l):")
    for k in range(v):
        for j in range(o):
            var = gf_var_map[(k, j)]
            print("  %s = %s" % (var, full[var]))