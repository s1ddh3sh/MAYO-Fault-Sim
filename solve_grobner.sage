
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
    else "mayo_equations_quadratic.txt"
)

BOOLEAN_EQ_FILE = (
    sys.argv[2]
    if len(sys.argv) > 2
    else "boolean_anf_equations_quadratic.txt"
)


# ============================================================
# GF(16)
#
# GF(16) = GF(2)[a] / (a^4 + a + 1)
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
# GF(16) POLYNOMIAL RING
#
# Unknowns are Oil matrix entries:
#
#   x_0_0 ... x_77_7
#
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
        gf_var_map[(k, j)] = (
            Rvars[k * o + j]
        )


print("========================================")
print("GF(16) PARTIAL-FAULT EQUATION CONVERTER")
print("========================================")

print("Input file       :", INPUT_FILE)
print("Output file      :", BOOLEAN_EQ_FILE)
print("GF(16) variables :", NUM_GF16_VARS)
print("GF(16) field     :", F)
print("Modulus          :", F.modulus())


# ============================================================
# INTEGER/NIBBLE -> GF(16)
# ============================================================

def int_to_gf16(c):

    c = int(c) & 0xF

    value = F(0)

    for bit in range(4):

        if (c >> bit) & 1:
            value += a**bit

    return value


# ============================================================
# PARSING GF(16) EQUATIONS
#
# Supports:
#
#   F(3)*x_0_0
#
#   F(3)*x_0_0*x_1_2
#
#   F(3)*x_0_0^2
#
#   x_0_0*x_1_2
#
#   F(7)
#
#   7
#
# Equations:
#
#   lhs == rhs
#
# or
#
#   lhs = rhs
#
# ============================================================

VAR_RE = re.compile(
    r"^x_(\d+)_(\d+)(?:\^(\d+))?$"
)

COEFF_RE = re.compile(
    r"^F\(\s*(\d+)\s*\)$"
)

EQ_PREFIX_RE = re.compile(
    r"^\s*Eq\s+\d+\s*:\s*",
    re.IGNORECASE
)


def parse_variable_factor(text):

    m = VAR_RE.fullmatch(text.strip())

    if not m:
        raise ValueError(
            "Invalid variable factor: %r"
            % text
        )

    k = int(m.group(1))
    j = int(m.group(2))

    exponent = (
        int(m.group(3))
        if m.group(3)
        else 1
    )

    if not (0 <= k < v):
        raise ValueError(
            "Invalid row index: %d"
            % k
        )

    if not (0 <= j < o):
        raise ValueError(
            "Invalid column index: %d"
            % j
        )

    if exponent < 1:
        raise ValueError(
            "Invalid exponent: %d"
            % exponent
        )

    return (
        gf_var_map[(k, j)],
        exponent
    )


def parse_term(term):
    """
    Examples:

        F(5)*x_0_0

        F(5)*x_0_0*x_3_2

        F(5)*x_0_0^2

        x_0_0*x_3_2

        F(5)

        5

    Returns a Sage polynomial in R.
    """

    term = term.strip()

    if not term:
        return R(0)

    factors = [
        x.strip()
        for x in term.split("*")
        if x.strip()
    ]

    coeff = F(1)

    variable_factors = []

    for factor in factors:

        # ----------------------------------------
        # F(c)
        # ----------------------------------------

        m = COEFF_RE.fullmatch(factor)

        if m:

            coeff *= int_to_gf16(
                int(m.group(1))
            )

            continue

        # ----------------------------------------
        # Bare integer
        # ----------------------------------------

        if factor.isdigit():

            coeff *= int_to_gf16(
                int(factor)
            )

            continue

        # ----------------------------------------
        # Variable
        # ----------------------------------------

        variable, exponent = (
            parse_variable_factor(factor)
        )

        variable_factors.append(
            (variable, exponent)
        )

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

    line = EQ_PREFIX_RE.sub(
        "",
        line,
        count=1
    )

    if "==" in line:

        lhs_text, rhs_text = (
            line.split("==", 1)
        )

    elif "=" in line:

        lhs_text, rhs_text = (
            line.split("=", 1)
        )

    else:

        raise ValueError(
            "Equation has no '=' or '==': %s"
            % line
        )

    lhs = parse_side(lhs_text)

    rhs = parse_side(rhs_text)

    # characteristic 2:
    #
    # lhs = rhs
    #
    # becomes
    #
    # lhs + rhs = 0

    return lhs + rhs


# ============================================================
# READ GF(16) SYSTEM
# ============================================================

gf16_system = []


with open(INPUT_FILE, "r") as fp:

    for line_no, line in enumerate(
        fp,
        start=1
    ):

        line = line.strip()

        if (
            not line
            or line.startswith("#")
        ):
            continue

        try:

            f = parse_equation(line)

        except Exception as e:

            print()
            print(
                "ERROR parsing line",
                line_no
            )

            print(line)

            print(
                "Reason:",
                e
            )

            raise

        gf16_system.append(f)


if not gf16_system:

    raise RuntimeError(
        "No GF(16) equations parsed"
    )


gf_degrees = [
    f.total_degree()
    for f in gf16_system
]


print()
print("Parsed GF(16) equations :", len(gf16_system))
print("Maximum GF(16) degree   :", max(gf_degrees))

print(
    "Linear equations        :",
    sum(d <= 1 for d in gf_degrees)
)

print(
    "Quadratic equations     :",
    sum(d == 2 for d in gf_degrees)
)

print(
    "Degree > 2 equations    :",
    sum(d > 2 for d in gf_degrees)
)


# ============================================================
# GF(16) -> GF(2) BIT BLAST
# ============================================================

def gf2m_system_to_gf2(system):

    R_ext = system[0].parent()

    F_ext = R_ext.base_ring()

    m = F_ext.degree()

    vars_ext = list(
        R_ext.gens()
    )

    n = len(vars_ext)

    modulus = F_ext.modulus()

    modulus_coeffs = list(
        modulus.list()
    )

    # ========================================================
    # BOOLEAN VARIABLES
    #
    # x_i_j ->
    #
    #   x_i_j_b0
    #   x_i_j_b1
    #   x_i_j_b2
    #   x_i_j_b3
    # ========================================================

    bool_names = []

    for var in vars_ext:

        for bit in range(m):

            bool_names.append(
                "%s_b%d"
                %
                (
                    str(var),
                    bit
                )
            )

    B = BooleanPolynomialRing(
        n * m,
        names=bool_names
    )

    Bvars = list(
        B.gens()
    )

    var_bits = {}

    index = 0

    for var in vars_ext:

        var_bits[var] = list(
            Bvars[
                index:
                index + m
            ]
        )

        index += m


    # ========================================================
    # FIELD ELEMENT -> BASIS BITS
    # ========================================================

    def field_element_to_bits(c):

        c = F_ext(c)

        poly = c.polynomial()

        coeffs = list(
            poly.list()
        )

        coeffs += (
            [0]
            *
            (m - len(coeffs))
        )

        return [

            B(
                int(coeffs[i])
            )

            for i in range(m)

        ]


    # ========================================================
    # SYMBOLIC GF(16) MULTIPLICATION
    #
    # This is what allows quadratic terms.
    #
    # If:
    #
    #   X = x0 + x1*a + x2*a² + x3*a³
    #
    #   Y = y0 + y1*a + y2*a² + y3*a³
    #
    # then field_mult() creates Boolean products:
    #
    #   xi*yj
    #
    # and reduces powers of a modulo:
    #
    #   a^4 + a + 1
    # ========================================================

    def field_mult(u_bits, v_bits):

        tmp = [

            B(0)

            for _ in range(
                2 * m - 1
            )

        ]

        # polynomial multiplication

        for i in range(m):

            for j in range(m):

                tmp[i + j] += (
                    u_bits[i]
                    *
                    v_bits[j]
                )

        # modular reduction

        for degree in range(
            2 * m - 2,
            m - 1,
            -1
        ):

            high = tmp[degree]

            if high == B(0):
                continue

            for i in range(m):

                if (
                    int(
                        modulus_coeffs[i]
                    )
                    & 1
                ):

                    tmp[
                        degree - m + i
                    ] += high

            tmp[degree] = B(0)

        return tmp[:m]


    def multiply_by_variable_power(
        term_bits,
        variable,
        exponent
    ):

        result = list(
            term_bits
        )

        bits = var_bits[
            variable
        ]

        for _ in range(exponent):

            result = field_mult(
                result,
                bits
            )

        return result


    # ========================================================
    # CONVERT SYSTEM
    # ========================================================

    boolean_eqs = []

    print()
    print("Bit-blasting GF(16) equations...")

    for eq_index, f in enumerate(system):

        result = [

            B(0)

            for _ in range(m)

        ]

        for mon in f.monomials():

            coeff = F_ext(
                f.monomial_coefficient(
                    mon
                )
            )

            term_bits = (
                field_element_to_bits(
                    coeff
                )
            )

            exponent_tuple = tuple(
                mon.exponents()[0]
            )

            for (
                var_index,
                exponent
            ) in enumerate(
                exponent_tuple
            ):

                exponent = int(exponent)

                if exponent == 0:
                    continue

                variable = (
                    vars_ext[var_index]
                )

                term_bits = (
                    multiply_by_variable_power(
                        term_bits,
                        variable,
                        exponent
                    )
                )

            for bit in range(m):

                result[bit] += (
                    term_bits[bit]
                )

        # One GF(16) equation ->
        # four GF(2) equations.

        for bit in range(m):

            boolean_eqs.append(
                B(result[bit])
            )

    return (
        B,
        var_bits,
        boolean_eqs
    )
# ============================================================
# GROBNER BASIS SOLVER
# ============================================================

print()
print("========================================")
print("Computing Grobner basis over GF(16)")
print("========================================")

t0 = time.time()

I = R.ideal(gf16_system)

print("Ideal created.")
print("Number of generators :", len(gf16_system))

# Change ordering if desired
# R = PolynomialRing(F, names=gf_names, order='lex')

G = I.groebner_basis()

elapsed = time.time() - t0

print()
print("Grobner basis computed.")
print("Basis size :", len(G))
print("Time       : %.3f s" % elapsed)

# --------------------------------------------------------
# Check consistency
# --------------------------------------------------------

if R(1) in G:
    print()
    print("========================================")
    print("UNSAT")
    print("Ideal contains 1")
    print("========================================")
    sys.exit(0)

print()
print("========================================")
print("SAT")
print("========================================")

# --------------------------------------------------------
# Try solving
# --------------------------------------------------------

print()
print("Computing solutions...")

t1 = time.time()

try:

    sols = I.variety()

    elapsed2 = time.time() - t1

    print("Solutions found :", len(sols))
    print("Enumeration time: %.3f s" % elapsed2)

    if len(sols):

        print()
        print("First solution:")

        sol = sols[0]

        for v in R.gens():

            print(v, "=", sol[v])

except Exception as e:

    print()
    print("Could not enumerate solutions.")
    print("Reason:", e)