
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

SUBSTITUTION_FILE = (
    sys.argv[3]
    if len(sys.argv) > 3
    else "mayo_substitutions.txt"
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
    names=gf_names
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

SUB_RE = re.compile(
    r"^\s*x_(\d+)_(\d+)\s*=\s*(.*)$"
)

def parse_gf16_poly(poly_text):
    """
    Converts strings like

        a^3+a+1
        a^2
        a
        1
        0

    into Sage GF(16) elements.
    """

    poly_text = poly_text.strip()

    if poly_text.startswith("(") and poly_text.endswith(")"):
        poly_text = poly_text[1:-1]

    poly_text = poly_text.replace(" ", "")

    if poly_text == "":
        return F(0)

    value = F(0)

    for term in poly_text.split("+"):

        term = term.strip()

        if term == "0":
            continue

        elif term == "1":
            value += F(1)

        elif term == "a":
            value += a

        elif term.startswith("a^"):

            exponent = int(term[2:])

            value += a**exponent

        else:

            raise ValueError(
                "Unknown GF16 coefficient: %s"
                % term
            )

    return value

def split_top_level_plus(expr):
    """
    Split an expression at '+' signs that are NOT inside parentheses.

    Example:

        (a^3+a+1)*x_0_0 + (a+1)*x_0_1 + x_0_2

    becomes

        [
            "(a^3+a+1)*x_0_0",
            "(a+1)*x_0_1",
            "x_0_2"
        ]
    """

    parts = []
    depth = 0
    current = []

    for ch in expr:

        if ch == "(":
            depth += 1
            current.append(ch)

        elif ch == ")":
            depth -= 1
            current.append(ch)

        elif ch == "+" and depth == 0:

            parts.append("".join(current).strip())
            current = []

        else:
            current.append(ch)

    if current:
        parts.append("".join(current).strip())

    return parts


def split_top_level_mul(expr):

    parts = []
    depth = 0
    current = []

    for ch in expr:

        if ch == "(":
            depth += 1
            current.append(ch)

        elif ch == ")":
            depth -= 1
            current.append(ch)

        elif ch == "*" and depth == 0:

            parts.append("".join(current).strip())
            current = []

        else:
            current.append(ch)

    if current:
        parts.append("".join(current).strip())

    return parts


def parse_substitution_rhs(text):

    result = R(0)

    for raw_term in split_top_level_plus(text):

        raw_term = raw_term.strip()

        if raw_term == "":
            continue

        factors = split_top_level_mul(raw_term)

        coeff = F(1)
        vars_part = []

        for f in factors:

            if f.startswith("("):

                coeff *= parse_gf16_poly(f)

            elif f == "a" or f.startswith("a^"):

                coeff *= parse_gf16_poly(f)

            elif f in ("0","1"):

                coeff *= parse_gf16_poly(f)

            elif f.startswith("x_"):

                var, exponent = parse_variable_factor(f)

                vars_part.append(
                    (var, exponent)
                )

            else:

                raise ValueError(
                    "Cannot parse factor: %s"
                    % f
                )

        term = R(coeff)

        for var, exp in vars_part:

            term *= var**exp

        result += term

    return result

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

substitutions = {}

if SUBSTITUTION_FILE is not None:

    print()
    print("Reading substitutions...")

    with open(SUBSTITUTION_FILE) as fp:

        for line in fp:

            line = line.strip()

            if line == "" or line.startswith("#"):
                continue

            m = SUB_RE.match(line)

            if m is None:
                continue

            row = int(m.group(1))
            col = int(m.group(2))

            rhs = m.group(3)

            variable = gf_var_map[(row,col)]

            substitutions[variable] = (
                parse_substitution_rhs(rhs)
            )

    print(
        "Loaded",
        len(substitutions),
        "substitutions."
    )


if substitutions:

    changed = True

    iteration = 0

    while changed:

        iteration += 1

        changed = False

        new_system = []

        for eq in gf16_system:

            new_eq = R(eq.subs(substitutions))

            if new_eq != eq:
                changed = True

            new_system.append(new_eq)

        gf16_system = new_system

        print(
            "Substitution pass",
            iteration,
            "changed =",
            changed
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
# RUN BIT BLAST
# ============================================================

t0 = time.time()

B, var_bits, boolean_eqs = (
    gf2m_system_to_gf2(
        gf16_system
    )
)

elapsed = time.time() - t0


bool_degrees = [

    f.degree()

    for f in boolean_eqs

    if f != 0

]


max_bool_degree = (
    max(bool_degrees)
    if bool_degrees
    else 0
)


num_linear_bool = sum(

    1

    for f in boolean_eqs

    if f != 0
    and f.degree() <= 1

)


num_quadratic_bool = sum(

    1

    for f in boolean_eqs

    if f != 0
    and f.degree() == 2

)


num_zero_bool = sum(

    1

    for f in boolean_eqs

    if f == 0

)


print()
print("========================================")
print("GF(16) -> GF(2) CONVERSION SUMMARY")
print("========================================")

print(
    "GF(16) equations       :",
    len(gf16_system)
)

print(
    "Boolean equations      :",
    len(boolean_eqs)
)

print(
    "Expected               :",
    4 * len(gf16_system)
)

print(
    "Boolean variables      :",
    B.ngens()
)

print(
    "Linear Boolean eqs     :",
    num_linear_bool
)

print(
    "Quadratic Boolean eqs  :",
    num_quadratic_bool
)

print(
    "Zero Boolean eqs       :",
    num_zero_bool
)

print(
    "Maximum Boolean degree :",
    max_bool_degree
)

print(
    "Bit-blast time         : %.3f s"
    % elapsed
)


if len(boolean_eqs) != (
    4 * len(gf16_system)
):

    raise RuntimeError(
        "Incorrect Boolean equation count"
    )


if max_bool_degree > 2:

    print()
    print(
        "WARNING: Boolean degree > 2."
    )

    print(
        "The CNF converter below supports "
        "arbitrary monomial degree, but this "
        "is unexpected for the intended "
        "partial-fault system."
    )


# ============================================================
# WRITE BOOLEAN ANF
# ============================================================

with open(
    BOOLEAN_EQ_FILE,
    "w"
) as fp:

    fp.write(
        "# Boolean ANF equations generated "
        "from partial-fault GF(16) system\n"
    )

    fp.write(
        "# GF(16) modulus: a^4 + a + 1\n"
    )

    fp.write(
        "# Each line is polynomial = 0 "
        "over GF(2)\n"
    )

    fp.write(
        "# Original Boolean variables: %d\n"
        % B.ngens()
    )

    fp.write(
        "# Boolean equations: %d\n"
        % len(boolean_eqs)
    )

    fp.write(
        "# Maximum Boolean degree: %d\n\n"
        % max_bool_degree
    )

    for i, f in enumerate(
        boolean_eqs,
        start=1
    ):

        fp.write(
            "Eq %5d: %s = 0\n"
            %
            (
                i,
                f
            )
        )


print()
print(
    "Written Boolean ANF system:",
    BOOLEAN_EQ_FILE
)