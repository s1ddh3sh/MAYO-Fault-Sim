#!/usr/bin/env python3

"""
Convert Boolean ANF equations containing linear and nonlinear
monomials into CryptoMiniSat DIMACS CNF+XOR.

Example input:

    Eq 1:
      x_0_0_b0*x_1_0_b2
      + x_3_0_b1
      + 1 = 0

Conceptually:

    (x_0_0_b0 AND x_1_0_b2)
        XOR x_3_0_b1
        = 1

A Tseitin variable is introduced:

    t1 <-> x_0_0_b0 AND x_1_0_b2

and the XOR becomes:

    t1 XOR x_3_0_b1 = 1

The AND relation is encoded with normal CNF clauses.

Supports arbitrary square-free Boolean monomial degree.

Usage:

    python3 anf_to_cnf.py \
        boolean_anf_equations_quadratic.txt \
        system_quadratic.cnf

Decode:

    python3 anf_to_cnf.py \
        --decode \
        solution.txt \
        system.cnf.varmap.txt \
        assignment.txt
"""

import re
import sys
from collections import Counter


# ============================================================
# REGEX
# ============================================================

VAR_RE = re.compile(
    r"[A-Za-z_][A-Za-z0-9_]*"
)

EQ_PREFIX_RE = re.compile(
    r"^\s*Eq\s+\d+\s*:\s*",
    re.IGNORECASE
)


# ============================================================
# PARSE ANF EQUATION
# ============================================================

def parse_line(line):

    line = line.strip()

    if (
        not line
        or line.startswith("#")
    ):
        return None

    # Remove:
    #
    #   Eq 123:

    line = EQ_PREFIX_RE.sub(
        "",
        line,
        count=1
    )

    # RHS

    m = re.search(
        r"(==|=)\s*([01])\s*$",
        line
    )

    if not m:

        raise ValueError(
            "Could not find '= 0' or '= 1': "
            + repr(line)
        )

    explicit_rhs = int(
        m.group(2)
    )

    lhs = line[
        :m.start()
    ].strip()

    if "F(" in lhs:

        raise ValueError(
            "GF(16) coefficient F(...) found. "
            "Input must already be bit-blasted."
        )

    # Each monomial represented canonically as
    # tuple(variable names).

    monomials = []

    const_parity = 0

    for raw_term in lhs.split("+"):

        term = raw_term.strip()

        if not term:
            continue

        if term == "0":
            continue

        if term == "1":

            const_parity ^= 1

            continue

        # ----------------------------------------
        # Parse product:
        #
        # x_a*x_b*x_c
        # ----------------------------------------

        raw_factors = [
            x.strip()
            for x in term.split("*")
            if x.strip()
        ]

        if not raw_factors:

            raise ValueError(
                "Empty monomial: %r"
                % term
            )

        factors = []

        for factor in raw_factors:

            if not VAR_RE.fullmatch(
                factor
            ):

                raise ValueError(
                    "Invalid Boolean variable %r "
                    "in monomial %r"
                    %
                    (
                        factor,
                        term
                    )
                )

            factors.append(factor)

        # ----------------------------------------
        # Boolean polynomial ring:
        #
        # x*x = x
        #
        # so remove repeated factors.
        # ----------------------------------------

        factors = sorted(
            set(factors)
        )

        monomials.append(
            tuple(factors)
        )

    # --------------------------------------------
    # ANF cancellation:
    #
    # m + m = 0
    #
    # --------------------------------------------

    counts = Counter(
        monomials
    )

    monomials = [

        mon

        for mon, count
        in counts.items()

        if count % 2 == 1

    ]

    # deterministic ordering

    monomials.sort()

    rhs = (
        explicit_rhs
        ^ const_parity
    )

    return (
        monomials,
        rhs
    )


# ============================================================
# CONVERSION
# ============================================================

def convert(
    input_path,
    output_cnf_path,
    verbose=True
):

    parsed_eqs = []

    original_vars = set()

    # ========================================================
    # PARSE EVERYTHING FIRST
    # ========================================================

    with open(
        input_path,
        "r"
    ) as f:

        for lineno, line in enumerate(
            f,
            start=1
        ):

            try:

                result = parse_line(
                    line
                )

            except ValueError as e:

                # Do NOT silently skip malformed equations.

                raise ValueError(
                    "Parse error at input line %d: %s"
                    %
                    (
                        lineno,
                        e
                    )
                ) from e

            if result is None:
                continue

            monomials, rhs = result

            for mon in monomials:

                for name in mon:

                    original_vars.add(
                        name
                    )

            parsed_eqs.append(
                (
                    monomials,
                    rhs
                )
            )


    # ========================================================
    # ORIGINAL VARIABLE NUMBERING
    #
    # Sort Oil bit variables naturally:
    #
    # x_0_0_b0
    # x_0_0_b1
    # ...
    # ========================================================

    def natural_key(name):

        return [

            int(x)
            if x.isdigit()
            else x

            for x in re.split(
                r"(\d+)",
                name
            )

        ]


    original_vars = sorted(
        original_vars,
        key=natural_key
    )


    var_map = {}

    for name in original_vars:

        var_map[name] = (
            len(var_map) + 1
        )


    n_original_vars = len(
        var_map
    )


    # ========================================================
    # TSEITIN VARIABLES
    #
    # Cache every product so the same Boolean monomial gets
    # exactly one auxiliary variable globally.
    #
    # Example:
    #
    #   x1*x7
    #
    # -> __and_1
    #
    # If x1*x7 occurs in 500 equations, the same Tseitin
    # variable is reused.
    # ========================================================

    product_map = {}

    tseitin_definitions = []

    next_var = (
        n_original_vars + 1
    )


    def get_product_var(
        monomial
    ):

        nonlocal next_var

        monomial = tuple(
            sorted(
                set(monomial)
            )
        )

        if len(monomial) == 1:

            return var_map[
                monomial[0]
            ]

        if monomial in product_map:

            return product_map[
                monomial
            ]

        z = next_var

        next_var += 1

        product_map[
            monomial
        ] = z

        # Direct n-input AND:
        #
        # z <-> x1 & x2 & ... & xn
        #
        # Clauses:
        #
        # (~z OR xi) for each i
        #
        # (z OR ~x1 OR ... OR ~xn)

        input_ids = [

            var_map[name]

            for name in monomial

        ]

        tseitin_definitions.append(
            (
                z,
                input_ids,
                monomial
            )
        )

        return z


    # ========================================================
    # BUILD XOR EQUATIONS
    # ========================================================

    encoded_xors = []

    constant_unsat = False


    for monomials, rhs in parsed_eqs:

        xor_vars = []

        for mon in monomials:

            if len(mon) == 1:

                xor_vars.append(
                    var_map[
                        mon[0]
                    ]
                )

            else:

                xor_vars.append(
                    get_product_var(
                        mon
                    )
                )

        # ----------------------------------------
        # Defensive XOR cancellation.
        # ----------------------------------------

        counts = Counter(
            xor_vars
        )

        xor_vars = [

            var_id

            for var_id, count
            in counts.items()

            if count % 2 == 1

        ]

        xor_vars.sort()

        if not xor_vars:

            # 0 = rhs

            if rhs == 1:

                constant_unsat = True

            # 0 = 0 contributes nothing.

            continue

        encoded_xors.append(
            (
                xor_vars,
                rhs
            )
        )


    # ========================================================
    # BUILD NORMAL CNF CLAUSES FOR AND VARIABLES
    # ========================================================

    cnf_clauses = []


    for (
        z,
        inputs,
        monomial
    ) in tseitin_definitions:

        # z -> xi
        #
        # (~z OR xi)

        for x in inputs:

            cnf_clauses.append(
                [
                    -z,
                    x
                ]
            )

        # x1 & x2 & ... -> z
        #
        # (z OR ~x1 OR ~x2 ...)

        cnf_clauses.append(

            [z]

            +

            [
                -x
                for x in inputs
            ]

        )


    if constant_unsat:

        # Empty clause = UNSAT

        cnf_clauses.append(
            []
        )


    # ========================================================
    # COUNTS
    #
    # In CryptoMiniSat's extended DIMACS, count both ordinary
    # CNF clauses and XOR constraints in the p cnf constraint
    # count.
    # ========================================================

    total_vars = (
        next_var - 1
    )

    total_constraints = (
        len(cnf_clauses)
        +
        len(encoded_xors)
    )


    # ========================================================
    # WRITE CNF+XOR
    # ========================================================

    with open(
        output_cnf_path,
        "w"
    ) as out:

        out.write(
            "c CNF+XOR generated from %s\n"
            % input_path
        )

        out.write(
            "c original variables: %d\n"
            % n_original_vars
        )

        out.write(
            "c Tseitin AND variables: %d\n"
            % len(product_map)
        )

        out.write(
            "c total SAT variables: %d\n"
            % total_vars
        )

        out.write(
            "c XOR constraints: %d\n"
            % len(encoded_xors)
        )

        out.write(
            "c ordinary CNF clauses: %d\n"
            % len(cnf_clauses)
        )

        out.write(
            "p cnf %d %d\n"
            %
            (
                total_vars,
                total_constraints
            )
        )


        # ----------------------------------------------------
        # AND CNF constraints
        # ----------------------------------------------------

        for clause in cnf_clauses:

            if clause:

                out.write(
                    " ".join(
                        str(x)
                        for x in clause
                    )
                    +
                    " 0\n"
                )

            else:

                out.write(
                    "0\n"
                )


        # ----------------------------------------------------
        # XOR constraints
        #
        # CryptoMiniSat convention used by your previous
        # working linear pipeline:
        #
        #   x1 2 3 0
        #
        # -> v1 XOR v2 XOR v3 = 1
        #
        #   x-1 2 3 0
        #
        # -> v1 XOR v2 XOR v3 = 0
        #
        # ----------------------------------------------------

        for lits, rhs in encoded_xors:

            encoded = list(
                lits
            )

            if rhs == 0:

                encoded[0] = (
                    -encoded[0]
                )

            out.write(

                "x"

                +

                " ".join(
                    str(x)
                    for x in encoded
                )

                +

                " 0\n"

            )


    # ========================================================
    # WRITE ORIGINAL VARIABLE MAP
    #
    # Keep this compatible with your old decoder.
    # Only original Oil Boolean variables go here.
    # ========================================================

    map_path = (
        output_cnf_path
        +
        ".varmap.txt"
    )


    with open(
        map_path,
        "w"
    ) as mf:

        for name in original_vars:

            mf.write(
                "%d %s\n"
                %
                (
                    var_map[name],
                    name
                )
            )


    # ========================================================
    # WRITE AUXILIARY MAP
    # ========================================================

    aux_path = (
        output_cnf_path
        +
        ".auxmap.txt"
    )


    with open(
        aux_path,
        "w"
    ) as af:

        for monomial, z in sorted(
            product_map.items(),
            key=lambda item: item[1]
        ):

            af.write(
                "%d %s\n"
                %
                (
                    z,
                    "*".join(
                        monomial
                    )
                )
            )


    if verbose:

        print()
        print("========================================")
        print("ANF -> CNF+XOR SUMMARY")
        print("========================================")

        print(
            "ANF equations          :",
            len(parsed_eqs)
        )

        print(
            "Original Boolean vars  :",
            n_original_vars
        )

        print(
            "Unique nonlinear terms :",
            len(product_map)
        )

        print(
            "Total SAT variables    :",
            total_vars
        )

        print(
            "XOR constraints        :",
            len(encoded_xors)
        )

        print(
            "AND CNF clauses        :",
            len(cnf_clauses)
        )

        print(
            "Total constraints      :",
            total_constraints
        )

        print()
        print(
            "Wrote:",
            output_cnf_path
        )

        print(
            "Variable map:",
            map_path
        )

        print(
            "Auxiliary map:",
            aux_path
        )

        print()
        print("Run:")

        print(
            "cryptominisat5 "
            "--maxmatrixrows 20000 "
            "--maxnummatrixes 10 "
            "%s > solution.txt"
            % output_cnf_path
        )


    return (
        var_map,
        product_map
    )


# ============================================================
# DECODE SOLUTION
# ============================================================

def decode_solution(
    solution_path,
    varmap_path,
    out_path=None
):

    idx_to_name = {}


    with open(
        varmap_path
    ) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            idx, name = (
                line.split(
                    None,
                    1
                )
            )

            idx_to_name[
                int(idx)
            ] = name


    lits = []

    sat = None


    with open(
        solution_path
    ) as f:

        for line in f:

            line = line.strip()

            if line.startswith(
                "s "
            ):

                if (
                    "UNSATISFIABLE"
                    in line
                ):

                    sat = False

                elif (
                    "SATISFIABLE"
                    in line
                ):

                    sat = True

            elif line.startswith(
                "v "
            ):

                lits.extend(

                    int(x)

                    for x
                    in line[2:].split()

                    if x != "0"

                )


    if sat is False:

        print("UNSAT")

        return None


    if sat is None:

        raise RuntimeError(
            "Could not find SAT/UNSAT status "
            "in solver output"
        )


    assignment = {}


    for lit in lits:

        idx = abs(lit)

        # Ignore Tseitin variables.
        #
        # Only IDs present in varmap.txt are
        # original Oil bits.

        if idx in idx_to_name:

            assignment[
                idx_to_name[idx]
            ] = (
                1
                if lit > 0
                else 0
            )


    if out_path:

        def natural_key(name):

            return [

                int(x)
                if x.isdigit()
                else x

                for x in re.split(
                    r"(\d+)",
                    name
                )

            ]


        with open(
            out_path,
            "w"
        ) as f:

            for name in sorted(
                assignment,
                key=natural_key
            ):

                f.write(
                    "%s %d\n"
                    %
                    (
                        name,
                        assignment[name]
                    )
                )


        print(
            "Wrote decoded assignment to",
            out_path
        )


    print(
        "Recovered original Boolean values:",
        len(assignment)
    )


    return assignment


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) < 3:

        print("Usage:")

        print(
            "  Convert:"
        )

        print(
            "    python3 anf_to_cnf.py "
            "<boolean_anf.txt> <system.cnf>"
        )

        print()

        print(
            "  Decode:"
        )

        print(
            "    python3 anf_to_cnf.py "
            "--decode <solution.txt> "
            "<system.cnf.varmap.txt> "
            "[assignment.txt]"
        )

        sys.exit(1)


    if sys.argv[1] == "--decode":

        if len(sys.argv) < 4:

            raise SystemExit(
                "Missing solution or varmap file"
            )

        solution_path = (
            sys.argv[2]
        )

        varmap_path = (
            sys.argv[3]
        )

        out_path = (

            sys.argv[4]

            if len(sys.argv) > 4

            else None

        )

        decode_solution(
            solution_path,
            varmap_path,
            out_path
        )

    else:

        convert(
            sys.argv[1],
            sys.argv[2]
        )