"""
Convert linear GF(2) ANF equations (from boolean_anf_equations.txt, format:
  F(...)-free lines like: x_0_0 + x_1_0 + ... + 1 = 0   /   ... = 0/1
or PolyBoRi-style "+"-joined monomials) into a DIMACS CNF file with native
CryptoMiniSat "x" XOR lines, so it can be fed directly to the cryptominisat5
CLI binary (which DOES expose --maxmatrixrows / --maxnummatrixes, unlike the
Python bindings).

CryptoMiniSat XOR-in-DIMACS convention (verified directly from CryptoMiniSat
source, src/cnf.h::clean_xor_vars_no_prop -- every "x" line starts at
rhs=True, then rhs ^= literal.sign() per literal):
  A line "x1 3 5 0"        means  var1 XOR var3 XOR var5 = 1   (0 negations)
  A line "x-1 3 5 0"       means  var1 XOR var3 XOR var5 = 0   (1 negation)
    (negating any ONE literal in the XOR line flips the RHS parity;
     by convention put the minus sign on the first literal only)

Usage:
    python3 anf_to_dimacs_xor.py boolean_anf_equations.txt system.cnf
    cryptominisat5 --maxmatrixrows 20000 --maxnummatrixes 5 system.cnf > out.txt
"""

import re
import sys


TERM_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)')

EQ_PREFIX_RE = re.compile(r'^\s*Eq\s+\d+\s*:\s*', re.IGNORECASE)

def parse_line(line):
    """
    Parse one GF(2) Boolean equation.

    Supported input examples:

        Eq     1: x_0_0_b1 + x_1_0_b0 = 0
        Eq     2: x_0_0_b1 + x_1_0_b0 + 1 = 0
        x_0_0_b1 + x_1_0_b0 = 1

    Returns:
        (list_of_var_names, rhs_bit)

    or None for blank/comment lines.

    All arithmetic is over GF(2).
    """
    line = line.strip()

    if not line or line.startswith('#'):
        return None

    # ------------------------------------------------------------
    # IMPORTANT:
    # Strip the human-readable equation-number prefix BEFORE
    # parsing the polynomial.
    #
    # Example:
    #
    #   Eq     3: x_0_0_b0 + x_1_0_b1 + 1 = 0
    #
    # becomes:
    #
    #   x_0_0_b0 + x_1_0_b1 + 1 = 0
    #
    # Without this step, "Eq" can be silently parsed as a variable.
    # ------------------------------------------------------------
    line = EQ_PREFIX_RE.sub('', line, count=1)

    # Split off RHS: = 0, = 1, == 0, == 1
    m = re.search(r'(==|=)\s*([01])\s*$', line)

    if not m:
        raise ValueError(
            f"Could not find trailing '= 0' or '= 1' in line: {line!r}"
        )

    explicit_rhs = int(m.group(2))
    lhs = line[:m.start()].strip()

    # This converter expects already bit-blasted GF(2) equations.
    if 'F(' in lhs:
        raise ValueError(
            "This line still has GF(16) 'F(k)*' coefficients -- "
            "pass the bit-blasted boolean_anf_equations.txt, "
            "not the raw GF(16) file: "
            f"{line!r}"
        )

    var_names = []
    const_parity = 0

    for term in lhs.split('+'):
        term = term.strip()

        if term == '':
            continue

        if term == '1':
            const_parity ^= 1

        elif term == '0':
            continue

        else:
            # IMPORTANT:
            # Use fullmatch(), NOT match().
            #
            # match() would incorrectly accept:
            #
            #   "Eq     1: x_0_0_b1"
            #
            # as variable "Eq".
            #
            # fullmatch() requires the complete term to be a valid
            # variable name.
            if not TERM_RE.fullmatch(term):
                raise ValueError(
                    f"Unrecognized term {term!r} in line: {line!r}"
                )

            var_names.append(term)

    # Move constant terms from LHS to RHS.
    #
    # Over GF(2):
    #
    #   x + y + 1 = 0
    #
    # is equivalent to:
    #
    #   x + y = 1
    #
    rhs = explicit_rhs ^ const_parity

    return var_names, rhs

def convert(input_path, output_cnf_path, verbose=True):
    var_map = {}          # name -> 1-indexed CNF var id
    parsed_eqs = []

    with open(input_path, 'r') as f:
        for lineno, line in enumerate(f, 1):
            try:
                result = parse_line(line)
            except ValueError as e:
                print(f"[line {lineno}] SKIPPED (parse error): {e}", file=sys.stderr)
                continue
            if result is None:
                continue
            var_names, rhs = result
            for name in var_names:
                if name not in var_map:
                    var_map[name] = len(var_map) + 1
            parsed_eqs.append((var_names, rhs))

    n_vars = len(var_map)
    n_eqs = len(parsed_eqs)

    with open(output_cnf_path, 'w') as out:
        out.write(f"c CNF+XOR system converted from {input_path}\n")
        out.write(f"c variables: {n_vars}  equations: {n_eqs}\n")
        out.write(f"p cnf {n_vars} {n_eqs}\n")
        for var_names, rhs in parsed_eqs:
            lits = [var_map[name] for name in var_names]
            if not lits:
                # degenerate equation like "1 = 1" (always true) or "1 = 0" (UNSAT)
                if rhs == 0:
                    continue  # trivially satisfied, no clause needed
                else:
                    # unsatisfiable constant equation -> force UNSAT via empty clause
                    out.write("0\n")
                    continue
            # CryptoMiniSat DIMACS "x" line convention (verified directly from
            # src/cnf.h clean_xor_vars_no_prop: rhs starts at True for every
            # "x" line, then rhs ^= literal.sign() for each literal):
            #   zero (or even) negated literals -> encoded equation == 1
            #   one  (or odd)  negated literal  -> encoded equation == 0
            # So to encode "sum of vars == 0" we must negate exactly one
            # literal; to encode "sum of vars == 1" we leave all positive.
            if rhs == 0:
                lits = [-lits[0]] + lits[1:]
            out.write("x" + " ".join(str(l) for l in lits) + " 0\n")

    # write the var_map alongside so you can decode the solver's output back
    # to variable names later
    map_path = output_cnf_path + ".varmap.txt"
    with open(map_path, 'w') as mf:
        for name, idx in sorted(var_map.items(), key=lambda kv: kv[1]):
            mf.write(f"{idx} {name}\n")

    if verbose:
        print(f"Wrote {output_cnf_path}")
        print(f"  variables : {n_vars}")
        print(f"  equations : {n_eqs}")
        print(f"Wrote variable map: {map_path}")
        print()
        print("Run with, e.g.:")
        print(f"  cryptominisat5 --maxmatrixrows 20000 --maxnummatrixes 10 "
              f"{output_cnf_path} > solution.txt")

    return var_map, parsed_eqs


def decode_solution(solution_path, varmap_path, out_path=None):
    """
    Parse cryptominisat5's stdout ('s SATISFIABLE' + 'v ...' lines) using the
    var map produced by convert(), and print/save {var_name: 0/1}.
    """
    idx_to_name = {}
    with open(varmap_path) as f:
        for line in f:
            idx, name = line.split()
            idx_to_name[int(idx)] = name

    lits = []
    sat = None
    with open(solution_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('s '):
                sat = 'SATISFIABLE' in line
            elif line.startswith('v '):
                lits.extend(int(x) for x in line[2:].split() if x != '0')

    if sat is False:
        print("UNSAT")
        return None

    assignment = {}
    for lit in lits:
        idx = abs(lit)
        if idx in idx_to_name:
            assignment[idx_to_name[idx]] = 1 if lit > 0 else 0

    if out_path:
        with open(out_path, 'w') as f:
            for name in sorted(assignment):
                f.write(f"{name} {assignment[name]}\n")
        print(f"Wrote decoded assignment to {out_path}")

    return assignment


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage:")
        print("  Convert:  python3 anf_to_dimacs_xor.py <in.txt> <out.cnf>")
        print("  Decode:   python3 anf_to_dimacs_xor.py --decode <cms_stdout.txt> "
              "<out.cnf.varmap.txt> [assignment_out.txt]")
        sys.exit(1)

    if sys.argv[1] == '--decode':
        sol_path = sys.argv[2]
        varmap_path = sys.argv[3]
        out_path = sys.argv[4] if len(sys.argv) > 4 else None
        assignment = decode_solution(sol_path, varmap_path, out_path)
        if assignment:
            print(f"Recovered {len(assignment)} variable values.")
    else:
        convert(sys.argv[1], sys.argv[2])