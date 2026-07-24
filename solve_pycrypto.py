"""
Fixed pycryptosat solver for linear GF(2) (XOR) systems.

Key fix vs. the broken version:
  - Uses Solver.add_xor_clause() natively instead of hand-expanding each
    XOR equation into plain CNF clauses. Manual expansion is where sign/
    parity bugs creep in for equations with many variables (common after
    bit-blasting GF(16) -> GF(2)).
  - Builds var_map from variable NAME (string), not from PolyBoRi
    generator object identity, so lookups can never silently miss and
    fall back to a wrong default.
  - Verifies pycryptosat's own solution against every boolean equation
    right after solving, before you spend time on GF(16) recovery.
  - Includes a "known-good assignment" check: if you already have a
    solution from Gaussian elimination, you can feed it in and confirm
    the XOR clauses you built are faithful to the original system.
"""

from pycryptosat import Solver


def parse_linear_anf_equation(poly):
    """
    Given a PolyBoRi/Sage BooleanPolynomial `poly` that is LINEAR
    (degree <= 1), return (list_of_variable_name_strings, rhs_bit).

    ANF linear equation looks like: x1 + x3 + x7 + 1 = 0
      -> variables = ['x1','x3','x7'], rhs = 1   (since x1+x3+x7 = 1)
    or: x2 + x5 = 0
      -> variables = ['x2','x5'], rhs = 0
    """
    if poly.degree() > 1:
        raise ValueError(f"Equation is not linear (degree {poly.degree()}): {poly}")

    variables = []
    rhs = 0
    for monomial in poly.terms():
        vars_in_term = monomial.variables()
        if len(vars_in_term) == 0:
            # constant term (the "1")
            rhs ^= 1
        elif len(vars_in_term) == 1:
            variables.append(str(vars_in_term[0]))
        else:
            # shouldn't happen if degree <= 1, but guard anyway
            raise ValueError(f"Unexpected higher-degree term in linear check: {monomial}")
    return variables, rhs


def build_var_map(boolean_eqs):
    """
    Collect every variable name used across all equations and assign
    each a stable 1-indexed pycryptosat variable id (pycryptosat/CryptoMiniSat
    variables must be indexed starting at 1, not 0).
    """
    names = set()
    for f in boolean_eqs:
        for v in f.variables():
            names.add(str(v))
    # sort for determinism/reproducibility across runs
    sorted_names = sorted(names)
    var_map = {name: i + 1 for i, name in enumerate(sorted_names)}  # 1-indexed
    return var_map


def solve_with_pycryptosat(boolean_eqs, verbose=True):
    """
    boolean_eqs: list of PolyBoRi/Sage BooleanPolynomial objects, all linear
                 (as confirmed by your bit-blast summary: "Boolean system linear: True").

    Returns: dict {variable_name_string: 0/1} or None if UNSAT.
    """
    var_map = build_var_map(boolean_eqs)
    n_vars = len(var_map)

    s = Solver()
    # NOTE: pycryptosat has no explicit "declare n variables" call in this
    # API version; variables are created on first reference inside
    # add_clause/add_xor_clause. Both add_clause and add_xor_clause use the
    # SAME 1-indexed convention (verified empirically: solution[i] lines up
    # with variable id i, with solution[0] == None as padding).

    parsed_eqs = []  # keep parsed form around for the post-solve self-check
    for f in boolean_eqs:
        var_names, rhs = parse_linear_anf_equation(f)
        lits = [var_map[name] for name in var_names]  # 1-indexed, matches solve()
        s.add_xor_clause(lits, rhs == 1)
        parsed_eqs.append((var_names, rhs))

    sat, solution = s.solve()

    if verbose:
        print("=" * 40)
        print("PYCRYPTOSAT (native XOR clauses)")
        print("=" * 40)
        print("Variables:", n_vars, "  Equations:", len(boolean_eqs))
        print("Result:", "SAT" if sat else "UNSAT")

    if not sat:
        return None

    # solution[0] is unused padding in pycryptosat; solution[i] corresponds
    # to variable i (1-indexed), i.e. same indexing as var_map values.
    assignment = {}
    for name, idx in var_map.items():
        val = solution[idx]
        assignment[name] = 1 if val else 0

    # --- self-check: does this assignment actually satisfy every equation? ---
    bad = 0
    for i, (var_names, rhs) in enumerate(parsed_eqs):
        parity = 0
        for name in var_names:
            parity ^= assignment[name]
        if parity != rhs:
            bad += 1
            if bad <= 5:
                print(f"  [self-check] VIOLATED eq {i}: vars={var_names} rhs={rhs} "
                      f"got_parity={parity}")

    if verbose:
        if bad == 0:
            print(f"Self-check: all {len(parsed_eqs)} equations satisfied. OK.")
        else:
            print(f"Self-check: {bad}/{len(parsed_eqs)} equations VIOLATED. "
                  f"Encoding is broken.")

    return assignment if bad == 0 else None


def check_known_assignment(boolean_eqs, known_assignment):
    """
    Sanity tool: verify a Gaussian-elimination-derived solution against the
    same boolean_eqs directly (no SAT solver involved at all). If this fails,
    the equations themselves (or your parsing of them) are the problem, not
    pycryptosat. If this PASSES but solve_with_pycryptosat's self-check
    FAILS, the bug is specifically in the encode/solve/decode path.

    known_assignment: dict {variable_name_string: 0/1}
    """
    bad = 0
    for i, f in enumerate(boolean_eqs):
        var_names, rhs = parse_linear_anf_equation(f)
        parity = 0
        for name in var_names:
            parity ^= known_assignment[name]
        if parity != rhs:
            bad += 1
            if bad <= 5:
                print(f"  [known-assignment check] VIOLATED eq {i}: "
                      f"vars={var_names} rhs={rhs} got_parity={parity}")
    print(f"Known assignment satisfies {len(boolean_eqs) - bad}/{len(boolean_eqs)} equations.")
    return bad == 0


if __name__ == "__main__":
    # ---- minimal synthetic self-test (no Sage required) proving the
    #      add_xor_clause encoding itself is correct ----
    # Simulate a tiny linear system by hand, without PolyBoRi:
    #   x1 + x2 + x3 = 1
    #   x2 + x3      = 0
    #   x1           = 1
    # Expected solution: x1=1, x2=x3 (any equal pair satisfying eq1 => x2=x3=0
    # since x1=1 -> x2+x3=0, consistent with eq2). So x1=1,x2=0,x3=0.

    class FakeVar:
        def __init__(self, name):
            self.name = name
        def __str__(self):
            return self.name

    class FakeMono:
        def __init__(self, vars_list):
            self._vars = vars_list
        def variables(self):
            return self._vars

    class FakePoly:
        def __init__(self, terms, degree):
            self._terms = terms
            self._degree = degree
        def degree(self):
            return self._degree
        def terms(self):
            return self._terms
        def variables(self):
            seen = {}
            for t in self._terms:
                for v in t.variables():
                    seen[str(v)] = v
            return list(seen.values())

    x1, x2, x3 = FakeVar("x1"), FakeVar("x2"), FakeVar("x3")
    eq1 = FakePoly([FakeMono([x1]), FakeMono([x2]), FakeMono([x3]), FakeMono([])], 1)  # x1+x2+x3+1=0
    eq2 = FakePoly([FakeMono([x2]), FakeMono([x3])], 1)                                # x2+x3=0
    eq3 = FakePoly([FakeMono([x1]), FakeMono([])], 1)                                  # x1+1=0

    result = solve_with_pycryptosat([eq1, eq2, eq3])
    print("Recovered assignment:", result)
    assert result == {"x1": 1, "x2": 0, "x3": 0}, "Self-test FAILED"
    print("Self-test PASSED: encoding is correct.")