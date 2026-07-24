

import re
import sys


# ============================================================
# GF(16)
# ============================================================

PR.<z> = PolynomialRing(GF(2))
F.<a> = GF(2^4, modulus=z^4 + z + 1)


# MAYO-1 dimensions
V = 78
O = 8


# ============================================================
# Read Boolean SAT assignment
# ============================================================

VAR_RE = re.compile(r'^x_(\d+)_(\d+)_b([0-3])$')


def read_assignment(filename):
    """
    Read:

        x_0_0_b0 1
        x_0_0_b1 0
        ...

    Returns:

        assignment[(row, col, bit)] = 0/1
    """

    assignment = {}

    with open(filename, "r") as f:
        for lineno, line in enumerate(f, 1):

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) != 2:
                raise ValueError(
                    f"{filename}:{lineno}: "
                    f"expected '<variable> <0/1>', got: {line!r}"
                )

            name, value = parts

            m = VAR_RE.fullmatch(name)

            if not m:
                raise ValueError(
                    f"{filename}:{lineno}: "
                    f"invalid variable name: {name!r}"
                )

            row = int(m.group(1))
            col = int(m.group(2))
            bit = int(m.group(3))

            value = int(value)

            if value not in (0, 1):
                raise ValueError(
                    f"{filename}:{lineno}: "
                    f"bit value must be 0 or 1, got {value}"
                )

            key = (row, col, bit)

            if key in assignment:
                raise ValueError(
                    f"{filename}:{lineno}: "
                    f"duplicate assignment for {name}"
                )

            assignment[key] = value

    return assignment


# ============================================================
# GF(2)^4 -> GF(16)
# ============================================================

def bits_to_gf16(b0, b1, b2, b3):
    """
    Convert polynomial-basis bits into GF(16):

        b0 + b1*a + b2*a^2 + b3*a^3

    Example:

        [1,0,1,0] -> 1 + a^2
    """

    return (
        F(b0)
        + F(b1) * a
        + F(b2) * a^2
        + F(b3) * a^3
    )


def reconstruct_matrix(assignment, rows=V, cols=O):
    """
    Reconstruct the complete rows x cols GF(16) matrix.
    """

    matrix = []

    missing = []

    for i in range(rows):

        row = []

        for j in range(cols):

            bits = []

            for b in range(4):

                key = (i, j, b)

                if key not in assignment:
                    missing.append(
                        f"x_{i}_{j}_b{b}"
                    )
                    bits.append(0)
                else:
                    bits.append(assignment[key])

            value = bits_to_gf16(
                bits[0],
                bits[1],
                bits[2],
                bits[3]
            )

            row.append(value)

        matrix.append(row)

    if missing:

        print()
        print("ERROR: Missing Boolean variables:")
        print()

        for name in missing[:20]:
            print("  ", name)

        if len(missing) > 20:
            print(
                f"  ... and {len(missing)-20} more"
            )

        raise RuntimeError(
            f"Missing {len(missing)} Boolean assignments"
        )

    return matrix


# ============================================================
# Parse expected O matrix
# ============================================================

O_LINE_RE = re.compile(
    r'^\s*O\[\s*(\d+)\s*\]\s*=\s*\[(.*)\]\s*$'
)


def parse_gf16_element(s):
    """
    Parse Sage-style GF(16) expressions such as:

        0
        1
        a
        a + 1
        a^2 + 1
        a^3 + a^2 + a + 1

    We explicitly parse the expression rather than using eval().
    """

    s = s.strip().replace(" ", "")

    if s == "0":
        return F(0)

    result = F(0)

    for term in s.split("+"):

        if term == "1":
            result += F(1)

        elif term == "a":
            result += a

        elif term == "a^2":
            result += a^2

        elif term == "a^3":
            result += a^3

        else:
            raise ValueError(
                f"Unknown GF(16) term: {term!r} "
                f"in expression {s!r}"
            )

    return result


def read_expected_o(filename, rows=V, cols=O):

    matrix = [None] * rows

    with open(filename, "r") as f:

        for lineno, line in enumerate(f, 1):

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            m = O_LINE_RE.fullmatch(line)

            if not m:
                raise ValueError(
                    f"{filename}:{lineno}: "
                    f"cannot parse line:\n{line}"
                )

            row_idx = int(m.group(1))

            if not (0 <= row_idx < rows):
                raise ValueError(
                    f"Invalid row index O[{row_idx}]"
                )

            entries = [
                x.strip()
                for x in m.group(2).split(",")
            ]

            if len(entries) != cols:
                raise ValueError(
                    f"O[{row_idx}] has {len(entries)} "
                    f"elements, expected {cols}"
                )

            matrix[row_idx] = [
                parse_gf16_element(x)
                for x in entries
            ]

    missing_rows = [
        i for i, row in enumerate(matrix)
        if row is None
    ]

    if missing_rows:
        raise RuntimeError(
            f"Missing expected O rows: {missing_rows}"
        )

    return matrix


# ============================================================
# Printing
# ============================================================

def print_matrix(M, name="O"):

    for i, row in enumerate(M):

        values = ", ".join(
            str(x) for x in row
        )

        print(
            f"{name}[{i:2d}] = [{values}]"
        )


# ============================================================
# Compare matrices
# ============================================================

def compare_matrices(recovered, expected):

    total = 0
    matches = 0
    mismatches = []

    for i in range(V):

        for j in range(O):

            total += 1

            if recovered[i][j] == expected[i][j]:

                matches += 1

            else:

                mismatches.append(
                    (
                        i,
                        j,
                        recovered[i][j],
                        expected[i][j]
                    )
                )

    print()
    print("=" * 70)
    print("COMPARISON RESULT")
    print("=" * 70)

    print(f"Total GF(16) elements : {total}")
    print(f"Matching elements     : {matches}")
    print(f"Mismatching elements  : {len(mismatches)}")

    percentage = 100.0 * matches / total

    print(
        f"Match percentage      : {percentage:.2f}%"
    )

    print()

    if not mismatches:

        print("========================================")
        print("SUCCESS: MATRICES MATCH EXACTLY")
        print("Recovered O == Expected O")
        print("========================================")

    else:

        print("========================================")
        print("MISMATCH: MATRICES ARE DIFFERENT")
        print("========================================")

        print()
        print("Mismatching elements:")
        print()

        for i, j, rec, exp in mismatches:

            print(
                f"O[{i:2d}][{j}]  "
                f"recovered = {str(rec):25s}  "
                f"expected = {exp}"
            )

    return mismatches


# ============================================================
# Detailed bit verification
# ============================================================

def show_first_elements(assignment, recovered, expected, count=10):
    """
    Useful sanity check for bit ordering.
    """

    print()
    print("=" * 70)
    print("FIRST ELEMENTS -- BIT RECONSTRUCTION CHECK")
    print("=" * 70)

    shown = 0

    for i in range(V):

        for j in range(O):

            bits = [
                assignment[(i, j, b)]
                for b in range(4)
            ]

            nibble = (
                bits[0]
                | (bits[1] << 1)
                | (bits[2] << 2)
                | (bits[3] << 3)
            )

            print(
                f"O[{i:2d}][{j}] : "
                f"bits b0..b3 = {bits}   "
                f"nibble = 0x{nibble:X}   "
                f"GF16 = {str(recovered[i][j]):25s}   "
                f"expected = {expected[i][j]}"
            )

            shown += 1

            if shown >= count:
                return


# ============================================================
# Main
# ============================================================

def main():

    if len(sys.argv) != 3:

        print(
            "Usage:"
        )

        print(
            "  sage compare_o.py "
            "assignment.txt expected_o.txt"
        )

        sys.exit(1)

    assignment_file = sys.argv[1]
    expected_file = sys.argv[2]

    print("=" * 70)
    print("GF(2) -> GF(16) O MATRIX RECONSTRUCTION")
    print("=" * 70)

    print()
    print("Field:")
    print("  GF(16) = GF(2)[a] / (a^4 + a + 1)")

    print()
    print("Matrix dimensions:")
    print(f"  rows = {V}")
    print(f"  cols = {O}")
    print(f"  GF(16) unknowns = {V * O}")
    print(f"  Boolean unknowns = {V * O * 4}")

    # --------------------------------------------------------
    # Read SAT assignment
    # --------------------------------------------------------

    assignment = read_assignment(
        assignment_file
    )

    print()
    print(
        f"Loaded Boolean assignments: "
        f"{len(assignment)}"
    )

    expected_boolean_vars = V * O * 4

    if len(assignment) != expected_boolean_vars:

        print(
            f"WARNING: expected "
            f"{expected_boolean_vars} assignments"
        )

    # --------------------------------------------------------
    # Reconstruct GF(16)
    # --------------------------------------------------------

    recovered = reconstruct_matrix(
        assignment
    )

    # --------------------------------------------------------
    # Read expected O
    # --------------------------------------------------------

    expected = read_expected_o(
        expected_file
    )

    # --------------------------------------------------------
    # Sanity check first few elements
    # --------------------------------------------------------

    show_first_elements(
        assignment,
        recovered,
        expected
    )

    # --------------------------------------------------------
    # Print recovered O
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("RECOVERED O MATRIX")
    print("=" * 70)
    print()

    print_matrix(
        recovered,
        "O_recovered"
    )

    # --------------------------------------------------------
    # Compare
    # --------------------------------------------------------

    mismatches = compare_matrices(
        recovered,
        expected
    )

    sys.exit(
        0 if not mismatches else 2
    )


if __name__ == "__main__":
    main()