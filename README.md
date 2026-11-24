# Matrix and Linear-Equation Solver (Java)

A command-line Java project that solves linear systems and computes determinants, inverses, ranks and reduced row echelon forms using Gaussian elimination. Every row operation can be printed step by step, and all arithmetic is **exact**: entries are fractions, so you see `1/2` instead of `0.5000000001`. A dependency-free test suite checks the code against known textbook answers and mathematical identities.

Requires Java 17 or newer. No dependencies.

## Build, test, run

```bash
mkdir out
javac -d out src/linalg/*.java test/linalg/*.java

java -cp out linalg.Tests                       # run the test suite
java -cp out linalg.Main solve "2 1 -1; -3 -1 2; -2 1 2" "8 -11 -3"
java -cp out linalg.Main det "6 1 1; 4 -2 5; 2 8 7"
java -cp out linalg.Main inverse "4 7; 2 6"
java -cp out linalg.Main rank "1 2 3; 4 5 6; 7 8 9"
java -cp out linalg.Main rref "1 2; 3 4"
```

Matrices are written with `;` between rows and spaces between entries. Entries can be integers (`3`), fractions (`-1/2`) or decimals (`0.25`, read exactly as 1/4). Add `--quiet` to hide the steps and `--decimal` to also show decimals next to fractions.

## What it does

| Command | Method |
|---|---|
| `solve A b` | Gauss-Jordan elimination on [A \| b]. Reports a unique solution, **no solution** (a row reads 0 = c), or **infinitely many** (`x = particular + t1*d1 + ...`). Works for square, over- and under-determined systems. |
| `det A` | Forward elimination without scaling rows. The determinant is (−1)^swaps × the product of the pivots, and 0 if a column has no pivot. |
| `inverse A` | Row-reduce [A \| I] to [I \| A⁻¹]. If the rank is below n, reports that A is singular. |
| `rank A` | Number of pivots found by elimination. |
| `rref A` | Reduced row echelon form. |

Each step prints the operation (`swap R1 <-> R2`, `R2 <- R2 + (3)*R1`, `R1 <- (1/2)*R1`) and the matrix after it.

## Example output

```
$ java -cp out linalg.Main det "0 2 1; 1 1 3; 2 0 1"
Step 1: swap R1 <-> R2
[ 1 1 3 ]
[ 0 2 1 ]
[ 2 0 1 ]

Step 2: R3 <- R3 - (2)*R1
[ 1  1  3 ]
[ 0  2  1 ]
[ 0 -2 -5 ]

Step 3: R3 <- R3 + (1)*R2
[ 1 1  3 ]
[ 0 2  1 ]
[ 0 0 -4 ]

det = -1 (from 1 swap) * 1 * 2 * -4 = 8
```

```
$ java -cp out linalg.Main solve "2 1 -1; -3 -1 2; -2 1 2" "8 -11 -3"
Augmented matrix [A | b]:
[  2  1 -1 |   8 ]
[ -3 -1  2 | -11 ]
[ -2  1  2 |  -3 ]

Step 1: R1 <- (1/2)*R1
[  1 1/2 -1/2 |   4 ]
...                          (nine steps in total)
Step 9: R2 <- R2 - (1)*R3
[ 1 0 0 |  2 ]
[ 0 1 0 |  3 ]
[ 0 0 1 | -1 ]

Unique solution:
  x1 = 2
  x2 = 3
  x3 = -1
```

```
$ java -cp out linalg.Main solve "1 2 1; 0 1 -1" "4 1" --quiet
Infinitely many solutions (rank 2, 1 free variable):
  x = (2, 1, 0) + t1*(-3, 1, 1)
```

## Testing

`test/linalg/Tests.java` runs 50 checks with a tiny built-in runner (exit code 1 on any failure), so it needs no JUnit. Expected values come from two sources, and none were produced by the code under test:

**Known answers**
- Fraction arithmetic: 1/2 + 1/3 = 5/6, `0.1 + 0.2 = 3/10` exactly, sign normalisation, parsing
- A classic 3 × 3 system with solution (2, 3, −1), a system that needs a row swap, an inconsistent system, an underdetermined one with the solution line verified by hand, and over-determined consistent and inconsistent systems
- Determinants: a textbook 3 × 3 (−306), the singular 1..9 matrix (0), a 4 × 4 cofactor-expansion example (30), identity, swap matrix, fractional entries
- Inverses: a 2 × 2 and a 3 × 3 textbook pair, the singular case, and the 4 × 4 **Hilbert matrix**, whose inverse is known to be all integers (a good test of exact arithmetic, since floating point loses digits on it)

**Identities on 300 random integer matrices (sizes 2 to 4, fixed seed)**
- A · A⁻¹ = A⁻¹ · A = I whenever the inverse exists
- det(A) ≠ 0 ⇔ inverse exists ⇔ rank = n
- det(AB) = det(A) · det(B) and det(Aᵀ) = det(A)
- For a right-hand side built from a known x, the solver never reports "no solution", and its answer satisfies A x = b

**Error handling:** ragged or non-numeric input, non-square determinant or inverse, wrong-length b, mismatched multiplication.

## Design notes

- **Exact arithmetic.** `Fraction` stores reduced `long` numerator and denominator and uses `Math.*Exact` operations, so overflow throws an `ArithmeticException` instead of silently giving wrong answers. That makes it a teaching tool for small and medium matrices; for large or ill-conditioned numerical work a floating-point version with partial pivoting is the right tool.
- **Pivoting.** With exact arithmetic any non-zero pivot is safe, so the code uses the first non-zero entry in the column and swaps rows only when needed. The swap count feeds the determinant's sign.
- **Cost.** Elimination is O(n³) operations on fractions, whose sizes can grow during elimination.
- **Separation of concerns.** `Fraction` and `Matrix` are immutable values. `GaussianElimination` works on a private mutable copy and reports each operation through a `Consumer<String>` log, so the same code gives silent results, printed steps or captured steps for tests.

## Files

```
src/linalg/
  Fraction.java            exact rational arithmetic
  Matrix.java              immutable matrix, parsing, multiply, transpose, display
  GaussianElimination.java rref, solve, determinant, inverse, rank (with step log)
  Solution.java            unique / none / infinite result with parametric form
  Main.java                command-line interface
test/linalg/Tests.java     test suite
```

## Ideas to extend

- Add a floating-point solver with partial pivoting and compare its output on Hilbert matrices of growing size
- LU decomposition (with row permutation) and solving many right-hand sides cheaply
- Least-squares solutions through the normal equations
- Use `BigInteger` fractions so large matrices cannot overflow
- Read matrices from a file or standard input
