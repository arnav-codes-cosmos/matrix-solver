package linalg;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

/** Gauss-Jordan elimination over exact fractions, with an optional step-by-step log. */
public final class GaussianElimination {
    public static final Consumer<String> SILENT = s -> {};

    private GaussianElimination() {}

    public record Reduced(Matrix matrix, int[] pivotColumns, int rank, int swaps) {}

    /** Mutable working copy that logs every row operation together with the resulting matrix. */
    private static final class Work {
        final Fraction[][] a; final int split; final Consumer<String> log; int step;

        Work(Matrix m, Consumer<String> log) { a = m.copyData(); split = m.split(); this.log = log; }

        void record(String op) { log.accept("Step " + (++step) + ": " + op + "\n" + new Matrix(a, split) + "\n"); }

        void swap(int i, int j) {
            Fraction[] t = a[i]; a[i] = a[j]; a[j] = t;
            record(String.format("swap R%d <-> R%d", i + 1, j + 1));
        }
        void scale(int i, Fraction k) {
            for (int j = 0; j < a[i].length; j++) a[i][j] = a[i][j].multiply(k);
            record(String.format("R%d <- (%s)*R%d", i + 1, k, i + 1));
        }
        /** R_target <- R_target + k * R_source */
        void addMultiple(int target, int source, Fraction k) {
            for (int j = 0; j < a[target].length; j++) a[target][j] = a[target][j].add(k.multiply(a[source][j]));
            record(k.signum() < 0 ? String.format("R%d <- R%d - (%s)*R%d", target + 1, target + 1, k.negate(), source + 1)
                                  : String.format("R%d <- R%d + (%s)*R%d", target + 1, target + 1, k, source + 1));
        }
        int findPivot(int fromRow, int col) {
            for (int i = fromRow; i < a.length; i++) if (!a[i][col].isZero()) return i;
            return -1;
        }
    }

    /** Reduced row echelon form, eliminating only in the first {@code coefCols} columns (the rest is an augmented part). */
    public static Reduced rref(Matrix m, int coefCols, Consumer<String> log) {
        Work w = new Work(m, log);
        List<Integer> pivots = new ArrayList<>();
        int r = 0, swaps = 0;
        for (int c = 0; c < coefCols && r < w.a.length; c++) {
            int p = w.findPivot(r, c);
            if (p < 0) continue;
            if (p != r) { w.swap(r, p); swaps++; }
            if (!w.a[r][c].isOne()) w.scale(r, Fraction.ONE.divide(w.a[r][c]));
            for (int i = 0; i < w.a.length; i++)
                if (i != r && !w.a[i][c].isZero()) w.addMultiple(i, r, w.a[i][c].negate());
            pivots.add(c); r++;
        }
        return new Reduced(new Matrix(w.a, w.split), pivots.stream().mapToInt(Integer::intValue).toArray(), r, swaps);
    }

    public static int rank(Matrix m, Consumer<String> log) { return rref(m, m.cols(), log).rank(); }

    public static Solution solve(Matrix a, Fraction[] b, Consumer<String> log) {
        if (b.length != a.rows()) throw new IllegalArgumentException("b has " + b.length + " entries but A has " + a.rows() + " rows");
        int n = a.cols();
        Matrix aug = a.augment(Matrix.column(b));
        log.accept("Augmented matrix [A | b]:\n" + aug + "\n");
        Reduced red = rref(aug, n, log);
        Matrix m = red.matrix();
        for (int i = red.rank(); i < m.rows(); i++) {
            if (!m.get(i, n).isZero()) {
                log.accept(String.format("Row %d reads 0 = %s, which is impossible.%n", i + 1, m.get(i, n)));
                return new Solution(Solution.Type.NONE, red.rank(), null, List.of());
            }
        }
        Fraction[] particular = new Fraction[n];
        java.util.Arrays.fill(particular, Fraction.ZERO);
        for (int i = 0; i < red.rank(); i++) particular[red.pivotColumns()[i]] = m.get(i, n);
        if (red.rank() == n) return new Solution(Solution.Type.UNIQUE, n, particular, List.of());

        List<Fraction[]> dirs = new ArrayList<>();
        boolean[] isPivot = new boolean[n];
        for (int p : red.pivotColumns()) isPivot[p] = true;
        for (int f = 0; f < n; f++) {
            if (isPivot[f]) continue;
            Fraction[] d = new Fraction[n];
            java.util.Arrays.fill(d, Fraction.ZERO);
            d[f] = Fraction.ONE;
            for (int i = 0; i < red.rank(); i++) d[red.pivotColumns()[i]] = m.get(i, f).negate();
            dirs.add(d);
        }
        log.accept(String.format("Rank %d < %d unknowns: the other variables are free parameters.%n", red.rank(), n));
        return new Solution(Solution.Type.INFINITE, red.rank(), particular, dirs);
    }

    /** Forward elimination without scaling rows: det = (-1)^swaps * product of the pivots. */
    public static Fraction determinant(Matrix m, Consumer<String> log) {
        if (m.rows() != m.cols()) throw new IllegalArgumentException("determinant needs a square matrix, got " + m.rows() + "x" + m.cols());
        int n = m.rows(), swaps = 0;
        Work w = new Work(m, log);
        log.accept("Forward elimination, with no row scaling, so det = (+/-) product of the pivots:\n" + m + "\n");
        for (int c = 0; c < n; c++) {
            int p = w.findPivot(c, c);
            if (p < 0) { log.accept(String.format("Column %d has no pivot, so det = 0.%n", c + 1)); return Fraction.ZERO; }
            if (p != c) { w.swap(c, p); swaps++; }
            for (int i = c + 1; i < n; i++)
                if (!w.a[i][c].isZero()) w.addMultiple(i, c, w.a[i][c].divide(w.a[c][c]).negate());
        }
        Fraction det = swaps % 2 == 0 ? Fraction.ONE : Fraction.of(-1);
        StringBuilder prod = new StringBuilder();
        for (int i = 0; i < n; i++) { det = det.multiply(w.a[i][i]); prod.append(i > 0 ? " * " : "").append(w.a[i][i]); }
        log.accept(String.format("det = %s%s = %s%n", swaps % 2 == 0 ? "" : "-1 (from " + swaps + " swap" + (swaps == 1 ? "" : "s") + ") * ", prod, det));
        return det;
    }

    /** Inverse by row-reducing [A | I]; empty if A is singular. */
    public static Optional<Matrix> inverse(Matrix m, Consumer<String> log) {
        if (m.rows() != m.cols()) throw new IllegalArgumentException("inverse needs a square matrix, got " + m.rows() + "x" + m.cols());
        int n = m.rows();
        Matrix aug = m.augment(Matrix.identity(n));
        log.accept("Start with [A | I]:\n" + aug + "\n");
        Reduced red = rref(aug, n, log);
        if (red.rank() < n) {
            log.accept(String.format("Rank %d < %d, so the matrix is singular and has no inverse.%n", red.rank(), n));
            return Optional.empty();
        }
        return Optional.of(red.matrix().columns(n, 2 * n));
    }
}
