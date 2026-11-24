package linalg;

import java.util.Arrays;

/** Immutable matrix of exact fractions. {@code split} marks where an augmented part starts (0 = none) and is only used for display. */
public final class Matrix {
    private final Fraction[][] a;
    private final int split;

    Matrix(Fraction[][] data, int split) {
        a = new Fraction[data.length][];
        for (int i = 0; i < data.length; i++) a[i] = data[i].clone();
        this.split = split;
    }

    /** Rows separated by ';' or newlines, entries by spaces or commas: "1 2; 3 4". Entries may be 3, -1/2 or 0.25. */
    public static Matrix parse(String text) {
        String[] lines = text.trim().split("\\s*[;\\n]\\s*");
        Fraction[][] d = new Fraction[lines.length][];
        for (int i = 0; i < lines.length; i++) {
            String[] parts = lines[i].trim().split("[\\s,]+");
            d[i] = Arrays.stream(parts).map(Fraction::parse).toArray(Fraction[]::new);
            if (d[i].length != d[0].length) throw new IllegalArgumentException("row " + (i + 1) + " has " + d[i].length + " entries, expected " + d[0].length);
        }
        if (d[0].length == 0) throw new IllegalArgumentException("empty matrix");
        return new Matrix(d, 0);
    }

    public static Matrix identity(int n) {
        Fraction[][] d = new Fraction[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d[i][j] = i == j ? Fraction.ONE : Fraction.ZERO;
        return new Matrix(d, 0);
    }

    public static Matrix column(Fraction[] v) {
        Fraction[][] d = new Fraction[v.length][1];
        for (int i = 0; i < v.length; i++) d[i][0] = v[i];
        return new Matrix(d, 0);
    }

    public int rows() { return a.length; }
    public int cols() { return a[0].length; }
    public int split() { return split; }
    public Fraction get(int i, int j) { return a[i][j]; }
    Fraction[][] copyData() { Fraction[][] d = new Fraction[a.length][]; for (int i = 0; i < a.length; i++) d[i] = a[i].clone(); return d; }

    /** All entries in row-major order, so "8 -11 -3" and "8; -11; -3" both give the same vector. */
    public Fraction[] toVector() {
        return Arrays.stream(a).flatMap(Arrays::stream).toArray(Fraction[]::new);
    }

    public Matrix augment(Matrix right) {
        if (right.rows() != rows()) throw new IllegalArgumentException("row counts differ");
        Fraction[][] d = new Fraction[rows()][cols() + right.cols()];
        for (int i = 0; i < rows(); i++) {
            System.arraycopy(a[i], 0, d[i], 0, cols());
            System.arraycopy(right.a[i], 0, d[i], cols(), right.cols());
        }
        return new Matrix(d, cols());
    }

    public Matrix columns(int from, int to) {
        Fraction[][] d = new Fraction[rows()][to - from];
        for (int i = 0; i < rows(); i++) System.arraycopy(a[i], from, d[i], 0, to - from);
        return new Matrix(d, 0);
    }

    public Matrix transpose() {
        Fraction[][] d = new Fraction[cols()][rows()];
        for (int i = 0; i < rows(); i++) for (int j = 0; j < cols(); j++) d[j][i] = a[i][j];
        return new Matrix(d, 0);
    }

    public Matrix multiply(Matrix o) {
        if (cols() != o.rows()) throw new IllegalArgumentException("cannot multiply " + rows() + "x" + cols() + " by " + o.rows() + "x" + o.cols());
        Fraction[][] d = new Fraction[rows()][o.cols()];
        for (int i = 0; i < rows(); i++) for (int j = 0; j < o.cols(); j++) {
            Fraction s = Fraction.ZERO;
            for (int k = 0; k < cols(); k++) s = s.add(a[i][k].multiply(o.a[k][j]));
            d[i][j] = s;
        }
        return new Matrix(d, 0);
    }

    @Override public boolean equals(Object o) { return o instanceof Matrix m && Arrays.deepEquals(a, m.a); }
    @Override public int hashCode() { return Arrays.deepHashCode(a); }

    @Override public String toString() {
        int[] w = new int[cols()];
        for (Fraction[] row : a) for (int j = 0; j < row.length; j++) w[j] = Math.max(w[j], row[j].toString().length());
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rows(); i++) {
            sb.append("[ ");
            for (int j = 0; j < cols(); j++) {
                if (split > 0 && j == split) sb.append("| ");
                sb.append(" ".repeat(w[j] - a[i][j].toString().length())).append(a[i][j]).append(' ');
            }
            sb.append("]");
            if (i < rows() - 1) sb.append('\n');
        }
        return sb.toString();
    }
}
