package linalg;

import java.util.List;
import java.util.Optional;
import java.util.Random;

/** Dependency-free tests: every expected value below is a known textbook answer or a mathematical identity. */
public final class Tests {
    private static int passed, failed;

    private static void check(String name, boolean ok) {
        if (ok) passed++; else { failed++; System.out.println("FAIL: " + name); }
    }
    private static Fraction f(String s) { return Fraction.parse(s); }
    private static Matrix m(String s) { return Matrix.parse(s); }
    private static Solution solve(String a, String b) { return GaussianElimination.solve(m(a), m(b).toVector(), GaussianElimination.SILENT); }
    private static Fraction det(String a) { return GaussianElimination.determinant(m(a), GaussianElimination.SILENT); }
    private static Optional<Matrix> inv(String a) { return GaussianElimination.inverse(m(a), GaussianElimination.SILENT); }
    private static boolean throwsIAE(Runnable r) { try { r.run(); return false; } catch (IllegalArgumentException e) { return true; } }

    public static void main(String[] args) {
        fractions(); solving(); determinants(); inverses(); ranks(); properties(); errors();
        System.out.printf("%d passed, %d failed%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    static void fractions() {
        check("1/2 + 1/3 = 5/6", f("1/2").add(f("1/3")).equals(f("5/6")));
        check("2/4 reduces to 1/2", f("2/4").equals(f("1/2")));
        check("negative denominator moves sign", f("3/-6").equals(f("-1/2")) && f("3/-6").toString().equals("-1/2"));
        check("decimal 0.75 = 3/4 exactly", f("0.75").equals(f("3/4")));
        check("(2/3) * (3/4) = 1/2", f("2/3").multiply(f("3/4")).equals(f("1/2")));
        check("(1/2) / (1/4) = 2", f("1/2").divide(f("1/4")).equals(Fraction.of(2)));
        check("x - x = 0", f("7/9").subtract(f("7/9")).isZero());
        check("0.1 + 0.2 = 3/10 (no binary rounding)", f("0.1").add(f("0.2")).equals(f("3/10")));
        check("compare 1/3 < 1/2", f("1/3").compareTo(f("1/2")) < 0);
        boolean threw = false;
        try { f("1").divide(Fraction.ZERO); } catch (ArithmeticException e) { threw = true; }
        check("division by zero throws", threw);
    }

    static void solving() {
        Solution s = solve("2 1 -1; -3 -1 2; -2 1 2", "8 -11 -3");   // classic: x=2, y=3, z=-1
        check("3x3 unique type", s.type() == Solution.Type.UNIQUE);
        check("3x3 unique values", s.particular()[0].equals(f("2")) && s.particular()[1].equals(f("3")) && s.particular()[2].equals(f("-1")));

        s = solve("1 1; 1 -1", "3 1");
        check("2x2 unique (2, 1)", s.particular()[0].equals(f("2")) && s.particular()[1].equals(f("1")));

        s = solve("0 1; 1 0", "2 3");   // zero pivot forces a row swap
        check("row swap needed: x=3, y=2", s.particular()[0].equals(f("3")) && s.particular()[1].equals(f("2")));

        s = solve("2 3; 4 9", "1 1");   // fractional answer (x = 1, y = -1/3), verified by substitution
        check("fractional solution satisfies A x = b",
              m("2 3; 4 9").multiply(Matrix.column(s.particular())).equals(Matrix.column(new Fraction[]{f("1"), f("1")})));

        s = solve("1 1; 1 1", "1 2");
        check("inconsistent system has no solution", s.type() == Solution.Type.NONE);

        s = solve("1 2 1; 0 1 -1", "4 1");   // x+2y+z=4, y-z=1
        check("underdetermined is INFINITE with rank 2", s.type() == Solution.Type.INFINITE && s.rank() == 2);
        check("particular solution (2, 1, 0)", s.particular()[0].equals(f("2")) && s.particular()[1].equals(f("1")) && s.particular()[2].equals(f("0")));
        check("direction (-3, 1, 1)", s.directions().size() == 1 && s.directions().get(0)[0].equals(f("-3")) && s.directions().get(0)[1].equals(f("1")) && s.directions().get(0)[2].equals(f("1")));

        s = solve("1 1; 1 -1; 2 1", "3 1 5");   // overdetermined but consistent
        check("overdetermined consistent -> (2, 1)", s.type() == Solution.Type.UNIQUE && s.particular()[0].equals(f("2")) && s.particular()[1].equals(f("1")));
        s = solve("1 1; 1 -1; 2 1", "3 1 6");
        check("overdetermined inconsistent", s.type() == Solution.Type.NONE);

        s = solve("0 0; 0 0", "0 0");
        check("zero system: everything solves it (2 free variables)", s.type() == Solution.Type.INFINITE && s.directions().size() == 2);
    }

    static void determinants() {
        check("diagonal 2,3 -> 6", det("2 0; 0 3").equals(f("6")));
        check("[[1,2],[3,4]] -> -2", det("1 2; 3 4").equals(f("-2")));
        check("3x3 textbook -> -306", det("6 1 1; 4 -2 5; 2 8 7").equals(f("-306")));
        check("singular 1..9 -> 0", det("1 2 3; 4 5 6; 7 8 9").isZero());
        check("identity -> 1", det("1 0 0; 0 1 0; 0 0 1").equals(f("1")));
        check("swap matrix -> -1", det("0 1; 1 0").equals(f("-1")));
        check("4x4 cofactor example -> 30", det("1 0 2 -1; 3 0 0 5; 2 1 4 -3; 1 0 5 0").equals(f("30")));
        check("1x1 -> the entry", det("7/3").equals(f("7/3")));
        check("fractional entries: [[1/2,1/3],[1/4,1/5]] -> 1/10 - 1/12 = 1/60", det("1/2 1/3; 1/4 1/5").equals(f("1/60")));
    }

    static void inverses() {
        check("[[4,7],[2,6]] inverse", inv("4 7; 2 6").get().equals(m("3/5 -7/10; -1/5 2/5")));
        check("3x3 classic inverse", inv("1 2 3; 0 1 4; 5 6 0").get().equals(m("-24 18 5; 20 -15 -4; -5 4 1")));
        check("singular has no inverse", inv("1 2; 2 4").isEmpty());
        check("1x1 inverse", inv("4").get().equals(m("1/4")));
        check("4x4 Hilbert inverse is all integers", inv("1 1/2 1/3 1/4; 1/2 1/3 1/4 1/5; 1/3 1/4 1/5 1/6; 1/4 1/5 1/6 1/7").get()
              .equals(m("16 -120 240 -140; -120 1200 -2700 1680; 240 -2700 6480 -4200; -140 1680 -4200 2800")));
    }

    static void ranks() {
        check("rank [[1,2],[2,4]] = 1", GaussianElimination.rank(m("1 2; 2 4"), GaussianElimination.SILENT) == 1);
        check("rank of 1..9 = 2", GaussianElimination.rank(m("1 2 3; 4 5 6; 7 8 9"), GaussianElimination.SILENT) == 2);
        check("rank of identity = 3", GaussianElimination.rank(Matrix.identity(3), GaussianElimination.SILENT) == 3);
        check("rank of 2x3 full = 2", GaussianElimination.rank(m("1 0 1; 0 1 1"), GaussianElimination.SILENT) == 2);
        check("RREF of [[1,2],[3,4]] is identity", GaussianElimination.rref(m("1 2; 3 4"), 2, GaussianElimination.SILENT).matrix().equals(Matrix.identity(2)));
    }

    /** Identities that must hold for any matrix, checked on 300 random integer matrices. */
    static void properties() {
        Random rnd = new Random(12345);
        int invertible = 0, bad = 0;
        for (int t = 0; t < 300; t++) {
            int n = 2 + rnd.nextInt(3);
            Matrix a = randomMatrix(rnd, n), b = randomMatrix(rnd, n);
            Fraction da = GaussianElimination.determinant(a, GaussianElimination.SILENT);
            Optional<Matrix> ia = GaussianElimination.inverse(a, GaussianElimination.SILENT);
            int rank = GaussianElimination.rank(a, GaussianElimination.SILENT);
            if (ia.isPresent() != !da.isZero() || (rank == n) != !da.isZero()) bad++;
            if (ia.isPresent()) { invertible++; if (!a.multiply(ia.get()).equals(Matrix.identity(n)) || !ia.get().multiply(a).equals(Matrix.identity(n))) bad++; }
            if (!GaussianElimination.determinant(a.multiply(b), GaussianElimination.SILENT).equals(da.multiply(GaussianElimination.determinant(b, GaussianElimination.SILENT)))) bad++;
            if (!GaussianElimination.determinant(a.transpose(), GaussianElimination.SILENT).equals(da)) bad++;
            Fraction[] x = new Fraction[n], rhs;
            for (int i = 0; i < n; i++) x[i] = Fraction.of(rnd.nextInt(9) - 4);
            rhs = a.multiply(Matrix.column(x)).toVector();
            Solution s = GaussianElimination.solve(a, rhs, GaussianElimination.SILENT);   // b built from a known x, so it must be solvable
            if (s.type() == Solution.Type.NONE) bad++;
            else if (!a.multiply(Matrix.column(s.particular())).equals(Matrix.column(rhs))) bad++;
        }
        check("300 random matrices satisfy all identities", bad == 0);
        check("random sample included invertible matrices", invertible > 100);
    }

    static Matrix randomMatrix(Random rnd, int n) {
        Fraction[][] d = new Fraction[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d[i][j] = Fraction.of(rnd.nextInt(11) - 5);
        return new Matrix(d, 0);
    }

    static void errors() {
        check("ragged matrix rejected", throwsIAE(() -> m("1 2; 3")));
        check("non-number rejected", throwsIAE(() -> m("1 x")));
        check("determinant of non-square rejected", throwsIAE(() -> det("1 2 3; 4 5 6")));
        check("inverse of non-square rejected", throwsIAE(() -> inv("1 2 3; 4 5 6")));
        check("b of wrong length rejected", throwsIAE(() -> solve("1 0; 0 1", "1 2 3")));
        check("multiply shape mismatch rejected", throwsIAE(() -> m("1 2").multiply(m("1 2"))));
        check("step log is produced", capture(() -> GaussianElimination.solve(m("2 1; 1 3"), m("3 5").toVector(), l -> sink.add(l))) > 2);
    }

    private static final List<String> sink = new java.util.ArrayList<>();
    private static int capture(Runnable r) { sink.clear(); r.run(); return sink.size(); }
}
