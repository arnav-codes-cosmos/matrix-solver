package linalg;

import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Optional;
import java.util.function.Consumer;

public final class Main {
    public static void main(String[] args) {
        List<String> rest = Arrays.stream(args).filter(a -> !a.startsWith("--")).toList();
        boolean quiet = Arrays.asList(args).contains("--quiet"), decimals = Arrays.asList(args).contains("--decimal");
        if (rest.isEmpty() || rest.get(0).equals("help")) { usage(); return; }
        Consumer<String> log = quiet ? GaussianElimination.SILENT : System.out::println;
        try {
            String cmd = rest.get(0);
            Matrix a = Matrix.parse(need(rest, 1));
            switch (cmd) {
                case "solve" -> {
                    Fraction[] b = Matrix.parse(need(rest, 2)).toVector();
                    System.out.println(GaussianElimination.solve(a, b, log).describe(decimals));
                }
                case "det" -> System.out.println("det = " + fmt(GaussianElimination.determinant(a, log), decimals));
                case "inverse" -> {
                    Optional<Matrix> inv = GaussianElimination.inverse(a, log);
                    System.out.println(inv.isPresent() ? "Inverse:\n" + inv.get() : "The matrix is singular: no inverse exists.");
                }
                case "rank" -> System.out.println("rank = " + GaussianElimination.rank(a, log));
                case "rref" -> System.out.println("Reduced row echelon form:\n" + GaussianElimination.rref(a, a.cols(), log).matrix());
                default -> { System.err.println("Unknown command: " + cmd); usage(); System.exit(1); }
            }
        } catch (IllegalArgumentException | ArithmeticException e) {
            System.err.println("Error: " + e.getMessage());
            System.exit(1);
        }
    }

    private static String need(List<String> a, int i) {
        if (i >= a.size()) throw new IllegalArgumentException("missing argument; run with 'help' for usage");
        return a.get(i);
    }

    private static String fmt(Fraction f, boolean decimals) {
        return decimals && f.toString().contains("/") ? String.format(Locale.ROOT, "%s (about %.6f)", f, f.toDouble()) : f.toString();
    }

    private static void usage() {
        System.out.println("""
                Matrix and linear-equation solver (exact fractions, step by step)
                Usage: java -cp out linalg.Main <command> "<matrix>" ["<b vector>"] [--quiet] [--decimal]

                Commands:
                  solve   "<A>" "<b>"   solve A x = b (unique, none, or infinitely many solutions)
                  det     "<A>"         determinant
                  inverse "<A>"         inverse, or report that A is singular
                  rank    "<A>"         rank
                  rref    "<A>"         reduced row echelon form

                Matrices: rows separated by ';', entries by spaces. Entries can be 3, -1/2 or 0.25.
                Example:  solve "2 1 -1; -3 -1 2; -2 1 2" "8 -11 -3"
                Options:  --quiet hides the steps, --decimal also shows decimals for fractions.""");
    }
}
