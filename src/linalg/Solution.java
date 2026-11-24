package linalg;

import java.util.List;
import java.util.Locale;
import java.util.StringJoiner;

/** Outcome of solving A x = b. For INFINITE, every solution is {@code particular + t1*d1 + t2*d2 + ...}. */
public record Solution(Type type, int rank, Fraction[] particular, List<Fraction[]> directions) {
    public enum Type { UNIQUE, NONE, INFINITE }

    public String describe(boolean decimals) {
        switch (type) {
            case NONE: return "No solution: the system is inconsistent.";
            case UNIQUE: {
                StringBuilder sb = new StringBuilder("Unique solution:");
                for (int i = 0; i < particular.length; i++) sb.append("\n  x").append(i + 1).append(" = ").append(show(particular[i], decimals));
                return sb.toString();
            }
            default: {
                StringBuilder sb = new StringBuilder(String.format("Infinitely many solutions (rank %d, %d free variable%s):\n  x = %s",
                        rank, directions.size(), directions.size() == 1 ? "" : "s", vec(particular)));
                for (int k = 0; k < directions.size(); k++) sb.append(" + t").append(k + 1).append("*").append(vec(directions.get(k)));
                return sb.toString();
            }
        }
    }

    private static String vec(Fraction[] v) {
        StringJoiner j = new StringJoiner(", ", "(", ")");
        for (Fraction f : v) j.add(f.toString());
        return j.toString();
    }

    private static String show(Fraction f, boolean decimals) {
        return decimals && f.toString().contains("/") ? String.format(Locale.ROOT, "%s (about %.6f)", f, f.toDouble()) : f.toString();
    }
}
