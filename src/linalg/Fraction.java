package linalg;

import java.math.BigDecimal;

/** Exact rational number stored in lowest terms with a positive denominator. Overflow throws ArithmeticException. */
public final class Fraction implements Comparable<Fraction> {
    public static final Fraction ZERO = new Fraction(0, 1);
    public static final Fraction ONE = new Fraction(1, 1);

    private final long num, den;

    private Fraction(long num, long den) { this.num = num; this.den = den; }

    public static Fraction of(long n) { return new Fraction(n, 1); }

    public static Fraction of(long n, long d) {
        if (d == 0) throw new ArithmeticException("division by zero");
        if (d < 0) { n = Math.negateExact(n); d = Math.negateExact(d); }
        long g = gcd(Math.absExact(n), d);
        return new Fraction(n / g, d / g);
    }

    /** Accepts "3", "-3/6", or a decimal such as "0.75" (read exactly, not as a binary double). */
    public static Fraction parse(String text) {
        String s = text.trim();
        try {
            int slash = s.indexOf('/');
            if (slash >= 0) return of(Long.parseLong(s.substring(0, slash).trim()), Long.parseLong(s.substring(slash + 1).trim()));
            BigDecimal bd = new BigDecimal(s);
            if (bd.scale() <= 0) return of(bd.longValueExact());
            if (bd.scale() > 18) throw new ArithmeticException("too many decimal places");
            return of(bd.unscaledValue().longValueExact(), Math.round(Math.pow(10, bd.scale())));
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("not a number: '" + text + "'");
        }
    }

    public Fraction add(Fraction o) {
        long g = gcd(den, o.den);
        long d = Math.multiplyExact(den / g, o.den);
        long n = Math.addExact(Math.multiplyExact(num, o.den / g), Math.multiplyExact(o.num, den / g));
        return of(n, d);
    }
    public Fraction subtract(Fraction o) { return add(o.negate()); }
    public Fraction negate() { return new Fraction(Math.negateExact(num), den); }

    public Fraction multiply(Fraction o) {
        long g1 = gcd(Math.absExact(num), o.den), g2 = gcd(Math.absExact(o.num), den);   // cross-reduce first to limit overflow
        return of(Math.multiplyExact(num / g1, o.num / g2), Math.multiplyExact(den / g2, o.den / g1));
    }
    public Fraction divide(Fraction o) {
        if (o.num == 0) throw new ArithmeticException("division by zero");
        return multiply(of(o.den, o.num));
    }

    public boolean isZero() { return num == 0; }
    public boolean isOne() { return num == 1 && den == 1; }
    public int signum() { return Long.signum(num); }
    public double toDouble() { return (double) num / den; }

    @Override public int compareTo(Fraction o) { return Long.compare(Math.multiplyExact(num, o.den), Math.multiplyExact(o.num, den)); }
    @Override public boolean equals(Object o) { return o instanceof Fraction f && f.num == num && f.den == den; }
    @Override public int hashCode() { return Long.hashCode(num) * 31 + Long.hashCode(den); }
    @Override public String toString() { return den == 1 ? Long.toString(num) : num + "/" + den; }

    private static long gcd(long a, long b) { while (b != 0) { long t = a % b; a = b; b = t; } return a == 0 ? 1 : a; }
}
