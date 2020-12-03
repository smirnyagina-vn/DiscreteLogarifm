import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Random;

public class DiscreteLogarithm {

    public static final double LOG_2 = Math.log(2.0);
    private static final int MAX_DIGITS_2 = 977;

    private static BigInteger minusOne = BigInteger.valueOf(-1);
    private static BigInteger zero = BigInteger.valueOf(0);
    private static BigInteger one = BigInteger.valueOf(1);
    private static BigInteger two = BigInteger.valueOf(2);

    private BigInteger aValue = new BigInteger("13");//13//10//13
    private BigInteger bValue = new BigInteger("13957");//3//64//13957
    private BigInteger pValue = new BigInteger("50091");//50091786122438801387//107//50091
    private BigInteger rValue = new BigInteger("8178");//25045893061219400693//53//8178

    long startTime;

    public DiscreteLogarithm()
    {
        startTime = System.currentTimeMillis();

        //BigInteger result = two.multiply(rValue.sqrt().add(BigInteger.valueOf((long)logBigInteger(rValue))) );
        //System.out.println(2.82*rValue.sqrt().doubleValue());

        P0PollardFunction(pValue,aValue,bValue,rValue);
        //Gelfond(pValue,aValue,bValue,rValue);

        System.out.println("Time: " + (double) (System.currentTimeMillis() - startTime)  + "\n");
    }

    public static double logBigInteger(BigInteger val) {
        if (val.signum() < 1)
            return val.signum() < 0 ? Double.NaN : Double.NEGATIVE_INFINITY;
        int blex = val.bitLength() - MAX_DIGITS_2; // any value in 60..1023 works here
        if (blex > 0)
            val = val.shiftRight(blex);
        double res = Math.log(val.doubleValue());
        return blex > 0 ? res + blex * LOG_2 : res;
    }

    public static void main(String[] args) {
        DiscreteLogarithm discreteLogarithm = new DiscreteLogarithm();
    }

    public BigInteger getRandomBigInteger(BigInteger maxLimit, BigInteger minLimit) {
        BigInteger bigInteger = maxLimit.subtract(minLimit);
        Random randNum = new Random();
        int len = maxLimit.bitLength();
        BigInteger resultBigInteger = new BigInteger(len, randNum);
        if (resultBigInteger.compareTo(minLimit) < 0)
            resultBigInteger = resultBigInteger.add(minLimit);
        if (resultBigInteger.compareTo(bigInteger) >= 0)
            resultBigInteger = resultBigInteger.mod(bigInteger).add(minLimit);
        return resultBigInteger;
    }

    public BigInteger P0SubFunction(BigInteger n)
    {
        BigInteger z = n.sqrt().add(one);
        BigInteger result = n;

        for (BigInteger counter = two; counter.compareTo(z) == -1; counter = counter.add(one))
        {
            if (n.mod(counter) == zero)
            {
                while (n.mod(counter) == zero)
                {
                    n = n.divide(counter);
                    result.subtract(result.divide(counter));
                }
            }
        }

        if (n.compareTo(one) == 1)
            result = result.subtract(result.divide(n));
        return result;
    }

    public BigInteger CompareSolve(BigInteger k1, BigInteger b1, BigInteger k2, BigInteger b2, BigInteger r)
    {
        BigInteger k = k1.subtract(k2);
        BigInteger b = b2.subtract(b1);
        if (k.compareTo(zero) == -1)
        {
            k = k.add(r);
        }
        if (b.compareTo(zero) == -1)
        {
            b = b.add(r);
        }
        BigInteger x = zero;

        while (true)
        {
            BigInteger tmp = k.multiply(x).mod(r);
            if (tmp.equals(b))
                break;
            x = x.add(one);
        }

        return x;
    }

    public void P0PollardFunction(BigInteger p, BigInteger a, BigInteger b, BigInteger r)
    {
        p = p.abs();
        a = a.abs();
        b = b.abs();
        r = r.abs();

        String column1Format = "%-5s";
        String column2Format = "%-10s";
        String formatInfo = column1Format + " | " + column2Format + " | " + column2Format + " | " + column2Format + " | " + column2Format;

        System.out.println("                        Pollard");
        System.out.println("-----------------------------------------------------");
        System.out.format(formatInfo, "i" , "c" , "log(a)c", "d ", "log(a)d");
        System.out.println();
        System.out.println("-----------------------------------------------------");

        //BigInteger u = getRandomBigInteger(p.subtract(one), one);
        //BigInteger v = getRandomBigInteger(p.subtract(one), one);
        BigInteger u = two;
        BigInteger v = two;
        BigInteger cu = u;
        BigInteger cv = v;
        BigInteger du = u;
        BigInteger dv = v;

        BigInteger c = a.pow(u.intValue()).multiply(b.pow(v.intValue())).mod(p);
        BigInteger d = c;

        BigInteger iterations = zero;

        while (true) {

            //if (iterations.mod(BigInteger.valueOf(1000000)) == zero) {
            //    System.out.println("Iteration: " + iterations + " Time: " + (double) (System.currentTimeMillis() - startTime)  + "\n");
            //}

            //System.out.println("  " + iterations + "  |  " + c + "  | "
            //        + cu + "+" + cv + "x |  " + d + "  | " + du + "+" + dv + "x");

            System.out.format(formatInfo,iterations, c, cu + "+" + cv + "x",d,du + "+" + dv + "x");
            System.out.println();

            if (c.compareTo(p.divide(two)) == -1)
            {
                c = a.multiply(c).mod(p);
                cu = cu.add(one);
            }
            else
            {
                c = b.multiply(c).mod(p);
                cv = cv.add(one);
            }

            for (int counter = 0; counter < 2; counter++)
            {
                if (d.compareTo(p.divide(two)) == -1)
                {
                    d = a.multiply(d).mod(p);
                    du = du.add(one);
                }
                else
                {
                    d = b.multiply(d).mod(p);
                    dv = dv.add(one);
                }
            }

            iterations = iterations.add(one);

            if (c.equals(d)) break;
        }

        System.out.format(formatInfo,iterations, c, cu + "+" + cv + "x",d,du + "+" + dv + "x");
        System.out.println();

        System.out.println(cu + "+" + cv + "x"+ " = " +du + "+" + dv + "x");

        a = cv.subtract(dv).mod(r);
        b = du.subtract(cu).mod(r);
        System.out.println("x = " + CompareSolve(cv, cu, dv, du, r) + " mod(" + r + ")");
    }

    public void Gelfond(BigInteger p, BigInteger a, BigInteger b, BigInteger r)
    {
        p = p.abs();
        a = a.abs();
        b = b.abs();
        r = r.abs();

        System.out.println("                Gelfond");
        System.out.println("----------------------------------------");
        ArrayList<BigInteger> sequenceA = new ArrayList<BigInteger>();//1, a, a^2, ..., a^(s-1) (mod p)
        ArrayList<BigInteger> sequenceBA = new ArrayList<BigInteger>();//b, ba^(-1*s), ba^(-2*s), ..., ba^-(s-1)s (mod p)

        BigInteger u = zero;
        BigInteger v = zero;
        BigInteger i = zero;
        BigInteger s = r.sqrt().add(one);
        System.out.println("s = " + s);
        BigInteger tmp = zero;

        System.out.println("Counting ba sequence...");
        while (i.compareTo(s) == -1)
        {
            tmp = b.multiply(a.modPow(i.multiply(r.subtract(s)), p)).mod(p);
            sequenceBA.add(tmp);
            if (i.mod(BigInteger.valueOf(1000000)) == zero) {
                //System.out.println(" " + i + " | " + tmp);
                System.out.println("Iteration: " + i + " Time: " + (double) (System.currentTimeMillis() - startTime)  + "\n");
            }
            i = i.add(one);
        }

        ArrayList<BigInteger> sortedCopySequenceBA = (ArrayList<BigInteger>) sequenceBA.clone();
        sortedCopySequenceBA.sort(BigInteger::compareTo);
        System.out.println("Sorted sequence b, ba^(-1*s), ba^(-2*s), ..., ba^-(s-1)s (mod p) :");
        for (BigInteger counter:sortedCopySequenceBA)
        {
            System.out.println("(" + counter + "," + sequenceBA.indexOf(counter) + ") ");
        }

        i = zero;

        while(i.compareTo(s) == -1)
        {
            sequenceA.add(a.modPow(i,p));
            for (BigInteger k = zero; k.compareTo(i) == -1; k = k.add(one))
            {
                if (sequenceBA.contains(sequenceA.get(k.intValue())))
                {
                    u = BigInteger.valueOf(sequenceBA.indexOf(sequenceA.get(k.intValue())));
                    v = k;
                    break;
                }
            }

            i = i.add(one);
        }

        System.out.println("x = " + u.multiply(s).add(v).mod(r));

    }

}
