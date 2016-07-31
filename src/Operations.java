import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.rational.ExtendedRational;

import java.util.Iterator;
import java.util.PriorityQueue;


public class Operations {
    public static SetInterval scalarMul(SetInterval[] X, SetInterval[] Y, SetIntervalContext ic) {
        SetInterval result = ic.numsToInterval(0, 0);
        for (int i = 0; i < X.length; i++) {
            result = ic.add(result, ic.mul(X[i], Y[i]));
        }
        return result;
    }

    public static Gradient scalarMul(Gradient[] X, Gradient[] Y, SetIntervalContext ic) {
        Gradient result = Gradient.num(ExtendedRational.valueOf(0));
        for (int i = 0; i < X.length; i++) {
            result = result.add(X[i].mul(Y[i]));
        }
        return result;
    }

    public static GradientSD scalarMul(GradientSD[] X, GradientSD[] Y, SetIntervalContext ic) {
        GradientSD result = GradientSD.num(ExtendedRational.valueOf(0));
        for (int i = 0; i < X.length; i++) {
            result = result.add(X[i].mul(Y[i]));
        }
        return result;
    }


    public static void rmUseless(ExtendedRational supMin, PriorityQueue<ListElem> wList) {
        Iterator<ListElem> iterator = wList.iterator();
        while (iterator.hasNext()) {
            if (iterator.next().getAssessment().ge(supMin))
                iterator.remove();
        }
    }

    public static String[] keysReader(String args) {
        String[] result = new String[args.length() - 5];
        String tmp;
        int keyConuter = 0;
        int i = 0;
        int j = 0;
        while (i < result.length) {
            if (args.charAt(i) == '-') {
                keyConuter++;
                j = i + 1;
                tmp = "-";
                while ((args.charAt(j) != ' ') || (args.charAt(j) != '\n') || (args.charAt(j) != '-')) { //regular expressions needed
                    tmp = tmp + args.charAt(j);
                    j++;
                }
                i = j;
                result[keyConuter] = tmp;
            }
            if (i != j)
                i++;
        }
        return result;
    }
}
