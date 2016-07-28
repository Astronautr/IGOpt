import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.*;


public class Gradient {
/**
 *   Gradient class for interval computing on Java.
 *   With this class you can compute partial derivatives and the range of values for some functions.
 *   Here was used the idea of automatic differentiation.
 *   JInterval library was used for all interval arithmetic and functions.
*/
    private static SetIntervalContext ic; //This field determine accuracy of computing.
    private static int dim; //Represents the number of variables.
    private SetInterval X; //Interval value of expression.
    private SetInterval[] dX; //Interval values of derivative.
    private static SetInterval ZERO = SetIntervalContexts.getExact().numsToInterval(0,0);
    private static SetInterval ONE = SetIntervalContexts.getExact().numsToInterval(1,1);

    public Gradient() {
        X = null;
        dX = null;
    }

    public Gradient(SetInterval box, SetIntervalContext ic) {
        dim = 1;
        dX = new SetInterval[1];
        dX[0] = ONE;
        X = box;
        Gradient.ic = ic;
    }

    public static Gradient[] init(SetInterval[] box, SetIntervalContext ic) {
        dim = box.length;
        Gradient.ic = ic;
        Gradient[] result = new Gradient[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = new Gradient();
            result[i].dX = new SetInterval[dim];
            result[i].X = box[i];
            for (int j = 0; j < i; j++) {
                result[i].dX[j] = ZERO;
            }
            result[i].dX[i] = ONE;
            for (int j = i + 1; j < dim; j++) {
                result[i].dX[j] = ZERO;
            }
        }
        return result;
    }

    public static Gradient num(double number) {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(number,number);
        return result;
    }

    public static Gradient num(Rational number) {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(number,number);
        return result;
    }

    public static Gradient num(ExtendedRational number) {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(number,number);
        return result;
    }

    public static Gradient nums(double lower,double upper) {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower,upper);
        return result;
    }

    public static Gradient nums(Rational lower,Rational upper) {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower,upper);
        return result;
    }

    public static Gradient nums(ExtendedRational lower,ExtendedRational upper) {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower,upper);
        return result;
    }

    public Gradient neg() {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(dX[i]);
        }
        result.X = ic.neg(X);
        return result;
    }

    public Gradient intersection(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i],Y.dX[i]);
        }
        return result;
    }

    public Gradient intersection(GradientMV Y) {
        Gradient result = new Gradient();
        result.X = ic.intersection(this.X, Y.getX());
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i],Y.getDX()[i]);
        }
        return result;
    }

    public Gradient intersectionX(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.dX = this.dX.clone();
        return result;
    }

    public Gradient intersectionX(SetInterval Y) {
        Gradient result = new Gradient();
        result.X = ic.intersection(this.X, Y);
        result.dX = new SetInterval[dim];
        result.dX = this.dX.clone();
        return result;
    }

    public Gradient intersectionDX(Gradient Y) {
        Gradient result = new Gradient();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public Gradient intersectionDX(SetInterval[] Y) {
        Gradient result = new Gradient();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y[i]);
        }
        return result;
    }

    public Gradient intersectionDXInd(Gradient Y, int ind) {
        Gradient result = new Gradient();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        result.dX = this.dX.clone();
        result.dX[ind] = ic.intersection(this.dX[ind],Y.dX[ind]);
        return result;
    }

    public Gradient intersectionDXInd(SetInterval Y, int ind) {
        Gradient result = new Gradient();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        result.dX = this.dX.clone();
        result.dX[ind] = ic.intersection(this.dX[ind],Y);
        return result;
    }

    public Gradient intersectionDXInd(SetInterval[] Y, int ind) {
        Gradient result = new Gradient();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        result.dX = this.dX.clone();
        result.dX[ind] = ic.intersection(this.dX[ind],Y[ind]);
        return result;
    }

    public Gradient hull(Gradient Y) {
        Gradient result = new Gradient();
        ExtendedRational min, max;
        if (this.X.inf().le(Y.X.inf()))
            min = this.X.inf();
        else
            min = Y.X.inf();
        if (this.X.sup().ge(Y.X.sup()))
            max = this.X.sup();
        else
            max = Y.X.sup();
        result.X = ic.numsToInterval(min, max);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            if (this.dX[i].inf().le(Y.dX[i].inf()))
                min = this.dX[i].inf();
            else
                min = Y.dX[i].inf();
            if (this.dX[i].sup().ge(Y.dX[i].sup()))
                max = this.dX[i].sup();
            else
                max = Y.dX[i].sup();
            result.dX[i] = ic.numsToInterval(min, max);
        }
        return result;
    }

    public Gradient max (Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.max(this.X,Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.max(this.dX[i],Y.dX[i]);
        }
        return result;
    }

    public Gradient abs() {
        Gradient result = new Gradient();
        result.X = ic.abs(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
//            result.dX[i] = ic.convexHull(ic.neg(this.dX[i]),this.dX[i]);
            result.dX[i] = ic.abs(this.dX[i]);
        }
        return result;
    }

    public Gradient add(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.add(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(this.dX[i],Y.dX[i]);
        }
        return result;
    }

    public Gradient sub(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.sub(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.sub(this.dX[i],Y.dX[i]);
        }
        return result;
    }

    public Gradient mul(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.mul(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(ic.mul(this.X, Y.dX[i]),ic.mul(this.dX[i],Y.X));
        }
        return result;
    }

    public Gradient div(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.div(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(ic.sub(ic.mul(this.dX[i], Y.X),ic.mul(this.X,Y.dX[i])),ic.sqr(Y.X));
        }
        return result;
    }

    public Gradient recip() {
        Gradient result = new Gradient();
        result.X = ic.recip(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.div(this.dX[i],ic.sqr(this.X)));
        }
        return result;
    }

    public Gradient pow(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.pow(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y.X, ic.pow(X,ic.sub(Y.X, ONE))),this.dX[i]);
/*            result.dX[i] = ic.add(ic.mul(ic.mul(ic.pow(this.X,Y.X),Y.dX[i]),ic.log(this.X)),
                    ic.mul(ic.mul(Y.X,this.dX[i]),ic.pow(this.X,ic.sub(Y.X,ONE))));*/
        }
        return result;
    }

    public Gradient pown(int n) {
        Gradient result = new Gradient();
        result.X = ic.pown(this.X, n);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n,n), ic.pown(X,n-1)),this.dX[i]);
        }
        return result;
    }

    public Gradient pown(long n) {
        Gradient result = new Gradient();
        result.X = ic.pown(this.X, n);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n,n), ic.pown(X,n-1)),this.dX[i]);
        }
        return result;
    }

    public Gradient rootn(int n) {
        Gradient result = new Gradient();
        SetInterval Y = ic.recip(ic.numsToInterval(n,n));
        result.X = ic.pow(this.X, Y);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y, ic.pow(X,ic.sub(Y, ONE))),this.dX[i]);
        }
        return result;
    }

    public Gradient rootn(long n) {
        Gradient result = new Gradient();
        SetInterval Y = ic.recip(ic.numsToInterval(n,n));
        result.X = ic.pow(this.X, Y);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y, ic.pow(this.X,ic.sub(Y, ONE))),this.dX[i]);
        }
        return result;
    }

    public Gradient sqr() {
        Gradient result = new Gradient();
        result.X = ic.sqr(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(2,2),X),this.dX[i]);
        }
        return result;
    }

    public Gradient sqrt() {
        Gradient result = new Gradient();
        result.X = ic.sqrt(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(dX[i], ic.mul(ic.numsToInterval(2,2),ic.sqrt(X)));
        }
        return result;
    }

    public Gradient sin() {
        Gradient result = new Gradient();
        result.X = ic.sin(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.cos(X),dX[i]);
        }
        return result;
    }

    public Gradient cos() {
        Gradient result = new Gradient();
        result.X = ic.cos(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.mul(ic.sin(X),dX[i]));
        }
        return result;
    }

    public Gradient tan() {
        Gradient result = new Gradient();
        result.X = ic.tan(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ONE,ic.sqr(ic.cos(X))),dX[i]);
        }
        return result;
    }

    public Gradient asin() {
        Gradient result = new Gradient();
        result.X = ic.asin(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ONE,ic.sqrt(ic.sub(ONE,ic.sqr(X)))),dX[i]);
        }
        return result;
    }

    public Gradient acos() {
        Gradient result = new Gradient();
        result.X = ic.acos(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.mul(ic.div(ONE,ic.sqrt(ic.sub(ONE,ic.sqr(X)))),dX[i]));
        }
        return result;
    }

    public Gradient atan() {
        Gradient result = new Gradient();
        result.X = ic.atan(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ONE,ic.add(ONE,ic.sqr(X))),dX[i]);
        }
        return result;
    }

    public Gradient exp() {
        Gradient result = new Gradient();
        result.X = ic.exp(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.exp(X),dX[i]);
        }
        return result;
    }

    public Gradient log() {
        Gradient result = new Gradient();
        result.X = ic.log(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.recip(this.X),this.dX[i]);
        }
        return result;
    }

    public Gradient log2() {
        Gradient result = new Gradient();
        result.X = ic.log2(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i],ic.mul(ic.log(ic.numsToInterval(2,2)),this.X));
        }
        return result;
    }

    public Gradient log10() {
        Gradient result = new Gradient();
        result.X = ic.log10(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i],ic.mul(ic.log(ic.numsToInterval(10,10)),this.X));
        }
        return result;
    }

    public SetInterval getX() {
        return X;
    }

    public void setX(SetInterval X) {
        this.X = X;
    }

    public SetInterval[] getDX() {
        return dX;
    }

    public void setDX(SetInterval[] dX) {
        for (int i = 0; i < dim; i++) {
            this.dX[i] = dX[i];
        }
    }

    public static double getDim() {return dim;}

    public void show() {
        System.out.println("[" + this.X.doubleInf() + ", " + this.X.doubleSup() + "]");
        System.out.print("(");
        for (int i = 0; i < dim; i++) {
            System.out.print(" [" + this.dX[i].doubleInf() + ", " + this.dX[i].doubleSup() + "]");
        }
        System.out.println(" )");
    }

}
