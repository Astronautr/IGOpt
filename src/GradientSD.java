import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.ExtendedRational;
import net.java.jinterval.rational.Rational;

import java.util.Arrays;


public class GradientSD {
    private static final SetInterval ZERO = SetIntervalContexts.getExact().numsToInterval(0, 0);
    private static final SetInterval ONE = SetIntervalContexts.getExact().numsToInterval(1, 1);
    static SetInterval[] box;
    static SetIntervalContext ic;
    static int dim;
    private SetInterval X;
    private SetInterval[] dX;
    private SetInterval[][] ddX;

    public static GradientSD[] init(SetInterval[] box, SetIntervalContext ic) {
        dim = box.length;
        GradientSD.ic = ic;
        GradientSD.box = box;
        GradientSD[] result = new GradientSD[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = new GradientSD();
            result[i].dX = new SetInterval[dim];
            result[i].ddX = new SetInterval[dim][dim];
            result[i].X = box[i];
            for (SetInterval[] row : result[i].ddX) {
                Arrays.fill(row, ZERO);
            }
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

    public static GradientSD num(double number) {
        GradientSD result = new GradientSD();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row : result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = SetIntervalContexts.getExact().numsToInterval(number, number);
        return result;
    }

    public static GradientSD num(Rational number) {
        GradientSD result = new GradientSD();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row : result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = SetIntervalContexts.getExact().numsToInterval(number, number);
        return result;
    }

    public static GradientSD num(ExtendedRational number) {
        GradientSD result = new GradientSD();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row : result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = SetIntervalContexts.getExact().numsToInterval(number, number);
        return result;
    }

    public static GradientSD nums(double lower, double upper) {
        GradientSD result = new GradientSD();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row : result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower, upper);
        return result;
    }

    public static GradientSD nums(Rational lower, Rational upper) {
        GradientSD result = new GradientSD();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row : result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower, upper);
        return result;
    }

    public static GradientSD nums(ExtendedRational lower, ExtendedRational upper) {
        GradientSD result = new GradientSD();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row : result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower, upper);
        return result;
    }

    private void hessianFill() {
        for (int i = 1; i < dim; i++)
            for (int j = i - 1; j >= 0; j--) {
                this.ddX[i][j] = this.ddX[j][i];
            }
    }

    public GradientSD neg() {
        GradientSD result = new GradientSD();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        result.X = ic.neg(X);
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(dX[i]);
            for (int j = 0; j < dim; j++) {
                result.ddX[i][j] = ic.neg(ddX[i][j]);
            }
        }
        return result;
    }

    public GradientSD max(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.max(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.max(this.dX[i], Y.dX[i]);
            for (int j = 0; j < dim; j++) {
                result.ddX[i][j] = ic.max(this.ddX[i][j], Y.ddX[i][j]);
            }
        }
        return result;
    }

    public GradientSD abs() {
        GradientSD result = new GradientSD();
        result.X = ic.abs(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
//            result.dX[i] = /*ic.convexHull(ic.neg(this.dX[i]),this.dX[i]);*/this.dX[i];
            result.dX[i] = ic.abs(this.dX[i]);
//            result.dxCentral[i] = /*ic.convexHull(ic.neg(this.dX[i]),this.dX[i]);*/this.dxCentral[i];
            for (int j = i; j < dim; j++) {
//                result.ddX[i][j] = this.ddX[i][j];
                result.ddX[i][j] = ic.abs(this.ddX[i][j]);
            }
        }
        result.hessianFill();
        return result;
    }

    public GradientSD add(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.add(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(this.dX[i], Y.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(this.ddX[i][j], Y.ddX[i][j]);
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD sub(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.sub(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.sub(this.dX[i], Y.dX[i]);

            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.sub(this.ddX[i][j], Y.ddX[i][j]);
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD mul(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.mul(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(ic.mul(this.X, Y.dX[i]), ic.mul(this.dX[i], Y.X));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.add(ic.mul(this.ddX[i][j], Y.X), ic.add(ic.mul(this.dX[i], Y.dX[j]), ic.mul(this.dX[j], Y.dX[i])))
                        , ic.mul(this.X, Y.ddX[i][j]));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD div(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.div(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.sub(ic.div(this.dX[i], Y.X), ic.div(ic.mul(this.X, Y.dX[i]), ic.sqr(Y.X)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.sub(ic.div(this.ddX[i][j], Y.X),
                        ic.div(
                                ic.add(ic.add(ic.mul(this.dX[i], Y.dX[j]), ic.mul(this.dX[j], Y.dX[i])), ic.mul(this.X, Y.ddX[i][j])), ic.sqr(Y.X))),
                        ic.mul(ic.numsToInterval(2, 2), ic.div(ic.mul(this.X, ic.mul(Y.dX[i], Y.dX[j])), ic.pown(Y.X, 3))));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD recip() {
        GradientSD result = new GradientSD();
        result.X = ic.recip(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.div(this.dX[i], ic.sqr(this.X)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.numsToInterval(2, 2), ic.pow(this.X, ic.numsToInterval(3, 3)));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD pow(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.pow(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y.X, ic.pow(X, ic.sub(Y.X, ic.numsToInterval(1, 1)))), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.mul(Y.X, ic.pow(this.X, ic.sub(Y.X, ic.numsToInterval(2, 2)))),
                        ic.add(ic.mul(ic.mul(ic.sub(Y.X, ic.numsToInterval(1, 1)), this.dX[i]), this.dX[j]),
                                ic.mul(this.X, this.ddX[i][j])));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD pown(int n) {
        GradientSD result = new GradientSD();
        result.X = ic.pown(this.X, n);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n, n), ic.pown(X, n - 1)), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.mul(ic.numsToInterval(n, n), ic.numsToInterval(n - 1, n - 1)), ic.pown(this.X, n - 2));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD pown(long n) {
        GradientSD result = new GradientSD();
        result.X = ic.pown(this.X, n);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n, n), ic.pown(X, n - 1)), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.mul(ic.numsToInterval(n, n), ic.numsToInterval(n - 1, n - 1)), ic.pown(this.X, n - 2));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD rootn(int n) {
        GradientSD result = new GradientSD();
        result.X = ic.rootn(this.X, n);
        SetInterval Y = ic.recip(ic.numsToInterval(n, n));
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y, ic.pow(X, ic.sub(Y, ONE))), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.mul(Y, ic.sub(Y, ONE)), ic.mul(ic.pow(this.X, ic.sub(Y, ic.numsToInterval(2, 2))), ic.sqr(this.dX[i]))),
                        ic.mul(Y, ic.mul(ic.pow(this.X, ic.sub(Y, ONE)), this.ddX[i][j])));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD rootn(long n) {
        GradientSD result = new GradientSD();
        result.X = ic.rootn(this.X, n);
        SetInterval Y = ic.recip(ic.numsToInterval(n, n));
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y, ic.pow(X, ic.sub(Y, ONE))), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.mul(Y, ic.sub(Y, ONE)), ic.mul(ic.pow(this.X, ic.sub(Y, ic.numsToInterval(2, 2))), ic.sqr(this.dX[i]))),
                        ic.mul(Y, ic.mul(ic.pow(this.X, ic.sub(Y, ONE)), this.ddX[i][j])));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD sqr() {
        GradientSD result = new GradientSD();
        result.X = ic.sqr(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(2, 2), X), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.numsToInterval(2, 2), ic.add(ic.mul(this.dX[j], this.dX[i]), ic.mul(this.X, this.ddX[i][j])));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD sqrt() {
        GradientSD result = new GradientSD();
        result.X = ic.sqrt(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(dX[i], ic.mul(ic.numsToInterval(2, 2), ic.sqrt(X)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(ic.mul(this.ddX[i][j], ic.numsToInterval(2, 2)), this.X), ic.mul(this.dX[i], this.dX[j])),
                        ic.mul(ic.numsToInterval(4, 4), ic.sqrt(ic.pown(this.X, 3))));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD sin() {
        GradientSD result = new GradientSD();
        result.X = ic.sin(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.cos(X), dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.mul(ic.neg(ic.sin(this.X)), this.dX[i]), this.dX[j]),
                        ic.mul(ic.cos(this.X), this.ddX[i][j]));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD cos() {
        GradientSD result = new GradientSD();
        result.X = ic.cos(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.mul(ic.sin(X), dX[i]));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.neg(ic.add(ic.mul(ic.mul((ic.cos(this.X)), this.dX[i]), this.dX[j]),
                        ic.mul(ic.sin(this.X), this.ddX[i][j])));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD tan() {
        GradientSD result = new GradientSD();
        result.X = ic.tan(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i], ic.sqr(ic.cos(this.X)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.sub(ic.div(this.ddX[i][j], ic.sqr(ic.cos(this.X))),
                        ic.div(ic.mul(ic.sin(ic.mul(ic.numsToInterval(2, 2), this.X)), ic.sqr(this.dX[i])), ic.pown(ic.cos(this.X), 4)));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD asin() {
        GradientSD result = new GradientSD();
        result.X = ic.asin(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.recip(ic.sqrt(ic.sub(ONE, ic.sqr(this.X)))), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.sub(ic.mul(ic.recip(ic.sqrt(ic.sub(ONE, ic.sqr(this.X)))), this.ddX[i][j]),
                        ic.mul(ic.div(this.X, ic.pow(ic.sub(ONE, ic.sqr(this.X)), ic.numsToInterval(1.5, 1.5))), ic.sqr(this.dX[i])));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD acos() {
        GradientSD result = new GradientSD();
        result.X = ic.acos(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.div(dX[i], ic.sqrt(ic.sub(ic.numsToInterval(1, 1), ic.sqr(X)))));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.neg(ic.div(ic.add(ic.mul(this.ddX[i][j], ic.sub(ic.numsToInterval(1, 1), ic.sqr(this.X))),
                        ic.mul(ic.mul(this.X, this.dX[i]), this.dX[j])),
                        ic.pow(ic.sub(ic.numsToInterval(1, 1), ic.sqr(this.X)), ic.numsToInterval(1.5, 1.5))));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD atan() {
        GradientSD result = new GradientSD();
        result.X = ic.atan(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1, 1), ic.add(ic.numsToInterval(1, 1), ic.sqr(X))), dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j], ic.add(ic.numsToInterval(1, 1), ic.sqr(this.X))),
                        ic.mul(ic.mul(ic.mul(ic.numsToInterval(2, 2), this.X), this.dX[i]), this.dX[j])),
                        ic.sqr(ic.add(ic.numsToInterval(1, 1), ic.sqr(this.X))));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD exp() {
        GradientSD result = new GradientSD();
        result.X = ic.exp(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.exp(this.X), this.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.exp(this.X), ic.sqr(this.dX[i])), ic.mul(ic.exp(this.X), this.ddX[i][j]));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD log() {
        GradientSD result = new GradientSD();
        result.X = ic.log(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i], this.X);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j], this.X), ic.sqr(this.dX[i])), ic.sqr(this.X));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD log2() {
        GradientSD result = new GradientSD();
        result.X = ic.log2(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i], ic.mul(this.X, ic.log(ic.numsToInterval(2, 2))));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j], ic.mul(this.X, ic.log(ic.numsToInterval(2, 2)))),
                        ic.mul(ic.sqr(this.dX[i]), ic.log(ic.numsToInterval(2, 2)))), ic.sqr(ic.mul(this.X, ic.log(ic.numsToInterval(2, 2)))));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD log10() {
        GradientSD result = new GradientSD();
        result.X = ic.log10(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i], ic.mul(this.X, ic.log(ic.numsToInterval(10, 10))));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j], ic.mul(this.X, ic.log(ic.numsToInterval(10, 10)))),
                        ic.mul(ic.sqr(this.dX[i]), ic.log(ic.numsToInterval(10, 10)))), ic.sqr(ic.mul(this.X, ic.log(ic.numsToInterval(10, 10)))));
            }
        }
        result.hessianFill();

        return result;
    }

    public GradientSD intersection(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.intersection(this.ddX[i][j], Y.ddX[i][j]);
            }
        }
        result.hessianFill();
        return result;
    }

    public GradientSD intersectionX(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        result.dX = this.dX.clone();
        result.ddX = this.ddX.clone();
        return result;
    }

    public GradientSD intersectionX(SetInterval Y) {
        GradientSD result = new GradientSD();
        result.X = ic.intersection(this.X, Y);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        result.dX = this.dX.clone();
        result.ddX = this.ddX.clone();
        return result;
    }

    public GradientSD intersectionDX(GradientSD Y) {
        GradientSD result = new GradientSD();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        result.ddX = this.ddX.clone();
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public GradientSD intersectionDX(SetInterval[] Y) {
        GradientSD result = new GradientSD();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        result.ddX = this.ddX.clone();
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y[i]);
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

    public SetInterval[][] getDDX() {
        return this.ddX;
    }

    public void setDDX(SetInterval[][] ddX) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                this.ddX[i][i] = ddX[i][j];
            }
        }
    }

    public void show() {
        System.out.println("[" + this.X.doubleInf() + ", " + this.X.doubleSup() + "]");
        System.out.print("(");
        for (int i = 0; i < dim; i++) {
            System.out.print(" [" + this.dX[i].doubleInf() + ", " + this.dX[i].doubleSup() + "]");
            //another formatting style
            // System.out.printf(" [%+10.4f, %+10.4f]\t", this.dX[i].doubleInf(), this.dX[i].doubleSup());
        }
        System.out.println(" )");
        for (int i = 0; i < dim; i++) {
            System.out.print("|");
            for (int j = 0; j < dim; j++) {
                System.out.printf(" [%+10.8f, %+10.8f]\t", this.ddX[i][j].doubleInf(), this.ddX[i][j].doubleSup());
            }
            System.out.println(" |");
        }
    }

}