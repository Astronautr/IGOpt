import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.ExtendedRational;
import net.java.jinterval.rational.Rational;

import java.util.Arrays;

public class GradientMV {
    static SetInterval[] box;
    static SetInterval[] boxBC;
    static SetIntervalContext ic;
    static int dim;
    private SetInterval xCentral;
    private SetInterval[] dxCentral;
    private SetInterval X;
    private SetInterval[] dX;
    private SetInterval[][] ddX;

    private static final SetInterval ZERO = SetIntervalContexts.getExact().numsToInterval(0,0);
    private static final SetInterval ONE = SetIntervalContexts.getExact().numsToInterval(1, 1);

    private void hessianFill() {
        for (int i = 1; i < dim; i++)
            for (int j = i-1; j >= 0; j--) {
                this.ddX[i][j] = this.ddX[j][i];
            }
    }

    private void improving() {
        for (int i = 0; i < dim; i++) {
            SetInterval dxiMV = ic.add(this.dxCentral[i],Operations.scalarMul(this.ddX[i],boxBC,ic));
            this.dX[i] = ic.intersection(this.dX[i], dxiMV);
        }
        SetInterval xMV = ic.add(this.xCentral,Operations.scalarMul(this.dX,boxBC,ic));
        this.X = ic.intersection(X,xMV);
    }

    public static GradientMV[] init(SetInterval[] box, SetIntervalContext ic) {
        dim = box.length;
        GradientMV.ic = ic;
        GradientMV.box = box;
        GradientMV[] result = new GradientMV[dim];
        GradientMV.boxBC = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = new GradientMV();
            result[i].xCentral = SetIntervalContexts.getExact().numsToInterval(box[i].mid(), box[i].mid());
            result[i].dxCentral = new SetInterval[dim];
            GradientMV.boxBC[i] = ic.sub(box[i],result[i].xCentral);
            result[i].dX = new SetInterval[dim];
            result[i].ddX = new SetInterval[dim][dim];
            result[i].X = box[i];
            for (SetInterval[] row: result[i].ddX) {
                Arrays.fill(row, ZERO);
            }
            for (int j = 0; j < i; j++) {
                result[i].dX[j] = ZERO;
                result[i].dxCentral[j] = ZERO;
            }
            result[i].dX[i] = ONE;
            result[i].dxCentral[i] = ONE;
            for (int j = i + 1; j < dim; j++) {
                result[i].dX[j] = ZERO;
                result[i].dxCentral[j] = ZERO;
            }
        }
        return result;
    }

    public static GradientMV num(double number) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.numsToInterval(number,number);
        result.dxCentral = new SetInterval[dim];
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row: result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
            result.dxCentral[i] = ZERO;
        }
        result.X = SetIntervalContexts.getExact().numsToInterval(number,number);
        return result;
    }

    public static GradientMV num(ExtendedRational number) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.numsToInterval(number,number);
        result.dxCentral = new SetInterval[dim];
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row: result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
            result.dxCentral[i] = ZERO;
        }
        result.X = SetIntervalContexts.getExact().numsToInterval(number,number);
        return result;
    }

    public static GradientMV nums(double lower,double upper) {
        GradientMV result = new GradientMV();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row: result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower,upper);
        return result;
    }

    public static GradientMV nums(Rational lower, Rational upper) {
        GradientMV result = new GradientMV();
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row: result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
        }
        result.X = ic.numsToInterval(lower,upper);
        return result;
    }

    public static GradientMV nums(ExtendedRational inf, ExtendedRational sup) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.numsToInterval(inf,sup);
        result.dxCentral = new SetInterval[dim];
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (SetInterval[] row: result.ddX) {
            Arrays.fill(row, ZERO);
        }
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ZERO;
            result.dxCentral[i] = ZERO;
        }
        result.X = SetIntervalContexts.getExact().numsToInterval(inf,sup);
        return result;
    }

    public GradientMV neg() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.neg(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        result.X = ic.neg(X);
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(dX[i]);
            result.dxCentral[i] = ic.neg(dxCentral[i]);
            for (int j = 0; j < dim; j++) {
                result.ddX[i][j]=ic.neg(ddX[i][j]);
            }
        }
        return result;
    }

    public GradientMV max(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.max(this.xCentral, Y.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.max(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        if (this.xCentral.strictPrecedes(Y.xCentral)) {
            for (int i = 0; i < dim; i++) {
                result.dxCentral[i] = Y.dxCentral[i];
            }
        } else if (Y.xCentral.strictPrecedes(this.xCentral)) {
            for (int i = 0; i < dim; i++) {
                result.dxCentral[i] = this.dxCentral[i];
            }
        } else {
            for (int i = 0; i < dim; i++) {
                result.dxCentral[i] = ic.convexHull(this.dxCentral[i], Y.dxCentral[i]);
            }
        }
        if (this.X.strictPrecedes(Y.X)) {
            for (int i = 0; i < dim; i++) {
                result.dX[i] = Y.dX[i];
                for (int j = i; j < dim; j++) {
                    result.ddX[i][j] = Y.ddX[i][j];
                }
            }
        } else if (Y.X.strictPrecedes(this.X)) {
            for (int i = 0; i < dim; i++) {
                result.dX[i] = this.dX[i];
                for (int j = i; j < dim; j++) {
                    result.ddX[i][j] = this.ddX[i][j];
                }
            }
        } else {
            for (int i = 0; i < dim; i++) {
                result.dX[i] = ic.convexHull(this.dX[i], Y.dX[i]);
                for (int j = i; j < dim; j++) {
                    result.ddX[i][j] = ic.entire();
                }
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV abs() {
        //NotIsRealAbs
/*        GradientMV result = new GradientMV();
        result.xCentral = this.xCentral;
        result.dxCentral = new SetInterval[dim];
        result.X = ic.intersection(this.X,ic.numsToInterval(1, Double.POSITIVE_INFINITY));
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = *//*ic.convexHull(ic.neg(this.dX[i]),this.dX[i]);*//*this.dX[i];
            result.dxCentral[i] = *//*ic.convexHull(ic.neg(this.dX[i]),this.dX[i]);*//*this.dxCentral[i];
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = this.ddX[i][j];
            }
        }
        result.hessianFill();
        return result;*/
        GradientMV result = new GradientMV();
        result.X = ic.abs(this.X);
        result.xCentral = ic.abs(this.xCentral);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
//            result.dX[i] = /*ic.convexHull(ic.neg(this.dX[i]),this.dX[i]);*/this.dX[i];
            result.dX[i] = ic.abs(this.dX[i]);
            result.dxCentral[i] = ic.abs(this.dxCentral[i]);
//            result.dxCentral[i] = /*ic.convexHull(ic.neg(this.dX[i]),this.dX[i]);*/this.dxCentral[i];
            for (int j = i; j < dim; j++) {
//                result.ddX[i][j] = this.ddX[i][j];
                result.ddX[i][j] = ic.abs(this.ddX[i][j]);
            }
        }
        result.hessianFill();
        return result;
    }

    public GradientMV add(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.add(this.xCentral,Y.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.add(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(this.dX[i],Y.dX[i]);
            result.dxCentral[i] = ic.add(this.dxCentral[i],Y.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(this.ddX[i][j],Y.ddX[i][j]);
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV sub(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.sub(this.xCentral,Y.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.sub(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.sub(this.dX[i],Y.dX[i]);
            result.dxCentral[i] = ic.sub(this.dxCentral[i],Y.dxCentral[i]);

            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.sub(this.ddX[i][j],Y.ddX[i][j]);
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV mul(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.mul(this.xCentral,Y.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.mul(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(ic.mul(this.X, Y.dX[i]),ic.mul(this.dX[i],Y.X));
            result.dxCentral[i] = ic.add(ic.mul(this.xCentral, Y.dxCentral[i]),ic.mul(this.dxCentral[i],Y.xCentral));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.add(ic.mul(this.ddX[i][j],Y.X),ic.add(ic.mul(this.dX[i],Y.dX[j]),ic.mul(this.dX[j],Y.dX[i])))
                        ,ic.mul(this.X,Y.ddX[i][j]));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV div(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.div(this.xCentral,Y.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.div(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.sub(ic.div(this.dX[i],Y.X),ic.div(ic.mul(this.X,Y.dX[i]),ic.sqr(Y.X)));
            result.dxCentral[i] = ic.sub(ic.div(this.dxCentral[i],Y.xCentral),ic.div(ic.mul(this.xCentral,Y.dxCentral[i]),ic.sqr(Y.xCentral)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] =ic.add(ic.sub(ic.div(this.ddX[i][j],Y.X),
                        ic.div(
                                ic.add(ic.add(ic.mul(this.dX[i],Y.dX[j]),ic.mul(this.dX[j],Y.dX[i])),ic.mul(this.X,Y.ddX[i][j])),ic.sqr(Y.X))),
                        ic.mul(ic.numsToInterval(2,2),ic.div(ic.mul(this.X,ic.mul(Y.dX[i],Y.dX[j])),ic.pown(Y.X,3))));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV recip() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.recip(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.recip(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.div(this.dX[i],ic.sqr(this.X)));
            result.dxCentral[i] = ic.neg(ic.div(this.dxCentral[i],ic.sqr(xCentral)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.sub(
                        ic.mul(ic.numsToInterval(2,2),ic.div(ic.mul(this.dX[i],this.dX[j]),ic.pown(this.X,3))),
                        ic.div(this.ddX[i][j],ic.sqr(this.X)));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV pow(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.xCentral = ic.pow(this.xCentral,Y.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.pow(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y.X, ic.pow(X,ic.sub(Y.X, ic.numsToInterval(1,1)))),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.mul(Y.xCentral, ic.pow(xCentral,ic.sub(Y.xCentral, ic.numsToInterval(1,1)))),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.mul(Y.X,ic.pow(this.X,ic.sub(Y.X,ic.numsToInterval(2,2)))),
                        ic.add(ic.mul(ic.mul(ic.sub(Y.X,ic.numsToInterval(1,1)),this.dX[i]),this.dX[j]),
                                ic.mul(this.X,this.ddX[i][j])));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV pown(int n) {
        GradientMV result = new GradientMV();
        result.X = ic.pown(this.X, n);
        result.xCentral = ic.pown(this.xCentral,n);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n,n), ic.pown(X,n-1)),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.mul(ic.numsToInterval(n,n), ic.pown(xCentral,n-1)),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.mul(ic.numsToInterval(n,n),ic.numsToInterval(n-1,n-1)),ic.pown(this.X,n-2));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV pown(long n) {
        GradientMV result = new GradientMV();
        result.X = ic.pown(this.X, n);
        result.xCentral = ic.pown(this.xCentral,n);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n,n), ic.pown(X,n-1)),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.mul(ic.numsToInterval(n,n), ic.pown(xCentral,n-1)),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.mul(ic.numsToInterval(n,n),ic.numsToInterval(n-1,n-1)),ic.pown(this.X,n-2));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV rootn(int n) {
        GradientMV result = new GradientMV();
        result.X = ic.rootn(this.X, n);
        SetInterval Y = ic.recip(ic.numsToInterval(n,n));
        result.xCentral = ic.rootn(this.xCentral,n);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y, ic.pow(X,ic.sub(Y, ONE))),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.mul(Y, ic.pow(this.xCentral,ic.sub(Y, ONE))),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.mul(Y,ic.sub(Y,ONE)),ic.mul(ic.pow(this.X,ic.sub(Y,ic.numsToInterval(2,2))),ic.sqr(this.dX[i]))),
                        ic.mul(Y,ic.mul(ic.pow(this.X,ic.sub(Y,ONE)),this.ddX[i][j])));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV rootn(long n) {
        GradientMV result = new GradientMV();
        result.X = ic.rootn(this.X, n);
        SetInterval Y = ic.recip(ic.numsToInterval(n,n));
        result.xCentral = ic.rootn(this.xCentral,n);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y, ic.pow(X,ic.sub(Y, ONE))),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.mul(Y, ic.pow(this.xCentral,ic.sub(Y, ONE))),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.mul(Y,ic.sub(Y,ONE)),ic.mul(ic.pow(this.X,ic.sub(Y,ic.numsToInterval(2,2))),ic.sqr(this.dX[i]))),
                        ic.mul(Y,ic.mul(ic.pow(this.X,ic.sub(Y,ONE)),this.ddX[i][j])));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV sqr() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.sqr(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.sqr(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(2,2),X),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.mul(ic.numsToInterval(2,2),xCentral),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.mul(ic.numsToInterval(2,2),ic.add(ic.mul(this.dX[j],this.dX[i]),ic.mul(this.X,this.ddX[i][j])));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV sqrt() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.sqrt(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.sqrt(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(dX[i], ic.mul(ic.numsToInterval(2,2),ic.sqrt(X)));
            result.dxCentral[i] = ic.div(dxCentral[i], ic.mul(ic.numsToInterval(2,2),ic.sqrt(xCentral)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(ic.mul(this.ddX[i][j],ic.numsToInterval(2,2)),this.X),ic.mul(this.dX[i],this.dX[j])),
                        ic.mul(ic.numsToInterval(4,4),ic.sqrt(ic.pown(this.X,3))));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV sin() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.sin(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.sin(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.cos(X),dX[i]);
            result.dxCentral[i] = ic.mul(ic.cos(xCentral),dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.mul(ic.neg(ic.sin(this.X)),this.dX[i]),this.dX[j]),
                        ic.mul(ic.cos(this.X),this.ddX[i][j]));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV cos() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.cos(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.cos(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.mul(ic.sin(X),dX[i]));
            result.dxCentral[i] = ic.neg(ic.mul(ic.sin(xCentral),dxCentral[i]));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.neg(ic.add(ic.mul(ic.mul((ic.cos(this.X)),this.dX[i]),this.dX[j]),
                        ic.mul(ic.sin(this.X),this.ddX[i][j])));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV tan() {
        GradientMV result = new GradientMV();
        result.X = ic.tan(this.X);
        result.xCentral = ic.tan(this.xCentral);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i],ic.sqr(ic.cos(this.X)));
            result.dxCentral[i] = ic.div(this.dxCentral[i],ic.sqr(ic.cos(this.xCentral)));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.sub(ic.div(this.ddX[i][j],ic.sqr(ic.cos(this.X))),
                        ic.div(ic.mul(ic.sin(ic.mul(ic.numsToInterval(2,2),this.X)),ic.sqr(this.dX[i])),ic.pown(ic.cos(this.X),4)));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV asin() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.asin(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.asin(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.recip(ic.sqrt(ic.sub(ONE,ic.sqr(this.X)))),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.recip(ic.sqrt(ic.sub(ONE,ic.sqr(this.xCentral)))),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.sub(ic.mul(ic.recip(ic.sqrt(ic.sub(ONE,ic.sqr(this.X)))),this.ddX[i][j]),
                        ic.mul(ic.div(this.X,ic.pow(ic.sub(ONE,ic.sqr(this.X)),ic.numsToInterval(1.5,1.5))),ic.sqr(this.dX[i])));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV acos() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.acos(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.acos(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.div(dX[i],ic.sqrt(ic.sub(ic.numsToInterval(1,1),ic.sqr(X)))));
            result.dxCentral[i] = ic.neg(ic.div(dxCentral[i],ic.sqrt(ic.sub(ic.numsToInterval(1,1),ic.sqr(xCentral)))));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.neg(ic.div(ic.add(ic.mul(this.ddX[i][j],ic.sub(ic.numsToInterval(1,1),ic.sqr(this.X))),
                        ic.mul(ic.mul(this.X,this.dX[i]),this.dX[j])),
                        ic.pow(ic.sub(ic.numsToInterval(1,1),ic.sqr(this.X)),ic.numsToInterval(1.5,1.5))));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV atan() {
        GradientMV result = new GradientMV();
        result.xCentral = ic.atan(this.xCentral);
        result.dxCentral = new SetInterval[dim];
        result.X = ic.atan(this.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1,1),ic.add(ic.numsToInterval(1,1),ic.sqr(X))),dX[i]);
            result.dxCentral[i] = ic.mul(ic.div(ic.numsToInterval(1,1),ic.add(ic.numsToInterval(1,1),ic.sqr(xCentral))),dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j],ic.add(ic.numsToInterval(1,1),ic.sqr(this.X))),
                        ic.mul(ic.mul(ic.mul(ic.numsToInterval(2,2),this.X),this.dX[i]),this.dX[j])),
                        ic.sqr(ic.add(ic.numsToInterval(1,1),ic.sqr(this.X))));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV exp() {
        GradientMV result = new GradientMV();
        result.X = ic.exp(this.X);
        result.xCentral = ic.exp(this.xCentral);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.exp(this.X),this.dX[i]);
            result.dxCentral[i] = ic.mul(ic.exp(this.xCentral),this.dxCentral[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.add(ic.mul(ic.exp(this.X),ic.sqr(this.dX[i])),ic.mul(ic.exp(this.X),this.ddX[i][j]));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV log() {
        GradientMV result = new GradientMV();
        result.X = ic.log(this.X);
        result.xCentral = ic.log(this.xCentral);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i],this.X);
            result.dxCentral[i] = ic.div(this.dxCentral[i],this.xCentral);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j],this.X),ic.sqr(this.dX[i])),ic.sqr(this.X));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV log2() {
        GradientMV result = new GradientMV();
        result.X = ic.log2(this.X);
        result.xCentral = ic.log2(this.xCentral);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i],ic.mul(this.X,ic.log(ic.numsToInterval(2,2))));
            result.dxCentral[i] = ic.div(this.dxCentral[i],ic.mul(this.xCentral,ic.log(ic.numsToInterval(2,2))));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j],ic.mul(this.X,ic.log(ic.numsToInterval(2,2)))),
                        ic.mul(ic.sqr(this.dX[i]),ic.log(ic.numsToInterval(2,2)))),ic.sqr(ic.mul(this.X,ic.log(ic.numsToInterval(2,2)))));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV log10() {
        GradientMV result = new GradientMV();
        result.X = ic.log10(this.X);
        result.xCentral = ic.log10(this.xCentral);
        result.dX = new SetInterval[dim];
        result.dxCentral = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(this.dX[i],ic.mul(this.X,ic.log(ic.numsToInterval(10,10))));
            result.dxCentral[i] = ic.div(this.dxCentral[i],ic.mul(this.xCentral,ic.log(ic.numsToInterval(10,10))));
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.div(ic.sub(ic.mul(this.ddX[i][j],ic.mul(this.X,ic.log(ic.numsToInterval(10,10)))),
                        ic.mul(ic.sqr(this.dX[i]),ic.log(ic.numsToInterval(10,10)))),ic.sqr(ic.mul(this.X,ic.log(ic.numsToInterval(10,10)))));
            }
        }
        result.hessianFill();
        result.improving();
        return result;
    }

    public GradientMV intersection(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i],Y.dX[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = ic.intersection(this.ddX[i][j],Y.ddX[i][j]);
            }
        }
        result.hessianFill();
        return result;
    }

    public GradientMV intersection(Gradient Y) {
        GradientMV result = new GradientMV();
        result.X = ic.intersection(this.X, Y.getX());
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i],Y.getDX()[i]);
            for (int j = i; j < dim; j++) {
                result.ddX[i][j] = this.ddX[i][j];
            }
        }
        result.hessianFill();
        return result;
    }

    public GradientMV intersectionX(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        result.dX = this.dX.clone();
        result.ddX = this.ddX.clone();
        return result;
    }

    public GradientMV intersectionX(SetInterval Y) {
        GradientMV result = new GradientMV();
        result.X = ic.intersection(this.X, Y);
        result.dX = new SetInterval[dim];
        result.ddX = new SetInterval[dim][dim];
        result.dX = this.dX.clone();
        result.ddX = this.ddX.clone();
        return result;
    }

    public GradientMV intersectionDX(GradientMV Y) {
        GradientMV result = new GradientMV();
        result.X = this.X;
        result.dX = new SetInterval[dim];
        result.ddX = this.ddX.clone();
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public GradientMV intersectionDX(SetInterval[] Y) {
        GradientMV result = new GradientMV();
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