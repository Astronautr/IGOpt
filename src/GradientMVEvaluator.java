import net.java.jinterval.expression.AbstractExpressionVisitor;
import net.java.jinterval.expression.CodeList;
import net.java.jinterval.expression.Expression;
import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.rational.Rational;


public class GradientMVEvaluator extends AbstractExpressionVisitor {
    final CodeList list;
    final Expression[] results;
    private final GradientMV[] values;
    private final SetIntervalContext ic;

    GradientMVEvaluator(SetIntervalContext ic, CodeList list, Expression... results) {
        this.list = list;
        this.results = results.clone();
        values = new GradientMV[list.getNumExprs()];
        this.ic = ic;
    }

    public GradientMV[] evaluate(GradientMV[] args) {
        if (args.length != list.getNumInps()) {
            throw new IllegalArgumentException();
        }
        System.arraycopy(args, 0, this.values, 0, list.getNumInps());
        list.acceptConstants(this);
        list.acceptForward(this);
        GradientMV[] retVal = new GradientMV[results.length];
        for (int i = 0; i < results.length; i++) {
            retVal[i] = this.values[results[i].getIndex()];
        }
        return retVal;
    }

    @Override
    public void visitInp(int r, String name) {
    }

    @Override
    public void visitLit(int r, String numerator, String denominator) {
        Rational rv = Rational.valueOf(numerator + "/" + denominator);
        SetInterval v = ic.numsToInterval(rv, rv);
        values[r] = GradientMV.nums(v.inf(), v.sup());
    }

    @Override
    public void visitLit(int r, String literal) {
        Rational rv = Rational.valueOf(literal);
        SetInterval v = ic.numsToInterval(rv, rv);
        values[r] = GradientMV.nums(v.inf(), v.sup());
    }

    @Override
    public void visitNum(int r, Number num) {
        Rational rv = Rational.valueOf(num);
        SetInterval v = ic.numsToInterval(rv, rv);
        values[r] = GradientMV.nums(v.inf(), v.sup());
    }

    @Override
    public void visitAbs(int r, int x) {
        values[r] = values[x].abs();
    }

    @Override
    public void visitPi(int r) {
        SetInterval v = ic.pi();
        values[r] = GradientMV.nums(v.inf(), v.sup());
    }

    @Override
    public void visitEuler(int r) {
        SetInterval v = ic.euler();
        values[r] = GradientMV.nums(v.inf(), v.sup());
    }

    @Override
    public void visitNeg(int r, int x) {
        values[r] = values[x].neg();
    }

    @Override
    public void visitAdd(int r, int x, int y) {
        values[r] = values[x].add(values[y]);
    }

    @Override
    public void visitSub(int r, int x, int y) {
        values[r] = values[x].sub(values[y]);
    }

    @Override
    public void visitMul(int r, int x, int y) {
        values[r] = values[x].mul(values[y]);
    }

    @Override
    public void visitDiv(int r, int x, int y) {
        values[r] = values[x].div(values[y]);
    }

    @Override
    public void visitRecip(int r, int x) {
        values[r] = values[x].recip();
    }

    @Override
    public void visitSqr(int r, int x) {
        values[r] = values[x].sqr();
    }

    @Override
    public void visitSqrt(int r, int x) {
        values[r] = values[x].sqrt();
    }

    @Override
    public void visitPow(int r, int x, int y) {
        if (!list.getExpr(y).isConstant()) {
            throw new UnsupportedOperationException();
        }
        values[r] = values[x].pow(values[y]);
    }

    @Override
    public void visitPown(int r, int x, int n) {
        values[r] = values[x].pown(n);
    }

    @Override
    public void visitPown(int r, int x, long n) {
        values[r] = values[x].pown(n);
    }

    @Override
    public void visitLog(int r, int x) {
        values[r] = values[x].log();
    }

    @Override
    public void visitLog2(int r, int x) {
        values[r] = values[x].log2();
    }

    @Override
    public void visitLog10(int r, int x) {
        values[r] = values[x].log10();
    }

    @Override
    public void visitExp(int r, int x) {
        values[r] = values[x].exp();
    }

    @Override
    public void visitSin(int r, int x) {
        values[r] = values[x].sin();
    }

    @Override
    public void visitCos(int r, int x) {
        values[r] = values[x].cos();
    }

    @Override
    public void visitTan(int r, int x) {
        values[r] = values[x].tan();
    }

    @Override
    public void visitAsin(int r, int x) {
        values[r] = values[x].asin();
    }

    @Override
    public void visitAcos(int r, int x) {
        values[r] = values[x].acos();
    }

    @Override
    public void visitAtan(int r, int x) {
        values[r] = values[x].atan();
    }

    @Override
    public void visitMax(int r, int x, int y) {
        values[r] = values[x].max(values[y]);
    }

    @Override
    public void visitRootn(int r, int x, int q) {
        values[r] = values[x].rootn(q);
    }

    @Override
    public void visitRootn(int r, int x, long q) {
        values[r] = values[x].rootn(q);
    }

/*    void setContext(SetIntervalContext ic) {
        this.ic = ic;
    }*/
}


