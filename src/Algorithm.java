import net.java.jinterval.expression.AbstractExpressionVisitor;
import net.java.jinterval.expression.CodeList;
import net.java.jinterval.expression.Expression;
import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.interval.set.SetIntervalEvaluator;
import net.java.jinterval.rational.BinaryValueSet;
import net.java.jinterval.rational.ExtendedRational;
import net.java.jinterval.rational.ExtendedRationalContext;
import net.java.jinterval.rational.ExtendedRationalContexts;

import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;

public class Algorithm {
    static PriorityQueue<ListElem> wList = new PriorityQueue<ListElem>(new AssessmentComp());
    private static SetIntervalContext ic;
    private static ExtendedRationalContext rc;
    private static CodeList list;
    private static SetIntervalEvaluator setEv;
    private static GradientEvaluator grEv;
    private static GradientSDEvaluator grSDEv;
    private static GradientMVEvaluator grMVEv;
    private static Expression obj;
    private static SetInterval[] initialBox;
    private static ExtendedRational tolerance;
    private static ExtendedRational supMin;
    Algorithm(SetInterval[] box, SetIntervalContext ic, ExtendedRational tolerance, String[] inps, int exprNum) {
        initialBox = box.clone();
        Algorithm.tolerance = tolerance;
        Functions func = new Functions(inps,exprNum);
        Algorithm.ic = ic;
        rc = ExtendedRationalContexts.mkNearest(BinaryValueSet.BINARY64);
        obj = func.getObjective();
        list = obj.getCodeList();
        supMin = ExtendedRational.POSITIVE_INFINITY;
        setEv = SetIntervalEvaluator.create(Algorithm.ic, list, obj);
        grEv = new GradientEvaluator(Algorithm.ic, list ,obj);
        grSDEv = new GradientSDEvaluator(Algorithm.ic, list, obj);
        grMVEv = new GradientMVEvaluator(Algorithm.ic, list, obj);
    }

    public Gradient meanValue(SetInterval[] box) {
        Gradient grObjective = grEv.evaluate(Gradient.init(box,ic))[0];
        SetInterval[] centralPoint = new SetInterval[box.length];
        SetInterval[] bias = new SetInterval[box.length];
        for (int i = 0; i < box.length; i++) {
            centralPoint[i] = ic.numsToInterval(box[i].mid(),box[i].mid());
            bias[i] = ic.sub(box[i], centralPoint[i]);
        }
        SetInterval centralVal = setEv.evaluate(centralPoint)[0];
        grObjective.setX(ic.intersection(ic.add(centralVal,Operations.scalarMul(bias,grObjective.getDX(),ic)),grObjective.getX()));
        return grObjective;
    }

    public GradientMV meanValueSO(SetInterval[] box) {
        SetInterval[] centralPoint = new SetInterval[box.length];
        SetInterval[] bias = new SetInterval[box.length];
        for (int i = 0; i < box.length; i++) {
            centralPoint[i] = ic.numsToInterval(box[i].mid(),box[i].mid());
            bias[i] = ic.sub(box[i], centralPoint[i]);
        }
        SetInterval centralVal = setEv.evaluate(centralPoint)[0];
        Gradient grObjective = grEv.evaluate(Gradient.init(box,ic))[0];
        GradientMV objectiveFunctionMV = grMVEv.evaluate(GradientMV.init(box, ic))[0];
        GradientMV result = objectiveFunctionMV.intersection(grObjective);
        result.setX(ic.intersection(ic.add(centralVal,Operations.scalarMul(result.getDX(),bias,ic)),result.getX()));
        return result;
    }

    public SetInterval[] krawczyk(SetInterval[] box, GradientMV grObjective) {
        boolean isSubPropBox = true;
        for (int i = 0; i < box.length; i++) {
            if (box[i].inf().eq(initialBox[i].inf()) || (box[i].sup().eq(initialBox[i].sup()))) {
                isSubPropBox = false;
                break;
            }
        }
        if (isSubPropBox) {
            for (int i = 0; i < box.length; i++) {
                if (ic.intersection(grObjective.getDX()[i], ic.numsToInterval(0, 0)).isEmpty()) {
                    return null;
                }
            }
            Matrix jacoby = new Matrix(grObjective.getDDX(), ic);
            Matrix lambda = jacoby.mid().inverse();
            SetInterval[] centralPoint = new SetInterval[box.length];
            SetInterval[] bias = new SetInterval[box.length];
            for (int i = 0; i < box.length; i++) {
                centralPoint[i] = ic.numsToInterval(box[i].mid(), box[i].mid());
                bias[i] = ic.sub(box[i], centralPoint[i]);
            }
            GradientMV centralVal = grMVEv.evaluate((GradientMV.init(centralPoint, ic)))[0];
            Matrix mXBiasCentral = new Matrix(bias, ic);
            Matrix mObjFuncCentrPointVal = new Matrix(centralVal.getDX(), ic);
            SetInterval[] krawczyk = ((mXBiasCentral.sub(lambda.mul(mObjFuncCentrPointVal))).
                    add((Matrix.identity(box.length).sub(lambda.mul(jacoby))).mul(mXBiasCentral))).getVector();
            for (int i = 0; i < box.length; i++) {
                box[i] = ic.intersection(box[i], krawczyk[i]);
                if (box[i].isEmpty())
                    return null;
            }
            return box;
        }
        return box;
    }

    public SetInterval[] monotonyCheck(SetInterval[] box, Gradient result) {
        for (int i = 0; i < box.length; i++) {
            if (result.getDX()[i].inf().ge(ExtendedRational.zero())) {
                box[i] = SetIntervalContexts.getExact().numsToInterval(box[i].sup(),box[i].sup());
            }
            else
            if (result.getDX()[i].sup().le(ExtendedRational.zero())) {
                box[i] = SetIntervalContexts.getExact().numsToInterval(box[i].inf(),box[i].inf());
            }
        }
        return box;
    }

    public SetInterval[] monotonyCheck(SetInterval[] box, GradientMV result) {
        for (int i = 0; i < box.length; i++) {
            if (result.getDX()[i].inf().ge(ExtendedRational.zero())) {
                box[i] = SetIntervalContexts.getExact().numsToInterval(box[i].sup(),box[i].sup());
            }
            else
            if (result.getDX()[i].sup().le(ExtendedRational.zero())) {
                box[i] = SetIntervalContexts.getExact().numsToInterval(box[i].inf(),box[i].inf());
            }
        }
        return box;
    }

    public void cleaning(SetInterval result) {
        if (wList.size() >1) {
            if (supMin.eq(ExtendedRational.POSITIVE_INFINITY) && (result.sup().eq(ExtendedRational.POSITIVE_INFINITY)))
                return;
            if (rc.abs(rc.sub(supMin, result.sup())).ge(tolerance)) {
                Operations.rmUseless(supMin, wList);
                supMin = result.sup();
            }
        }
    }

    public int basicChooseRule(SetInterval[] box) {
        int max = 0;
        for (int i = 0; i < box.length; i++) {
            if (box[i].wid().gt(box[max].wid()))
                max = i;
        }
        return max;
    }

    public int altChooseRule(SetInterval[] box, SetInterval[] partials) {
        if (partials != null) {
            int pMax = 0;
            for (int i = 0; i < partials.length; i++) {
                if (partials[i].wid().gt(partials[pMax].wid()))
                    pMax = i;
            }
            if (partials[pMax].wid().lt(ExtendedRational.POSITIVE_INFINITY)) {
                for (int i = 0; i < box.length; i++) {
                    if ((rc.mul(ic.abs(partials[i]).sup(),box[i].wid())).gt(rc.mul(ic.abs(partials[pMax]).sup(),box[pMax].wid()))) {
                        pMax = i;
                    }
                }
                return pMax;
            }
            else
                return basicChooseRule(box);
        }
        else return basicChooseRule(box);
    }

    void optStep(SetInterval[] box, String[] keys) {
        if (keys[0].equals("-s")) {
            SetInterval result = setEv.evaluate(box)[0];
            wList.add(new ListElem(box,result.inf(),result.wid(),null));
        }

        if (keys[0].equals("-b")) {
            Gradient result = meanValue(box);
            box = monotonyCheck(box,result).clone();
            wList.add(new ListElem(box,result.getX().inf(),result.getX().wid(),result.getDX()));
        }

        if (keys[0].equals("-a")) {
            GradientMV objectiveMV = meanValueSO(box);
            box = krawczyk(box,objectiveMV);
            if (box == null) {
                return;
            }
            box = monotonyCheck(box,objectiveMV).clone();
            cleaning(objectiveMV.getX());
            wList.add(new ListElem(box,objectiveMV.getX().inf(),objectiveMV.getX().wid(),objectiveMV.getDX()));
        }
    }

    public ExtendedRational start(String[] keys) {
        long steps = 0;
        supMin = ExtendedRational.POSITIVE_INFINITY;
        optStep(initialBox, keys);
        steps++;
        while (wList.peek().getWid().ge(tolerance)) {
            int max;
            SetInterval[] box = wList.peek().getData();
            SetInterval[] partials = wList.poll().getPartials();
            if ((keys[0].equals("-a") || (keys[0].equals("-b"))))
                max = altChooseRule(box,partials);
            else
                max = basicChooseRule(box);
            SetInterval[] fDescendant = box.clone();
            SetInterval[] sDescendant = box.clone();
            fDescendant[max] = ic.numsToInterval(fDescendant[max].inf(),fDescendant[max].mid());
            sDescendant[max] = ic.numsToInterval(sDescendant[max].mid(),sDescendant[max].sup());
            optStep(fDescendant, keys);
            optStep(sDescendant, keys);
            steps++;
//            System.out.println("Assessment = " + wList.peek().getAssessment().doubleValue());
        }
        System.out.println("Number of algorithm steps = " + steps);
        System.out.println("Assessment = " + wList.peek().getAssessment().doubleValue());
        return wList.peek().getAssessment();
    }
}

class AssessmentComp implements Comparator<ListElem> {
    public int compare(ListElem el1, ListElem el2) {
        return ExtendedRational.valueOf(el1.getAssessment()).compareTo(el2.getAssessment());
    }
}

class ListElem implements Comparable<ListElem> {

    private SetInterval[] data;
    private ExtendedRational assessment;
    private ExtendedRational wid;
    private SetInterval[] partials;

    ListElem() {
        data = null;
        assessment = ExtendedRational.valueOf(0.0);
        wid = ExtendedRational.valueOf(0.0);
        Arrays.fill(partials, SetIntervalContexts.getExact().nai());
    }

    ListElem(SetInterval[] data, ExtendedRational assessment, ExtendedRational wid, SetInterval[] partials) {
        this.data = data;
        this.assessment = assessment;
        this.wid = wid;
        this.partials = partials;
    }

    public ExtendedRational getAssessment() {
        return assessment;
    }

    public ExtendedRational getWid() {
        return wid;
    }

    public SetInterval[] getData() {
        return data;
    }

    public SetInterval[] getPartials() {
        return partials;
    }

    public int compareTo(ListElem el) {
        return ExtendedRational.valueOf(this.assessment).compareTo(ExtendedRational.valueOf(el.assessment));
    }
}