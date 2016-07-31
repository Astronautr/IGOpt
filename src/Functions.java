import net.java.jinterval.expression.CodeList;
import net.java.jinterval.expression.Expression;


public class Functions {
    private CodeList list;
    private Expression[] expr;
    private Expression objectiveFunction;

    Functions(String[] inps, int exprNum) {
        list = CodeList.create(inps);
        expr = new Expression[exprNum];
        expr[0] = list.getInp(0).mul(list.lit("2")).div(list.getInp(1)); //a
        expr[1] = (list.lit("200").div(list.getInp(1)));//b
        expr[2] = ((expr[0].sqr().add((expr[1].add(list.lit("1"))).sqr())).
                sub(list.lit("2").mul(expr[0].mul(expr[1].add(list.lit("1"))).mul(list.getInp(2).sin())))).abs().sqrt(); //A
        expr[3] = (expr[0].sqr().add((expr[1].sub(list.lit("1"))).sqr()).
                sub(list.lit("2").mul(expr[0].mul(expr[1].sub(list.lit("1"))).mul(list.getInp(2).sin())))).abs().sqrt(); //B
        expr[4] = (list.lit("1").add((expr[1].sqr().sub(list.lit("1"))).mul(list.getInp(2).cos().sqr()))).sqrt(); //C
        expr[5] = ((list.lit("200").sub(list.getInp(1))).div(list.lit("200").add(list.getInp(1)))).sqrt(); //D
        expr[6] = (list.getInp(0).mul(list.getInp(2).cos())).div(list.lit("100").sub(list.getInp(0).mul(list.getInp(2).sin()))); //E
        expr[7] = (expr[1].sqr().sub(list.lit("1"))).sqrt();//F
        expr[8] = (list.lit("1").div(list.pi())). //Fh
                mul(
                (list.lit("1").div(expr[5])).atan().
                        add((list.getInp(2).sin().div(expr[4])).
                                mul(
                                        ((expr[0].mul(expr[1]).sub(expr[7].sqr().mul(list.getInp(2).sin()))).div(expr[7].mul(expr[4]))).atan().
                                                add((expr[7].mul(list.getInp(2).sin()).div(expr[4])).atan()
                                                )
                                )
                        ).
                        sub(
                                (expr[0].sqr().add((expr[1].add(list.lit("1"))).sqr()).sub(list.lit("2").
                                        mul(expr[1].add(list.lit("1").add(expr[0].mul(expr[1]).mul(list.getInp(2).sin())))))
                                ).div(expr[2].mul(expr[3])).mul((expr[2].mul(expr[5]).div(expr[3])).atan())
                        )
        );
        expr[9] = (list.lit("1").div(list.pi())). //Fv
                mul(expr[6].neg().mul(expr[5].atan()).
                add(expr[6].
                        mul(expr[0].sqr().add((expr[1].add(list.lit("1"))).sqr()).sub(list.lit("2").mul(expr[1]).
                                mul(list.lit("1").add(expr[0].mul(list.getInp(2).sin()))))
                        ).div(expr[2].mul(expr[3])).mul((expr[2].mul(expr[5]).div(expr[3])).atan()).
                        add(list.getInp(2).cos().div(expr[4]).
                                mul(((expr[0].mul(expr[1]).sub(expr[7].sqr().mul(list.getInp(2).sin()))).div(expr[7].mul(expr[4]))).atan().
                                        add((expr[7].mul(list.getInp(2).sin()).div(expr[4])).atan()
                                        )
                                )
                        )
                )
        );
        objectiveFunction = (expr[8].sqr().add(expr[9].sqr()));
    }

    public Expression getObjective() {
        return this.objectiveFunction;
    }
}