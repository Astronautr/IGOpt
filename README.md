# IGOpt
IGOpt (Interval Global Optimization) is an application for finding global maximum or minimum of objective function
using interval algorithms. There were used JInterval library for all interval operations
and Gradient class for automatic differentiation.

#Usage
To run this program you need to set dependent and objective functions placed (for this moment)
in Functions.java and some input data placed in txt file (you can enter this data from console):
* Name of independent variables;
* Tolerance (error);
* The range of independent variables;
* Number of dependent expressions.  


After that you can run IGOpt.sh with input file as script argument, computation modes key and one of keys which set combination
of algorithm modifications.


Computation modes provided by JInterval (full name used in command line input):
* -pl (plain);
* -ac (accur);
* -t64 (tightest64);
* -ex (exact);
* -dn (dnearest).


Algorithm modifications:
* -s for simplest interval algorithm of global optimization;
* -b for algorithm with some effective additions for one or more times differentiable functions;
* -a for algorithm with more modifications like Krawczyk operator, second order operations, etc.; it is applied
for two or more times differentiable functions.

#Links
Introduction to JInterval in "Reliable Computing Journal": http://goo.gl/IcQSDP.  
Source of JInterval library: https://java.net/projects/jinterval.  
Source of Gradient class: https://github.com/Astronautr/gradient.
