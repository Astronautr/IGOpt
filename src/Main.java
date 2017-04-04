import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.ExtendedRational;

import java.io.*;


public class Main {
    public static void run(String mode, SetIntervalContext ic, String[] keys, Boolean inpFromFile, BufferedReader bf) throws IOException {
        System.out.println("Running with " + mode + " mode");
        if (inpFromFile) {
            System.out.println("input from file");
            String[] inps = bf.readLine().split(" "); //reading input variables
            float tolerance = Float.parseFloat(bf.readLine()); //reading the tolerance
            String strIntv = bf.readLine(); //reading input intervals
            strIntv = strIntv.replaceAll("\\s+", "");
            String[] strIntval = strIntv.split("(?<=])");
            SetInterval[] box = new SetInterval[strIntval.length];
            for (int i = 0; i < box.length; i++) {
                box[i] = ic.textToInterval(strIntval[i]);
            }
            int exprNum = Integer.parseInt(bf.readLine()); //reading number of dependent expressions

            ExtendedRational goptInf, goptSup;
            Algorithm gopt = new Algorithm(box, ic, ExtendedRational.valueOf(tolerance), inps, exprNum, false);
            long startTime = System.currentTimeMillis();
            goptInf = gopt.start(keys);
            gopt = new Algorithm(box, ic, ExtendedRational.valueOf(tolerance), inps, exprNum, true);
            goptSup = gopt.start(keys);
            System.out.println("Computation time = " + (System.currentTimeMillis() - startTime) / 1E+3 + " sec");
            System.out.println("Assessment = [" + goptInf.doubleValue() + ", " + (-goptSup.doubleValue()) + "]." );
        } else {
            System.out.println("input from console");
            System.out.println("Enter the name of independent variables");
            String[] inps = bf.readLine().split(" ");
            System.out.println("Enter the tolerance");
            float tolerance = Float.parseFloat(bf.readLine());
            System.out.println("Enter input intervals");
            String strIntv = bf.readLine();
            strIntv = strIntv.replaceAll("\\s+", "");
            String[] strIntval = strIntv.split("(?<=])");
            SetInterval[] box = new SetInterval[strIntval.length];
            for (int i = 0; i < box.length; i++) {
                box[i] = ic.textToInterval(strIntval[i]);
            }
            System.out.println("Enter the number of dependent expressions");
            int exprNum = Integer.parseInt(bf.readLine());

            ExtendedRational goptInf, goptSup;
            Algorithm gopt = new Algorithm(box, ic, ExtendedRational.valueOf(tolerance), inps, exprNum, false);
            long startTime = System.currentTimeMillis();
            goptInf = gopt.start(keys);
            gopt = new Algorithm(box, ic, ExtendedRational.valueOf(tolerance), inps, exprNum, true);
            goptSup = gopt.start(keys);
            System.out.println("Computation time = " + (System.currentTimeMillis() - startTime) / 1E+3 + " sec");
            System.out.println("Assessment = [" + goptInf.doubleValue() + ", " + (-goptSup.doubleValue()) + "]." );
        }
    }

    public static void main(String[] args) throws IOException {
        boolean argsEmpty = false;
        try {
            args[0].isEmpty();
        } catch (ArrayIndexOutOfBoundsException indOut) {
            argsEmpty = true;
            System.out.println("Please, enter input data from console.");
            System.out.println();
            BufferedReader consoleReader = new BufferedReader(new InputStreamReader(System.in));
            System.out.println("Enter computing mode");
            String cMode = consoleReader.readLine();
            String[] keys;
            System.out.println("Enter keys for computing mode");
            String tmp = consoleReader.readLine();
            keys = tmp.split("(?=-)");
            switch (cMode) {
                case "plain":
                    run("plain", SetIntervalContexts.getPlain(), keys, false, consoleReader);
                    break;
                case "accur":
                    run("accur", SetIntervalContexts.getAccur64(), keys, false, consoleReader);
                    break;
                case "tightest64":
                    run("tightest64", SetIntervalContexts.getTightest64(), keys, false, consoleReader);
                    break;
                case "exact":
                    run("exact", SetIntervalContexts.getExact(), keys, false, consoleReader);
                    break;
                case "dnearest":
                    run("double nearest", SetIntervalContexts.getDoubleNearest(), keys, false, consoleReader);
                    break;
            }
        }
        if ((!argsEmpty) && (args[0].contains(".txt"))) {
            argsEmpty = true;
            try {
                args[1].isEmpty();
                args[2].isEmpty();
            }
            catch (ArrayIndexOutOfBoundsException indOut) {
                System.err.print("Missing some input arguments.");
                return;
            }
            try {
                BufferedReader fileReader = new BufferedReader(new FileReader(args[0]));
                String cMode = args[1];
                StringBuilder sb = new StringBuilder();
                for (String s : args) {
                    sb.append(s);
                }
                sb.delete(0, args[0].length() + args[1].length());
                String[] keys = sb.toString().split("(?=-)");
                switch (cMode) {
                    case "-pl":
                        run("plain", SetIntervalContexts.getPlain(), keys, true, fileReader);
                        break;
                    case "-ac":
                        run("accur", SetIntervalContexts.getAccur64(), keys, true, fileReader);
                        break;
                    case "-t64":
                        run("tightest64", SetIntervalContexts.getTightest64(), keys, true, fileReader);
                        break;
                    case "-ex":
                        run("exact", SetIntervalContexts.getExact(), keys, true, fileReader);
                        break;
                    case "-dn":
                        run("double nearest", SetIntervalContexts.getDoubleNearest(), keys, true, fileReader);
                        break;
                }
            } catch (FileNotFoundException fnf) {
                System.out.println("No such file in directory. Please, enter input data from command line");
                BufferedReader consoleReader = new BufferedReader(new InputStreamReader(System.in));
                System.out.println("Enter computing mode");
                String cMode = consoleReader.readLine();
                String[] keys;
                if (args.length == 1) {
                    System.out.println("Enter keys for computing mode");
                    String tmp = consoleReader.readLine();
                    keys = tmp.split("(?=-)");
                } else {
                    StringBuilder sb = new StringBuilder();
                    for (String s : args) {
                        sb.append(s);
                    }
                    sb.delete(0, args[0].length());
                    keys = sb.toString().split("(?=-)");
                }
                switch (cMode) {
                    case "plain":
                        run("plain", SetIntervalContexts.getPlain(), keys, false, consoleReader);
                        break;
                    case "accur":
                        run("accur", SetIntervalContexts.getAccur64(), keys, false, consoleReader);
                        break;
                    case "tightest64":
                        run("tightest64", SetIntervalContexts.getTightest64(), keys, false, consoleReader);
                        break;
                    case "exact":
                        run("exact", SetIntervalContexts.getExact(), keys, false, consoleReader);
                        break;
                    case "dnearest":
                        run("double nearest", SetIntervalContexts.getDoubleNearest(), keys, false, consoleReader);
                        break;
                }
            }
        }
    }
}
