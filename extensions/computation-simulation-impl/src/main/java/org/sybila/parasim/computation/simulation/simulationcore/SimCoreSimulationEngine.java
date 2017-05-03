package org.sybila.parasim.computation.simulation.simulationcore;

import org.apache.commons.math.ode.DerivativeException;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.validator.ModelOverdeterminedException;
import org.simulator.math.odes.*;
import org.simulator.sbml.SBMLinterpreter;
import org.sybila.parasim.computation.simulation.api.PrecisionConfiguration;
import org.sybila.parasim.computation.simulation.cpu.SimulationEngine;
import org.sybila.parasim.model.math.Parameter;
import org.sybila.parasim.model.math.Variable;
import org.sybila.parasim.model.ode.OdeSystem;
import org.sybila.parasim.model.trajectory.*;


/**
 * @author Vojtech Bruza
 */
public class SimCoreSimulationEngine implements SimulationEngine {


    //TODO Vojta - add getStepLimit()? octave limitations or limitations of parasim?

    @Override
    public void close() {
        //TODO dont have to close this one - remove from SimulationEngine interface
    }

    @Override
    public Trajectory simulate(Point point, OdeSystem odeSystem, double timeLimit, PrecisionConfiguration precision) {
        long basicSettingStartTime = System.nanoTime();
        //todo create own Model instance instead of modifying the current one
        Model model = odeSystem.getOriginalModel(); //todo return copy of original model? not modify the previous one?
//        System.out.println("PARAMETER VALUES: " + odeSystem.getAvailableParameters().keySet() +" VARIABLES: " + odeSystem.getVariables().keySet());
//        // DONE Vojta - how to recognize if the index is same as in my Model (corresponding parameters and variables)
        //SETTING VARIABLES
        //DONE!! SET VARIABLES
        for(Variable variable : odeSystem.getVariables().values()) {
//            System.out.println(variable.getName() + " - VALUE: " + variable.evaluate(point));
//            System.out.println("INITIAL: " + point.getValue(odeSystem.getInitialVariableValue(variable).getExpression().getIndex()));
            if (!variable.isSubstituted()) { //TODO redundant condition? variable is never substituted
                //set species (variables) values in model
                double currentValue = point.getValue(variable.getIndex());
                model.getSpecies(variable.getName()).setValue(currentValue);
                //TODO throw exceptions if no variables with this name are found?
                //TODO Vojta - how does initial conditions setting work (or perturbating over variables not parameters) ->need to debug on linux? //I think it works like expected (get set values from the current point)
            }
        }
        //TODO - pritority - check difference when perturbating over subset of parameters or variables
//        System.out.println("Parameter count\nModel: " + model.getParameterCount() + "\nOde system: " + odeSystem.getAvailableParameters().size()); //2
//        System.out.println("Model variable count: " + model.getVariableCount()); //10
//        System.out.println("Model number of variables: " + model.getNumVariables()); 10
//        System.out.println("Ode system variable count: " + odeSystem.getVariables().size()); /3

//        //odeSystem.getAvailableParameters() returns only parameters, use odeSystem.getVariables() to get variables
//        //DONE findOut if parameters are also variables (if it is synonym)

//        //SETTING PARAMETERS
//        //DONE Vojta - set parameter value in model according to ode system parameter value
//        //DONE Vojta - create substituted odeSystem (model)

        //DONE!! SET PARAMETERS
        for(Parameter parameter : odeSystem.getAvailableParameters().values()){
//                System.out.println(parameter.getName() + " - VALUE: " + point.getValue(parameter.getIndex()));
            if (!parameter.isSubstituted()) { //substituted are those that are set at the beginning (not perturbating over them)
                //set parameters values in model
                model.getParameter(parameter.getName()).setValue(point.getValue(parameter.getIndex()));//what is the difference?? "parameter.evaluate(point)" - doesnt work if parameter is not substituted vs "point.getValue(parameter.getIndex())" - works
            }
        }


        //TODO!! SET SIMULATION TIME (NUM OF ITERATIONS) - START, END, NUM OF ITERATION, TIME STEP
        long numOfIterations = Math.round(Math.ceil((1.05 * timeLimit - point.getTime()) / precision.getTimeStep())); //magical constant 1.05 taked from Jan Papousek's implementation => guess it is for rendering reasons (to simulate a bit more data than a user excpects)
//        if (numOfIterations > getStepLimit()) { //TODO max num iterations limit is Integer.MAX_VALUE because I have to change the type to integer
//            throw new IllegalStateException("Can't simulate the trajectory because the number of iterations <" + numOfIterations + "> is higher than the given limit <" + getStepLimit() + ">.");
//        }

        long basicSettingTime = System.nanoTime() - basicSettingStartTime;

        long timeStartTime = System.nanoTime();
        //TIME - need to be float //TODO use only double? ask safranek
        double[] times = new double[(int) numOfIterations];
        float[] timesFloat = new float[(int) numOfIterations]; //TODO make time computation static to improve performance (check if other things can be static)
        float time = point.getTime(); //TODO use primary double or float?
        for (int j = 0; j < times.length; j++) {
            time += precision.getTimeStep();
            times[j] = (double) time;
            timesFloat[j] = time;
        }
        long timeTime = System.nanoTime() - timeStartTime;

        //DONE!! SIMULATION
        //SIMULATION
        long simulationSettingStartTime = System.nanoTime();
        SBMLinterpreter interpreter = null;
        try {
            interpreter = new SBMLinterpreter(model);
        } catch (ModelOverdeterminedException e) {
            e.printStackTrace();
            return null; //good?
        }
        if(interpreter == null){
            //TODO throw exception
        }

//        //TODO find out what is 'size' parameter -> use this constructor
//        AbstractDESSolver solver = new RosenbrockSolver(0, precision.getTimeStep());
//        in case of useage of the second constructor: solution = solver.solve(interpreter, interpreter.getInitialValues(), point.getTime(), timeLimit);

//        //DONE Vojta - where to get time array or start time and end time

        AdaptiveStepsizeIntegrator solver = new RosenbrockSolver(); //only +-10x slower than octave on lotkav model
//        AdaptiveStepsizeIntegrator solver = new AdamsBashforthSolver(); //Slower +-100x than octave on lotkav model
//        AdaptiveStepsizeIntegrator solver = new AdamsMoultonSolver(); //Slower more than +-100x than octave on lotkav model

        MultiTable solution = null;

        //DONE!! SET MAX RELATIVE
        solver.setRelTol(precision.getMaxRelativeError());
        //TODO SET MAX ABSOLUTE ERROR
//        float maxAbsoluteError = precision.getMaxAbsoluteError(0);
//        for (int i = 0; i < precision.getDimension(); i++){
//            if(precision.getMaxAbsoluteError(i) < maxAbsoluteError){
//                maxAbsoluteError = precision.getMaxAbsoluteError(i); //taking the smallest error
//            }
//        }
//        solver.setAbsTol(maxAbsoluteError); //TODO correct dimension setting???
        double[] intialValues = interpreter.getInitialValues();
        long simulationSettingTime = (System.nanoTime() - simulationSettingStartTime);
        long simulationStartTime = System.nanoTime();
        try {
            solution = solver.solve(interpreter, intialValues, times);
        } catch (DerivativeException e) {
            e.printStackTrace();
            return null; //good?
        }
        long simulationTime = System.nanoTime() - simulationStartTime;
//        if (solution == null) {
////            //TODO throw exception?
//        }
        //DONE!! DATA FROM MULTITABLE TO FLOAT ARRAY
//        //PARSING DATA TO TRAJECTORY
        //help printing
//        for (int i = 0; i < solution.getBlockCount(); i++) {
//            System.out.println("block" + i + " name: " + solution.getBlock(i).getName());
//            for (int c = 1; c < solution.getBlock(i).getColumnCount(); c++) {
//                System.out.println("\tcolumn" + c + " name: " + solution.getBlock(i).getColumn(c).getColumnName());
//            }
//        }
//        System.out.println("COLUMNS:");
//        for (Parameter p : odeSystem.getAvailableParameters().values()){
//            System.out.println("parameter: " + p.getName() + " column: " + solution.getColumn(p.getName()).getColumnName());
//        }
//        for (Variable v : odeSystem.getVariables().values()){
//            System.out.println("variable: " + v.getName() + " column: " + solution.getColumn(v.getName()).getColumnName());
//        }
//        System.out.println("Start time: " + point.getTime() + " End Time: " + timeLimit);
//        System.out.println("Number of iterations: " + numOfIterations);

        long parsingStartTime = System.nanoTime();
        int numberOfSteps = 0;
        try{
            numberOfSteps = solution.getRowCount();
        } catch (NullPointerException e) {
            e.printStackTrace(); //if solver failed and solution is null
            return null; //good? throw exception instead?
        }
        float[] simulatedData = new float[numberOfSteps*odeSystem.getVariables().size()];
        for (int currentVariable = 0; currentVariable < odeSystem.getVariables().size(); currentVariable++) { //simulating only variables, not parameters
//            System.out.println(odeSystem.getVariable(currentVariable).getName());
            for (int i = 0; i < solution.getRowCount(); i++) {
//                System.out.println("row: " + i + " time: " + times[i] + " value: " + solution.getColumn(odeSystem.getVariable(currentVariable).getName()).getValue(i));
                simulatedData[currentVariable + i * odeSystem.getVariables().size()] = (float) solution.getColumn(odeSystem.getVariable(currentVariable).getName()).getValue(i);
            }
        }

       long parsingTime = System.nanoTime() - parsingStartTime;

        //TODO - top priority - verification of sizes of arrays, dimensions and so
        //DONE!! OUTPUT TRAJECTORY
        //DONE Vojta - how to create new trajectory from multitable - correct data parsing
        System.out.println("TIME");
        System.out.printf("Setting time: %.9f ms\n", (simulationSettingTime + basicSettingTime) / 1000000000.0);
        System.out.printf("Simulation time: %.9f ms\n", simulationTime / 1000000000.0);
        System.out.printf("Time time: %.9f ms\n", timeTime / 1000000000.0);
        System.out.printf("Parsing time: %.9f ms\n", parsingTime / 1000000000.0);


        if (odeSystem.getAvailableParameters().isEmpty()) {
            return new ArrayTrajectory(simulatedData, timesFloat, point.getDimension()); //TODO howcome this line never runs?
        } else {
            return new ArrayTrajectory(point, simulatedData, timesFloat, odeSystem.getVariables().size());
        }
    }
}
