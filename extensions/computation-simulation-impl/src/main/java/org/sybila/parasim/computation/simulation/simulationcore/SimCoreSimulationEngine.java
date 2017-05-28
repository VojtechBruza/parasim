/**
 * Copyright 2011-2016, Sybila, Systems Biology Laboratory and individual
 * contributors by the @authors tag.
 *
 * This file is part of Parasim.
 *
 * Parasim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.sybila.parasim.computation.simulation.simulationcore;

import org.apache.commons.math.ode.DerivativeException;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.validator.ModelOverdeterminedException;
import org.simulator.math.odes.AdaptiveStepsizeIntegrator;
import org.simulator.math.odes.MultiTable;
import org.simulator.math.odes.RosenbrockSolver;
import org.simulator.sbml.SBMLinterpreter;
import org.sybila.parasim.computation.simulation.api.PrecisionConfiguration;
import org.sybila.parasim.computation.simulation.cpu.SimulationEngine;
import org.sybila.parasim.model.math.Parameter;
import org.sybila.parasim.model.math.Variable;
import org.sybila.parasim.model.ode.OdeSystem;
import org.sybila.parasim.model.trajectory.*;

/**
 * @author <a href="mailto:433392@fi.muni.cz">Vojtech Bruza</a>
 */
public class SimCoreSimulationEngine implements SimulationEngine {

    @Override
    public void close() {
    }

    /**
     * array for caching time of simulations for Parasim
     */
    private float[] l_times;
    /**
     * array for caching time of simulations for Simulatin Core Library
     */
    private double [] l_timesDouble;
    /**
     * start of last simulation
     */
    private float l_startTime;
    /**
     * end of last simulation
     */
    private double l_endTime;
    /**
     * time step of previous simulation
     */
    private float timeStep;

    @Override
    public Trajectory simulate(Point point, OdeSystem odeSystem, double timeLimit, PrecisionConfiguration precision) {
        Model model = odeSystem.getOriginalModel();
        if(model == null){
            throw new IllegalStateException("Can't simulate the trajectory due to missing SBML model.");
        }
//        //SETTING VARIABLES
        for(Variable variable : odeSystem.getVariables().values()){
            if (!variable.isSubstituted()) {
                //set species (variables) values in model
                model.getSpecies(variable.getName()).setValue(point.getValue(variable.getIndex()));
            }
        }
        //SETTING PARAMETERS
        for(Parameter parameter : odeSystem.getAvailableParameters().values()){
            if (!parameter.isSubstituted()) { //substituted are those that are set at the beginning (not perturbing over them)
                model.getParameter(parameter.getName()).setValue(point.getValue(parameter.getIndex()));
            }
        }
        //SET SIMULATION TIME (NUM OF ITERATIONS) - START, END, NUM OF ITERATION, TIME STEP
        long numOfIterations = Math.round(Math.ceil((1.05 * timeLimit - point.getTime()) / precision.getTimeStep()));
        if (numOfIterations > Integer.MAX_VALUE) { // step limit is max integer value, because it needs to be retyped from long later on
            throw new IllegalStateException("Can't simulate the trajectory because the number of iterations <" + numOfIterations + "> has exceeded its maximum <" + Integer.MAX_VALUE + ">.");
        }
        //TIME - needs to be float due to communication with other modules of Parasim and double due to SimCore
        if(l_startTime != point.getTime() || l_endTime != timeLimit || timeStep != precision.getTimeStep()){
            l_startTime = point.getTime();
            l_endTime = timeLimit;
            timeStep = precision.getTimeStep();
            l_times = new float[(int) numOfIterations];
            float timeFloat = point.getTime();
            l_timesDouble = new double[(int) numOfIterations];
            for (int j = 0; j < l_times.length; j++) {
                timeFloat += precision.getTimeStep();
                l_times[j] = timeFloat;
                l_timesDouble[j] = timeFloat;
            }
        }
        //SET SIMULATION METHOD
        SBMLinterpreter interpreter = null;
        try {
            interpreter = new SBMLinterpreter(model);
        } catch (ModelOverdeterminedException e) {
            e.printStackTrace();
        }
        AdaptiveStepsizeIntegrator solver = new RosenbrockSolver();
        MultiTable solution = null;
        //TOLERANCE SETTING
        //relative error
        solver.setRelTol(precision.getMaxRelativeError());
        //absolute error for each variable not supported => using minimum of all variables as the tolerance for each of them
        float maxAbsoluteError = precision.getMaxAbsoluteError(0);
        for (int i = 0; i < odeSystem.getVariables().size(); i++){
            if(precision.getMaxAbsoluteError(i) < maxAbsoluteError){
                maxAbsoluteError = precision.getMaxAbsoluteError(i);
            }
        }
        solver.setAbsTol(maxAbsoluteError);
        //SIMULATION
        try {
            solution = solver.solve(interpreter, interpreter.getInitialValues(), l_timesDouble);
        } catch (DerivativeException | NullPointerException e) {
            e.printStackTrace();
        }
        int numberOfSteps = 0;
        try{
            numberOfSteps = solution.getRowCount();
        } catch (NullPointerException e) {
            e.printStackTrace(); //if solver failed and solution is null
        }
        float[] simulatedData = new float[numberOfSteps*odeSystem.getVariables().size()];
        for (int currentVariable = 0; currentVariable < odeSystem.getVariables().size(); currentVariable++) {
            for (int i = 0; i < solution.getRowCount(); i++) {
                simulatedData[currentVariable + i * odeSystem.getVariables().size()] = (float) solution.getColumn(odeSystem.getVariable(currentVariable).getName()).getValue(i);
            }
        }
        //OUTPUT TRAJECTORY
        if (odeSystem.getAvailableParameters().isEmpty()) {
            return new ArrayTrajectory(simulatedData, l_times, point.getDimension());
        } else {
            return new ArrayTrajectory(point, simulatedData, l_times, odeSystem.getVariables().size());
        }
    }
}
