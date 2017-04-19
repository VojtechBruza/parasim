package org.sybila.parasim.computation.simulation.libsbmlsim;

import org.sybila.parasim.computation.simulation.api.PrecisionConfiguration;
import org.sybila.parasim.computation.simulation.cpu.SimulationEngine;
import org.sybila.parasim.model.ode.OdeSystem;
import org.sybila.parasim.model.trajectory.Point;
import org.sybila.parasim.model.trajectory.Trajectory;
import jp.ac.keio.bio.fun.libsbmlsim.*;

/**
 * @author Vojta_2
 */
public class LibSBMLSimSimulationEngine implements SimulationEngine {

    @Override
    public void close() {

    }

    static {
        System.loadLibrary("sbmlsimj");
    }

    @Override
    public Trajectory simulate(Point point, OdeSystem odeSystem, double timeLimit, PrecisionConfiguration configuration) {
        libsbmlsim.simulateSBMLFromString(odeSystem.getOriginalModel().toString(),);
        return null;
    }
}
