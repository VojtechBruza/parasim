package org.sybila.parasim.computation.simulation;

import org.sybila.parasim.model.trajectory.Trajectory;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface Simulator<Conf extends Configuration, Out extends SimulatedDataBlock> {

    Out simulate(Conf configuration, org.sybila.parasim.computation.DataBlock<Trajectory> trajectories);
   
}