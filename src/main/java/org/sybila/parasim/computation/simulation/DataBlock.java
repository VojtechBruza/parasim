package org.sybila.parasim.computation.simulation;

import org.sybila.parasim.model.trajectory.Trajectory;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface DataBlock<T extends Trajectory> extends org.sybila.parasim.computation.DataBlock<T> {
   
    /**
     * Returns trajectory status
     * 
     * @param index number from interval [0, number of trajectories)
     * @return trajectory status
     */
    Status getStatus(int index);
    
}
