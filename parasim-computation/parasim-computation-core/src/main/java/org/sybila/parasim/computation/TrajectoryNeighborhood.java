package org.sybila.parasim.computation;

import org.sybila.parasim.computation.DataBlock;
import org.sybila.parasim.model.trajectory.Trajectory;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface TrajectoryNeighborhood<T extends Trajectory> {
    
    /**
     * Returns a trajectory neighborhood
     * 
     * @param trajectory
     * @return trajectories in neighborhood
     */
    DataBlock<T> getNeighbors(Trajectory trajectory);
    
    /**
     * Sets a trajectory neighborhood
     * 
     * @param trajectory
     * @param neighborhood 
     */
    void setNeighbors(Trajectory trajectory, DataBlock<T> neighborhood);
    
}