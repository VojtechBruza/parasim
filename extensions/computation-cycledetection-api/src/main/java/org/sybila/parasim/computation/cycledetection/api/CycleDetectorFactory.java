/**
 * Copyright 2011 - 2013, Sybila, Systems Biology Laboratory and individual
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
package org.sybila.parasim.computation.cycledetection.api;

import org.sybila.parasim.model.trajectory.DataBlock;
import org.sybila.parasim.model.trajectory.Trajectory;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface CycleDetectorFactory {

    /**
     * Performs cycle detection on the trajectory
     *
     * @param trajectory Trajectory on which to carry out cycle detection.
     * @return cycle detector
     */
    CycleDetector detect(Trajectory trajectory);

    /**
     * Performs cycle detection on the trajectory data block
     *
     * @param trajectories trajectories on which to carry out cycle detection.
     * @return cycle detected data block
     */
    <T extends Trajectory> CycleDetectedDataBlock<T> detect(DataBlock<T> trajectories);

}
