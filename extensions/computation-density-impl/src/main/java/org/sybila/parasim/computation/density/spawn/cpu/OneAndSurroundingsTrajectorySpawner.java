/**
 * Copyright 2011 - 2012, Sybila, Systems Biology Laboratory and individual
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
package org.sybila.parasim.computation.density.spawn.cpu;

import java.lang.ref.ReferenceQueue;
import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.sybila.parasim.computation.density.api.Configuration;
import org.sybila.parasim.model.trajectory.Distance;
import org.sybila.parasim.model.trajectory.ListDataBlock;
import org.sybila.parasim.model.trajectory.Point;
import org.sybila.parasim.model.trajectory.PointTrajectory;
import org.sybila.parasim.model.trajectory.Trajectory;
import org.sybila.parasim.model.trajectory.TrajectoryWithNeighborhood;
import org.sybila.parasim.model.trajectory.TrajectoryWithNeighborhoodWrapper;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class OneAndSurroundingsTrajectorySpawner extends AbstractTrajectorySpawner {

    private final Map<Point, WeakReference<Trajectory>> alreadySpawnedCollisionTrajectories = new HashMap<>();

    private final Map<Point, WeakReference<Trajectory>> alreadySpawnedPrimaryTrajectories = new HashMap<>();

    @Override
    protected SpawnedResult spawnTrajectories(Configuration configuration, Trajectory trajectory, Trajectory neighbor, Distance distance) {
        // mark secondary trajectory as a primary one
        Trajectory newPrimary = neighbor;
        // find dimension which the first points of the given trajectories differ in
        int diffDimension = -1;
        for (int dim = 0; dim < trajectory.getDimension(); dim++) {
            if (trajectory.getFirstPoint().getValue(dim) != neighbor.getFirstPoint().getValue(dim)) {
                if (diffDimension != -1) {
                    throw new IllegalStateException("The first points of the given trajectories differ in more than one dimension.");
                }
                diffDimension = dim;
            }
        }
        // compute half of their distance
        float radius = Math.abs(trajectory.getFirstPoint().getValue(diffDimension) - neighbor.getFirstPoint().getValue(diffDimension)) / 2;
        // memory for spawned trajectories
        List<Trajectory> neighborTrajectories = new ArrayList<>();
        List<Trajectory> spawnedSecondaryTrajectories = new ArrayList<>();
        // create neighbor trajectories which can have collision
        for (int dim = 0; dim < trajectory.getDimension(); dim++) {
            for (int sign = -1; sign <= 1; sign += 2) {
                float[] newPointData = newPrimary.getFirstPoint().toArrayCopy();
                newPointData[dim] += sign * radius;
                Trajectory newTrajectory = new PointTrajectory(trajectory.getFirstPoint().getTime(), newPointData);
                if (!configuration.getInitialSpace().isIn(newTrajectory.getFirstPoint())) {
                    continue;
                }
                synchronized(alreadySpawnedCollisionTrajectories) {
                    if (alreadySpawnedCollisionTrajectories.containsKey(newTrajectory.getFirstPoint())) {
                        WeakReference<Trajectory> reference = alreadySpawnedCollisionTrajectories.get(newTrajectory.getFirstPoint());
                        if (reference.get() != null) {
                            newTrajectory = reference.get() instanceof TrajectoryWithNeighborhood ? ((TrajectoryWithNeighborhood) reference.get()).withoutNeighbors() : reference.get();
                        } else {
                            alreadySpawnedCollisionTrajectories.put(newTrajectory.getFirstPoint(), new WeakReferenceWithEquals<>(newTrajectory));
                        }
                    } else {
                        alreadySpawnedCollisionTrajectories.put(newTrajectory.getFirstPoint(), new WeakReferenceWithEquals<>(newTrajectory));
                        spawnedSecondaryTrajectories.add(newTrajectory);
                    }
                }
                neighborTrajectories.add(newTrajectory);
            }
        }
        // reorganize

        Collection<TrajectoryWithNeighborhood> spawnedCol = new ArrayList<>();
        synchronized(alreadySpawnedPrimaryTrajectories) {
            if (!alreadySpawnedPrimaryTrajectories.containsKey(newPrimary.getFirstPoint())) {
                spawnedCol.add(TrajectoryWithNeighborhoodWrapper.createAndUpdateReference(newPrimary, new ListDataBlock<>(neighborTrajectories)));
                alreadySpawnedPrimaryTrajectories.put(newPrimary.getFirstPoint(), new WeakReferenceWithEquals<>(newPrimary));
            }
        }
        return new SpawnedResult(spawnedCol, spawnedCol.isEmpty() ? Collections.EMPTY_LIST : spawnedSecondaryTrajectories);
    }

    private static class WeakReferenceWithEquals<T> extends WeakReference<T> {

        public WeakReferenceWithEquals(T referent) {
            super(referent);
        }

        public WeakReferenceWithEquals(T referent, ReferenceQueue<? super T> q) {
            super(referent, q);
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof WeakReferenceWithEquals)) {
                return false;
            }
            WeakReferenceWithEquals other = (WeakReferenceWithEquals) obj;
            if (this.get() == null && other.get() == null) {
                return true;
            }
            if (this.get() == null || other.get() == null) {
                return false;
            }
            return this.get().equals(other.get());
        }

        @Override
        public int hashCode() {
            if (this.get() == null) {
                return 113 * WeakReferenceWithEquals.class.hashCode();
            }
            return 23 * this.get().hashCode();
        }

    }
}
