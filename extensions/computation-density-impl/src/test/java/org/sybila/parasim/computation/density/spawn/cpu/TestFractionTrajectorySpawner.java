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
package org.sybila.parasim.computation.density.spawn.cpu;

import org.sybila.parasim.computation.density.spawn.api.TrajectorySpawner;
import org.testng.annotations.Test;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class TestFractionTrajectorySpawner extends AbstractTrajectorySpawnerTest {

    @Test
    public void testNumberOfTrajectoriesAfterSpawn() {
        for (int dim = 1; dim <= 5; dim++) {
            super.testNumberOfTrajectoriesAfterInitialSpawn(dim);
        }
    }

    @Test
    public void testNumberOfTrajectoriesInNeighborhoodAfterInitialSpawn() {
        for (int dim = 1; dim <= 5; dim++) {
            super.testNumberOfTrajectoriesInNeighborhoodAfterInitialSpawn(dim);
        }
    }
    @Override
    protected TrajectorySpawner createTrajectorySpawner() {
        return new FractionTrajectorySpawner(new SpawnedTrajectoriesCacheImpl(), new SpawnedTrajectoriesCacheImpl());
    }

}
