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
package org.sybila.parasim.computation.simulation.cpu;

import java.util.concurrent.ConcurrentHashMap;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface SimulationEngineFactory<E extends SimulationEngine>{

    boolean isAvailable();

    /**
     * To get simulation engine for the thread which this method is called from.
     * New SimulationeEngine instance is constructed, if none exists according to THREAD_SIMULATION_ENGINE_MAP.
     * @param stepLimit max number of simulation steps of the simulation engine
     * @return SimulationEngine
     */
    E simulationEngine(long stepLimit);

    /**
     * Hash map that allows closing of all SimulationEngine instances after the simulation of the whole space finishes
     */
    ConcurrentHashMap<Thread,SimulationEngine> THREAD_SIMULATION_ENGINE_MAP = new ConcurrentHashMap<>();

}
