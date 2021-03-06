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
package org.sybila.parasim.computation.lifecycle.api;

/**
 * This service provides available computation instances to the executor.
 * It's mainly internal service and shouldn't be used outside of computation
 * life cycle.
 *
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface Offerer extends Emitter, ProgressListener {

    /**
     * Retrieves and removes the head of this queue, or null if this offerer
     * is empty. Computation instance is returned for the execution purpose.
     */
    Computation poll();

    /**
     * Retrieves and removes the head of this queue, or null if this offerer
     * is empty. Computation instance is returned for the balancing purpose.
     */
    Computation balance();

    /**
     * Returns the number of elements in this offerer.
     */
    int size();

}
