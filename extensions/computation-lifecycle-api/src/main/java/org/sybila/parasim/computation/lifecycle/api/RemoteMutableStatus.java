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
package org.sybila.parasim.computation.lifecycle.api;

import java.rmi.Remote;
import java.rmi.RemoteException;
import java.util.UUID;
import org.sybila.parasim.model.Mergeable;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface RemoteMutableStatus extends Remote {

    /**
     * Number of done computation instances.
     */
    long getDone() throws RemoteException;

    /**
     * Number of currently computing computation instances.
     */
    long getComputing() throws RemoteException;

    /**
     * Number of remaining (not done) computation instances.
     */
    long getRemaining() throws RemoteException;

    /**
     * Checks whether the computation is finished.
     */
    boolean isFinished() throws RemoteException;

    void compute(UUID node, java.util.concurrent.Future event) throws RemoteException;

    void done(UUID node, Mergeable event) throws RemoteException;

    void emit(UUID node, Computation computation) throws RemoteException;

    void reschedule(UUID node, Computation computation) throws RemoteException;
}
