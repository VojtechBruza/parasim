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

import java.util.UUID;
import org.sybila.parasim.model.Mergeable;

/**
 * The mutable status is used to control the computation execution. See
 * documentation of {@link ProgressListener} about the contract between
 * {@link ProgressListener} and {@link MutableStatus}.
 *
 * @see {@link ProgressListener}
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface MutableStatus extends Status {

    /**
     * Add a progress listener.
     */
    void addProgressListerner(ProgressListener progressListener);

    /**
     * Increment the number of computed computation instances and invoke
     * {@link ProgressListener#computing(java.util.concurrent.Future) }.
     *
     * @see #getComputing()
     */
    void compute(UUID node, java.util.concurrent.Future event);

    /**
     * Invoke {@link ProgressListener#done(org.sybila.parasim.model.Mergeable) }.
     * Increment the number of done computation instances and decrements number of
     * computed computation instances and number of remaining computation instances.
     * If the the computation is finished, invoke {@link ProgressListener#finished(org.sybila.parasim.model.Mergeable) }
     * and mark the computation as finished.
     *
     * @see #getDone()
     * @see #isFinished()
     */
    void done(UUID node, Mergeable event);

    /**
     * Increment the number of remaining computation instances and invoke
     * {@link ProgressListener#emitted(org.sybila.parasim.computation.lifecycle.api.Computation) }.
     */
    void emit(UUID node, Computation event);

    /**
     * Invoke {@link ProgressListener#rescheduled(org.sybila.parasim.computation.lifecycle.api.Computation) }.
     * This method doesn't increment any counter.
     */
    void balance(UUID node, Computation event);

}
