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
package org.sybila.parasim.computation.lifecycle.impl.common;

import java.util.ArrayList;
import java.util.List;
import java.util.UUID;
import net.jcip.annotations.GuardedBy;
import org.apache.commons.lang3.Validate;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.sybila.parasim.computation.lifecycle.api.Computation;
import org.sybila.parasim.computation.lifecycle.api.MutableStatus;
import org.sybila.parasim.computation.lifecycle.api.Offerer;
import org.sybila.parasim.computation.lifecycle.api.ProgressAdapter;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class SimpleOfferer extends ProgressAdapter implements Offerer {

    private final static Logger LOGGER = LoggerFactory.getLogger(SimpleOfferer.class);

    @GuardedBy("this")
    private final List<Computation> computations = new ArrayList<>();
    private final MutableStatus status;
    private final UUID node;

    public SimpleOfferer(UUID node, MutableStatus status) {
        Validate.notNull(node);
        Validate.notNull(status);
        this.status = status;
        this.node = node;
    }

    @Override
    public synchronized Computation poll() {
        if (computations.isEmpty()) {
            return null;
        }
        return computations.remove(0);
    }

    @Override
    public synchronized Computation reschedule() {
        if (computations.isEmpty()) {
            return null;
        } else {
            return computations.remove(0);
        }
    }

    @Override
    public void reschedule(Computation computation) {
        synchronized(this) {
            computations.add(computation);
            LOGGER.debug("rescheduling " + computation + " => " + computations);
        }
        status.reschedule(node, computation);
    }

    @Override
    public synchronized int size() {
        return computations.size();
    }

    @Override
    public void emit(Computation computation) {
        status.emit(node, computation);
    }

    @Override
    public void emitted(UUID node, Computation event) {
        synchronized(this) {
            computations.add(event);
        }
    }

    protected final UUID getNode() {
        return node;
    }


}