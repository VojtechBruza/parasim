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
package org.sybila.parasim.computation.lifecycle.impl.shared;

import java.util.concurrent.ExecutorService;
import org.sybila.parasim.computation.lifecycle.api.Computation;
import org.sybila.parasim.computation.lifecycle.api.Emitter;
import org.sybila.parasim.computation.lifecycle.api.Future;
import org.sybila.parasim.computation.lifecycle.api.MutableStatus;
import org.sybila.parasim.computation.lifecycle.api.Offerer;
import org.sybila.parasim.computation.lifecycle.api.SharedMemoryExecutor;
import org.sybila.parasim.computation.lifecycle.api.Status;
import org.sybila.parasim.computation.lifecycle.api.annotations.ComputationScope;
import org.sybila.parasim.computation.lifecycle.impl.common.CallableFactory;
import org.sybila.parasim.computation.lifecycle.impl.common.AbstractExecutor;
import org.sybila.parasim.core.annotation.Default;
import org.sybila.parasim.core.api.Binder;
import org.sybila.parasim.core.api.Context;
import org.sybila.parasim.core.api.enrichment.Enrichment;
import org.sybila.parasim.model.Mergeable;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class SharedMemoryExecutorImpl extends AbstractExecutor implements SharedMemoryExecutor {

    public SharedMemoryExecutorImpl(Enrichment enrichment, Context context) {
        super(enrichment, context);
    }

    @Override
    public <M extends Mergeable<M>> Future<M> submit(Computation<M> computation) {
        // init context
        Context context = getContext().context(ComputationScope.class);
        // prepare services
        MutableStatus status = new SimpleStatus();
        Offerer offerer = new SimpleOfferer(status);
        ExecutorService executorService = context.resolve(ExecutorService.class, Default.class);
        CallableFactory callableFactory = new CallableFactory(context, status);
        SharedMemoryFuture<M> future = new SharedMemoryFuture<>(context, status);
        SharedMemoryMucker mucker = new SharedMemoryMucker(status, executorService, offerer, Long.MAX_VALUE, callableFactory);
        // bind services
        Binder binder = context.resolve(Binder.class, Default.class);
        binder.bind(Emitter.class, Default.class, offerer);
        binder.bind(Status.class, Default.class, status);
        // register progress listeners
        status.addProgressListerner(future);
        status.addProgressListerner(mucker);
        // emit computation
        offerer.emit(computation);
        // return future
        return future;
        // TODO: destroy computation context
    }



}
