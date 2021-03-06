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
package org.sybila.parasim.computation.lifecycle.impl.common;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.Validate;
import org.sybila.parasim.computation.lifecycle.api.Computation;
import org.sybila.parasim.computation.lifecycle.api.ComputationContainer;
import org.sybila.parasim.computation.lifecycle.api.Executor;
import org.sybila.parasim.computation.lifecycle.api.Future;
import org.sybila.parasim.computation.lifecycle.api.annotations.RunWith;
import org.sybila.parasim.core.annotation.Default;
import org.sybila.parasim.core.api.Resolver;
import org.sybila.parasim.model.Mergeable;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class DefaultComputationContainer implements ComputationContainer {

    private Resolver resolver;

    public DefaultComputationContainer(Resolver resolver) {
        Validate.notNull(resolver);
        this.resolver = resolver;
    }

    @Override
    public <Result extends Mergeable<Result>> Future<Result> compute(Computation<Result> computation) {
        Validate.notNull(computation);
        RunWith runWith = computation.getClass().getAnnotation(RunWith.class);
        if (runWith == null) {
            try {
                runWith = computation.getClass().getMethod("call").getAnnotation(RunWith.class);
            } catch (NoSuchMethodException ex) {
                Logger.getLogger(DefaultComputationContainer.class.getName()).log(Level.SEVERE, null, ex);
            } catch (SecurityException ex) {
                Logger.getLogger(DefaultComputationContainer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        Class<? extends Executor> executorClass = runWith == null ? Executor.class : runWith.executor();
        Executor executor = resolver.resolve(executorClass, Default.class);
        return executor.submit(computation);
    }
}