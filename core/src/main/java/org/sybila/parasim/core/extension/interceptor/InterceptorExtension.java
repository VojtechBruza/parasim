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
package org.sybila.parasim.core.extension.interceptor;

import org.sybila.parasim.core.Instance;
import org.sybila.parasim.core.annotations.Inject;
import org.sybila.parasim.core.annotations.Observes;
import org.sybila.parasim.core.event.ManagerProcessing;
import org.sybila.parasim.core.extension.interceptor.api.InterceptorRegistry;
import org.sybila.parasim.core.extension.interceptor.api.Managed;
import org.sybila.parasim.core.extension.interceptor.api.Standalone;
import org.sybila.parasim.core.extension.interceptor.impl.ManagedInterceptorRegistry;
import org.sybila.parasim.core.extension.interceptor.impl.StandaloneInterceptorRegistry;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class InterceptorExtension {

    @Inject
    private Instance<InterceptorRegistry> defaultRegistry;
    @Inject
    @Managed
    private Instance<InterceptorRegistry> managedRegistry;
    @Inject
    @Standalone
    private Instance<InterceptorRegistry> standaloneRegistry;

    public void registerDefaultRegistry(@Observes ManagerProcessing event) {
        defaultRegistry.set(new StandaloneInterceptorRegistry());
    }

    public void registerManagedRegistry(@Observes ManagerProcessing event) {
        managedRegistry.set(new ManagedInterceptorRegistry(event.getManager()));
    }

    public void registerStandaloneRegistry(@Observes ManagerProcessing event) {
        standaloneRegistry.set(new StandaloneInterceptorRegistry());
    }

}