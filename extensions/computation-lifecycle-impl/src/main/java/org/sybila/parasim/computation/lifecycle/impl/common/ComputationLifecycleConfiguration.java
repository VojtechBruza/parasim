/**
 * Copyright 2011 - 2013, Sybila, Systems Biology Laboratory and individual
 * contributors by the
 *
 * @authors tag.
 *
 * This file is part of Parasim.
 *
 * Parasim is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
package org.sybila.parasim.computation.lifecycle.impl.common;

import java.util.concurrent.TimeUnit;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class ComputationLifecycleConfiguration {
    private int queueSize = Runtime.getRuntime().availableProcessors() * 1000;
    private int corePoolSize = Runtime.getRuntime().availableProcessors();
    private int maxPoolSize = corePoolSize * 5;
    private long keepAliveTimeAmount = 5;
    private TimeUnit keepAliveTimeUnit = TimeUnit.SECONDS;
    private long shutdownTimeAmout = 2;
    private TimeUnit shutdownTimeUnit = TimeUnit.SECONDS;

    public int getQueueSize() {
        return queueSize;
    }

    public int getCorePoolSize() {
        return corePoolSize;
    }

    public int getMaxPoolSize() {
        return maxPoolSize;
    }

    public long getKeepAliveTimeAmount() {
        return keepAliveTimeAmount;
    }

    public TimeUnit getKeepAliveTimeUnit() {
        return keepAliveTimeUnit;
    }

    public long getShutdownTimeAmout() {
        return shutdownTimeAmout;
    }

    public TimeUnit getShutdownTimeUnit() {
        return shutdownTimeUnit;
    }
}