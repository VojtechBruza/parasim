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
package org.sybila.parasim.visualisation.plot.impl;

import org.sybila.parasim.model.trajectory.Point;
import org.sybila.parasim.visualisation.plot.api.MouseOnResultListener;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class SimpleResultEvent implements MouseOnResultListener.ResultEvent {

    private Point point;
    private Float robustness;

    public SimpleResultEvent(Point point, Float robustness) {
        this.point = point;
        this.robustness = robustness;
    }

    @Override
    public Point getPoint() {
        return point;
    }

    @Override
    public Float getRobustness() {
        return robustness;
    }

}
