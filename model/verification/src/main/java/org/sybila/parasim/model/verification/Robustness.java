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
package org.sybila.parasim.model.verification;

import org.sybila.parasim.model.trajectory.LimitedDistance;
import org.sybila.parasim.model.trajectory.LimitedPointDistanceMetric;
import org.sybila.parasim.model.trajectory.Point;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface Robustness extends LimitedPointDistanceMetric {

    public static final Robustness UNDEFINED = new Robustness() {

        @Override
        public float getTime() {
            return Float.NaN;
        }

        @Override
        public float getValue() {
            return Float.NaN;
        }

        @Override
        public Robustness invert() {
            return UNDEFINED;
        }

        @Override
        public LimitedDistance distance(float[] first, float[] second) {
            throw new UnsupportedOperationException("Undefined");
        }

        @Override
        public LimitedDistance distance(Point first, Point second) {
            throw new UnsupportedOperationException("Undefined");
        }

        @Override
        public String toString() {
            return "Not defined.";
        }
    };

    float getTime();

    float getValue();

    Robustness invert();
}
