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
package org.sybila.parasim.computation.verification;

import org.sybila.parasim.computation.verification.api.VerifiedDataBlock;
import org.sybila.parasim.model.verification.stl.UntilFormula;
import org.sybila.parasim.model.verification.stl.OrFormula;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.sybila.parasim.computation.verification.api.STLVerifier;
import org.sybila.parasim.core.annotations.Default;
import org.sybila.parasim.model.ode.PointVariableMapping;
import org.sybila.parasim.model.trajectory.ArrayDataBlock;
import org.sybila.parasim.model.trajectory.ArrayPoint;
import org.sybila.parasim.model.trajectory.ListTrajectory;
import org.sybila.parasim.model.trajectory.Point;
import org.sybila.parasim.model.trajectory.Trajectory;
import org.sybila.parasim.model.verification.stl.AndFormula;
import org.sybila.parasim.model.verification.stl.Formula;
import org.sybila.parasim.model.verification.stl.FutureFormula;
import org.sybila.parasim.model.verification.stl.GloballyFormula;
import org.sybila.parasim.model.verification.stl.IntervalBoundaryType;
import org.sybila.parasim.model.verification.stl.LinearPredicate;
import org.sybila.parasim.model.verification.stl.LinearPredicate.Type;
import org.sybila.parasim.model.verification.stl.Predicate;
import org.sybila.parasim.model.verification.stl.TimeInterval;
import org.testng.annotations.Test;
import static org.testng.Assert.*;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class TestSTLFormulaVerification extends AbstractVerificationTest {

    @Test
    public void testAnd() {
        Formula formula = new AndFormula(
                createPredicate(createMapping(), Type.GREATER, 1),
                createPredicate(createMapping(), Type.GREATER, 10));
        STLVerifier verifier = getManager().resolve(STLVerifier.class, Default.class, getManager().getRootContext());
        VerifiedDataBlock verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(1, 2, 3) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), -9f);
        verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(11, 10, 9) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), 1f);
    }

    @Test
    public void testFuture() {
        Formula formula = new FutureFormula(createPredicate(createMapping(), Type.GREATER, 10), new TimeInterval(2, 4, IntervalBoundaryType.CLOSED));
        STLVerifier verifier = getManager().resolve(STLVerifier.class, Default.class, getManager().getRootContext());
        VerifiedDataBlock verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(1, 2, 3, 11, 12, 13) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), 2f);
        verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(11, 10, 9, 8, 7, 6, 5) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), -1f);
    }

    @Test
    public void testGlobally() {
        Formula formula = new GloballyFormula(createPredicate(createMapping(), Type.GREATER, 10), new TimeInterval(2, 4, IntervalBoundaryType.CLOSED));
        STLVerifier verifier = getManager().resolve(STLVerifier.class, Default.class, getManager().getRootContext());
        VerifiedDataBlock verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(1, 2, 3, 10, 11, 12) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), -7f);
        verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(11, 11, 12, 24, 13, 6, 5) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), 2f);
    }

    @Test
    public void testOr() {
        Formula formula = new OrFormula(
                createPredicate(createMapping(), Type.GREATER, 1),
                createPredicate(createMapping(), Type.GREATER, 10));
        STLVerifier verifier = getManager().resolve(STLVerifier.class, Default.class, getManager().getRootContext());
        VerifiedDataBlock verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(2, 3, 4) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), 1f);
        verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(0, 1, 2) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), -1f);
    }

    @Test
    public void testUntil() {
        Formula formula = new UntilFormula(
                createPredicate(createMapping(), Type.GREATER, 9),
                createPredicate(createMapping(), Type.LESSER, 15),
                new TimeInterval(2, 4, IntervalBoundaryType.CLOSED)
                );
        STLVerifier verifier = getManager().resolve(STLVerifier.class, Default.class, getManager().getRootContext());
        VerifiedDataBlock verified = verifier.verify(
                new ArrayDataBlock<Trajectory>(new Trajectory[] { createTrajectory(10, 10, 16, 12, 15, 9) }),
                formula);
        assertEquals(verified.getRobustness(0).getValue(), 1f);
    }

    private Trajectory createTrajectory(float...oneDimensionPoints) {
        List<Point> points = new ArrayList<Point>();
        for (int i=0; i<oneDimensionPoints.length; i++) {
            points.add(new ArrayPoint(i, oneDimensionPoints[i]));
        }
        return new ListTrajectory(points);
    }

    private Predicate createPredicate(PointVariableMapping mapping, Type type, float rightSide) {
        Map<Integer, Float> left = new HashMap<Integer, Float>();
        left.put(0, 1f);
        return new LinearPredicate(left, rightSide, type, mapping);
    }

    private PointVariableMapping createMapping() {
        return new PointVariableMapping() {
            public int getDimension() {
                return 1;
            }
            public Integer getKey(String variableName) {
                return 0;
            }
            public String getName(Integer variableKey) {
                return "X";
            }
        };
    }
}
