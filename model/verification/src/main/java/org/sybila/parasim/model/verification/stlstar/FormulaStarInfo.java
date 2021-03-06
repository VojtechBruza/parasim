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
package org.sybila.parasim.model.verification.stlstar;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang3.Validate;
import org.sybila.parasim.model.verification.stl.Formula;
import org.sybila.parasim.model.verification.stl.FormulaType;
import org.sybila.parasim.model.verification.stl.FreezeFormula;
import org.sybila.parasim.model.verification.stl.Predicate;

/**
 * Contains information about frozen time indices in all subformulas of a
 * formula.
 *
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 */
public class FormulaStarInfo {

    private final Formula root;
    private final Map<Formula, Set<Integer>> unfrozen = new HashMap<>();
    private final int starNumber;

    private int analyze(Formula target) {
        // if we are at predicate, finish //
        if (target.getType().equals(FormulaType.PREDICATE)) {
            Predicate targetPredicate = (Predicate) target;
            unfrozen.put(target, new HashSet<>(targetPredicate.getStars()));
            return targetPredicate.getStarNumber();
        }

        // analyse all subformulae //
        int max = 0;
        Set<Integer> current = new HashSet<>();
        for (int i = 0; i < target.getArity(); i++) {
            Formula subformula = target.getSubformula(i);
            // analyze subformula //
            int analyzed = analyze(subformula);

            // update star number //
            if (analyzed > max) {
                max = analyzed;
            }

            // update unfrozen indices -- subformula already analyzed //
            current.addAll(getUnfrozenIndices(subformula));
        }

        // solve freeze operator //
        if (target.getType().equals(FormulaType.FREEZE)) {
            current.remove(((FreezeFormula) target).getFreezeIndex());
        }

        // put current //
        unfrozen.put(target, current);

        return max;
    }

    /**
     * Compute information about a formula.
     *
     * @param root Analyzed formula.
     */
    public FormulaStarInfo(Formula root) {
        Validate.notNull(root);
        this.root = root;
        starNumber = analyze(root);
    }

    /**
     * Access the root of the formula.
     *
     * @return The analyzed formula.
     */
    public Formula getRoot() {
        return root;
    }

    /**
     * Return the collection of indices which are not frozen at the subformula
     * (i.e. their associated frozen-time values have to be distinguished).
     *
     * @param subformula Subformula whose unfrozen indices should be returned.
     * @return Set of unfrozen indices.
     */
    public Collection<Integer> getUnfrozenIndices(Formula subformula) {
        Set<Integer> target = unfrozen.get(subformula);
        if (target == null) {
            throw new IllegalArgumentException("Given formula is not a subformula of the root formula.");
        }
        return Collections.unmodifiableCollection(target);
    }

    /**
     * Get the number of the highest frozen-time index used in formula.
     *
     * @return Highest frozen-time index.
     */
    public int getStarNumber() {
        return starNumber;
    }

    /**
     * Return maximum number of unfrozen indices within one formula.
     *
     * @return Number of unfrozen indices.
     */
    public int getMaxUnfrozenIndices() {
        int result = 0;
        for (Set<Integer> indices : unfrozen.values()) {
            if (indices.size() > result) {
                result = indices.size();
            }
        }
        return result;
    }
}
