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
package org.sybila.parasim.model.math;

import java.util.Collection;
import java.util.Objects;
import org.sybila.parasim.model.trajectory.Point;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public final class Parameter implements Expression<Parameter>, Indexable {

    private final int index;
    private final int originalIndex;
    private final String name;
    private final SubstitutionValue<Parameter> substitution;

    public Parameter(String name, int index) {
        if (name == null) {
            throw new IllegalArgumentException("The paramater [name] is null.");
        }
        if (index < 0) {
            throw new IllegalArgumentException("The paramater [index] has to be a non negative number.");
        }
        this.name = name;
        this.index = index;
        this.originalIndex = index;
        this.substitution = null;
    }

    private Parameter(String name, int index, int originalIndex) {
        if (name == null) {
            throw new IllegalArgumentException("The paramater [name] is null.");
        }
        if (index < 0) {
            throw new IllegalArgumentException("The paramater [index] has to be a non negative number.");
        }
        if (originalIndex < 0) {
            throw new IllegalArgumentException("The paramater [originalIndex] has to be a non negative number.");
        }
        this.name = name;
        this.index = index;
        this.originalIndex = originalIndex;
        this.substitution = null;
    }

    private Parameter(String name, int index, int originalIndex, SubstitutionValue<Parameter> substitution) {
        if (name == null) {
            throw new IllegalArgumentException("The paramater [name] is null.");
        }
        if (index < 0) {
            throw new IllegalArgumentException("The paramater [index] has to be a non negative number.");
        }
        if (substitution == null) {
            throw new IllegalArgumentException("The paramater [substitution] is null.");
        }
        this.name = name;
        this.index = index;
        this.originalIndex = originalIndex;
        this.substitution = substitution;
    }

    @Override
    public float evaluate(Point point) {
        if (substitution == null) {
            throw new IllegalStateException("Can't evaliuate the parameter, because it hasn't been substituted.");
        }
        return substitution.getValue();
    }

    @Override
    public float evaluate(float[] point) {
        if (substitution == null) {
            throw new IllegalStateException("Can't evaliuate the parameter, because it hasn't been substituted.");
        }
        return substitution.getValue();
    }

    @Override
    public int getIndex() {
        return index;
    }

    public String getName() {
        return name;
    }

    public boolean isSubstituted() {
        return substitution != null;
    }

    @Override
    public Parameter release(Expression... expressions) {
        int indexBefore = 0;
        boolean release = false;
        for (Expression e: expressions) {
            if (!(e instanceof Parameter)) {
                if (e instanceof Indexable && ((Indexable)e).getIndex() < originalIndex) {
                    indexBefore++;
                }
            }
            if (e.equals(this)) {
                release = true;
            }
        }
        if ((release && isSubstituted()) || indexBefore != 0) {
            return new Parameter(name, index+indexBefore, originalIndex);
        } else {
            return this;
        }
    }

    @Override
    public Parameter release(Collection<Expression> expressions) {
        int indexBefore = 0;
        boolean release = false;
        for (Expression e: expressions) {
            if (!(e instanceof Parameter)) {
                if (e instanceof Indexable && ((Indexable)e).getIndex() < originalIndex) {
                    indexBefore++;
                }
            }
            if (e.equals(this)) {
                release = true;
            }
        }
        if ((release && isSubstituted()) || indexBefore != 0) {
            return new Parameter(name, index+indexBefore, originalIndex);
        } else {
            return this;
        }
    }

    @Override
    public Parameter substitute(SubstitutionValue... substitutionValues) {
        int indexBefore = 0;
        for (SubstitutionValue v: substitutionValues) {
            if (!(v instanceof ParameterValue)) {
                if (v.getExpression() instanceof Indexable && ((Indexable) v.getExpression()).getIndex() < index) {
                    indexBefore++;
                }
                continue;
            }
            ParameterValue paramValue = (ParameterValue) v;
            if (paramValue.getExpression().equals(this)) {
                return new Parameter(name, index, originalIndex, paramValue);
            }
            if (paramValue.getExpression().getIndex() < index) {
                indexBefore++;
            }
        }
        if (indexBefore != 0) {
            return new Parameter(name, index - indexBefore, originalIndex);
        } else {
            return this;
        }
    }

    @Override
    public Parameter substitute(Collection<SubstitutionValue> substitutionValues) {
        int indexBefore = 0;
        for (SubstitutionValue v: substitutionValues) {
            if (!(v instanceof ParameterValue)) {
                if (v.getExpression() instanceof Indexable && ((Indexable) v.getExpression()).getIndex() < index) {
                    indexBefore++;
                }
                continue;
            }
            ParameterValue paramValue = (ParameterValue) v;
            if (paramValue.getExpression().equals(this)) {
                return new Parameter(name, index, originalIndex, paramValue);
            }
            if (paramValue.getExpression().getIndex() < index) {
                indexBefore++;
            }
        }
        if (indexBefore != 0) {
            return new Parameter(name, index - indexBefore, originalIndex);
        } else {
            return this;
        }
    }

    @Override
    public String toFormula() {
        return substitution == null ? name : Float.toString(substitution.getValue());
    }

    @Override
    public String toFormula(VariableRenderer renderer) {
        return substitution == null ? name : Float.toString(substitution.getValue());
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 23 * hash + Objects.hashCode(this.name);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Parameter other = (Parameter) obj;
        if (!Objects.equals(this.name, other.name)) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return toFormula();
    }

    @Override
    public <T> T traverse(TraverseFunction<T> function) {
        return function.apply(this);
    }
}