package org.sybila.parasim.model.verification.stl;

/**
 * Represents an STL (Signal Temporal Logic) formula.
 *
 * @author <a href="mailto:sven@mail.muni.cz">Sven Dra�an</a>
 */
public interface Formula
{
    /**
     * Returns the arity of this formula, that is how many subformulas it consists of.
     *
     * @return Number of subformulas.
     */
    int getArity();

    /**
     * If getArity() > 0 then each subformula may be obtained by this method.
     *
     * @param index Index of subformula to return.
     * @return Subformula with given index.
     */
    Formula getSubformula(int index);

    /**
     * Formulas such as Until, Future and Globaly have as a parameter
     * an interval. To evaluate the satisfaction of a formula on a trajectory
     * a minimal length is needed which can be computed from the structure
     * of the formula.
     * 
     * This method returns the time needed to evaluate this formula and all
     * its subformulas. If a formula is evaluated in a single time point 0
     * will be returned.
     *
     * @return Time or trajectory length needed to evaluate formula.
     */
    float getTimeNeeded();

    @Override
    String toString();

    /**
     * Returns the type of the formula.
     * 
     * @return Type of formula.
     */
    FormulaType getType();

    /**
     * Enables the comparison of formulas. Two formulas are equal if they
     * represent the same predicate or all their subformulas are equal.
     *
     * @param object Formula to compare to.
     * @return True if formula is equal to this object, false otherwise.
     */
    boolean equals(Formula formula);
}