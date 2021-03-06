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
package org.sybila.parasim.model.verification.stl;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import org.sybila.parasim.model.ode.PointVariableMapping;
import org.sybila.parasim.model.xml.FloatFactory;
import org.sybila.parasim.model.xml.XMLFormatException;
import org.sybila.parasim.model.xml.XMLRepresentableFactory;
import org.sybila.parasim.util.Pair;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * Factory creating {@link LinearPredicate} objects from XML.
 *
 * Needs mapping between variable names and indices.
 *
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 *
 */
public class LinearPredicateFactory implements
        XMLRepresentableFactory<LinearPredicate> {

    public static final String PREDICATE_NAME = "predicate";
    public static final String VARIABLE_NAME = "variable";
    public static final String MULTIPLIER_ATTRIBUTE = "multiplier";
    public static final String FREEZE_ATTRIBUTE = "frozen";
    public static final String VALUE_NAME = "value";
    private PointVariableMapping mapping;

    /**
     * Creates a new factory with designated variable name-to-index mapping.
     *
     * @param mapping Mapping between variable names and indices.
     */
    public LinearPredicateFactory(PointVariableMapping mapping) {
        if (mapping == null) {
            throw new IllegalArgumentException(
                    "Must receive non-null variable mapping.");
        }
        this.mapping = mapping;
    }

    @Override
    public LinearPredicate getObject(Node source) throws XMLFormatException {
        if (!source.getNodeName().equals(PREDICATE_NAME)) {
            throw new XMLFormatException(
                    "Predicate should have <predicate> node as a root (found `" + source.getNodeName() + "').");
        }
        Float value = null;
        Map<Pair<Integer, Integer>, Float> multipliers = new HashMap<>();
        LinearPredicate.Type type = null;

        NodeList children = source.getChildNodes();
        for (int index = 0; index < children.getLength(); index++) {
            Node child = children.item(index);
            String name = child.getNodeName();
            // variable and multiplier //
            switch (name) {
                case VARIABLE_NAME:
                    String varName = child.getFirstChild().getNodeValue();
                    Float multiplier;
                    try {
                        multiplier = new Float(child.getAttributes().getNamedItem(MULTIPLIER_ATTRIBUTE).getNodeValue());
                    } catch (NumberFormatException nfe) {
                        throw new XMLFormatException("Illegible number.", nfe);
                    }
                    Integer frozen = new Integer(child.getAttributes().getNamedItem(FREEZE_ATTRIBUTE).getNodeValue());
                    Integer var = mapping.getKey(varName);
                    if (var == null) {
                        throw new XMLFormatException("Variable `" + varName
                                + "' is not a part of the model.");
                    }
                    Pair<Integer, Integer> key = new Pair<>(var, frozen);
                    if (multipliers.get(key) != null) {
                        throw new XMLFormatException(
                                "Two occurences of a variable of the same name ("
                                + varName + ").");
                    }
                    multipliers.put(key, multiplier);

                    // value //
                    break;
                case VALUE_NAME:
                    if (value == null) {
                        value = new Float(FloatFactory.getObject(child.getFirstChild()));
                    } else {
                        throw new XMLFormatException(
                                "Mulitple value nodes encountered. Expected only one.");
                    }
                    break;
                default:
                    // type //
                    String typeString = child.getNodeName().toUpperCase(
                            Locale.ENGLISH);
                    try {
                        LinearPredicate.Type childType = LinearPredicate.Type.valueOf(typeString);
                        if (type != null) {
                            throw new XMLFormatException(
                                    "Two occurences of predicate type.");
                        }
                        type = childType;
                    } catch (IllegalArgumentException iae) {
                        throw new XMLFormatException("Unknown element: "
                                + child.getNodeName());
                    }
                    break;
            }
        }

        // got everything? //
        if (value == null) {
            throw new XMLFormatException(
                    "Predicate has to contain a right-side value.");
        }
        if (multipliers.isEmpty()) {
            throw new XMLFormatException(
                    "Predicate has to contain at least one variable.");
        }
        if (type == null) {
            throw new XMLFormatException(
                    "Predicate has to be of a defined type.");
        }

        return new LinearPredicate(multipliers, value, type, mapping);
    }
}
