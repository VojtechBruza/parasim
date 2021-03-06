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

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;
import org.sybila.parasim.model.ode.PointVariableMapping;
import org.sybila.parasim.model.verification.stl.Formula;
import org.sybila.parasim.model.verification.stl.FormulaResource;
import org.sybila.parasim.model.xml.XMLException;
import org.testng.Assert;

/**
 *
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 */
public class TestFormulae {

    private final Formula phi, psi;
    private static final PointVariableMapping MAPPING = new PointVariableMapping() {

        @Override
        public int getDimension() {
            return 3;
        }

        @Override
        public Integer getKey(String variableName) {
            switch (variableName) {
                case "x":
                    return 0;
                case "y":
                    return 1;
                case "z":
                    return 2;
                default:
                    return null;
            }
        }

        @Override
        public String getName(Integer variableKey) {
            switch (variableKey) {
                case 0:
                    return "x";
                case 1:
                    return "y";
                case 2:
                    return "z";
                default:
                    return null;
            }
        }
    };

    private File getFormulaFile(String resource) {

        URL res = getClass().getClassLoader().getResource(resource);
        try {
            return new File(res.toURI());
        } catch (URISyntaxException urise) {
            Assert.fail("Wrong formula resource.", urise);
        }
        return null;
    }

    public TestFormulae(String phiResource, String psiResource) {
        FormulaResource resource = new FormulaResource(getFormulaFile(phiResource));
        resource.setVariableMapping(MAPPING);
        try {
            resource.load();
        } catch (XMLException xmle) {
            Assert.fail("Cannot load first formula.", xmle);
        }
        phi = resource.getRoot();

        resource = new FormulaResource(getFormulaFile(psiResource));
        resource.setVariableMapping(MAPPING);
        try {
            resource.load();
        } catch (XMLException xmle) {
            Assert.fail("Cannot load second formula.", xmle);
        }
        psi = resource.getRoot();
    }

    public Formula phi() {
        return phi;
    }

    public Formula psi() {
        return psi;
    }
}
