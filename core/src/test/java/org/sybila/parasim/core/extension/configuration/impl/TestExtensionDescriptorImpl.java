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
package org.sybila.parasim.core.extension.configuration.impl;

import org.sybila.parasim.core.extension.configuration.api.ParasimDescriptor;
import org.sybila.parasim.core.extension.configuration.api.ExtensionDescriptor;
import java.io.IOException;
import java.util.Collection;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import org.xml.sax.SAXException;
import static org.testng.Assert.*;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class TestExtensionDescriptorImpl {

    private ParasimDescriptor descriptor;

    @BeforeMethod
    public void prepareDescriptor() throws IOException, SAXException {
        descriptor = ParasimDescriptorImpl.fromXMLFile("something", "src/test/resources/org/sybila/parasim/core/extension/configuration/parasim.xml");
    }

    @Test
    public void testExtensions() {
       Collection<ExtensionDescriptor> extensions = descriptor.getExtensionDescriptors();
       assertEquals(extensions.size(), 1);
       ExtensionDescriptor extDescriptor = descriptor.getExtensionDescriptor("test");
       assertEquals(extDescriptor.getProperty("test-property"), "test-value");
       assertFalse(extDescriptor.isPropertyArray("test-property"));
       assertTrue(extDescriptor.isPropertyArray("test-array-property"));
       String[] arrayValue = new String[] {"1", "2", "3"};
       assertEquals(extDescriptor.getPropertyAsArray("test-array-property"), arrayValue);
    }

}