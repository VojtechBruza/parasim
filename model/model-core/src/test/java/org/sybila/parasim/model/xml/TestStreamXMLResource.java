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
package org.sybila.parasim.model.xml;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import java.net.URL;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

public class TestStreamXMLResource extends StreamXMLResource<Target> {
    private static final String RESULT = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
            + "<test xmlns=\"http://www.sybila.org/parasim/TestSchema\">This is test target.</test>\n";
    private class TestInput extends ByteArrayInputStream {
        private boolean closed = false;

        public TestInput() {
            super(TestStreamXMLResource.RESULT.getBytes());
        }

        @Override
        public void close() throws IOException {
            super.close();
            closed = true;
        }

        public boolean isClosed() {
            return closed;
        }
    }
    private class TestOutput extends ByteArrayOutputStream {
        private boolean closed = false;

        @Override
        public void close() throws IOException {
            super.close();
            closed = true;
        }

        public boolean isClosed() {
            return closed;
        }
    }
    private TestInput in;
    private TestOutput out;
    
    @Override
    protected String getNamespace() {
        return "http://www.sybila.org/parasim/TestSchema";
    }

    @Override
    protected URL getXMLSchema() {
        return getClass().getClassLoader().getResource("org/sybila/parasim/model/xml/testSchema.xsd");
    }

    @Override
    protected XMLRepresentableFactory<Target> getFactory() {
        return new TargetFactory();
    }
    
    @Override
    protected InputStream openInputStream() throws XMLException {
        in = new TestInput();
        return in;
    }
    
    @Override
    protected OutputStream openOutputStream() throws XMLException {
        out = new TestOutput();
        return out;
    }
    
    @BeforeTest
    public void clearBuffers() {
        in = null;
        out = null;
    }
    
    /**
     * Creates a new contained object (using {@link #setRoot(Target)}).
     * Tests whether {@link #getRoot()} returns the same object.
     */
    @Test
    public void testRootGetterReturnsTheSame() {
        Target testTarget = new Target();
        setRoot(testTarget);
        Assert.assertSame(testTarget, getRoot(), "getRoot() does not return the object set by setRoot().");
    }
    
    /**
     * Creates a new contained object (using {@link #setRoot(Target)}).
     * Tests whether after using {@link #store()}, {@link #getRoot()}
     * returns the same object.
     */
    @Test
    public void testStoreDoesNotChangeRoot() {
        Target testTarget = new Target();
        setRoot(testTarget);
        try {
            store();
        } catch (XMLException xmle) {
            xmle.printStackTrace();
        } finally {
            Assert.assertSame(testTarget, getRoot(), "getRoot() does not return the object set by setRoot() (aftre store() was invoked)");
        }
    }
    
    /**
     * Creates a new contained object, stores it and tests whether it is stored correctly.
     */
    @Test
    public void testStore() {
        setRoot(new Target());
        try {
            store();
        } catch (XMLException xmle) {
            xmle.printStackTrace();
        }
        Assert.assertTrue(out.isClosed(), "Output stream was not closed");
        Assert.assertEquals(RESULT, out.toString(), "XML was not stored correctly");
    }
    
    /**
     * Tries to load an object from an input and tests whether it is loaded correctly.
     */
    @Test
    public void testLoad() {
        try {
            load();
        } catch (XMLException xmle) {
            xmle.printStackTrace();
        }
        Assert.assertTrue(in.isClosed(), "Input stream was not closed.");
        Assert.assertEquals(new Target(), getRoot(), "XML was not loaded correctly");
    }

}

class Target implements XMLRepresentable {
    private String text = null;
    public Target() {
        text = "This is test target.";
    }
    public Target(String text) {
        this.text = text;
    }
    @Override
    public Element toXML(Document doc) {
        Element target = doc.createElement("test");
        target.appendChild(doc.createTextNode(text));
        return target;
    }
    @Override
    public boolean equals(Object obj) {
        if (obj == this)
            return true;
        if (!(obj instanceof Target))
            return false;
        return text.equals(((Target) obj).text);
    }
    @Override
    public int hashCode() {
        return text.hashCode();
    }
}

class TargetFactory implements XMLRepresentableFactory<Target> {
    @Override
    public Target getObject(Node source) throws XMLFormatException {
        String text = source.getFirstChild().getNodeValue();
        return new Target(text);
    }
}
