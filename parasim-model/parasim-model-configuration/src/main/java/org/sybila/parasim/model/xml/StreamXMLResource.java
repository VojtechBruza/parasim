package org.sybila.parasim.model.xml;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMResult;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;
import javax.xml.validation.Validator;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

/**
 * Stores and loads objects from streams (i.e. InputStream and OutputStream).
 * 
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 *
 * @param <T> Type of object being stored/loaded.
 */
public abstract class StreamXMLResource<T extends XMLRepresentable> implements XMLResource<T> {

    private static final String SCHEMA_XMLNS = "http://www.w3.org/2001/XMLSchema";
    private T root = null;

    @Override
    public void setRoot(T target) {
        root = target;
    }

    @Override
    public T getRoot() {
        return root;
    }

    /**
     * Specifies way of obtaining input stream.
     * @return a new open input stream.
     * @throws XMLException when stream could not be opened.
     */
    protected abstract InputStream openInputStream() throws XMLException;

    /**
     * Specifies way of obtaining output stream.
     * @return a new open output stream.
     * @throws XMLException when stream could not be opened.
     */
    protected abstract OutputStream openOutputStream() throws XMLException;

    /**
     * Specifies way of obtaining factory for contained objects.
     * @return Factory used to transform XML into contained objects.
     */
    protected abstract XMLRepresentableFactory<T> getFactory();

    /**
     * Specifies way of obtaining XLM schema for validation.
     * @return URL of XML schema used to validate input.
     */
    protected abstract URL getXMLSchema();

    /**
     * Specifies way of obtaining name of XML namespace
     * @return XML namespace of contained objects.
     */
    protected abstract String getNamespace();

    @Override
    public void load() throws XMLException {
        /* get schema */
        SchemaFactory schFact = SchemaFactory.newInstance(SCHEMA_XMLNS);
        Schema schema;
        try {
            schema = schFact.newSchema(getXMLSchema());
        } catch (SAXException saxe) {
            throw new XMLException("Schema could not be lodaded.", saxe);
        }
        Validator valid = schema.newValidator();

        /* get document builder (input parser) */
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        dbf.setNamespaceAware(true);
        DocumentBuilder docBuild;
        try {
            docBuild = dbf.newDocumentBuilder();
        } catch (ParserConfigurationException pce) {
            throw new XMLException("Parser could not be configured.", pce);
        }

        /* parsing */
        Document doc;
        InputStream is = openInputStream();
        try {
            doc = docBuild.parse(is);
        } catch (SAXException saxe) {
            throw new XMLException("Parse error during document parsing.", saxe);
        } catch (IOException ioe) {
            throw new XMLException("IO error during document parsing.", ioe);
        }

        /* validation; also enhances result (default attribute values, ...) */
        DOMResult result = new DOMResult();
        try {
            valid.validate(new DOMSource(doc), result);
        } catch (SAXException saxe) {
            throw new XMLException("Parse error during document validation.", saxe);
        } catch (IOException ioe) {
            throw new XMLException("IO error during document validation.", ioe);
        }
        Document input = (Document) result.getNode();
        input.normalize();

        /* get target */
        setRoot(getFactory().getObject(input.getDocumentElement()));

        /* close stream */
        try {
            is.close();
        } catch (IOException ioe) {
            throw new XMLException("Input could not be closed", ioe);
        }
    }

    @Override
    public void store() throws XMLException {
        if (root == null) {
            throw new IllegalStateException("Nothing to store."); //TODO
        }
        OutputStream os = openOutputStream();

        /* get empty [!] document */
        DocumentBuilderFactory docFact = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuild;
        try {
            docBuild = docFact.newDocumentBuilder();
        } catch (ParserConfigurationException pce) {
            throw new XMLException("Transfer could not be configured.", pce);
        }
        Document doc = docBuild.newDocument();

        /* transform to xml */
        Element target = root.toXML(doc);
        setNamespace(target);
        doc.appendChild(target);

        /* get transformer (output parser) */
        TransformerFactory transFact = TransformerFactory.newInstance();
        Transformer trans;
        try {
            trans = transFact.newTransformer();
        } catch (TransformerConfigurationException tce) {
            throw new XMLException("Output parser could not be configured.");
        }
        trans.setOutputProperty(OutputKeys.INDENT, "yes"); //spreads output over multiple lines

        /* print XML */
        try {
            trans.transform(new DOMSource(doc), new StreamResult(os));
        } catch (TransformerException te) {
            throw new XMLException("Error during XML output.", te);
        }

        /* close stream */
        try {
            os.close();
        } catch (IOException ioe) {
            throw new XMLException("Output could not be closed properly.", ioe);
        }
    }

    private void setNamespace(Element target) {
        target.setAttribute("xmlns", getNamespace());
    }
}
