package org.sybila.parasim.model.verification.result;

import org.sybila.parasim.model.space.OrthogonalSpace;
import org.sybila.parasim.model.trajectory.ArrayPoint;
import org.sybila.parasim.model.trajectory.Point;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Implements some auxiliary methods of {@link VerificationResult}
 * according with use of its interface methods.
 *
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 */
public abstract class AbstractVerificationResult implements VerificationResult {
    
    /**
     * Return smallest space including all points comprising this result.
     * Linear to the number of points.
     * 
     * Might be elevated into a factory in the future.
     * 
     * @return Smallest space including all points.
     */
    public OrthogonalSpace getEncompassingSpace() {
        if (size() == 0) {
            Point empty = new ArrayPoint(0, new float [0], 0, 0);
            return new OrthogonalSpace(empty, empty);
        }
        int dims = getPoint(0).getDimension();

        float [] mins = new float[dims];
        float [] maxs = new float[dims];
        
        for (int dim = 0; dim < dims; dim++) {
            mins[dim] = getPoint(0).getValue(dim);
            maxs[dim] = getPoint(0).getValue(dim);
        }
        
        for (int pointIndex = 1; pointIndex < size(); pointIndex++) {
            for (int dim = 0; dim < dims; dim++) {
                float val = getPoint(pointIndex).getValue(dim);
                if (val < mins[dim]) {
                    mins[dim] = val;
                } else if (val > maxs[dim]) {
                    maxs[dim] = val;
                }
            }
        }
        
        return new OrthogonalSpace(new ArrayPoint(0, mins, 0, dims), new ArrayPoint(0, maxs, 0, dims));
    }

    @Override
    public Element toXML(Document doc) {
        Element target = doc.createElement(VerificationResultFactory.RESULT_NAME);
        for (int i = 0; i < size(); i++) {
            target.appendChild(pointToXML(doc, i));
        }

        return target;
    }

    private Element pointToXML(Document doc, int index) {
        Point p = getPoint(index);
        float r = getRobustness(index);

        Element target = doc.createElement(VerificationResultFactory.POINT_NAME);
        for (int i = 0; i < p.getDimension(); i++) {
            Element dim = doc.createElement(VerificationResultFactory.DIMENSION_NAME);
            dim.appendChild(doc.createTextNode(Float.toString(p.getValue(i))));
            target.appendChild(dim);
        }
        Element rob = doc.createElement(VerificationResultFactory.ROBUSTNESS_NAME);
        rob.appendChild(doc.createTextNode(Float.toString(r)));
        target.appendChild(rob);

        return target;
    }

    @Override
    public boolean equals(Object obj) {
        //if (obj == this) return true;
        if (!(obj instanceof VerificationResult)) {
            return false;
        }
        VerificationResult target = (VerificationResult) obj;
        if (size() != target.size()) {
            return false;
        }
        for (int i = 0; i < size(); i++) {
            if (!getPoint(i).equals(target.getPoint(i))) {
                return false;
            }
            if (Float.floatToIntBits(getRobustness(i)) != Float.floatToIntBits(target.getRobustness(i))) {
                return false;
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        final int prime = 41;
        int result = size();
        for (int i = 0; i < size(); i++) {
            result *= prime;
            result += getPoint(i).hashCode();
            result *= prime;
            result += Float.floatToIntBits(getRobustness(i));
        }
        return result;
    }
}
