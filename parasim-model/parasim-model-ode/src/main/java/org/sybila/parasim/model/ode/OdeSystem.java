package org.sybila.parasim.model.ode;

import org.sybila.parasim.model.trajectory.Point;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface OdeSystem {
 
    int dimension();
    
    String octaveName();
    
    String octaveString();
    
    float value(Point point, int dimension);
    
    float value(float[] point, int dimension);
}
