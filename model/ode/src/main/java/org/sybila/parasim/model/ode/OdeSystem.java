package org.sybila.parasim.model.ode;

import org.sybila.parasim.model.trajectory.Point;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public interface OdeSystem {
 
    int dimension();

    Variable getVariable(int dimension);    

    OdeSystemEncoding encoding();
    
    float value(Point point, int dimension);
    
    float value(float[] point, int dimension);

}