package org.sybila.parasim.visualisation.plot.impl;

/**
 *
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 */
public class ResultPlotterConfiguration {
    private int plotterWindowHeight = 500;
    private int plotterWindowWidth = 750;
    private float minimumDifference = 0.001f;

    public int getPlotterWindowHeight() {
        return plotterWindowHeight;
    }

    public int getPlotterWindowWidth() {
        return plotterWindowWidth;
    }

    public float getMinimumDifference() {
        return minimumDifference;
    }

}