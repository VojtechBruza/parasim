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
package org.sybila.parasim.visualisation.plot.impl.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.Set;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import org.sybila.parasim.model.ode.OdeSystem;
import org.sybila.parasim.model.ode.OdeVariableMapping;
import org.sybila.parasim.model.ode.PointVariableMapping;
import org.sybila.parasim.model.space.OrthogonalSpace;
import org.sybila.parasim.model.trajectory.ArrayPoint;
import org.sybila.parasim.model.verification.Robustness;
import org.sybila.parasim.model.verification.result.AbstractVerificationResult;
import org.sybila.parasim.model.verification.result.VerificationResult;
import org.sybila.parasim.visualisation.plot.api.Plotter;
import org.sybila.parasim.visualisation.plot.api.PlotterWindowListener;
import org.sybila.parasim.visualisation.plot.api.MouseOnResultListener;
import org.sybila.parasim.visualisation.plot.impl.LayerFactory;
import org.sybila.parasim.visualisation.plot.impl.LayerMetaFactory;
import org.sybila.parasim.visualisation.plot.impl.ResultPlotterConfiguration;
import org.sybila.parasim.visualisation.plot.impl.SimpleResultEvent;
import org.sybila.parasim.visualisation.plot.impl.SpaceUtils;

/**
 * Plots a 2D projection of a generally n-D space in a window. Exactly two axes
 * are taken for the 2D space base, the other coordinates are left out. For the
 * points not to overlap, left out axes are divided into intervals which may be
 * chosen using sliders.
 *
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 */
public class ProjectionPlotter extends JFrame implements Plotter {

    private static final int INSET = 10;
    private static final Insets PADDING = new Insets(INSET, INSET, INSET, INSET);
    //components//
    private JPanel sliders, axes;
    private Canvas canvas;
    private AxisChooser xAxis, yAxis;
    private AxisSlider[] axisSliders;
    private Rule hRule, vRule;
    private StatusBar status;
    private CanvasPane canvasPane;
    //variables//
    private PointVariableMapping names;
    private OrthogonalSpace extent;
    private int dimension;
    private LayerMetaFactory metaLayers;
    private LayerFactory layers;
    private Set<MouseOnResultListener> mouseOnResultListeners = new HashSet<>();
    private Set<PlotterWindowListener> plotterWindowListeners = new HashSet<>();

    /**
     * Creates new plotter on a given verification result with specified axes
     * labels, algorithm of projection into 2D and point appearance.
     *
     * @param result Result of verification. Is not necessarily rendered.
     * @param names Labels of axes.
     * @param pointSource Specifies manner of projection into 2D and contains
     * rendered points.
     * @param pointAppearance Specifies point appearance.
     */
    public ProjectionPlotter(ResultPlotterConfiguration conf, VerificationResult result, OdeSystem odeSystem, LayerMetaFactory pointSource, PointRenderer pointAppearance) {
        dimension = result.getPoint(0).getDimension();
        this.names = new OdeVariableMapping(odeSystem);
        extent = (new SpaceUtils(conf)).provideWithPadding(AbstractVerificationResult.getEncompassingSpace(result, odeSystem));

        init(conf, pointAppearance, conf.getPlotterWindowWidth(), conf.getPlotterWindowHeight());
        initRobustnessLabel(result.getGlobalRobustness(), conf);

        metaLayers = pointSource;
        //initially, (0,1) are chosen//
        layers = metaLayers.getLayerFactory(0, 1);
        //updating sliders//
        for (int i = 0; i < dimension; i++) {
            axisSliders[i].update(layers.ticks(i), 0);
        }
        updateRules(0, 1);
        updateView();
    }

    @Override
    public void addMouseOnResultListener(MouseOnResultListener listener) {
        mouseOnResultListeners.add(listener);
    }

    @Override
    public void addPlotterWindowListener(PlotterWindowListener listener) {
        plotterWindowListeners.add(listener);
    }

    @Override
    public void removeMouseOnResultListener(MouseOnResultListener listener) {
        mouseOnResultListeners.remove(listener);
    }

    @Override
    public void removePlotterWindowListener(PlotterWindowListener listener) {
        plotterWindowListeners.remove(listener);
    }

    private void init(ResultPlotterConfiguration conf, PointRenderer appearance, int width, int height) {
        ResourceBundle strings = ResourceBundle.getBundle(getClass().getSimpleName());
        setTitle(strings.getString("title"));

        setSize(width, height);
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        //setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);

        setLayout(new BorderLayout());

        initSliders();
        initAxes(conf, strings.getString("x_axis"), strings.getString("y_axis"));
        initCanvas(appearance, conf);
        addWindowListener(new WindowListener() {

            @Override
            public void windowOpened(WindowEvent e) {
            }

            @Override
            public void windowClosing(WindowEvent e) {
            }

            @Override
            public void windowClosed(WindowEvent e) {
                for (PlotterWindowListener listener : plotterWindowListeners) {
                    listener.windowClosed(PlotterWindowListener.PlotterWindowEvent.CLOSED);
                }
            }

            @Override
            public void windowIconified(WindowEvent e) {
            }

            @Override
            public void windowDeiconified(WindowEvent e) {
            }

            @Override
            public void windowActivated(WindowEvent e) {
            }

            @Override
            public void windowDeactivated(WindowEvent e) {
            }
        });
    }

    private void initRobustnessLabel(Robustness robustness, ResultPlotterConfiguration conf) {
        JPanel robustnessPanel = new JPanel();
        robustnessPanel.setLayout(new BoxLayout(robustnessPanel, BoxLayout.LINE_AXIS));

        JLabel label = new JLabel("Global robustness: ");
        JLabel robustnessLabel = new JLabel(robustness.toString());
        Font font = label.getFont().deriveFont(15f);
        label.setFont(font);
        robustnessLabel.setFont(font);
        Color robustnessColor = (robustness.getValue() > 0) ? conf.getPointColorValid() : conf.getPointColorInvalid();
        robustnessLabel.setForeground(robustnessColor.darker());

        robustnessPanel.add(label);
        robustnessPanel.add(robustnessLabel);
        robustnessPanel.setBorder(new EmptyBorder(5, 15, 0, 5));
        add(robustnessPanel, BorderLayout.PAGE_START);
    }

    private void initCanvas(PointRenderer appearance, ResultPlotterConfiguration conf) {
        canvas = new Canvas(appearance);
        canvasPane = new CanvasPane(conf, canvas, new CanvasPane.PositionChangeListener() {

            @Override
            public void updatePosition(float x, float y) {
                status.setValue(xAxis.getSelected(), x);
                status.setValue(yAxis.getSelected(), y);
            }
        });
        canvasPane.addMouseListener(new MouseAdapter() {

            @Override
            public void mouseClicked(MouseEvent e) {
                // call 'click on result' listeners
                float[] pointData = new float[extent.getDimension()];
                for (int dim = 0; dim < axisSliders.length; dim++) {
                    pointData[dim] = status.getValue(dim);
                }
                MouseOnResultListener.ResultEvent event = new SimpleResultEvent(new ArrayPoint(0, pointData), null);
                for (MouseOnResultListener listener : mouseOnResultListeners) {
                    listener.click(event);
                }
            }
        });
        hRule = new Rule(conf, Rule.Orientation.HORIZONTAL);
        vRule = new Rule(conf, Rule.Orientation.VERTICAL);

        JPanel canvasPanel = new JPanel(new GridBagLayout());
        canvasPanel.setBorder(new EmptyBorder(PADDING));

        GridBagConstraints pos = getDefaultConstraints();
        pos.weightx = 1;
        pos.weighty = 1;
        pos.gridx = 1;
        pos.gridy = 0;
        canvasPanel.add(canvasPane, pos);

        pos = getDefaultConstraints();
        pos.gridx = 0;
        pos.gridy = 0;
        canvasPanel.add(vRule, pos);

        pos = getDefaultConstraints();
        pos.gridx = 1;
        pos.gridy = 1;
        canvasPanel.add(hRule, pos);

        add(canvasPanel, BorderLayout.CENTER);
    }

    private GridBagConstraints getDefaultConstraints() {
        GridBagConstraints result = new GridBagConstraints();
        result.gridwidth = 1;
        result.gridheight = 1;
        result.fill = GridBagConstraints.BOTH;
        result.weightx = 0;
        result.weighty = 0;
        return result;
    }

    private void initSliders() {
        sliders = new JPanel();
        sliders.setLayout(new BoxLayout(sliders, BoxLayout.LINE_AXIS));
        sliders.setBorder(new EmptyBorder(PADDING));

        ChangeListener changed = new ChangeListener() {

            @Override
            public void stateChanged(ChangeEvent ce) {
                updateView();
            }
        };
        axisSliders = new AxisSlider[dimension];
        for (int i = 0; i < dimension; i++) {
            axisSliders[i] = new AxisSlider(names.getName(i), changed, extent.getMinBounds().getValue(i), extent.getMaxBounds().getValue(i));
            sliders.add(axisSliders[i]);
        }
        axisSliders[0].setActive(false);
        axisSliders[1].setActive(false);

        add(sliders, BorderLayout.LINE_END);
    }

    private void initAxes(ResultPlotterConfiguration conf, String xAxe, String yAxe) {
        AxisChooser[] axisChoosers = AxisChooser.getPairedAxes(dimension, names, new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                updateAxes();
            }
        });
        xAxis = axisChoosers[0];
        yAxis = axisChoosers[1];

        JPanel bottom = new JPanel();
        bottom.setLayout(new BoxLayout(bottom, BoxLayout.PAGE_AXIS));

        axes = new JPanel();
        axes.setLayout(new BoxLayout(axes, BoxLayout.LINE_AXIS));
        axes.setBorder(new EmptyBorder(PADDING));
        axes.add(new JLabel(xAxe));
        axes.add(Box.createRigidArea(new Dimension(INSET, INSET)));
        axes.add(xAxis);
        axes.add(Box.createRigidArea(new Dimension(2 * INSET, INSET)));
        axes.add(new JLabel(yAxe));
        axes.add(Box.createRigidArea(new Dimension(INSET, INSET)));
        axes.add(yAxis);
        axes.add(Box.createHorizontalGlue());
        bottom.add(axes);

        status = new StatusBar(conf, dimension, names);
        bottom.add(status);

        add(bottom, BorderLayout.PAGE_END);
    }

    /**
     * Called when an axis is selected by AxisChoosers.
     */
    private void updateAxes() {
        int xSelected = xAxis.getSelected();
        int ySelected = yAxis.getSelected();

        // change sliders //
        for (AxisSlider as : axisSliders) {
            as.setActive(true);
        }
        axisSliders[xSelected].setActive(false);
        axisSliders[ySelected].setActive(false);

        float[] values = new float[dimension];
        //get real values from axisSliders//
        for (int i = 0; i < dimension; i++) {
            values[i] = layers.getValue(i, axisSliders[i].getValue());
        }
        //create new LayerFactory//
        layers = metaLayers.getLayerFactory(xSelected, ySelected);
        //update values and maximums of axissliders//
        for (int i = 0; i < dimension; i++) {
            axisSliders[i].update(layers.ticks(i), layers.getTicks(i, values[i]));
        }
        updateRules(xSelected, ySelected);
        updateView();
    }

    /**
     * Called when an axis slider is moved and view has to be changed
     * accordingly.
     */
    private void updateView() {
        Map<Integer, Integer> projections = new HashMap<Integer, Integer>();
        for (int i = 0; i < dimension; i++) {
            int ticks = axisSliders[i].getValue();
            projections.put(i, ticks);
            status.setValue(i, layers.getValue(i, ticks));
        }
        canvas.setPoints(layers.getLayer(projections));
    }

    private void updateRules(int selHor, int selVer) {
        hRule.update(extent.getMinBounds().getValue(selHor), extent.getMaxBounds().getValue(selHor));
        vRule.update(extent.getMinBounds().getValue(selVer), extent.getMaxBounds().getValue(selVer));
    }

    @Override
    public void plot() {
        setVisible(true);
    }
}
