package org.sybila.parasim.extension.projectmanager.names;

import java.util.Properties;
import org.sybila.parasim.extension.projectmanager.project.ResourceException;

/**
 *
 * @author <a href="mailto:xvejpust@fi.muni.cz">Tomáš Vejpustek</a>
 */
public class ExperimentNames {

    private String formulaName, initialSpaceName, simulationSpaceName, precisionConfigurationName, initialSamplingName, verificationResultName, modelName, annotation;
    private int iterationLimit = 0;
    private long timeout = 0;

    public ExperimentNames() {
        formulaName = null;
        initialSpaceName = null;
        simulationSpaceName = null;
        precisionConfigurationName = null;
        initialSamplingName = null;
        verificationResultName = null;
        modelName = null;
        annotation = null;
    }

    public ExperimentNames(ExperimentNames source) {
        formulaName = source.formulaName;
        initialSpaceName = source.initialSpaceName;
        simulationSpaceName = source.simulationSpaceName;
        precisionConfigurationName = source.precisionConfigurationName;
        initialSamplingName = source.initialSamplingName;
        verificationResultName = source.verificationResultName;
        modelName = source.modelName;
        annotation = source.annotation;
        iterationLimit = source.iterationLimit;
        timeout = source.timeout;
    }

    public String getModelName() {
        return modelName;
    }

    public void setModelName(String modelName) {
        this.modelName = modelName;
    }

    public String getFormulaName() {
        return formulaName;
    }

    public void setFormulaName(String formulaName) {
        this.formulaName = formulaName;
    }

    public String getInitialSamplingName() {
        return initialSamplingName;
    }

    public void setInitialSamplingName(String initialSamplingName) {
        this.initialSamplingName = initialSamplingName;
    }

    public String getInitialSpaceName() {
        return initialSpaceName;
    }

    public void setInitialSpaceName(String initialSpaceName) {
        this.initialSpaceName = initialSpaceName;
    }

    public int getIterationLimit() {
        return iterationLimit;
    }

    public void setIterationLimit(int iterationLimit) {
        this.iterationLimit = iterationLimit;
    }

    public String getPrecisionConfigurationName() {
        return precisionConfigurationName;
    }

    public void setPrecisionConfigurationName(String precisionConfigurationName) {
        this.precisionConfigurationName = precisionConfigurationName;
    }

    public String getSimulationSpaceName() {
        return simulationSpaceName;
    }

    public void setSimulationSpaceName(String simulationSpaceName) {
        this.simulationSpaceName = simulationSpaceName;
    }

    public long getTimeout() {
        return timeout;
    }

    public void setTimeout(long timeout) {
        this.timeout = timeout;
    }

    public String getVerificationResultName() {
        return verificationResultName;
    }

    public void setVerificationResultName(String verificationResultName) {
        this.verificationResultName = verificationResultName;
    }

    public String getAnnotation() {
        return annotation;
    }

    public void setAnnotation(String annotation) {
        this.annotation = annotation;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof ExperimentNames)) {
            return false;
        }
        ExperimentNames target = (ExperimentNames) obj;
        if (!getModelName().equals(target.getModelName())) {
            return false;
        }
        if (!getFormulaName().equals(target.getFormulaName())) {
            return false;
        }
        if (!getInitialSpaceName().equals(target.getInitialSpaceName())) {
            return false;
        }
        if (!getSimulationSpaceName().equals(target.getSimulationSpaceName())) {
            return false;
        }
        if (!getPrecisionConfigurationName().equals(target.getPrecisionConfigurationName())) {
            return false;
        }
        if (!getInitialSamplingName().equals(target.getInitialSamplingName())) {
            return false;
        }
        if (!getVerificationResultName().equals(target.getVerificationResultName())) {
            return false;
        }
        if (getTimeout() != target.getTimeout()) {
            return false;
        }
        if (getIterationLimit() != target.getIterationLimit()) {
            return false;
        }
        if (!getAnnotation().equals(target.getAnnotation())) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        final int prime = 83;
        int result = getModelName().hashCode();
        result = result * prime + getFormulaName().hashCode();
        result = result * prime + getInitialSpaceName().hashCode();
        result = result * prime + getSimulationSpaceName().hashCode();
        result = result * prime + getPrecisionConfigurationName().hashCode();
        result = result * prime + getInitialSamplingName().hashCode();
        result = result * prime + getVerificationResultName().hashCode();
        result = result * prime + Long.valueOf(getTimeout()).hashCode();
        result = result * prime + getIterationLimit();
        result = result * prime + getAnnotation().hashCode();
        return result;
    }

    public static enum Factory {

        INSTANCE;

        private static String getAndValidate(Properties src, String property) throws ResourceException {
            String result = src.getProperty(property);
            if (result == null) {
                throw new ResourceException("`" + property + "' property missing.");
            }
            return result;
        }

        private static String getAndRemoveSuffix(Properties src, String property, ExperimentSuffixes suff) throws ResourceException {
            String fileName = getAndValidate(src, property);
            String result = suff.remove(fileName);
            if (result == null) {
                throw new ResourceException("Wrong file name format: " + fileName + " (missing `" + suff.getSuffix() + "' suffix)");
            }
            return result;
        }

        public ExperimentNames getExperimentNames(Properties properties) throws ResourceException {
            ExperimentNames result = new ExperimentNames();
            result.setModelName(getAndRemoveSuffix(properties, MODEL_PRP, ExperimentSuffixes.MODEL));
            result.setFormulaName(getAndRemoveSuffix(properties, FORMULA_PRP, ExperimentSuffixes.FORMULA));
            result.setInitialSpaceName(getAndRemoveSuffix(properties, INITSPACE_PRP, ExperimentSuffixes.INITIAL_SPACE));
            result.setSimulationSpaceName(getAndRemoveSuffix(properties, SIMSPACE_PRP, ExperimentSuffixes.SIMULATION_SPACE));
            result.setPrecisionConfigurationName(getAndRemoveSuffix(properties, PRECISION_PRP, ExperimentSuffixes.PRECISION_CONFIGURATION));
            result.setInitialSamplingName(getAndRemoveSuffix(properties, SAMPLING_PRP, ExperimentSuffixes.INITIAL_SAMPLING));
            result.setVerificationResultName(getAndRemoveSuffix(properties, RESULT_PRP, ExperimentSuffixes.VERIFICATION_RESULT));
            String property = null;
            try {
                property = getAndValidate(properties, TIMEOUT_PRP);
                result.setTimeout(Long.valueOf(property));
            } catch (NumberFormatException nfe) {
                throw new ResourceException("Wrong number format of timeout: " + property, nfe);
            }
            try {
                property = getAndValidate(properties, ITERATION_PRP);
                result.setIterationLimit(Integer.valueOf(property));
            } catch (NumberFormatException nfe) {
                throw new ResourceException("Wrong number format of iteration limit: " + property, nfe);
            }
            result.setAnnotation(properties.getProperty(ANNOTATION_PRP));
            return result;
        }

        public Properties getProperties(ExperimentNames names) {
            Properties result = new Properties();
            result.setProperty(MODEL_PRP, ExperimentSuffixes.MODEL.add(names.getModelName()));
            result.setProperty(FORMULA_PRP, ExperimentSuffixes.FORMULA.add(names.getFormulaName()));
            result.setProperty(INITSPACE_PRP, ExperimentSuffixes.INITIAL_SPACE.add(names.getInitialSpaceName()));
            result.setProperty(SIMSPACE_PRP, ExperimentSuffixes.SIMULATION_SPACE.add(names.getSimulationSpaceName()));
            result.setProperty(PRECISION_PRP, ExperimentSuffixes.PRECISION_CONFIGURATION.add(names.getPrecisionConfigurationName()));
            result.setProperty(SAMPLING_PRP, ExperimentSuffixes.INITIAL_SAMPLING.add(names.getInitialSamplingName()));
            result.setProperty(RESULT_PRP, ExperimentSuffixes.VERIFICATION_RESULT.add(names.getVerificationResultName()));
            result.setProperty(TIMEOUT_PRP, Long.toString(names.getTimeout()));
            result.setProperty(ITERATION_PRP, Integer.toString(names.getIterationLimit()));
            if (names.getAnnotation() != null && !names.getAnnotation().isEmpty()) {
                result.setProperty(ANNOTATION_PRP, names.getAnnotation());
            }
            return result;
        }

        public static Factory getInstance() {
            return INSTANCE;
        }
    }
    private static final String MODEL_PRP = "sbml.file";
    private static final String FORMULA_PRP = "stl.file";
    private static final String INITSPACE_PRP = "space.initial.file";
    private static final String SIMSPACE_PRP = "space.simulation.file";
    private static final String PRECISION_PRP = "simulation.precision.file";
    private static final String SAMPLING_PRP = "density.sampling.file";
    private static final String RESULT_PRP = "result.output.file";
    private static final String TIMEOUT_PRP = "timeout";
    private static final String ITERATION_PRP = "iteration.limit";
    private static final String ANNOTATION_PRP = "annotation";
}
