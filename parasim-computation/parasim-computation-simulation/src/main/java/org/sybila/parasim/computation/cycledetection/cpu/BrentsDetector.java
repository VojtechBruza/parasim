package org.sybila.parasim.computation.cycledetection.cpu;

import org.sybila.parasim.computation.cycledetection.Detector;
import org.sybila.parasim.computation.cycledetection.CycleDetectionStatus;
import org.sybila.parasim.computation.cycledetection.ArrayCycleDetectionDataBlock;
import org.sybila.parasim.model.trajectory.DataBlock;
import org.sybila.parasim.model.trajectory.ListMutableDataBlock;
import org.sybila.parasim.model.trajectory.Trajectory;
import org.sybila.parasim.model.trajectory.PointComparator;
import java.util.Iterator;
import java.util.Arrays;

/**
 * Detects cycles on data blocks of trajectories using Brent's cycle detector.
 *
 * @author <a href="mailto:sven@mail.muni.cz">Sven Dra�an</a>
 */
public class BrentsDetector implements Detector<PointComparator, BrentsCycleDetector, ArrayCycleDetectionDataBlock>
{
    /**
     * Detects cycles on given trajectories using the comparator with a maximum
     * of stepLimit points processed from every trajectory.
     *
     * @param comparator Point comparator to for similarity test between points.
     * @param trajectories Set of trajectories on which to detect cycles.
     * @param stepLimit Maximum number of points to process from each trajectory.
     * @return Cycle detection data block with results of the computation.
     */
    @Override
    public ArrayCycleDetectionDataBlock detect(PointComparator comparator,
                   DataBlock<Trajectory> trajectories, int stepLimit)
    {
        if (trajectories == null)
        {
            throw new IllegalArgumentException("The parameter trajectories is null.");
        }
        if (trajectories.size() <= 0)
        {
            throw new IllegalArgumentException("Trajectories must contain at least one trajectory.");
        }
        BrentsCycleDetector[] cycleDetectors = new BrentsCycleDetector[trajectories.size()];
        for (int i=0; i<trajectories.size(); i++)
        {
            cycleDetectors[i] = new BrentsCycleDetector(comparator);
        }
        return this.detect(comparator, trajectories, cycleDetectors, stepLimit);
    }

    /**
     * Detects cycles on given trajectories using the comparator with a maximum
     * of stepLimit points processed from every trajectory. Cycle detection
     * carried out in previous iterations is conserved inside the detectors array.
     *
     * @param comparator Point comparator for similarity tests between points.
     * @param trajectories Set of trajectories on which to detect cycles.
     * @param detectors Brent's cycle detectors used in previous iterations on
     *                  the given trajectories.
     * @param stepLimit Maximum number of points to process from each trajectory.
     * @return Cycle detection data block with results of the computation.
     */
    @Override
    public ArrayCycleDetectionDataBlock detect(PointComparator comparator,
                   DataBlock<Trajectory> trajectories,
                   BrentsCycleDetector[] detectors,
                   int stepLimit)
    {
        if (comparator == null)
        {
            throw new IllegalArgumentException("The parameter comparator is null.");
        }
        if (trajectories == null)
        {
            throw new IllegalArgumentException("The parameter trajectories is null.");
        }
        if (trajectories.size() <= 0)
        {
            throw new IllegalArgumentException("Trajectories must contain at least one trajectory.");
        }
        if (detectors == null)
        {
            throw new IllegalArgumentException("The parameter detecotrs is null.");
        }
        if (trajectories.size() != detectors.length)
        {
            throw new IllegalArgumentException("Number of trajectories must match number of comparators.");
        }
        if (stepLimit < 0)
        {
            throw new IllegalArgumentException("The parameter stepLimit must be none negative.");
        }       
        CycleDetectionStatus[] statuses = new CycleDetectionStatus[trajectories.size()];
        Arrays.fill(statuses, CycleDetectionStatus.COMPUTING);

        Iterator<Trajectory> tIterator = trajectories.iterator();
        int index = 0;
        while (tIterator.hasNext())
        {
            Trajectory t = tIterator.next();
            BrentsCycleDetector detector = detectors[index];
            int stepsUsed = detector.detectCycle(t, stepLimit);
            if (detector.cycleDetected())
            {
                statuses[index] = CycleDetectionStatus.CYCLE;
            }
            else if (stepLimit == 0 || stepLimit > stepsUsed)
            {
                statuses[index] = CycleDetectionStatus.NOCYCLE;
            }
            index++;
        }
        return new ArrayCycleDetectionDataBlock(trajectories, detectors, statuses);
    }

    /**
     * Detects cycles on two blocks of trajectories merging the results.
     * The first set given by oldTrajectories is expected to have already
     * some cycle detection performed on them and the array oldDetectors belongs
     * to these. The second set given by newTrajectories is expected to have
     * no cycle detection yet performed.
     *
     * @param comparator Point comparator for similarity tests between points.
     * @param oldTrajectories Set of partialy processed trajectories.
     * @param oldDetectors Array of cycle detectors of all oldTrajectories.
     * @param newTrajectories Set of new trajectories on which to detect cycles.
     * @param stepLimit Maximum number of points to process from each trajectory.
     * @return Cycle detection data block with results of the computation for both
     *         DataBlocks, old first then new ones.
     */
    public ArrayCycleDetectionDataBlock detect(PointComparator comparator,
                                               DataBlock<Trajectory> oldTrajectories,
                                               BrentsCycleDetector[] oldDetectors,
                                               DataBlock<Trajectory> newTrajectories,
                                               int stepLimit)
    {
        if (oldTrajectories == null)
        {
            throw new IllegalArgumentException("The parameter oldTrajectories is null.");
        }
        if (newTrajectories == null || newTrajectories.size() == 0)
        {
            return detect(comparator, oldTrajectories, oldDetectors, stepLimit);
        }

        BrentsCycleDetector[] cycleDetectors = new BrentsCycleDetector[oldDetectors.length + newTrajectories.size()];
        System.arraycopy(oldDetectors, 0, cycleDetectors, 0, oldDetectors.length);
        for (int i=oldDetectors.length; i<cycleDetectors.length; i++)
        {
            cycleDetectors[i] = new BrentsCycleDetector(comparator);
        }

        ListMutableDataBlock<Trajectory> allTrajectories = new ListMutableDataBlock(oldTrajectories);
        allTrajectories.merge(newTrajectories);

        return this.detect(comparator, allTrajectories, cycleDetectors, stepLimit);        
    }

}
