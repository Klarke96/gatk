package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadEnds;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

import java.util.Map;

/**
 * Class representing a single read fragment at a particular start location without a mapped mate.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
public class Fragment extends PairedEnds {
    protected transient ReadsKey key;

    private final int firstStartPosition;
    private final boolean R1R;

    protected final int score;

    public Fragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap) {
        super(partitionIndex, first.getName());

        this.firstStartPosition = first.getAssignedStart();
        this.score = scoringStrategy.score(first);
        this.R1R = first.isReverseStrand();
        this.key = ReadsKey.getKeyForFragment(ReadUtils.getStrandedUnclippedStart(first),
                isRead1ReverseStrand(),
                (short)ReadUtils.getReferenceIndex(first, header),
                headerLibraryMap.get(ReadUtils.getLibrary(first, header, LibraryIdGenerator.UNKNOWN_LIBRARY)));
    }

    @Override
    public Type getType() {
      return Type.FRAGMENT;
    }
    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public ReadsKey key() {
        return key;
    }
    @Override
    public int getScore() {
      return score;
    }
    @Override
    public int getFirstStartPosition() {
      return firstStartPosition;
    }
    @Override
    public boolean isRead1ReverseStrand() {
      return R1R;
    }
    @Override
    public byte getOrientationForPCRDuplicates() {
        return (R1R)? ReadEnds.R : ReadEnds.F;
    }
    @Override
    public String toString() {
        return "fragment: " + name + " " + firstStartPosition;
    }
}
