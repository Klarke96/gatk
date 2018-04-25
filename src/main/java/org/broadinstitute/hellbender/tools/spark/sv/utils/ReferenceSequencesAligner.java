package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;

import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public final class ReferenceSequencesAligner implements AutoCloseable {

    private final Path fasta;
    private final Path image;
    private final BwaMemAligner aligner;
    private final BwaMemIndex index;
    private final List<String> refNames;

    public static final class DescribedRefContig {
        private final String name;
        private final String description;
        private final byte[] bases;

        public DescribedRefContig(final String name, final String description, final byte[] bases) {
            this.name = name;
            this.description = description;
            this.bases = bases;
        }

        public String getName() {
            return name;
        }

        public String getDescription() {
            return description;
        }

        public byte[] getBases() {
            return bases;
        }
    }

    public ReferenceSequencesAligner(final String name, final byte[] bases) {
        try {
            Utils.nonNull(bases);
            refNames = Collections.singletonList( Utils.nonNull(name) );
            fasta = Files.createTempFile("ssvh-temp", ".fasta");
            fasta.toFile().deleteOnExit();
            image = Files.createTempFile(fasta.getParent().toString(), fasta.toString().replace(".fasta", ".img"));
            image.toFile().deleteOnExit();
            FastaReferenceWriter.writeSingleSequenceReference(fasta.toAbsolutePath(), false, false, name, null, bases);
            BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), image.toString());
            index = new BwaMemIndex(image.toString());
            aligner = new BwaMemAligner(index);
        } catch (final IOException ex) {
            throw new GATKException("could not create index files", ex);
        }
    }

    /**
     * Note: when {@code pathToSaveFasta} is {@code null},
     * the fasta is saved in a temporary directory that will be deleted when the alignment step is finished.
     */
    public ReferenceSequencesAligner(final List<String> refNames, final List<DescribedRefContig> refContigs,
                                     final boolean makeIndex, final boolean makeDict,
                                     final int basesPerLine, final String pathToSaveFasta) {

        this.refNames = Utils.nonNull(refNames);

        try {
            Path fastaLocalVal;
            if (pathToSaveFasta == null) {
                fastaLocalVal = Files.createTempFile("ssvh-temp", ".fasta");
            } else {
                try {
                    fastaLocalVal = Files.createFile(IOUtils.getPath(pathToSaveFasta));
                } catch (final FileAlreadyExistsException faeex) {
                    // if fasta already exists, assume outdated and overwrite both fasta and fai (assumed to be "*.fasta" & "*.fasta.fai")
                    Files.delete(IOUtils.getPath(pathToSaveFasta));
                    try {Files.delete(IOUtils.getPath(pathToSaveFasta + ".fai"));} catch (final NoSuchFileException ex) {}
                    fastaLocalVal = Files.createFile(IOUtils.getPath(pathToSaveFasta));
                }
            }
            fasta = fastaLocalVal;
            if (pathToSaveFasta == null) fasta.toFile().deleteOnExit();

            final Path imgPath = IOUtils.getPath(fasta.toString().replace(".fasta", ".img"));
            Path imgLocalVal;
            try {
                imgLocalVal = Files.createFile( imgPath );
            } catch (final FileAlreadyExistsException faeex) {
                // if image file already exists, assume outdated and overwrite
                Files.delete( imgPath );
                imgLocalVal = Files.createFile( imgPath );
            }
            image = imgLocalVal;
            image.toFile().deleteOnExit();
        } catch (final IOException ioex) {
            throw new GATKException("could not create index files", ioex);
        }

        try (final FastaReferenceWriter fastaReferenceWriter =
                     makeIndex ? new FastaReferenceWriter(fasta.toAbsolutePath(), basesPerLine, makeIndex, makeDict)
                               : new FastaReferenceWriter(fasta.toAbsolutePath(), makeIndex, makeDict)
        ) {
            for (final DescribedRefContig contig : refContigs) {
                fastaReferenceWriter.appendSequence(contig.name, contig.description, contig.bases);
            }
            BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), image.toString());
            index = new BwaMemIndex(image.toString());
            aligner = new BwaMemAligner(index);
        } catch (final IOException ioex) {
            throw new GATKException("could not create index files", ioex);
        }
    }

    public final List<List<AlignmentInterval>> align(final List<byte[]> seqs) {
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs);
        final List<List<AlignmentInterval>> result = new ArrayList<>(alignments.size());
        for (int i = 0; i < alignments.size(); i++) {
            final int queryLength = seqs.get(i).length;
            final List<AlignmentInterval> intervals = alignments.get(i).stream()
                    .filter(bwa -> bwa.getRefId() >= 0)
                    .filter(bwa -> SAMFlag.SECONDARY_ALIGNMENT.isUnset(bwa.getSamFlag()))
                    .map(bma -> new AlignmentInterval(bma, refNames, queryLength))
                    .collect(Collectors.toList()); // ignore secondary alignments.
            result.add(intervals);
        }
        return result;
    }

    public final List<List<SAMRecord>> align(final List<SVFastqUtils.FastqRead> reads,
                                             final SAMFileHeader header,
                                             final boolean pairedReadsInterleaved) {

        final int readCount = reads.size();

        final List<String> readNames = new ArrayList<>(readCount);
        final List<byte[]> bases = new ArrayList<>(readCount);
        reads.forEach(fastqRead -> {
            readNames.add(fastqRead.getName());
            bases.add(fastqRead.getBases());
        });

        if (pairedReadsInterleaved) aligner.alignPairs();
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(bases);

        final List<List<SAMRecord>> result = new ArrayList<>(readCount);
        for (int i = 0; i < readNames.size(); ++i) {
            final String readName = readNames.get(i);
            final List<BwaMemAlignment> bwaMemAlignments = alignments.get(i);
            final SAMReadGroupRecord readGroup = getReadGroup(readName);
            final List<SAMRecord> samRecords =
                    BwaMemAlignmentUtils.toSAMStreamForRead(readName, bases.get(i), bwaMemAlignments,
                            header, refNames, readGroup).collect(Collectors.toList());
            result.add(samRecords);
        }
        return result;
    }

    // TODO: 5/8/18 how?
    private static SAMReadGroupRecord getReadGroup(final String readName) {
        return null;
    }

    public SAMSequenceDictionary getDict() {
        return SAMSequenceDictionaryExtractor
                .extractDictionary(IOUtils.getPath(fasta.toString().replace(".fasta", ".dict")));
    }

    @Override
    public void close() throws IOException {
        aligner.close();
        index.close();
    }
}