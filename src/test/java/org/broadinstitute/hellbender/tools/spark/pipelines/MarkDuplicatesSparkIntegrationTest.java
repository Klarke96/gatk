package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.SparkException;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.markduplicates.MarkDuplicatesGATKIntegrationTest;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.GATKDuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesSparkTester;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.markduplicates.AbstractMarkDuplicatesCommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.testers.AbstractMarkDuplicatesTester;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@Test(groups = "spark")
public class MarkDuplicatesSparkIntegrationTest extends AbstractMarkDuplicatesCommandLineProgramTest {

    @Override
    protected AbstractMarkDuplicatesTester getTester() {
        MarkDuplicatesSparkTester markDuplicatesSparkTester = new MarkDuplicatesSparkTester();
        markDuplicatesSparkTester.addArg("--"+ MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME);
        return markDuplicatesSparkTester;
    }

    @Override
    protected CommandLineProgram getCommandLineProgramInstance() {
        return new MarkDuplicatesSpark();
    }

    @Override
    protected boolean markSecondaryAndSupplementaryRecordsLikeTheCanonical() { return true; }

    @Test(dataProvider = "testMDdata", groups = "spark")
    @Override
    public void testMDOrder(final File input, final File expectedOutput) throws Exception {
        // Override this test case to provide a --sharded-output false argument, so that we write a single, sorted
        // bam (since sharded output is not sorted, and this test case is sensitive to order).
        testMDOrderImpl(input, expectedOutput, "--" + GATKSparkTool.SHARDED_OUTPUT_LONG_NAME +" false");
    }

    @DataProvider(name = "md")
    public Object[][] md(){
        return new Object[][]{
            // The first two values are total reads and duplicate reads. The list is an encoding of the metrics
            // file output by this bam file. These metrics files all match the outputs of picard mark duplicates.

             //Note: in each of those cases, we'd really want to pass null as the last parameter (not 0L) but IntelliJ
             // does not like it and skips the test (rendering issue) - so we pass 0L and account for it at test time
             // (see comment in testMarkDuplicatesSparkIntegrationTestLocal)
            {new File(MarkDuplicatesGATKIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam"), 20, 0,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16406", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L))},
            {new File(MarkDuplicatesGATKIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.bam"), 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File(MarkDuplicatesGATKIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.markedDups.bam"), 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File(MarkDuplicatesGATKIntegrationTest.TEST_DATA_DIR, "optical_dupes.bam"), 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
            {new File(MarkDuplicatesGATKIntegrationTest.TEST_DATA_DIR, "optical_dupes_casava.bam"), 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
        };
    }

    @Test( dataProvider = "md")
    public void testMarkDuplicatesSparkIntegrationTestLocal(
        final File input, final long totalExpected, final long dupsExpected,
        Map<String, List<String>> metricsExpected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(input.getPath());
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.add(outputFile.getAbsolutePath());

        args.add("--"+StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.add(metricsFile.getAbsolutePath());

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        int totalReads = 0;
        int duplicateReads = 0;
        try ( final ReadsDataSource outputReads = new ReadsDataSource(outputFile.toPath()) ) {
            for ( GATKRead read : outputReads ) {
                ++totalReads;

                if ( read.isDuplicate() ) {
                    ++duplicateReads;
                }
            }
        }

        Assert.assertEquals(totalReads, totalExpected, "Wrong number of reads in output BAM");
        Assert.assertEquals(duplicateReads, dupsExpected, "Wrong number of duplicate reads in output BAM");

        final MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<GATKDuplicationMetrics> nonEmptyMetrics = metricsOutput.getMetrics().stream().filter(
                metric ->
                    metric.UNPAIRED_READS_EXAMINED != 0L ||
                    metric.READ_PAIRS_EXAMINED != 0L ||
                    metric.UNMAPPED_READS != 0L ||
                    metric.UNPAIRED_READ_DUPLICATES != 0L ||
                    metric.READ_PAIR_DUPLICATES != 0L ||
                    metric.READ_PAIR_OPTICAL_DUPLICATES != 0L ||
                    (metric.PERCENT_DUPLICATION != null && metric.PERCENT_DUPLICATION != 0.0 && !Double.isNaN(metric.PERCENT_DUPLICATION)) ||
                    (metric.ESTIMATED_LIBRARY_SIZE != null && metric.ESTIMATED_LIBRARY_SIZE != 0L)
        ).collect(Collectors.toList());

        Assert.assertEquals(nonEmptyMetrics.size(), metricsExpected.size(),
                            "Wrong number of metrics with non-zero fields.");
        for (int i = 0; i < nonEmptyMetrics.size(); i++ ){
            final GATKDuplicationMetrics observedMetrics = nonEmptyMetrics.get(i);
            List<?> expectedList = metricsExpected.get(observedMetrics.LIBRARY);
            Assert.assertNotNull(expectedList, "Unexpected library found: " + observedMetrics.LIBRARY);
            Assert.assertEquals(observedMetrics.UNPAIRED_READS_EXAMINED, expectedList.get(0));
            Assert.assertEquals(observedMetrics.READ_PAIRS_EXAMINED, expectedList.get(1));
            Assert.assertEquals(observedMetrics.UNMAPPED_READS, expectedList.get(2));
            Assert.assertEquals(observedMetrics.UNPAIRED_READ_DUPLICATES, expectedList.get(3));
            Assert.assertEquals(observedMetrics.READ_PAIR_DUPLICATES, expectedList.get(4));
            Assert.assertEquals(observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, expectedList.get(5));
            Assert.assertEquals(observedMetrics.PERCENT_DUPLICATION, expectedList.get(6));

            //Note: IntelliJ does not like it when a parameter for a test is null (can't print it and skips the test)
            //so we work around it by passing in an 'expected 0L' and only comparing to it if the actual value is non-null
            if (observedMetrics.ESTIMATED_LIBRARY_SIZE != null && (Long)expectedList.get(7) != 0L)  {
                Assert.assertEquals(observedMetrics.ESTIMATED_LIBRARY_SIZE, expectedList.get(7));
            }
        }
    }

    @Test
    public void testSupplementaryReadUnmappedMate() {
        File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(getTestFile("supplementaryReadUnmappedmate.bam"));
        runCommandLine(args);

        try ( final ReadsDataSource outputReadsSource = new ReadsDataSource(output.toPath()) ) {
            final List<GATKRead> actualReads = new ArrayList<>();
            for ( final GATKRead read : outputReadsSource ) {
                Assert.assertFalse(read.isDuplicate());
                actualReads.add(read);
            }

            Assert.assertEquals(actualReads.size(), 3, "Wrong number of reads output");
        }
    }

    @Test
    public void testHashCollisionHandling() {
        // This test asserts that the handling of two read pairs with the same start positions but on different in such a way
        // that they might cause hash collisions are handled properly.
        final File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(getTestFile("hashCollisionedReads.bam"));
        runCommandLine(args);

        try ( final ReadsDataSource outputReadsSource = new ReadsDataSource(output.toPath()) ) {
            final List<GATKRead> actualReads = new ArrayList<>();
            for ( final GATKRead read : outputReadsSource ) {
                Assert.assertFalse(read.isDuplicate());
                actualReads.add(read);
            }

            Assert.assertEquals(actualReads.size(), 4, "Wrong number of reads output");
        }
    }

    // Tests asserting that without --do-not-mark-unmapped-mates argument that unmapped mates are still duplicate marked with their partner
    @Test
    public void testMappedPairAndMappedFragmentAndMatePairSecondUnmapped() {
        final AbstractMarkDuplicatesTester tester = new MarkDuplicatesSparkTester(true);
        tester.addMatePair(1, 10040, 10040, false, true, true, true, "76M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, ELIGIBLE_BASE_QUALITY); // mapped OK
        tester.addMappedFragment(1, 10040, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testNonExistantReadGroupInRead() {
        final AbstractMarkDuplicatesTester tester = new MarkDuplicatesSparkTester(true);
        tester.addMatePair("RUNID:7:1203:2886:82292",  19, 19, 485253, 485253, false, false, true, true, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, "NotADuplicateGroup");
        try {
            tester.runTest();
            Assert.fail("Should have thrown an exception");
        } catch (Exception e){
           Assert.assertTrue(e instanceof SparkException);
           Assert.assertTrue(e.getCause() instanceof UserException.HeaderMissingReadGroup);
        }
    }

    @Test
    public void testNoReadGroupInRead() {
        final AbstractMarkDuplicatesTester tester = new MarkDuplicatesSparkTester(true);
        tester.addMatePair("RUNID:7:1203:2886:82292",  19, 19, 485253, 485253, false, false, true, true, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, null);

        try {
            tester.runTest();
            Assert.fail("Should have thrown an exception");
        } catch (Exception e){
            Assert.assertTrue(e instanceof SparkException);
            Assert.assertTrue(e.getCause() instanceof UserException.ReadMissingReadGroup);
        }
    }
}
