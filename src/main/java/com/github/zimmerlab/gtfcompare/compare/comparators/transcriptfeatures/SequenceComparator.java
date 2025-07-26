package com.github.zimmerlab.gtfcompare.compare.comparators.transcriptfeatures;

import com.github.zimmerlab.gtfcompare.compare.CDSComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SequenceComparator implements ComparisonFeature, CDSComparisonFeature {

    private final static Logger LOG = LogManager.getLogger(SequenceComparator.class);
    @Override
    public String getName() {
        return Constants.SEQUENCE_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targetSequenceExtractor = ctx.getTargetExtractor();
        var querySequenceExtractor = ctx.getQueryExtractor();

        if(targetSequenceExtractor.isEmpty() || querySequenceExtractor.isEmpty()) {
            LOG.warn("One of the sequence extractors is empty. Cannot compare sequences.");
            return false;
        }

        var targetBaseData = ctx.getTargetFeature().getBaseData();
        var queryBaseData = ctx.getQueryFeature().getBaseData();
        try{
            var targetSequence = targetSequenceExtractor.get().getSequence(targetBaseData.getContig(), targetBaseData.getStart(), targetBaseData.getEnd());
            var querySequence = targetSequenceExtractor.get().getSequence(queryBaseData.getContig(), queryBaseData.getStart(), queryBaseData.getEnd());
            return !targetSequence.equals(querySequence);
        } catch (Exception e){
            LOG.error("Error while extracting sequences: {}", e.getMessage());
            return false;
        }


    }
}
