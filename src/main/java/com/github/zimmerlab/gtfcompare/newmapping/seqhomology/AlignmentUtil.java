package com.github.zimmerlab.gtfcompare.newmapping.seqhomology;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;

public class AlignmentUtil implements SeqHomologyUtil {

    private final static Logger logger = LogManager.getLogger(AlignmentUtil.class);

    public double calculate(String targetSeq, String querySeq){
        try{
            var targetDnaSeq = new DNASequence(targetSeq);
            var queryDnaSeq = new DNASequence(querySeq);

            var matrix = SubstitutionMatrixHelper.getNuc4_4();

            var gapPenalty = new SimpleGapPenalty();
            gapPenalty.setOpenPenalty(10);
            gapPenalty.setExtensionPenalty(1);

            logger.debug("Target Length: {}, Query Length: {}", targetDnaSeq.getLength(), queryDnaSeq.getLength());

            var pair = Alignments.getPairwiseAlignment(
                    targetDnaSeq,
                    queryDnaSeq,
                    Alignments.PairwiseSequenceAlignerType.GLOBAL,
                    gapPenalty,
                    matrix
            );

            return pair.getPercentageOfIdentity(true);
        } catch (CompoundNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
}
