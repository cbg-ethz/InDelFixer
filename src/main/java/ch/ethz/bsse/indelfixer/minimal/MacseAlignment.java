//package ch.ethz.bsse.indelfixer.minimal;
//
//import align.BasicProfile;
//import align.CodingMSA;
//import align.ElementaryCost;
//import bioObject.CodingDnaSeq;
//import codesInterfaces.MacsEparamCode;
//import java.util.ArrayList;
//import java.util.Hashtable;
//import utils.AlignmentParameterWrappers;
//
//public class MacseAlignment {
//
//    private static AlignmentParameterWrappers createAlignmentParameterWrappers(Hashtable<String, String> sequences2add, String refSequence, String refSequenceName) {
//        float gapE = -1f;
//        float gapO = -7f;
//        float gapC = 0f;
//        float stop = -100f;
//        float fs = -30f;
//        float begEndGapFact = .95f;
//        float optPessFact = 1.0f;
//        float fsLR = -10f;
//        float stopLR = -10f;
//
//        ElementaryCost cost = new ElementaryCost(fs, gapE, gapO, gapC, stop, begEndGapFact, optPessFact);
//        ElementaryCost costLR = new ElementaryCost(fsLR, gapE, gapO, gapC, stopLR, begEndGapFact, optPessFact);
//
//        ArrayList<CodingDnaSeq> sequences = new ArrayList<CodingDnaSeq>();
//        sequences.add(new CodingDnaSeq(refSequenceName, refSequence, true, cost));
//
//        ArrayList<CodingDnaSeq> lessReliableSequences = new ArrayList<CodingDnaSeq>();
//        for (String seqNmae : sequences2add.keySet()) {
//            lessReliableSequences.add(new CodingDnaSeq(seqNmae, sequences2add.get(seqNmae), true, costLR));
//        }
//
//        return new AlignmentParameterWrappers(cost, costLR, sequences, lessReliableSequences, MacsEparamCode.default_GC);
//    }
//
//    public static BasicProfile alignSequences(Hashtable<String, String> sequences2add, String refSequence, String refSequenceName) throws Exception {
//        AlignmentParameterWrappers parameters = createAlignmentParameterWrappers(sequences2add, refSequence, refSequenceName);
//        return CodingMSA.run(parameters);
//    }
//
//    public static void main(String[] args) throws Exception {
//        Hashtable<String, String> seq2add = new Hashtable<String, String>();
//        seq2add.put("seq1", "GATGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGC---AAGAAATGGAGCCAGTAGATCCTAGACTAGAGCC");
//        seq2add.put("seq2", "GATGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGC---AAGAAATGGAGCCAGTAGATCCTAGACTAGAGCC");
//        seq2add.put("seq3", "GATGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGC---AAGAAATGGAGCCAGTAGATCCTAGACTAGAGCC");
//        String seqRef = "AAGAAATGGAGCCAGTAGATCCTAACCTAGAGCC";
//        String seqRefName = "seqRef";
//
//        BasicProfile alignedSeq = alignSequences(seq2add, seqRef, seqRefName);
//
//        //accessing to result
//        //System.out.println(alignedSeq.toString());
//
//        //or by accessing to individual sequence
//        System.out.println("\n\nDNA alignement");
//        for (int i = 0; i < alignedSeq.nbSeq(); i++) {
//            CodingDnaSeq seq = alignedSeq.getSeq(i);
//            System.out.println(">" + seq.getRealFullName() + "\n" + seq.getSeq()); // similar to seq.toFasta()
//        }
//        System.out.println("\n\nAA alignement");
//        for (int i = 0; i < alignedSeq.nbSeq(); i++) {
//            CodingDnaSeq seq = alignedSeq.getSeq(i);
//            System.out.println(">" + seq.getRealFullName() + "\n" + seq.getAAtranslation(1)); // similar to seq.toAAfasta(1)
//        }
//
//
//    }
//}
