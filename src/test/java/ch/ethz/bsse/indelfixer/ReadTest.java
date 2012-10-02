/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.indelfixer;

import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.Read;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import org.biojava3.alignment.NeedlemanWunsch;
import org.junit.*;

/**
 *
 * @author toepfera
 */
public class ReadTest {
    
    public ReadTest() {
        Globals.KMER_LENGTH = 5;
    }

    /**
     * Test of getKmers method, of class Read.
     */
    @Test
    public void testGetKmers() throws MatrixLoaderException {
        String genome = "CCTCTGATCACTCTTTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTGTTAGAAGAAATGAATTTACCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACCCATAGAAGTTTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGGAATCTGTTGACCGAGATTGGTTGCACTTTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCGGGAATGGATGGCCCAAAGGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACTGAAATGGAAAAGGAGGGGAAAATTTCAAAAATTGGACCTGAAAATCCATACAATACTCCAATATTTGCCATAAAGAAAAAAGATGGGACTAAATGGAGAAAATTAGTGGATTTCAGAGAACTTAATAAAAGAACTCAAGACTTCTGGGAAGTTCAGTTAGGAATACCCCATCCCGCAGGGTTTAAAAAGAAAAAATCAATAACAGTACTGGATGTGGGTGATGCGTATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTCACCATACCTAGTATAAACAATGAGACACCAGGAATTAGATATCAGTACAATGTGCTTCCACAAGGACGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAGTCCAGACATAGTTATCTATCAGTACATGGATGACTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAATAAAAATAGAAGAACTGAGGCAGCATCTGTTGCAGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGAGCTGCCAGAAAAAGATAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGAAAATTAAATTGGGCAAGTCAAATTTATCCAGGGATTAAAGTAAGGCAATTATGTAAGCCCCTTAGGGGAACCAAAGCACTAACAGAAATAATACCACTAACAGAAGAAGCGGAGCTAGAACTAGCAGAAAACAGGGAGATTTTAAAAGAATCAGTACATGGA";
        Alignment align = SmithWatermanGotoh.align(
                                                  new Sequence(genome.substring(0,1000), "", "", Sequence.NUCLEIC),
                                                  new Sequence("TGATGCGTATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTCACCATACCTAGTATAAACAATGAGACACCAGGAATTAGATATCAGTACAATGTGCTTCCACAAGGACGGAAAGGATCACCAGCAATATTCCAAAG", "", "", Sequence.NUCLEIC),
//                                                  new Sequence("ACTCTTTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTGTTAGAAGAAATGAATTTACCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTT", "", "", Sequence.NUCLEIC),
                                                  MatrixLoader.load("EDNAFULL"), 10, 1);
        System.out.println(align.getStart1());
        //        System.out.println("split reads to kmers");
        //        Read instance = new Read("ACGT");
        //        String[] expResult = new String[]{"AC","CG","GT"};
        //        String[] result = instance.getKmers();
        //        System.out.println(Arrays.toString(result));
        //        Alignment align = SmithWatermanGotoh.align(
        //                new Sequence("ACTT", "", "", Sequence.NUCLEIC),
        //                            new Sequence("CGACCGACGTTACTAAA", "", "", Sequence.NUCLEIC),
        //
        //                            MatrixLoader.load("EDNAFULL"), 10, 4);
        //        System.out.println("");
        //        Read r = new Read("ACTGGACTCATATACTTCT");
    }
}
