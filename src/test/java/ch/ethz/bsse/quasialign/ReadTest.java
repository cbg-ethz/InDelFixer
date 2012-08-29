/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.InDelFixer;

import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.Read;
import jaligner.matrix.MatrixLoaderException;
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
