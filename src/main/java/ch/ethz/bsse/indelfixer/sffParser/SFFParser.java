package ch.ethz.bsse.indelfixer.sffParser;

import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This class keeps track of one or multiple sff files from 454 sequencing
 * machines. Always call this class to initiate parsing.
 *
 * TODO: create an interator, in whole read pool. how to realize that? make use
 * of mft index structur ArrayList for addtional information not best solution,
 * better way?
 *
 *
 * @author Charles Imbusch < charles at imbusch dot com >
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SFFParser {

    private SFFFile sffFile = null;

    /**
     * Constructor for one single sff file
     *
     * @param filename
     * @throws FileNotFoundException
     * @throws IOException
     */
    public SFFParser(String filename) throws FileNotFoundException, IOException {
        sffFile = new SFFFile(filename);
    }

    /**
     * Returns a DNASequence object. Implies that read identifier is unique in
     * pool of sff files specified.
     *
     * @param identifier
     * @return DNASequence object.
     * @throws IOException
     */
    public List<SFFRead> getReads() {
        List<SFFRead> reads = new ArrayList<>();
        int i = 0;
        for (String s : this.sffFile.getReadNames()) {
            SFFRead r = null;
            try {
                r = sffFile.getRead(s);
            } catch (IOException ex) {
                Logger.getLogger(SFFParser.class.getName()).log(Level.SEVERE, null, ex);
            }
            Globals.printPercentageExtractingReads();
            if (r.getRead().length() >= 150) {
                reads.add(r);
            }
        }
        
        StatusUpdate.println("Extracting reads:\t"+reads.size());
        return reads;
    }

    /**
     * @return Number of reads in readpool.
     */
    public int getNumberOfReads() {
        return sffFile.getNumberOfReads();
    }
}
