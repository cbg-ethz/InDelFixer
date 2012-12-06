/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.indelfixer.sffWriter;

import ch.ethz.bsse.indelfixer.stored.Read;
import java.io.File;
import java.util.LinkedList;
import java.util.List;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.BAMRecord;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.DefaultSAMRecordFactory;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

/**
 *
 * @author toepfera
 */
public class SFFWriter {

    /**
     *
     */
    public SFFWriter() {
        SAMFileHeader sfh = new SAMFileHeader();
        List<SAMSequenceRecord> ssrlist = new LinkedList<>();
        SAMSequenceRecord ssr = new SAMSequenceRecord("refi", 20);
        ssrlist.add(ssr);
        SAMSequenceDictionary samSequenceDictionary = new SAMSequenceDictionary(ssrlist);
        sfh.setSequenceDictionary(samSequenceDictionary);
        Read r = new Read(0, 3, "ACGT");
        BAMRecord bamr = new DefaultSAMRecordFactory().createBAMRecord(sfh, 0, r.getBegin(), (short) 0, (short) 0, 0, 1, 0, r.getEnd() - r.getBegin(), 0, 0, 0, new byte[1]);
        List<CigarElement> clist = new LinkedList<>();
        clist.add(new CigarElement(r.getEnd() - r.getBegin(), CigarOperator.M));
        bamr.setCigar(new Cigar(clist));
        byte[] seq = new byte[r.getAlignedRead().length()];
        int i = 0;
        for (char c : r.getAlignedRead().toCharArray()) {
            seq[i++] = (byte) c;
        }
        bamr.setReadBases(seq);
//        bamr.
        final BAMFileWriter outputSam = new BAMFileWriter(new File("/Users/XLR/Desktop/test.bam"));
        outputSam.setHeader(sfh);
        outputSam.addAlignment(bamr);
        outputSam.close();
//        outputSam.
    
    }
}
