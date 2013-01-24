/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.indelfixer.stored;

/**
 *
 * @author xlr
 */
public class SequenceEntry {

    public String sequence;
    public String tag;
    public int pairedNumber;
    public String quality;

    /**
     *
     * @param sequence
     * @param tag
     * @param pairedNumber
     * @param quality
     */
    public SequenceEntry(String sequence, String tag, int pairedNumber, String quality) {
        this.sequence = sequence;
        this.tag = tag;
        this.pairedNumber = pairedNumber;
        this.quality = quality;
    }

    /**
     *
     * @param sequence
     */
    public SequenceEntry(String sequence) {
        this.sequence = sequence;
    }

    public SequenceEntry(String sequence, String quality) {
        this.sequence = sequence;
        this.quality = quality;
    }
    
    public SequenceEntry(SimpleRead sr) {
        this.sequence = sr.read;
        this.quality = sr.quality;
    }
}