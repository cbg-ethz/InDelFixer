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
    public String header;
    public String quality;
    public int startOffset;
    public int endOffset;

    public SequenceEntry(String sequence, String header, String quality, int startOffset, int endOffset) {
        this.sequence = sequence;
        this.header = header;
        this.quality = quality;
        this.startOffset = startOffset;
        this.endOffset = endOffset;
    }

    public SequenceEntry(String sequence, String header, String quality) {
        this.sequence = sequence;
        this.header = header;
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