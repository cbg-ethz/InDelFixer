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

    public SequenceEntry(String sequence, String tag, int pairedNumber, String quality) {
        this.sequence = sequence;
        this.tag = tag;
        this.pairedNumber = pairedNumber;
        this.quality = quality;
    }

    public SequenceEntry(String sequence) {
        this.sequence = sequence;
    }
}