/**
 * Copyright (c) 2011-2013 Armin Töpfer
 *
 * This file is part of InDelFixer.
 *
 * InDelFixer is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * InDelFixer is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * InDelFixer. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.indelfixer.stored;

import java.io.Serializable;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Kmer implements Serializable {

    private String sequence;
//    private int start = -1;
//    private int stop = -1;
    private int number;
//    private String alignedSequence;

    /**
     *
     * @return
     */
    public boolean isMapped() {
        return mapped;
    }

    /**
     *
     * @param mapped
     */
    public void setMapped(boolean mapped) {
        this.mapped = mapped;
    }
    private boolean mapped = false;

    /**
     *
     * @param sequence
     * @param number
     */
    public Kmer(String sequence, int number) {
        this.sequence = sequence;
        this.number = number;
    }

    /**
     *
     * @return
     */
    public int getNumber() {
        return number;
    }

    /**
     *
     * @return
     */
    public String getSequence() {
        return sequence;
    }

//    public int getStart() {
//        return start;
//    }
//
//    public int getStop() {
//        return stop;
//    }
//
//    public void setStart(int start) {
//        this.start = start;
//    }
//
//    public void setStop(int stop) {
//        this.stop = stop;
//    }
//
//    public String getAlignedSequence() {
//        return alignedSequence == null ? sequence : alignedSequence;
//    }
//
//    public void setAlignedSequence(String alignedSequence) {
//        this.alignedSequence = alignedSequence;
//    }
}
