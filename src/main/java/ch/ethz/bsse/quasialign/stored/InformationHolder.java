/**
 * Copyright (c) 2011-2012 Armin Töpfer
 *
 * This file is part of QuasiAlign.
 *
 * QuasiAlign is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * QuasiAlign is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * QuasiAlign. If not, see <http://www.gnu.org/licenses/>.
 */

package ch.ethz.bsse.quasialign.stored;

import java.io.Serializable;
import java.util.List;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class InformationHolder implements Serializable {
    private List<Read> reads;
    private int genomeLength;

    public int getGenomeLength() {
        return genomeLength;
    }

    public void setGenomeLength(int genomeLength) {
        this.genomeLength = genomeLength;
    }

    public List<Read> getReads() {
        return reads;
    }

    public void setReads(List<Read> reads) {
        this.reads = reads;
    }

    public InformationHolder(List<Read> reads, int genomeLength) {
        this.reads = reads;
        this.genomeLength = genomeLength;
    }

}