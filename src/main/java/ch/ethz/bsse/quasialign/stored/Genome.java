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

import java.util.HashMap;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Genome {

    private String sequence;
    private Map<String,Integer> kmerMap = new HashMap<>();
    
    public Genome(String sequence) throws IllegalStateException {
        this.sequence = sequence;
        this.split();
    }

    private void split() throws IllegalStateException {
        for (int i = Globals.KMER_LENGTH; i <= this.sequence.length(); i++) {
            String tmp = this.sequence.substring(i-Globals.KMER_LENGTH, i);
            if (kmerMap.containsKey(tmp)) {
                throw new IllegalStateException("Non unique k-mer");
            } else {
                kmerMap.put(tmp, i-Globals.KMER_LENGTH+1);
            }
        }
    }

    public Map<String, Integer> getKmerMap() {
        return kmerMap;
    }

    public String getSequence() {
        return sequence;
    }
    
}
