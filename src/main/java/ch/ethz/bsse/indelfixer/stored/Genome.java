/**
 * Copyright (c) 2011-2013 Armin Töpfer
 *
 * This file is part of InDelFixer.
 *
 * InDelFixer is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * InDelFixer is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * InDelFixer. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.indelfixer.stored;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Genome {

    private String sequence;
    private String header;
    private Map<String, List<Integer>> kmerMap = new HashMap<>();

    /**
     *
     * @param sequence
     * @throws IllegalStateException
     */
    public Genome(String sequence) throws IllegalStateException {
        this.sequence = sequence;
        this.split();
    }

    /**
     *
     * @param hap
     * @throws IllegalStateException
     */
    public Genome(Map.Entry<String, String> hap) throws IllegalStateException {
        this.sequence = hap.getKey();
        this.header = hap.getValue().replaceAll(">", "").replace(" ", "");
    }

    public void split() throws IllegalStateException {
        for (int i = Globals.KMER_LENGTH; i <= this.sequence.length(); i++) {
            String tmp = this.sequence.substring(i - Globals.KMER_LENGTH, i);
            if (kmerMap.containsKey(tmp)) {
                kmerMap.get(tmp).add(i - Globals.KMER_LENGTH + 1);
            } else {
                List<Integer> list = new ArrayList<>();
                list.add(i - Globals.KMER_LENGTH + 1);
                kmerMap.put(tmp, list);
            }
        }
    }

    /**
     *
     * @return
     */
    public Map<String, List<Integer>> getKmerMap() {
        return Collections.unmodifiableMap(kmerMap);
    }

    /**
     *
     * @return
     */
    public String getSequence() {
        return sequence;
    }

    /**
     *
     * @return
     */
    public String getHeader() {
        return header;
    }
}
