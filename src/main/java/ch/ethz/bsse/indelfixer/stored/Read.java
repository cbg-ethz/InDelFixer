/**
 * Copyright (c) 2011-2012 Armin Töpfer
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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Read implements Serializable {

    private String read;
    private List<Kmer> kmers;
    private Hits[] hits;
    private int begin;
    private int end;
    private String alignedRead;
    private int number;
    private boolean reverse;
//    private int[] quality;
    private String description;
    private int matePair = -1;
    private int bestGenomeIndex = -1;

    public Read(int begin, int end, String alignedRead) {
        this.begin = begin;
        this.end = end;
        this.alignedRead = alignedRead;
        System.out.println("#########");
//        this.quality = quality;
//        this.kmers = new ArrayList<>();
//        this.split();
    }

    public Read(String read, int number, boolean reverse) {
        this.read = read;
        this.number = number;
        this.reverse = reverse;
//        this.quality = quality;
        this.kmers = new ArrayList<>();
        this.init();
        this.split();
    }

    private void init() {
        this.hits = new Hits[Globals.GENOME_COUNT];
        for (int i = 0; i < Globals.GENOME_COUNT; i++) {
            this.hits[i] = new Hits();
        }
    }

    private void split() {
        int location = 0;
        for (int i = Globals.KMER_LENGTH; i <= read.length(); i += Globals.KMER_OVERLAP) {
            this.kmers.add(new Kmer(this.read.substring(i - Globals.KMER_LENGTH, i), location++));
        }
    }

    public void findMaxHitRegion() {
        for (int x = 0; x < hits.length; x++) {
            Hits h = hits[x];
            int currentHits = 0;
            for (int i = 0; i < Globals.GENOME_SEQUENCES[x].length(); i++) {
                currentHits += getHit(i, x);
                if (i >= (h.min + read.length())) {
                    if (currentHits > h.maximumHits) {
                        h.region[0] = i - read.length();
                        h.region[1] = i;
                        h.maximumHits = currentHits;
                    }
                    currentHits -= getHit(i - read.length(), x);
                }
            }
            if (currentHits > h.maximumHits) {
                h.region[0] = Globals.GENOME_SEQUENCES[x].length() - read.length();
                h.region[1] = Globals.GENOME_SEQUENCES[x].length();
                h.maximumHits = currentHits;
            }
        }
    }

    public void addHit(int from, int to, int genome) {
        Hits h = hits[genome];
        h.min = Math.min(from, h.min);
        h.max = Math.max(to, h.max);
        for (int i = from; i < to; i++) {
            if (h.hitMap.containsKey(i)) {
                h.hitMap.put(i, h.hitMap.get(i) + i);
            } else {
                h.hitMap.put(i, 1);
            }
        }
    }

    public int getHit(int i, int genome) {
        Hits h = hits[genome];
        if (h.hitMap.containsKey(i)) {
            return h.hitMap.get(i);
        } else {
            return 0;
        }
    }

    public List<Kmer> getKmers() {
        return kmers;
    }

    public String getRead() {
        return read;
    }

    public int getBegin() {
        return begin;
    }

    public void setBegin(int begin) {
        this.begin = begin;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getBestFittingGenome() {
        if (this.bestGenomeIndex == -1) {
            int maximumHits = 0;
            for (int i = 0; i < hits.length; i++) {
                Hits h = hits[i];
                if (maximumHits < h.maximumHits) {
                    maximumHits = h.maximumHits;
                    bestGenomeIndex = i;
                }
            }
        }
//        if (this.bestGenomeIndex == -1) {
//            System.out.println("");
//        }
        return bestGenomeIndex;
    }

    public int getMaximumHits() {
        if (this.getBestFittingGenome() == -1) {
            return -1;
        }
        return hits[this.getBestFittingGenome()].maximumHits;
    }

    public int[] getRegion() {
        return hits[this.getBestFittingGenome()].region;
    }

    public int[] getRegion(int genome) {
        return hits[genome].region;
    }

    public String getAlignedRead() {
        return alignedRead;
    }

    public boolean isAligned() {
        return alignedRead != null;
    }

    public void setAlignedRead(String alignedRead) {
        this.alignedRead = alignedRead;
    }

    public boolean isReverse() {
        return reverse;
    }

    public int getNumber() {
        return number;
    }

//    public int[] getQuality() {
//        return quality;
//    }
//
//    public void setQuality(int[] quality) {
//        this.quality = quality;
//    }
    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public int getMatePair() {
        return matePair;
    }

    public void setMatePair(int matePair) {
        this.matePair = matePair;
    }
}

class Hits {

    Map<Integer, Integer> hitMap = new HashMap<>();
    int min = Integer.MAX_VALUE;
    int max = Integer.MIN_VALUE;
    //from,to
    int[] region = new int[2];
    int maximumHits = Integer.MIN_VALUE;
}