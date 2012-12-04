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
package ch.ethz.bsse.indelfixer.parallel;

import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.Kmer;
import ch.ethz.bsse.indelfixer.stored.Read;
import ch.ethz.bsse.indelfixer.utils.Utils;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.javatuples.Pair;
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class CallableWorkerPairedCounter implements Callable<List<Map<Integer, Map<Integer, Integer>>>> {

    private Triplet<String, String, Integer> watsonTriple;
    private Triplet<String, String, Integer> crickTriple;
    private int L;
    private Genome[] genome;
    private Matrix matrix;
    private int number;
    private Map<Integer, Map<Integer, Integer>> substitutionsForward = new HashMap<>();
    private Map<Integer, Map<Integer, Integer>> substitutionsBackward = new HashMap<>();

    public CallableWorkerPairedCounter(Triplet<String, String, Integer> watson, Triplet<String, String, Integer> crick, int L, Genome[] genome, Matrix matrix, int number) {
        this.watsonTriple = watson;
        this.crickTriple = crick;
        this.L = L;
        this.genome = genome;
        this.matrix = matrix;
        this.number = number;
        this.initSubs();
    }

    @Override
    public List<Map<Integer, Map<Integer, Integer>>> call() throws Exception {
        Read watsonRead = map(createRead(watsonTriple, false));
        Read watsonRevRead = map(createRead(watsonTriple, true));
        Read crickRead = map(createRead(crickTriple, false));
        Read crickRevRead = map(createRead(crickTriple, true));

        Read watson = align(watsonRead.getMaximumHits() > watsonRevRead.getMaximumHits() ? watsonRead : watsonRevRead, this.substitutionsForward);
        Read crick = align(crickRead.getMaximumHits() > crickRevRead.getMaximumHits() ? crickRead : crickRevRead, this.substitutionsBackward);



        List<Map<Integer, Map<Integer, Integer>>> subList = new ArrayList<>();
        if (watson != null && watson.isAligned()) {
            subList.add(substitutionsForward);
            if (crick != null && crick.isAligned()) {
                subList.add(substitutionsBackward);
            }
        } else if (crick != null && crick.isAligned()) {
            subList.add(substitutionsBackward);
        }
        return subList;
    }

    private void initSubs() {
        for (int v = 0; v < 6; v++) {
            substitutionsForward.put(v, new HashMap<Integer, Integer>());
            substitutionsBackward.put(v, new HashMap<Integer, Integer>());
            for (int b = 0; b < 6; b++) {
                substitutionsForward.get(v).put(b, 0);
                substitutionsBackward.get(v).put(b, 0);
            }
        }
    }

    private Read createRead(Triplet<String, String, Integer> triple, boolean reverse) {
        Read r;
        if (reverse) {
            r = new Read(Utils.reverseComplement(triple.getValue0()), 0, true);
        } else {
            r = new Read(triple.getValue0(), 0, false);
        }

        r.setDescription(triple.getValue1());
        r.setMatePair(triple.getValue2());

        return r;
    }

    private Read map(Read r) {
        List<Kmer> kmers = r.getKmers();
        for (int x = 0; x < this.genome.length; x++) {
            for (Kmer kmer : kmers) {
                if (this.genome[x].getKmerMap().containsKey(kmer.getSequence())) {
                    for (int pos : this.genome[x].getKmerMap().get(kmer.getSequence())) {
                        r.addHit(pos, pos + Globals.KMER_LENGTH, x);
                    }
                }
            }
        }

        r.findMaxHitRegion();
        int y = r.getBestFittingGenome();
        if (y != -1) {
            if (r.getRegion(y)[0] - r.getRead().length() < 0) {
                r.setBegin(0);
            } else {
                r.setBegin(r.getRegion(y)[0] - r.getRead().length());
            }
            if (r.getRegion(y)[1] + r.getRead().length() > genome[y].getSequence().length()) {
                r.setEnd(genome[y].getSequence().length());
            } else {
                r.setEnd(r.getRegion(y)[1] + r.getRead().length());
            }
        }
        return r;
    }

    private Read align(Read r, Map<Integer, Map<Integer, Integer>> sub) {
        Alignment align;
        if (r.getBestFittingGenome() == -1) {
            return null;
//                    System.out.println(r.getRead());
        }
        if (r.getEnd() < 0) {
            return null;
//            Globals.printPercentageAligningReads();
//                    align = SmithWatermanGotoh.align(
//                            new Sequence(genome, "", "", Sequence.NUCLEIC),
//                            new Sequence(reads[i].getRead(), "", "", Sequence.NUCLEIC),
//                            matrix, 10, 1);
        } else {
            if (r.getBegin() <= 0) {
                align = SmithWatermanGotoh.align(
                        new Sequence(Globals.GENOME_SEQUENCES[r.getBestFittingGenome()], "", "", Sequence.NUCLEIC),
                        new Sequence(r.getRead(), "", "", Sequence.NUCLEIC),
                        matrix, 10, 1);
                r.setBegin(r.getBegin() + 1);
            } else {
                int readEnd = r.getEnd() >= Globals.GENOME_SEQUENCES[r.getBestFittingGenome()].length() ? Globals.GENOME_SEQUENCES[r.getBestFittingGenome()].length() : r.getEnd();
                align = SmithWatermanGotoh.align(
                        new Sequence(Globals.GENOME_SEQUENCES[r.getBestFittingGenome()].substring(r.getBegin() - 1, readEnd), "", "", Sequence.NUCLEIC),
                        new Sequence(r.getRead(), "", "", Sequence.NUCLEIC),
                        matrix, 10, 1);

            }
        }
        StringBuilder sb = new StringBuilder();
        r.setBegin(align.getStart1() + r.getBegin());
        char[] m = align.getMarkupLine();
        char[] c = align.getSequence2();
        char[] g = align.getSequence1();
        int L = m.length;
        if (m.length != c.length) {
            System.err.println("DIFFERENT LENGTHS");
        }
        int ins = 0;
        int del = 0;

        for (int j = 0; j < L; j++) {
            sub.get(convert(c[j])).put(convert(g[j]), sub.get(convert(c[j])).get(convert(g[j])) + 1);
            if (m[j] == '|') {
                sb.append(c[j]);
            } else if (m[j] == ' ') {
                if (isGAP(c[j]) && isGAP(g[j])) {
                    sb.append("-");
                } else if (isGAP(c[j])) {
                    del++;
                    if (Globals.FILL) {
                        sb.append(g[j]);
                    } else {
                        sb.append("-");
                    }
                } else if (isGAP(g[j])) {
                    ins++;
                }
            } else if (m[j] == '.') {
                if (isGAP(g[j])) {
                    sb.append(c[j]);
                } else if (isGAP(c[j])) {
                    sb.append(g[j]);
                    if (c[j] != 'N') {
                        System.err.println("strange");
                        System.err.println(". " + c[j] + " " + m[j] + " " + g[j]);
                    }
                } else {
                    sb.append(c[j]);
                }
            } else {
                System.err.println("E " + c[j] + " " + m[j] + " " + g[j]);
            }
        }
//        if (del < (L * 0.01d) && ins < (L * 0.01d)) 
        r.setAlignedRead(sb.toString());
        r.setEnd(r.getBegin() + sb.length());
        return r;
    }

    private boolean isGAP(char c) {
        return c == '-' || c == 'N';
    }

    private int convert(char c) {
        switch (c) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            case '-':
                return 4;
            default:
                return 5;
        }
    }

    private StringBuilder toString(Read read) {
        StringBuilder sb = new StringBuilder();
        sb.append(">READ").append(this.number).append("_").append(read.getBegin()).append("-").append(read.getEnd());
        sb.append("|").append(read.getDescription()).append("/").append(read.getMatePair());
        sb.append("\n");
        sb.append(read.getAlignedRead()).append("\n");
        return sb;
    }
}
