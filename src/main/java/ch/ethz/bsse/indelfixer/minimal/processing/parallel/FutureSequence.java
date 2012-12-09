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
package ch.ethz.bsse.indelfixer.minimal.processing.parallel;

import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.Kmer;
import ch.ethz.bsse.indelfixer.stored.Read;
import ch.ethz.bsse.indelfixer.stored.SequenceEntry;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import ch.ethz.bsse.indelfixer.utils.Utils;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class FutureSequence implements Callable<Pair<String, Map<Integer, Map<Integer, Integer>>>> {

    private SequenceEntry watsonTriple;
    private Genome[] genome;
    private Matrix matrix;
    private int number;
    private Map<Integer, Map<Integer, Integer>> substitutionsForward = new HashMap<>();

    public FutureSequence(SequenceEntry watson, int number) {
        this.watsonTriple = watson;
        this.genome = Globals.GENOMES;
        this.matrix = Globals.MATRIX;
        this.number = number;
        this.initSubs();
    }

    @Override
    public Pair<String, Map<Integer, Map<Integer, Integer>>> call() {
        StringBuilder sb = new StringBuilder();
        if (this.watsonTriple != null) {
            Read watsonRead = map(createRead(watsonTriple, false));
            Read watsonRevRead = map(createRead(watsonTriple, true));
            try {
                Read align = align(watsonRead.getMaximumHits() > watsonRevRead.getMaximumHits() ? watsonRead : watsonRevRead, this.substitutionsForward);
                Read watson = cut(align);
                if (watson != null) {
                    sb.append(toString(watson));
                    StatusUpdate.processReads();
                    return Pair.with(sb.toString(), substitutionsForward);
                }
            } catch (Exception e) {
                System.err.println(e);
            }
        }
        StatusUpdate.processReads();
        return null;
    }

    private void initSubs() {
        for (int v = 0; v < 6; v++) {
            substitutionsForward.put(v, new HashMap<Integer, Integer>());
            for (int b = 0; b < 6; b++) {
                substitutionsForward.get(v).put(b, 0);
            }
        }
    }

    private Read createRead(SequenceEntry entry, boolean reverse) {
        Read r;
        if (reverse) {
            r = new Read(Utils.reverseComplement(entry.sequence), 0, true);
        } else {
            r = new Read(entry.sequence, 0, false);
        }
        if (entry.tag != null) {
            r.setDescription(entry.tag);
        }
        if (entry.pairedNumber != 0) {
            r.setMatePair(entry.pairedNumber);
        }
        if (entry.quality != null) {
            r.setQuality(entry.quality);
        }
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
        double L = m.length;
        if (m.length != c.length) {
        }
        int ins = 0;
        int del = 0;
        int subs = 0;

        for (int j = 0; j < L; j++) {
            char currentConsensus = '*';
            if (m[j] == '|') {
                currentConsensus = c[j];
            } else if (m[j] == ' ') {
                if (isGAP(c[j]) && isGAP(g[j])) {
                    currentConsensus = '-';
                } else if (isGAP(c[j])) {
                    if (c[j] != 'N') {
                        del++;
                    }
                    if (Globals.FILL) {
                        currentConsensus = g[j];
                    } else {
                        currentConsensus = '-';
                    }
                } else if (isGAP(g[j])) {
                    ins++;
                }
            } else if (m[j] == '.') {
                subs++;
                if (isGAP(g[j])) {
                    currentConsensus = c[j];
                } else if (isGAP(c[j])) {
                    currentConsensus = g[j];
                    if (c[j] != 'N') {
                        System.err.println("strange");
                    }
                } else {
                    currentConsensus = c[j];
                }
            }
            if (currentConsensus != '*') {
                sb.append(currentConsensus);
                sub.get(convert(currentConsensus)).put(convert(g[j]), sub.get(convert(currentConsensus)).get(convert(g[j])) + 1);
            }
        }
        if (((del / L) > Globals.MAX_DEL) || ((ins / L) > Globals.MAX_INS) || ((subs / L) > Globals.MAX_SUB)) {
            return null;
        }
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
            case 'N':
                return 5;
            default:
                return 6;
        }
    }

    private StringBuilder toString(Read read) {
        StringBuilder sb = new StringBuilder();
        sb.append(">READ").append(this.number).append("_").append(read.getBegin()).append("-").append(read.getEnd());
        if (read.getDescription() != null) {
            sb.append("|").append(read.getDescription()).append("/").append(read.getMatePair());
        }
        sb.append("\n");
        sb.append(read.getAlignedRead()).append("\n");
        return sb;
    }

    private Read cut(Read read) {
        if (read == null || !read.isAligned()) {
            return null;
        }
        if (Globals.RS == null) {
            return read;
        }
        int[] rs = Globals.RS;
        if (read.getEnd() > rs[0] && read.getBegin() < rs[1] && read.isAligned()) {
            if (read.getBegin() < rs[0]) {
                if (read.getEnd() > rs[1]) {
                    read.cut(rs[0] - read.getBegin(), rs[1] - read.getBegin());
                    read.setBegin(rs[0]);
                    read.setEnd(rs[1]);
                } else if (read.getEnd() <= rs[1]) {
                    read.cut(rs[0] - read.getBegin());
                    read.setBegin(rs[0]);
                }
            } else if (read.getBegin() >= rs[0]) {
                if (read.getEnd() > rs[1]) {
                    read.cut(0, rs[1] - read.getBegin());
                    read.setEnd(rs[1]);
                }
            }
            return read;
        }
        return null;
    }
}
