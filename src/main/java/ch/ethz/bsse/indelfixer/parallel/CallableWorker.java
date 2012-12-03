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
import java.util.List;
import java.util.concurrent.Callable;
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class CallableWorker implements Callable<StringBuilder> {

    private Triplet<String, String, Integer> watsonTriple;
    private Triplet<String, String, Integer> crickTriple;
    private int L;
    private Genome[] genome;
    private Matrix matrix;
    private int number;

    public CallableWorker(Triplet<String, String, Integer> watson, Triplet<String, String, Integer> crick, int L, Genome[] genome, Matrix matrix, int number) {
        this.watsonTriple = watson;
        this.crickTriple = crick;
        this.L = L;
        this.genome = genome;
        this.matrix = matrix;
        this.number = number;
    }

    @Override
    public StringBuilder call() throws Exception {
        Read watsonRead = map(createRead(watsonTriple, false));
        Read watsonRevRead = map(createRead(watsonTriple, true));
        Read crickRead = map(createRead(crickTriple, false));
        Read crickRevRead = map(createRead(crickTriple, true));

        Read watson = align(watsonRead.getMaximumHits() > watsonRevRead.getMaximumHits() ? watsonRead : watsonRevRead);
        Read crick = align(crickRead.getMaximumHits() > crickRevRead.getMaximumHits() ? crickRead : crickRevRead);

        StringBuilder sb = new StringBuilder();
        if (watson != null && watson.isAligned()) {
            sb.append(toString(watson));
        }
        if (crick != null && crick.isAligned()) {
            sb.append(toString(crick));
        }
        
        return sb;
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

    private Read align(Read r) {
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
        boolean[] miss = new boolean[align.getSequence2().length];
        int sum10 = 0;
        boolean good = true;
        char[] m = align.getMarkupLine();
        char[] c = align.getSequence2();
        char[] g = align.getSequence1();
        int L = m.length;
        if (m.length != c.length) {
            System.err.println("DIFFERENT LENGTHS");
        }
        int maxMisses = 0;
        int consecutiveMisses = 0;
        for (int j = 0; j < L; j++) {
            if (g[j] != '-') {
                switch (c[j]) {
                    case '-':
                    case 'N':
                        if (Globals.FILL) {
                            sb.append(g[j]);
                        } else {
                            sb.append(c[j]);
                        }
                        if (miss[j - 1] && miss[j - 2]) {
                            miss[j - 1] = false;
                            miss[j - 2] = false;
                            sum10 -= 2;
                        }
                        miss[j] = true;
                        consecutiveMisses++;
                        break;
                    default:
                        maxMisses = Math.max(maxMisses, consecutiveMisses);
                        consecutiveMisses = 0;
                        sb.append(c[j]);
                        break;
                }
            } else {
                switch (c[j]) {
                    case '-':
                    case 'N':
                        sb.append("-");
                        break;
                    default:
                        sb.append(c[j]);
                        break;
                }
            }
            sum10 += miss[j] ? 1 : 0;
            if (j > 10) {
                if (sum10 >= 3) {
                    good = false;
                    break;
                }
                sum10 -= miss[j - 10] ? 1 : 0;
            }
            if (maxMisses > 6) {
                good = false;
                break;
            }
        }
        if (good) {
            r.setAlignedRead(sb.toString());
            r.setEnd(r.getBegin() + sb.length());
        } else {
//                    System.out.println("x");
        }
        return r;
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
