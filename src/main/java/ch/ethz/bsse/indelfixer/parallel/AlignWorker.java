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

import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.Read;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class AlignWorker extends RecursiveTask<Void> {

    private int start;
    private int end;
    private Matrix matrix;

    public AlignWorker(Matrix matrix, int start, int end) {
        this.start = start;
        this.end = end;
        this.matrix = matrix;
    }

    @Override
    protected Void compute() {
        if (end - start <= Globals.STEPSIZE) {
//            final Map<String, Alignment> map = new HashMap<>();
            for (int i = start; i < end; i++) {
                Alignment align;
                final Read r = Globals.READS.get(i);
                if (r.getBestFittingGenome() == -1) {
//                    System.out.println(r.getRead());
                    continue;
                }
                if (r.getEnd() < 0) {
                    Globals.printPercentageAligningReads();
                    continue;
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
                Globals.printPercentageAligningReads();
            }
        } else {
            final int mid = start + (end - start) / 2;
            final AlignWorker left = new AlignWorker(matrix, start, mid);
            final AlignWorker right = new AlignWorker(matrix, mid, end);
            left.fork();
            right.compute();
            left.join();
        }
        return null;
    }
}
