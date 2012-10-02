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
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class AlignWorker extends RecursiveTask<Void> {

    private int start;
    private int end;
    private Read[] reads;
    private String genome;
    private Matrix matrix;

    public AlignWorker(int start, int end, Read[] reads, String genome, Matrix matrix) {
        this.start = start;
        this.end = end;
        this.reads = reads;
        this.genome = genome;
        this.matrix = matrix;
    }

    @Override
    protected Void compute() {
        if (end - start <= Globals.STEPSIZE) {
//            final Map<String, Alignment> map = new HashMap<>();
            for (int i = start; i < end; i++) {
                Alignment align;
                final Read r = reads[i];
                if (r.getDescription().equals("generator-0_746_1245_0")) {
                    System.out.println("");
                }
                if (r.getEnd() < 0) {
                    Globals.printPercentageAligningReads();
                    continue;
//                    align = SmithWatermanGotoh.align(
//                            new Sequence(genome, "", "", Sequence.NUCLEIC),
//                            new Sequence(reads[i].getRead(), "", "", Sequence.NUCLEIC),
//                            matrix, 10, 1);
                } else {
                    if (r.getDescription().equals("generator-0_10_509_0")) {
                        System.out.println("");
                    }
                    if (r.getBegin() <= 0) {
                        align = SmithWatermanGotoh.align(
                                new Sequence(genome, "", "", Sequence.NUCLEIC),
                                new Sequence(r.getRead(), "", "", Sequence.NUCLEIC),
                                matrix, 10, 1);
                        r.setBegin(r.getBegin()+1);
                    } else {
                        int readEnd = r.getEnd() >= genome.length() ? genome.length() : r.getEnd();
                        align = SmithWatermanGotoh.align(
                                new Sequence(genome.substring(r.getBegin() - 1, readEnd), "", "", Sequence.NUCLEIC),
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
                for (int j = 0; j < L; j++) {
                    if (g[j] != '-') {
                        switch (c[j]) {
                            case '-':
                            case 'N':
                                sb.append(g[j]);
                                miss[j] = true;
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
                }
                if (good) {
                    r.setAlignedRead(sb.toString());
                    r.setEnd(r.getBegin() + sb.length());
                } else {
                    System.out.println("x");
                }
                Globals.printPercentageAligningReads();
            }
        } else {
            final int mid = start + (end - start) / 2;
            final AlignWorker left = new AlignWorker(start, mid, reads, genome, matrix);
            final AlignWorker right = new AlignWorker(mid, end, reads, genome, matrix);
            left.fork();
            right.compute();
            left.join();
        }
        return null;
    }
}
