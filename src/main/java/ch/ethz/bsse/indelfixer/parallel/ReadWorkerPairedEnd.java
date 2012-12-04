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
import ch.ethz.bsse.indelfixer.utils.Utils;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.RecursiveTask;
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadWorkerPairedEnd extends RecursiveTask<List<Read>> {

    private int start;
    private int end;
    private List<Triplet<String, String, Integer>> reads;

    public ReadWorkerPairedEnd(int start, int end, List<Triplet<String, String, Integer>> reads) {
        this.start = start;
        this.end = end;
        this.reads = reads;
    }

    @Override
    protected List<Read> compute() {
        if (end - start <= Globals.STEPSIZE) {
            final List<Read> list = new ArrayList<>();
            for (int i = start; i < end; i++) {
//                if (reads.get(i).getValue2() == 1) {
                final Read r = new Read(reads.get(i).getValue0(), i, false);
                r.setDescription(reads.get(i).getValue1());
                r.setMatePair(reads.get(i).getValue2());
                list.add(r);
//                } else if (reads.get(i).getValue2() == 2) {
                final Read rR = new Read(Utils.reverseComplement(reads.get(i).getValue0()), i, true);
                rR.setDescription(reads.get(i).getValue1());
                rR.setMatePair(reads.get(i).getValue2());
                Globals.printPercentageProcessingReads();
                list.add(rR);
//                }
//                int l = reads[i].getValue1().length;
//                int[] qualityR = new int[l];
//                for (int j = 0; j < l; j++) {
//                    qualityR[l-j-1] = reads[i].getValue1()[j];
//                }
            }
            return list;
        } else {
            final int mid = start + (end - start) / 2;
            final ReadWorkerPairedEnd left = new ReadWorkerPairedEnd(start, mid, reads);
            final ReadWorkerPairedEnd right = new ReadWorkerPairedEnd(mid, end, reads);
            left.fork();
            final List<Read> list = right.compute();
            list.addAll(left.join());
            return list;
        }
    }
}
