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
package ch.ethz.bsse.quasialign.parallel;

import ch.ethz.bsse.quasialign.stored.Globals;
import ch.ethz.bsse.quasialign.stored.Read;
import ch.ethz.bsse.quasialign.utils.Utils;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.RecursiveTask;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadWorker extends RecursiveTask<List<Read>> {

    private int start;
    private int end;
    private String[] reads;

    public ReadWorker(int start, int end, String[] reads) {
        this.start = start;
        this.end = end;
        this.reads = reads;
    }

    @Override
    protected List<Read> compute() {
        if (end - start <= Globals.STEPSIZE) {
            final List<Read> list = new ArrayList<>();
            for (int i = start; i < end; i++) {
                final Read r = new Read(reads[i],i,false);
//                int l = reads[i].getValue1().length;
//                int[] qualityR = new int[l];
//                for (int j = 0; j < l; j++) {
//                    qualityR[l-j-1] = reads[i].getValue1()[j];
//                }
                final Read rR = new Read(Utils.reverseComplement(reads[i]),i,true);
                Globals.printPercentageProcessingReads();
                list.add(r);
                list.add(rR);
            }
            return list;
        } else {
            final int mid = start + (end - start) / 2;
            final ReadWorker left = new ReadWorker(start, mid, reads);
            final ReadWorker right = new ReadWorker(mid, end, reads);
            left.fork();
            final List<Read> list = new LinkedList<>();
            list.addAll(right.compute());
            list.addAll(left.join());
            return list;
        }
    }
}
