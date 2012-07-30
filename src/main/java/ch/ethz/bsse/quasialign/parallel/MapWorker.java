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
public class MapWorker extends RecursiveTask<Map<String, Alignment>> {

    private int start;
    private int end;
    private String[] kmers;
    private String genome;
    private Matrix matrix;

    public MapWorker(int start, int end, String[] reads, String genome, Matrix matrix) {
        this.start = start;
        this.end = end;
        this.kmers = reads;
        this.genome = genome;
        this.matrix = matrix;
    }

    @Override
    protected Map<String, Alignment> compute() {
        if (end - start <= Globals.STEPSIZE) {
            final Map<String, Alignment> map = new HashMap<>();
            for (int i = start; i < end; i++) {
                Alignment align = SmithWatermanGotoh.align(
                        new Sequence(genome, "", "", Sequence.NUCLEIC),
                        new Sequence(kmers[i], "", "", Sequence.NUCLEIC),
                        matrix, 10, 1);
                map.put(kmers[i], align);
                Globals.printPercentageMappingReads();
            }
            return map;
        } else {
            final int mid = start + (end - start) / 2;
            final MapWorker left = new MapWorker(start, mid, kmers,genome,matrix);
            final MapWorker right = new MapWorker(mid, end, kmers,genome,matrix);
            left.fork();
            final Map<String, Alignment> map = new HashMap<>();
            map.putAll(right.compute());
            map.putAll(left.join());
            return map;
        }
    }
}
