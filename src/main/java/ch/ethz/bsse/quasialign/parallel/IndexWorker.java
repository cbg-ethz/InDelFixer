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
import ch.ethz.bsse.quasialign.stored.Genome;
import ch.ethz.bsse.quasialign.stored.Kmer;
import ch.ethz.bsse.quasialign.stored.Read;
import java.util.List;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class IndexWorker extends RecursiveTask<Void> {

    private int start;
    private int end;
    private Read[] reads;
    private Genome genome;

    public IndexWorker(int start, int end, Read[] reads, Genome genome) {
        this.start = start;
        this.end = end;
        this.reads = reads;
        this.genome = genome;
    }

    @Override
    protected Void compute() {
        if (end - start <= Globals.STEPSIZE) {
            for (int i = start; i < end; i++) {
                Read r = reads[i];
                List<Kmer> kmers = r.getKmers();
                for (Kmer kmer : kmers) {
                    if (this.genome.getKmerMap().containsKey(kmer.getSequence())) {
                        r.addHit(this.genome.getKmerMap().get(kmer.getSequence()), this.genome.getKmerMap().get(kmer.getSequence()) + Globals.KMER_LENGTH);
                    }
                }
                r.findMaxHitRegion();
                r.setBegin(r.getRegion()[0] - r.getRead().length());
                r.setEnd(r.getRegion()[1] + r.getRead().length());
                Globals.printPercentageIndexReads();
            }
        } else {
            final int mid = start + (end - start) / 2;
            final IndexWorker left = new IndexWorker(start, mid, reads, genome);
            final IndexWorker right = new IndexWorker(mid, end, reads, genome);
            left.fork();
            right.compute();
            left.join();
        }
        return null;
    }
}
