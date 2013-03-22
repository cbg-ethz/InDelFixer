/**
 * Copyright (c) 2011-2013 Armin Töpfer
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
package ch.ethz.bsse.indelfixer.minimal.processing;

import ch.ethz.bsse.indelfixer.minimal.processing.parallel.FutureSequence;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.SequenceEntry;
import ch.ethz.bsse.indelfixer.stored.SimpleRead;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessingSFFSingle extends ProcessingGeneral {

    private SimpleRead[] sequences;

    /**
     * Constructor.
     *
     * @param sequences
     */
    public ProcessingSFFSingle(SimpleRead[] sequences) {
        this.sequences = sequences;
        try {
            this.start();
        } catch (IOException | InterruptedException | ExecutionException e) {
        }
    }

    /**
     * Queues sequences for threading and initiates processing of results.
     */
    private void start() throws FileNotFoundException, IOException, InterruptedException, ExecutionException {
        for (int i = 0; i < this.sequences.length; i++) {
            SequenceEntry watsonS = new SequenceEntry(sequences[i]);
            if (watsonS.sequence.length() >= Globals.MIN_LENGTH) {
                List<SequenceEntry> l = new LinkedList<>();
                l.add(watsonS);
                synchronized (results) {
                    results.add(executor.submit(new FutureSequence(l)));
                }
            } else {
                StatusUpdate.processLength();
            }
            if (i % 10000 == 0) {
                this.processResults();
            }
        }

        this.processResults();
        this.printMatrix();
        this.saveConsensus();
        executor.shutdown();
    }
}