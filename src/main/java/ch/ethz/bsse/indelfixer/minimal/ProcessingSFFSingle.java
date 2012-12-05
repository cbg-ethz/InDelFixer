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
package ch.ethz.bsse.indelfixer.minimal;

import ch.ethz.bsse.indelfixer.stored.Globals;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessingSFFSingle extends ProcessingGeneral {

    private String[] sequences;

    public ProcessingSFFSingle(String[] sequences) {
        this.sequences = sequences;
        try {
            this.start();
        } catch (IOException | InterruptedException | ExecutionException e) {
            System.err.println("SFF problem " + e.getLocalizedMessage());
        }
    }

    private void start() throws FileNotFoundException, IOException, InterruptedException, ExecutionException {
        for (int i = 0; i < this.sequences.length; i++) {
            SequenceEntry watsonS = new SequenceEntry(sequences[i]);
            if (watsonS.sequence.length() >= Globals.MIN_LENGTH) {
                results.add(executor.submit(new FutureSequence(watsonS, i)));
            }
        }

        this.processResults();
        executor.shutdown();
    }
}