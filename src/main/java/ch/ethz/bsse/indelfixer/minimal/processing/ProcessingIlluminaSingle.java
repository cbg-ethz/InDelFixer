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

import static ch.ethz.bsse.indelfixer.minimal.processing.ProcessingFastaSingle.parseFastaEntry;
import ch.ethz.bsse.indelfixer.minimal.processing.parallel.FutureSequence;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.SequenceEntry;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessingIlluminaSingle extends ProcessingGeneral {

    private String inputWatson;

    /**
     * Constructor.
     *
     * @param inputWatson Path to multiple fastq file for forward.
     */
    public ProcessingIlluminaSingle(String inputWatson) {
        this.inputWatson = inputWatson;
        try {
            this.start();
        } catch (IOException | InterruptedException | ExecutionException e) {
        }
    }

    /**
     * Queues sequences for threading and initiates processing of results.
     */
    private void start() throws FileNotFoundException, IOException, InterruptedException, ExecutionException {
        BufferedReader brWatson;
        if (Globals.FLAT) {
            flattenFastq(this.inputWatson);
            brWatson = new BufferedReader(new FileReader(new File(this.inputWatson + "_flat")));
        } else {
            brWatson = new BufferedReader(new FileReader(new File(this.inputWatson)));
        }

        List<SequenceEntry> list = new LinkedList<>();
        for (int i = 0;;) {
            try {
                for (int j = 0; j < Globals.STEPS; j++) {
                    i++;
                    SequenceEntry watsonQ = parseFastq(brWatson);
                    if (watsonQ != null && watsonQ.sequence.length() >= Globals.MIN_LENGTH) {
                        list.add(watsonQ);
                    } else {
                        StatusUpdate.processLength();
                    }
                    if (i % 10000 == 0) {
                        this.processResults();
                    }
                }
                synchronized (results) {
                    results.add(executor.submit(new FutureSequence(list)));
                }
            } catch (IllegalAccessError e) {
                if (!list.isEmpty()) {
                    synchronized (results) {
                        results.add(executor.submit(new FutureSequence(list)));
                    }
                }
                // used to halt in case of EOF
                break;
            }
        }

        this.processResults();
        this.printMatrix();
        this.saveConsensus();
        executor.shutdown();
    }

    /**
     * Parses a single fastq block of given BufferedReader.
     *
     * @param br BufferedReader
     * @return Single fastq block of type SequenceEntry
     * @throws IllegalAccessError Thrown if BufferedReader has reached EOF.
     */
    private SequenceEntry parseFastq(BufferedReader br) throws IOException, IllegalAccessError {
        //head
        String header = br.readLine();
        if (header == null) {
            throw new IllegalAccessError();
        }
        //sequence
        String seq = br.readLine();
        char[] c = seq.toCharArray();
        //description
        br.readLine();
        //quality
        String qualityString = br.readLine();
        char[] quality = qualityString.toCharArray();
        int[] p = new int[quality.length];

        double qualitySum = 0;
        int begin = -1;
        boolean started = false;
        int end = -1;

        int Ns = 0;

        for (int i = 0; i < quality.length; i++) {
            p[i] = ((int) quality[i]) - 33;
            if (!started && p[i] >= Globals.THRESHOLD) {
                started = true;
                begin = i;
                qualitySum += p[i];
            }
            if (started) {
                qualitySum += p[i];
            }
        }
        boolean foundEnd = false;
        for (int i = quality.length - 1; i >= begin; i--) {
            qualitySum -= p[i];
            if (!foundEnd && p[i] >= Globals.THRESHOLD) {
                end = i;
                foundEnd = true;
            }
            if (foundEnd) {
                if (c[i] == 'N') {
                    Ns++;
                } else {
                    if (Ns >= Globals.maxN) {
                        break;
                    } else {
                        Ns = 0;
                    }
                }
            }
        }
        if (Ns >= Globals.maxN) {
            return null;
        }
        if (begin == -1 || end == -1) {
            return null;
        }
        if (begin < Globals.CUT) {
            begin = Globals.CUT;
        }
        if (begin > end) {
            return null;
        }
        qualitySum /= end - begin - 1;
        return new SequenceEntry(seq.substring(begin, end + 1),
                header,
                qualityString.substring(begin, end + 1));
    }
}
