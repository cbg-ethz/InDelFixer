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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessingFastaSingle extends ProcessingGeneral {

    private String inputWatson;

    /**
     * Constructor.
     * 
     * @param inputWatson Path to multiple fasta file
     */
    public ProcessingFastaSingle(String inputWatson) {
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
        BufferedReader brWatson = new BufferedReader(new FileReader(new File(this.inputWatson)));

        for (int i = 0;; i++) {
            try {
                SequenceEntry watsonS = new SequenceEntry(parseFastaEntry(brWatson));
                if (watsonS.sequence.length() >= Globals.MIN_LENGTH) {
                    results.add(executor.submit(new FutureSequence(watsonS, i)));
                }
            } catch (IllegalAccessError e) {
                // used to halt in case of EOF
                break;
            }
            if (i % 10000 == 0) {
                this.processResults();
            }
        }
        this.processResults();
        this.printMatrix();
        executor.shutdown();
    }

    /**
     * Parses a single sequence of given BufferedReader.
     * 
     * @param br BufferedReader
     * @return Single sequence
     * @throws IllegalAccessError Thrown if BufferedReader has reached EOF.
     */
    public static String parseFastaEntry(BufferedReader br) throws IOException, IllegalAccessError {
        boolean untouched = true;
        boolean cut = false;
        boolean first = true;

        StringBuilder sb = new StringBuilder();
        String strLine;
        while ((strLine = br.readLine()) != null) {
            untouched = false;
            if (strLine.startsWith(">")) {
                if (sb.length() > 0) {
                    return sb.toString();
                }
            } else {
                if (first) {
                    boolean upperCase = false;
                    boolean lowerCase = false;
                    for (char c : strLine.toCharArray()) {
                        if (!lowerCase) {
                            if (Character.isLowerCase(c)) {
                                lowerCase = true;
                                continue;
                            }
                        }
                        if (!upperCase) {
                            if (Character.isUpperCase(c)) {
                                upperCase = true;
                                continue;
                            }
                        }
                        if (upperCase && lowerCase) {
                            cut = true;
                            break;
                        }
                    }
                    first = false;
                }
                if (cut) {
                    for (char c : strLine.toCharArray()) {
                        if (Character.isUpperCase(c)) {
                            sb.append(c);
                        }
                    }
                } else {
                    sb.append(strLine);
                }
            }
        }
        if (untouched) {
            throw new IllegalAccessError("done");
        }
        return sb.toString();
    }
}