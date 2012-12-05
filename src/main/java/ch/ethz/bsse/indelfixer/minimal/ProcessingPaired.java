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
import ch.ethz.bsse.indelfixer.utils.Utils;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessingPaired extends ProcessingGeneral {

    private String inputWatson;
    private String inputCrick;

    public ProcessingPaired(String inputWatson, String inputCrick) {
        this.inputWatson = inputWatson;
        this.inputCrick = inputCrick;
        try {
            this.start();
        } catch (IOException | InterruptedException | ExecutionException e) {
            System.err.println("Fastq problem " + e.getLocalizedMessage());
        }
    }

    private void start() throws FileNotFoundException, IOException, InterruptedException, ExecutionException {
        BufferedReader brWatson = new BufferedReader(new FileReader(new File(this.inputWatson)));
        BufferedReader brCrick = new BufferedReader(new FileReader(new File(this.inputCrick)));

        List<Future<Pair<String, List<Map<Integer, Map<Integer, Integer>>>>>> results = new LinkedList<>();
        for (int i = 0;; i++) {
            try {
                FastqEntry watsonQ = parseFastq(brWatson);
                FastqEntry crickQ = parseFastq(brCrick);
                results.add(executor.submit(new FuturePaired(watsonQ, crickQ, i)));
            } catch (IllegalAccessError e) {
                // used to halt in case of EOF
                break;
            }
        }

        StringBuilder sb = new StringBuilder();
        for (Future<Pair<String, List<Map<Integer, Map<Integer, Integer>>>>> future : results) {
            Pair<String, List<Map<Integer, Map<Integer, Integer>>>> result = future.get();
            if (result != null) {
                sb.append(result.getValue0());
                this.updateMatrix(result);
            }
        }
        executor.shutdown();
        Utils.saveFile(Globals.output + "reads.fasta", sb.toString());
        this.printMatrix();
    }

    private FastqEntry parseFastq(BufferedReader br) throws IOException, IllegalAccessError {
        //head
        String header = br.readLine();
        if (header == null) {
            throw new IllegalAccessError();
        }
        String tag = header.split(" ")[0];
        int pairedNumber = Integer.parseInt(header.split(" ")[1].split(":")[0]);
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
                if (c[i] == 'N') {
                    Ns++;
                } else {
                    if (Ns >= 3) {
                        break;
                    } else {
                        Ns = 0;
                    }
                }
            }
        }
        if (Ns >= 3) {
            return null;
        }
        for (int i = quality.length - 1; i >= 0; i--) {
            qualitySum -= p[i];
            if (p[i] >= Globals.THRESHOLD) {
                end = i;
                break;
            }
        }
        if (begin == -1) {
            return null;
        }
        qualitySum /= end - begin - 1;
        return new FastqEntry(seq.substring(begin, end + 1),
                tag,
                pairedNumber,
                qualityString.substring(begin, end + 1));
    }
}

class FastqEntry {

    String sequence;
    String tag;
    int pairedNumber;
    String quality;

    public FastqEntry(String sequence, String tag, int pairedNumber, String quality) {
        this.sequence = sequence;
        this.tag = tag;
        this.pairedNumber = pairedNumber;
        this.quality = quality;
    }
}