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

import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.utils.Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.javatuples.Pair;
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessingGeneral {

    private boolean virgin = true;
    protected BlockingQueue<Runnable> blockingQueue = new ArrayBlockingQueue<>(1000 * Runtime.getRuntime().availableProcessors());
    protected RejectedExecutionHandler rejectedExecutionHandler = new ThreadPoolExecutor.CallerRunsPolicy();
    protected ExecutorService executor = new ThreadPoolExecutor(2 * Runtime.getRuntime().availableProcessors(), 2 * Runtime.getRuntime().availableProcessors(), 0L, TimeUnit.MILLISECONDS, blockingQueue, rejectedExecutionHandler);
    protected Map<Integer, Map<Integer, Integer>> substitutions = initSubs();
    protected List<Future<Pair<String, Map<Integer, Map<Integer, Integer>>>>> results = new LinkedList<>();
    protected List<Triplet<Double, Double, Double>> inDelSubsList = new LinkedList<>();

    private Map<Integer, Map<Integer, Integer>> initSubs() {
        Map<Integer, Map<Integer, Integer>> map = new HashMap<>();
        for (int v = 0; v < 6; v++) {
            map.put(v, new HashMap<Integer, Integer>());
            for (int b = 0; b < 6; b++) {
                map.get(v).put(b, 0);
            }
        }
        return map;
    }

    protected void flattenFastq(String input) {
        final String newline = System.getProperty("line.separator");
        StringBuilder sb = new StringBuilder();
        try (BufferedReader br = Files.newBufferedReader(FileSystems.getDefault().getPath(input), Charset.defaultCharset())) {
            String line = br.readLine();
            do {
                if (!line.startsWith("@")) {
                    Utils.error();
                }
                sb.append(line.trim());
                sb.append(newline);
                int lines = 0;
                for (;;) {
                    line = br.readLine();
                    if (line.startsWith("+")) {
                        break;
                    }
                    lines++;
                    sb.append(line.trim());
                }
                sb.append(newline);
                sb.append(line.trim());
                sb.append(newline);
                for (int i = 0; i < lines; i++) {
                    line = br.readLine();
                    sb.append(line.trim());
                }
                sb.append(newline);
                line = br.readLine();
            } while (line != null && !line.isEmpty());
        } catch (IOException e) {
            System.err.println("Error parsing file " + input);
        }
        try {
            Files.deleteIfExists(FileSystems.getDefault().getPath(input + "_flat"));
        } catch (IOException ex) {
            Logger.getLogger(ProcessingIlluminaSingle.class.getName()).log(Level.SEVERE, null, ex);
        }
        try (BufferedWriter bw = Files.newBufferedWriter(FileSystems.getDefault().getPath(input + "_flat"), Charset.defaultCharset(), StandardOpenOption.CREATE)) {
            bw.write(sb.toString());
            bw.flush();
        } catch (IOException e) {
            System.err.println("Error parsing file " + input);
        }
    }

    /**
     *
     * @param result
     */
    protected void updateMatrix(Pair<String, Map<Integer, Map<Integer, Integer>>> result) {
        for (int v = 0; v < 6; v++) {
            for (int b = 0; b < 6; b++) {
                substitutions.get(v).put(b, substitutions.get(v).get(b) + result.getValue1().get(v).get(b));
            }
        }
    }

    protected Triplet<Double, Double, Double> getInDelSubRates(Map<Integer, Map<Integer, Integer>> map) {
        double sub = 0;
        double ins = 0;
        double del = 0;
        double sum = 0;
        for (int v = 0; v < 6; v++) {
            for (int b = 0; b < 6; b++) {
                sum += map.get(v).get(b);
            }
        }
        for (int v = 0; v < 6; v++) {
            for (int b = 0; b < 6; b++) {
                double tmp = map.get(v).get(b) / sum;
                if (v == 4 && b != 4) {
                    del += tmp;
                } else if (v != 4 && b == 4) {
                    ins += tmp;
                } else if (v != b) {
                    sub += tmp;
                }
            }
        }
        return Triplet.with(sub, del, ins);
    }

    /**
     *
     */
    protected void printMatrix() {
        System.out.print("\n\nSubstitution matrix:\nr/g");
        for (int v = 0; v < 6; v++) {
            System.out.print("\t" + convert(v));
        }
        System.out.println("");
        double sub = 0;
        double ins = 0;
        double del = 0;
        double sum = 0;
        for (int v = 0; v < 6; v++) {
            for (int b = 0; b < 6; b++) {
                sum += substitutions.get(v).get(b);
            }
        }
        for (int v = 0; v < 6; v++) {
            System.out.print(convert(v));
            for (int b = 0; b < 6; b++) {
                double tmp = substitutions.get(v).get(b) / sum;
                System.out.print("\t" + shorten(tmp));
                if (v == 4 && b != 4) {
                    del += tmp;
                } else if (v != 4 && b == 4) {
                    ins += tmp;
                } else if (v != b) {
                    sub += tmp;
                }
            }
            System.out.println("");
        }

        System.out.println("SUBSTITUTIONS: " + shorten(sub));
        System.out.println("DELETIONS:     " + shorten(del));
        System.out.println("INSERTIONS:    " + shorten(ins));
    }

    private String convert(int c) {
        switch (c) {
            case 0:
                return "A";
            case 1:
                return "C";
            case 2:
                return "G";
            case 3:
                return "T";
            case 4:
                return "-";
            case 5:
                return "N";
            default:
                return "";
        }
    }

    /**
     *
     * @param value
     * @return
     */
    public static String shorten(double value) {
        String s;
//        if (value < 1e-20) {
//            s = "0      ";
//        } else 
        if (value == 1.0) {
            s = "1      ";
        } else {
            String t = "" + value;
            String r;
            if (t.length() > 7) {
                r = t.substring(0, 7);
                if (t.contains("E")) {
                    r = r.substring(0, 4);
                    r += "E" + t.split("E")[1];
                }
                s = r;
            } else {
                s = String.valueOf(value);
            }
        }
        return s;
    }

    protected void processResults() throws InterruptedException, ExecutionException {
        StringBuilder samSB = new StringBuilder();
        if (virgin) {
            samSB.append("@HD\tVN:1.0\tSO:unsorted\n");
            for (Genome g : Globals.GENOMES) {
                samSB.append("@SQ\tSN:").append(g.getHeader()).append("\tLN:").append(g.getSequence().length()).append("\n");
            }
            samSB.append("@PG\tID:InDelFixer\tPN:InDelFixer\tVN:0.6\n");
            virgin = false;
        }
        for (Future<Pair<String, Map<Integer, Map<Integer, Integer>>>> future : results) {
            Pair<String, Map<Integer, Map<Integer, Integer>>> result = future.get();
            if (result != null) {
                samSB.append(result.getValue0());
                this.updateMatrix(result);
            }
        }
        results.clear();
        Utils.appendFile(Globals.output + "reads.sam", samSB.toString());
    }
}