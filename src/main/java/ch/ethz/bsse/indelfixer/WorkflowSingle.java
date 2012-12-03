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
package ch.ethz.bsse.indelfixer;

import ch.ethz.bsse.indelfixer.parallel.CallableWorkerSingle;
import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.Read;
import ch.ethz.bsse.indelfixer.utils.FastaParser;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import ch.ethz.bsse.indelfixer.utils.Utils;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.javatuples.Pair;

/**
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class WorkflowSingle implements WorkflowI {

    private Genome[] genome;
    private Map<Character, Map<Character, Integer>> substitutions = new HashMap<>();

    public WorkflowSingle(Genome[] genome, String[] forward, int[][] rs) {
        for (char v = 0; v < 5; v++) {
            substitutions.put(v, new HashMap<Character, Integer>());
            for (char b = 0; b < 5; b++) {
                substitutions.get(v).put(b, 0);
            }
        }
        this.genome = genome;
        Globals.N = forward.length;
        WorkflowSingle.removeLogging();

        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1);
        List<FutureTask<Pair<Read, Map<Character, Map<Character, Integer>>>>> taskList = new ArrayList<>();

        int i = 0;
        Matrix matrix = loadMatrix();
        for (String f : forward) {
            FutureTask<Pair<Read, Map<Character, Map<Character, Integer>>>> futureTask_1 = new FutureTask<>(new CallableWorkerSingle(f, genome[0].getSequence().length(), genome, matrix, i++));
            taskList.add(futureTask_1);
            executor.execute(futureTask_1);
        }
        StringBuilder sb = new StringBuilder();
        List<Read> reads = new ArrayList<>();
        for (int j = 0; j < taskList.size(); j++) {
            FutureTask<Pair<Read, Map<Character, Map<Character, Integer>>>> futureTask = taskList.get(j);
            try {
                Pair<Read, Map<Character, Map<Character, Integer>>> pair = futureTask.get();
                if (pair != null && pair.getValue0() != null) {
                    Read read = pair.getValue0();
                    boolean hit = false;
                    if (rs != null) {
                        for (int x = 0; x < rs.length; x++) {
                            if (read.getEnd() > rs[x][0] && read.getBegin() < rs[x][1] && read.isAligned()) {
                                String s = null;
                                if (read.getBegin() < rs[x][0]) {
                                    if (read.getEnd() > rs[x][1]) {
                                        s = read.getAlignedRead().substring(rs[x][0] - read.getBegin(), rs[x][1] - read.getBegin());
                                        read.setBegin(rs[x][0]);
                                        read.setEnd(rs[x][1]);
                                    } else if (read.getEnd() <= rs[x][1]) {
                                        s = read.getAlignedRead().substring(rs[x][0] - read.getBegin());
                                        read.setBegin(rs[x][0]);
                                    }
                                } else if (read.getBegin() >= rs[x][0]) {
                                    if (read.getEnd() > rs[x][1]) {
                                        s = read.getAlignedRead().substring(0, rs[x][1] - read.getBegin());
                                        read.setEnd(rs[x][1]);
                                    } else {
                                        s = read.getAlignedRead();
                                    }
                                }
                                if (s.length() > rs[x][1] - rs[x][0]) {
                                    System.out.println("");
                                }
                                read.setAlignedRead(s);

                                hit = true;
                                break;
                            }
                        }
                    } else {
                        hit = true;
                    }

                    if (hit) {
                        reads.add(read);
                        sb.append(toString(read));
                        for (char v = 0; v < 5; v++) {
                            for (char b = 0; b < 5; b++) {
                                substitutions.get(v).put(b, substitutions.get(v).get(b) + pair.getValue1().get(v).get(b));
                            }
                        }
                    }
                }
            } catch (InterruptedException | ExecutionException ex) {
                System.out.println(j + "\tERROR");
                System.err.println(ex.getLocalizedMessage());
//                Logger.getLogger(WorkflowPaired.class.getName()).log(Level.SEVERE, null, ex);
            }

            StatusUpdate.print(
                    "Processing reads:\t" + (Math.round((1000d * j) / Globals.N / 10d)) + "%");
        }
        StatusUpdate.println("Processing reads:\t100%");
        Utils.saveFile(Globals.output + "region.fasta", sb.toString());
        printAlignment(reads);
        System.out.print("\nSUBS");
        for (int v = 0; v < 5; v++) {
            System.out.print("\t" + v);
        }
        System.out.println("");
        int a = 0;
        int sub = 0;
        int ins = 0;
        int del = 0;
        double sum = 0;
        for (char v = 0; v < 5; v++) {
            System.out.print(a++);
            for (char b = 0; b < 5; b++) {
                System.out.print("\t" + substitutions.get(v).get(b));
                if (v != 4 && b == 4) {
                    del += substitutions.get(v).get(b);
                } else if (v == 4 && b != 4) {
                    ins += substitutions.get(v).get(b);
                } else if (v != b) {
                    sub += substitutions.get(v).get(b);
                }
                sum += substitutions.get(v).get(b);
            }
            System.out.println("");
        }

        System.out.println("SUBSTITUTIONS: " + sub / sum);
        System.out.println("DELETIONS:     " + del / sum);
        System.out.println("INSETIONS:     " + ins / sum);
        System.exit(0);
    }
    int runningNumber = 0;

    private StringBuilder toString(Read read) {
        StringBuilder sb = new StringBuilder();
        sb.append(">READ").append(runningNumber++).append("_").append(read.getBegin()).append("-").append(read.getEnd());
        sb.append("\n");
        sb.append(read.getAlignedRead()).append("\n");
        return sb;
    }

    public void printAlignment(List<Read> reads) {
        Map<Integer, List<Read>> readMap = new HashMap<>();
        for (Read r : reads) {
            if (!readMap.containsKey(r.getBegin())) {
                readMap.put(r.getBegin(), new ArrayList<Read>());
            }
            readMap.get(r.getBegin()).add(r);
        }
        Integer[] order = readMap.keySet().toArray(new Integer[readMap.size()]);
        Arrays.sort(order);
        int min = order[0];
        StringBuilder sb = new StringBuilder();
        sb.append(min).append("\n");
        for (int item : order) {
            List<Read> currentReads = readMap.get(item);
            for (Read currentRead : currentReads) {
                for (int i = min; i < item; i++) {
                    sb.append(" ");
                }
                sb.append(currentRead.getAlignedRead()).append("\n");
            }
        }
        Utils.saveFile(Globals.output + "alignment.txt", sb.toString());
    }

    private Genome parseGenomeRead(String genomePath) {
        Genome g = new Genome(genomePath);
        String out = "";
//        for (int i = 50;; i++) {
//            out = "Searching for smallest kmer:\t";
//            try {
//                Globals.KMER_LENGTH = i;
//                for (int j = 0; j < i; j++) {
//                    out += ("|");
//                }
//                g = new Genome(genomePath);
//                break;
//            } catch (IllegalStateException e) {
//            }
//            StatusUpdate.print(out);
//            Globals.KMER_LENGTH = i;
//        }
        StatusUpdate.println(out + " (" + Globals.KMER_LENGTH + ")");
        return g;
    }

    private Genome[] parseGenome(String genomePath) {
        String[] g = FastaParser.parseFarFile(genomePath);
        Genome[] gs = new Genome[g.length];
        for (int i = 0; i < g.length; i++) {
            gs[i] = parseGenomeRead(g[i]);
        }
        Globals.GENOMES = gs;
        Globals.GENOME_COUNT = gs.length;
        Globals.GENOME_SEQUENCES = g;
        Globals.GENOME_LENGTH = g[0].length();
        return gs;
    }

    public void saveCoverage() {
    }

    public static void removeLogging() throws SecurityException {
        Logger rootLogger = Logger.getLogger("");
        Handler[] handlers = rootLogger.getHandlers();
        if (handlers.length > 0) {
            rootLogger.removeHandler(handlers[0]);
        }
    }

    public static Matrix loadMatrix() {
        Matrix matrix = null;
        try {
            matrix = MatrixLoader.load("EDNAFULL");


        } catch (MatrixLoaderException ex) {
            Logger.getLogger(WorkflowSingle.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        return matrix;
    }
}
