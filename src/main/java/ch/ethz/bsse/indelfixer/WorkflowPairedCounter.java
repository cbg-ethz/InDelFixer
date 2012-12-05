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

import ch.ethz.bsse.indelfixer.parallel.CallableWorkerPairedCounter;
import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.utils.FastaParser;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import java.util.ArrayList;
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
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class WorkflowPairedCounter implements WorkflowI {

    private Genome[] genome;
    private Map<Integer, Map<Integer, Integer>> substitutions = new HashMap<>();

    public WorkflowPairedCounter(String genomePath, List<Triplet<String, String, Integer>> forward, List<Triplet<String, String, Integer>> backward, int[][] rs) {
        for (int v = 0; v < 6; v++) {
            substitutions.put(v, new HashMap<Integer, Integer>());
            for (int b = 0; b < 6; b++) {
                substitutions.get(v).put(b, 0);
            }
        }

        this.genome = parseGenome(genomePath);
        Globals.N = forward.size();
        WorkflowPairedCounter.removeLogging();
        Map<String, Triplet<String, String, Integer>> backMap = new HashMap<>();
        for (Triplet<String, String, Integer> b : backward) {
            backMap.put(b.getValue1(), b);
        }


        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1);
        List<FutureTask<List<Map<Integer, Map<Integer, Integer>>>>> taskList = new ArrayList<>();

        int i = 0;
        Matrix matrix = loadMatrix();
        for (Triplet<String, String, Integer> f : forward) {
            FutureTask<List<Map<Integer, Map<Integer, Integer>>>> futureTask_1 = new FutureTask<>(new CallableWorkerPairedCounter(f, backMap.get(f.getValue1()), genome[0].getSequence().length(), genome, matrix, i++));
            taskList.add(futureTask_1);
            executor.execute(futureTask_1);
        }
        StringBuilder sb = new StringBuilder();
//        List<Read> reads = new ArrayList<>();
        for (int j = 0; j < taskList.size(); j++) {
            FutureTask<List<Map<Integer, Map<Integer, Integer>>>> futureTask = taskList.get(j);
            if (futureTask != null) {
                try {
                    List<Map<Integer, Map<Integer, Integer>>> result = futureTask.get();
                    for (Map<Integer, Map<Integer, Integer>> map : result) {
                        for (int v = 0; v < 6; v++) {
                            for (int b = 0; b < 6; b++) {
                                substitutions.get(v).put(b, substitutions.get(v).get(b) + map.get(v).get(b));
                            }
                        }
                    }
                } catch (InterruptedException | ExecutionException ex) {
//                    System.out.println(j + "\tERROR");
//                    System.err.println(ex.getLocalizedMessage());
//                Logger.getLogger(WorkflowPaired.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            StatusUpdate.print("Processing reads:\t" + (Math.round((1000d * j) / Globals.N / 10d)) + "%");
        }
        StatusUpdate.println("Processing reads:\t100%");

//        printAlignment(readss);
        System.out.print("\nr/g");
        for (int v = 0; v < 6; v++) {
            System.out.print("\t" + v);
        }
        System.out.println("");
        int a = 0;
        int sub = 0;
        int ins = 0;
        int del = 0;
        double sum = 0;
        for (int v = 0; v < 6; v++) {
            System.out.print(a++);
            for (int b = 0; b < 6; b++) {
                System.out.print("\t" + substitutions.get(v).get(b));
                if (v == 4 && b != 4) {
                    del += substitutions.get(v).get(b);
                } else if (v != 4 && b == 4) {
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
        System.out.println("INSERTIONS:     " + ins / sum);
        System.exit(0);
    }

    private Genome parseGenomeRead(String genomePath) {
        Genome g = new Genome(genomePath);
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

    @Override
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
            Logger.getLogger(WorkflowPairedCounter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return matrix;
    }
}
