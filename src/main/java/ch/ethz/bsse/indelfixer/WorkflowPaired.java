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

import ch.ethz.bsse.indelfixer.parallel.CallableWorker;
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
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class WorkflowPaired implements WorkflowI {

    private Genome[] genome;

    public WorkflowPaired(String genomePath, List<Triplet<String, String, Integer>> forward, List<Triplet<String, String, Integer>> backward) {
        this.genome = parseGenome(genomePath);
        Globals.N = forward.size();
        WorkflowPaired.removeLogging();
        Map<String, Triplet<String, String, Integer>> backMap = new HashMap<>();
        for (Triplet<String, String, Integer> b : backward) {
            backMap.put(b.getValue1(), b);
        }


        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1);
        List<FutureTask<StringBuilder>> taskList = new ArrayList<>();

        int i = 0;
        Matrix matrix = loadMatrix();
        for (Triplet<String, String, Integer> f : forward) {
            FutureTask<StringBuilder> futureTask_1 = new FutureTask<>(new CallableWorker(f, backMap.get(f.getValue1()), genome[0].getSequence().length(), genome, matrix, i++));
            taskList.add(futureTask_1);
            executor.execute(futureTask_1);
        }
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < taskList.size(); j++) {
            FutureTask<StringBuilder> futureTask = taskList.get(j);
            try {
                StringBuilder sbb = futureTask.get();
                if (sbb != null) {
                    sb.append(futureTask.get());
                }
            } catch (InterruptedException | ExecutionException ex) {
                System.out.println(j + "\tERROR");
                System.err.println(ex.getLocalizedMessage());
//                Logger.getLogger(WorkflowPaired.class.getName()).log(Level.SEVERE, null, ex);
            }
            StatusUpdate.print("Processing reads:\t" + (Math.round((1000d * j) / Globals.N / 10d)) + "%");
        }
        StatusUpdate.println("Processing reads:\t100%");
        Utils.saveFile(Globals.output + "region.fasta", sb.toString());
        System.exit(0);
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
            Logger.getLogger(WorkflowPaired.class.getName()).log(Level.SEVERE, null, ex);
        }
        return matrix;
    }
}
