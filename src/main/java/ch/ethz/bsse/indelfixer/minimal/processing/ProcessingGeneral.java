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
package ch.ethz.bsse.indelfixer.minimal.processing;

import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.utils.Utils;
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
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessingGeneral {

    protected BlockingQueue<Runnable> blockingQueue = new ArrayBlockingQueue<>(100 * Runtime.getRuntime().availableProcessors());
    protected RejectedExecutionHandler rejectedExecutionHandler = new ThreadPoolExecutor.CallerRunsPolicy();
    protected ExecutorService executor = new ThreadPoolExecutor(Runtime.getRuntime().availableProcessors() - 1, Runtime.getRuntime().availableProcessors() - 1, 0L, TimeUnit.MILLISECONDS, blockingQueue, rejectedExecutionHandler);
    protected Map<Integer, Map<Integer, Integer>> substitutions = initSubs();
    protected List<Future<Pair<String, List<Map<Integer, Map<Integer, Integer>>>>>> results = new LinkedList<>();

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

    /**
     *
     * @param result
     */
    protected void updateMatrix(Pair<String, List<Map<Integer, Map<Integer, Integer>>>> result) {
        for (Map<Integer, Map<Integer, Integer>> sub : result.getValue1()) {
            for (int v = 0; v < 6; v++) {
                for (int b = 0; b < 6; b++) {
                    substitutions.get(v).put(b, substitutions.get(v).get(b) + sub.get(v).get(b));
                }
            }
        }
    }

    /**
     *
     */
    protected void printMatrix() {
        System.out.print("\n\nr/g");
        for (int v = 0; v < 6; v++) {
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
                if (v == 4 && b != 4) {
                    del += tmp;
                } else if (v != 4 && b == 4) {
                    ins += tmp;
                } else if (v != b) {
                    sub += tmp;
                }
            }
        }

        System.out.println("SUBSTITUTIONS: " + shorten(sub));
        System.out.println("DELETIONS:     " + shorten(del));
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

    /**
     *
     * @throws InterruptedException
     * @throws ExecutionException
     */
    protected void processResults() throws InterruptedException, ExecutionException {
        StringBuilder sb = new StringBuilder();
        for (Future<Pair<String, List<Map<Integer, Map<Integer, Integer>>>>> future : results) {
            Pair<String, List<Map<Integer, Map<Integer, Integer>>>> result = future.get();
            if (result != null) {
                sb.append(result.getValue0());
                this.updateMatrix(result);
            }
        }
        Utils.saveFile(Globals.output + "reads.fasta", sb.toString());
        this.printMatrix();
    }
}