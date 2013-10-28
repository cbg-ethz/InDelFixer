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
package ch.ethz.bsse.indelfixer.stored;

import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import java.io.File;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.ForkJoinPool;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    public static boolean NO_HASHING;
    public static boolean SENSITIVE;
    public static boolean CONSENSUS;
    public static int READCOUNTER = 0;
    public static int STEPS;
    public static List<SequenceEntry> TRASH_READS = new CopyOnWriteArrayList<>();
    public static boolean RM_DEL;
    public static boolean REFINE;
    public static boolean FILTER;
    public static boolean FLAT;
    public static float GOP;
    public static float GEX;
    public static double MAX_DEL;
    public static double MAX_INS;
    public static double MAX_SUB;
    public static boolean KEEP;
    public static List<Read> READS;
    public static int MAX_CONSECUTIVE_DEL;
    public static int MIN_LENGTH_ALIGNED;
    public static int MIN_LENGTH;
    public static int GENOME_LENGTH;
    public static int THRESHOLD;
    public static String[] GENOME_SEQUENCES;
    public static boolean ADJUST;
    public static boolean SAVE;
    public static int GENOME_COUNT;
    public static Genome[] GENOMES;
    public static int CUT = 0;
    public static int KMER_OVERLAP = 10;
    public static int KMER_LENGTH = 1;
    public static final ForkJoinPool fjPool = new ForkJoinPool();
    public static int N = 0;
    public static int maxN = 0;
    public static double PERCENTAGE_PROCESSING_READS = 0;
    public static double PERCENTAGE_EXTRACTING_READS = 0;
    public static double UNIDENTICAL = 0;
    public static String output = System.getProperty("user.dir") + File.separator;
    public static Matrix MATRIX = loadMatrix();
    public static int[] RS;
    public static String[] ADAPTER;
    public static double PLURALITY;

    /**
     *
     * @return
     */
    public static Matrix loadMatrix() {
        Matrix m = null;
        try {
            m = MatrixLoader.load("EDNAFULL");
        } catch (MatrixLoaderException ex) {
            System.err.println("Matrix error loading");
        }
        return m;
    }

    public static void printPercentageExtractingReads() {
        PERCENTAGE_EXTRACTING_READS += 100d / N;
        StatusUpdate.print("Extracting reads:\t", PERCENTAGE_EXTRACTING_READS);
    }

    public static void printPercentageProcessingReads() {
        PERCENTAGE_PROCESSING_READS += 100d / N;
        StatusUpdate.print("Processing reads:\t", PERCENTAGE_PROCESSING_READS);
    }
}
