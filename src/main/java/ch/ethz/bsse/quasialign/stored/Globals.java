/**
 * Copyright (c) 2011-2012 Armin Töpfer
 *
 * This file is part of QuasiAlign.
 *
 * QuasiAlign is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * QuasiAlign is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * QuasiAlign. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.quasialign.stored;

import ch.ethz.bsse.quasialign.stored.InformationHolder;
import ch.ethz.bsse.quasialign.utils.StatusUpdate;
import java.io.File;
import java.util.concurrent.ForkJoinPool;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    public static boolean FILL;
    public static boolean SAVE;
    public static boolean HEURISTIC;
    public static int STEPSIZE = 1000;
    public static int KMER_OVERLAP = 10;
    public static int KMER_LENGTH = 50;
    public static final ForkJoinPool fjPool = new ForkJoinPool();
    public static int N = 0;
    public static double PERCENTAGE_PROCESSING_READS = 0;
    public static double PERCENTAGE_EXTRACTING_READS = 0;
    public static double UNIDENTICAL = 0;
    public static InformationHolder INFO_HOLDER;
    public static String output = System.getProperty("user.dir") + File.separator;

    public static void printPercentageExtractingReads() {
        PERCENTAGE_EXTRACTING_READS += 100d / N;
        StatusUpdate.print("Extracting reads:\t\t", PERCENTAGE_EXTRACTING_READS);
    }
    public static void printPercentageProcessingReads() {
        PERCENTAGE_PROCESSING_READS += 100d / N;
        StatusUpdate.print("Processing reads:\t\t", PERCENTAGE_PROCESSING_READS);
    }
    public static double PERCENTAGE_ALIGNING_READS = 0;

    public static void printPercentageAligningReads() {
        PERCENTAGE_ALIGNING_READS += 100d / UNIDENTICAL;
        if (PERCENTAGE_ALIGNING_READS > 100) {
            System.out.println("");
        }
        StatusUpdate.print("Aligning reads:\t\t", PERCENTAGE_ALIGNING_READS);
    }
    public static double PERCENTAGE_MAPPING_READS = 0;

    public static void printPercentageMappingReads() {
        PERCENTAGE_MAPPING_READS += 100d / UNIDENTICAL;
        StatusUpdate.print("Mapping non-identical reads:\t", PERCENTAGE_MAPPING_READS);
    }
    public static double PERCENTAGE_INDEX_READS = 0;

    public static void printPercentageIndexReads() {
        PERCENTAGE_INDEX_READS += 100d / N;
        StatusUpdate.print("Finding indizes for kmers:\t", PERCENTAGE_INDEX_READS);
    }
}
