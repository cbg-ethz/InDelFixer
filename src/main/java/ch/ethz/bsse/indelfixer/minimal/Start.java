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
package ch.ethz.bsse.indelfixer.minimal;

import ch.ethz.bsse.indelfixer.minimal.processing.ProcessingFastaSingle;
import ch.ethz.bsse.indelfixer.minimal.processing.ProcessingIlluminaPaired;
import ch.ethz.bsse.indelfixer.minimal.processing.ProcessingIlluminaSingle;
import ch.ethz.bsse.indelfixer.minimal.processing.ProcessingSFFSingle;
import ch.ethz.bsse.indelfixer.sffParser.SFFParsing;
import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.utils.FastaParser;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import ch.ethz.bsse.indelfixer.utils.Utils;
import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Handler;
import java.util.logging.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * Entry point.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Start {

    @Option(name = "-i")
    private String input;
    @Option(name = "-ir")
    private String inputReverse;
    @Option(name = "-g")
    private String genome;
    @Option(name = "-o", usage = "Path to the output directory (default: current directory)", metaVar = "PATH")
    private String output;
    @Option(name = "-v")
    private int overlap = 2;
    @Option(name = "-k")
    private int kmerLength = 10;
    @Option(name = "-adjust")
    private boolean adjust;
    @Option(name = "-t")
    private int threshold = 30;
    @Option(name = "-step")
    private int step = 10;
    @Option(name = "-r")
    private String regions;
    @Option(name = "-flat")
    private boolean flat;
    @Option(name = "-l")
    private int minlength = 0;
    @Option(name = "-la")
    private int minlengthAligned = 0;
    @Option(name = "-sub")
    private double sub = 1;
    @Option(name = "-del")
    private double del = 1;
    @Option(name = "-ins")
    private double ins = 1;
    @Option(name = "-gop")
    private float gop = 30;
    @Option(name = "-gex")
    private float gex = 3;
//    @Option(name = "-pacbio")
//    private boolean pacbio;
//    @Option(name = "-illumina")
//    private boolean illumina;
//    @Option(name = "-454")
//    private boolean roche;
    @Option(name = "-cut")
    private int cut;
    @Option(name = "-header")
    private String header;
    @Option(name = "--filter")
    private boolean filter;
    @Option(name = "--freq")
    private boolean freq;
    @Option(name = "-refine")
    private int refine;
    @Option(name = "--version")
    private boolean version;
    @Option(name = "-rmDel")
    private boolean rmDel;
    @Option(name = "-consensus")
    private boolean consensus;
    @Option(name = "-N")
    private int N = 3;
    @Option(name = "-sensitive")
    private boolean sensitive;
    @Option(name = "-maxDel")
    private int maxDel = -1;
    @Option(name = "-noHashing")
    private boolean noHashing;
    @Option(name = "-fix")
    private boolean fix;

    /**
     * Remove logging of jaligner.
     */
    {
        Logger rootLogger = Logger.getLogger("");
        Handler[] handlers = rootLogger.getHandlers();
        if (handlers.length > 0) {
            rootLogger.removeHandler(handlers[0]);
        }
    }

    /**
     * Forwards command-line parameters.
     *
     * @param args command-line parameters
     * @throws IOException If parameters are wrong
     */
    public static void main(String[] args) throws IOException {
        new Start().doMain(args);
        System.exit(0);
    }

    /**
     * Main.
     *
     * @param args command-line parameters
     */
    public void doMain(String[] args) {
        try {
            CmdLineParser parser = new CmdLineParser(this);

            parser.setUsageWidth(80);
            try {
                parser.parseArgument(args);
                if (this.version) {
                    System.out.println("InDelFixer version: " + Start.class.getPackage().getImplementationVersion());
                    return;
                }
                this.checkOutput();

                if (this.input == null && this.genome == null) {
                    throw new CmdLineException("");
                }
                if (this.fix) {
                    if (this.refine == 0) {
                        refine = 1;
                    }
                }
                this.setGlobals();
                Genome[] genomes = parseGenome(this.genome);
                if (this.regions != null) {
                    Globals.RS = this.splitRegion();
                    genomes = parseGenome(this.cutGenomes(genomes));
                }
                compute(genomes);
                if (this.fix) {
                    Globals.ADJUST = true;
                }
                if (this.refine > 0) {
                    this.genome = this.output + "consensus.fasta";
                    genomes = parseGenome(this.genome);
                    for (int i = 0; i < this.refine; i++) {
                        StatusUpdate.readCount = 0;
                        StatusUpdate.unmappedCount = 0;
                        StatusUpdate.tooSmallCount = 0;
                        StatusUpdate.alignCount1 = 0;
                        StatusUpdate.alignCount2 = 0;
                        StatusUpdate.alignCount3 = 0;
                        compute(genomes);
                    }
                }
            } catch (CmdLineException e) {
                System.err.println(e.getMessage());
                System.err.println("InDelFixer version: " + Start.class.getPackage().getImplementationVersion());
                System.err.println("Get latest version from http://bit.ly/indelfixer");
                System.err.println("");
                System.err.println("USAGE: java -jar InDelFixer.jar options...\n");
                System.err.println(" ------------------------");
                System.err.println(" === GENERAL options ===");
                System.err.println("  -o PATH\t\t: Path to the output directory (default: current directory)");
                System.err.println("  -i PATH\t\t: Path to the NGS input file (FASTA, FASTQ or SFF format) [REQUIRED]");
                System.err.println("  -ir PATH\t\t: Path to the second paired end file (FASTQ) [ONLY REQUIRED if first file is also fastq]");
                System.err.println("  -g PATH\t\t: Path to the reference genomes file (FASTA format) [REQUIRED]");
                System.err.println("  -r interval\t\t: Region on the reference genome (i.e. 342-944)");
                System.err.println("  -k INT\t\t: Kmer size (default 10)");
                System.err.println("  -v INT\t\t: Kmer offset (default 2)");
                System.err.println("  -cut INT\t\t: Cut given number of bases (primer) from beginning of the read (default 0)");
                System.err.println("  -refine INT\t\t: Computes a consensus sequence from alignment and re-aligns against that.");
                System.err.println("\t\t\t  Refinement is repeated as many times as specified.");
                System.err.println("  -rmDel\t\t: Removes conserved gaps from consensus sequence during refinement");
                System.err.println("  -sensitive\t\t: More sensitive but slower alignment");
                System.err.println("  -fix\t\t\t: Fill frame-shift causing deletions with consensus sequence.");
                System.err.println("  -noHashing\t\t: No fast kmer-matching to find approximate mapping region. Please use with PacBio data!");
                System.err.println("");
                System.err.println(" === FILTER ===");
                System.err.println("  -l INT\t\t: Minimal read-length prior alignment (default 0)");
                System.err.println("  -la INT\t\t: Minimal read-length after alignment (default 0)");
                System.err.println("  -ins INT\t\t: The maximum precentage of insertions [range 0.0 - 1.0]\n");
                System.err.println("  -del INT\t\t: The maximum precentage of deletions [range 0.0 - 1.0]\n");
                System.err.println("  -sub INT\t\t: The maximum precentage of insertions [range 0.0 - 1.0]\n");
                System.err.println("  -maxDel INT\t\t: The maximum number of consecutive deletions allowed");
                System.err.println("");
                System.err.println(" === GAP costs ===");
                System.err.println("  -gop\t\t\t: Gap opening costs for Smith-Waterman (default 30)");
                System.err.println("  -gex\t\t\t: Gap extension costs for Smith-Waterman (default 3)");
//                System.err.println("");
//                System.err.println(" === GAP costs predefined ===");
//                System.err.println("  -454\t\t\t: 10 open / 1 extend");
//                System.err.println("  -illumina\t\t: 30 open / 3 extend");
//                System.err.println("  -pacbio\t\t: 10 open / 10 extend");
                System.err.println("");
                System.err.println(" ------------------------");
                System.err.println(" === EXAMPLES ===");
                System.err.println("  454/Roche\t\t: java -jar InDelFixer.jar -i libCase102.fastq -g referenceGenomes.fasta");
                System.err.println("  PacBio\t\t: java -jar InDelFixer.jar -i libCase102.ccs.fastq -g referenceGenomes.fasta");
                System.err.println("  Illumina\t\t: java -jar InDelFixer.jar -i libCase102_R1.fastq -ir libCase102_R2.fastq -g referenceGenomes.fasta");
                System.err.println(" ------------------------");
            }
        } catch (OutOfMemoryError e) {
            Utils.error();
            System.err.println("Please increase the heap space.");
        }
    }

    /**
     * Creates output directory if not existing.
     */
    private void checkOutput() {
        if (this.output == null) {
            this.output = System.getProperty("user.dir") + File.separator;
        } else {
            Globals.output = this.output;
        }
        if (!new File(this.output).exists()) {
            new File(this.output).mkdirs();
        }
    }

    /**
     * Splits input region in format x-y
     *
     * @return
     */
    private int[] splitRegion() {
        String[] split = regions.split("-");
        return new int[]{Integer.parseInt(split[0]), Integer.parseInt(split[1])};
    }

    /**
     * Cuts genomes into defined region.
     *
     * @param genomes of type Genome
     */
    private String cutGenomes(Genome[] genomes) {
        int[] rs = Globals.RS;
        StringBuilder sb = new StringBuilder();
        for (Genome g : genomes) {
            try {
                sb.append(">").append(g.getHeader()).append("\n");
                sb.append(g.getSequence().substring(rs[0] - 1, rs[1] - 1)).append("\n");
            } catch (Exception e) {
                System.err.println(e);
                Utils.error();
            }
        }
        String output = Globals.output + "ref_" + (rs[0]) + "-" + (rs[1] - 1) + ".fasta";
        Utils.saveFile(output, sb.toString());
        return output;
    }

    /**
     * Parses genomes from multiple fasta file.
     *
     * @param genomePath multiple fasta file path
     * @return Genome sequences wrapped into Genome class
     */
    private Genome[] parseGenome(String genomePath) {
        Map<String, String> haps = FastaParser.parseHaplotypeFile(genomePath);
        List<Genome> genomeList = new LinkedList<>();
        for (Map.Entry<String, String> hap : haps.entrySet()) {
            if (header == null || hap.getValue().startsWith(this.header)) {
                genomeList.add(new Genome(hap));
            }
        }
        Genome[] gs = genomeList.toArray(new Genome[genomeList.size()]);
        Globals.GENOMES = gs;
        Globals.GENOME_COUNT = gs.length;
        Globals.GENOME_SEQUENCES = haps.keySet().toArray(new String[gs.length]);
        Globals.GENOME_LENGTH = Globals.GENOME_SEQUENCES[0].length();
        return gs;
    }

    /**
     * Set global variables from command-line parameters
     */
    private void setGlobals() {
//        if (this.pacbio) {
//            Globals.GOP = 10;
//            Globals.GEX = 10;
//        } else if (this.illumina) {
//            Globals.GOP = 30;
//            Globals.GEX = 3;
//        } else if (this.roche) {
//            Globals.GOP = 10;
//            Globals.GEX = 1;
//        } else {
        Globals.GOP = this.gop;
        Globals.GEX = this.gex;
//        }
        Globals.MIN_LENGTH_ALIGNED = minlengthAligned;
        Globals.MIN_LENGTH = minlength;
        Globals.ADJUST = this.adjust;
        Globals.STEPS = this.step;
        Globals.THRESHOLD = this.threshold;
        Globals.KMER_OVERLAP = this.overlap;
        Globals.KMER_LENGTH = this.kmerLength;
        Globals.MAX_DEL = this.del;
        Globals.MAX_INS = this.ins;
        Globals.MAX_SUB = this.sub;
        Globals.FLAT = this.flat;
        Globals.CUT = this.cut;
        Globals.FILTER = this.filter;
        Globals.REFINE = this.refine > 0;
        Globals.RM_DEL = this.rmDel;
        Globals.maxN = this.N;
        Globals.CONSENSUS = this.consensus;
        Globals.SENSITIVE = this.sensitive;
        Globals.MAX_CONSECUTIVE_DEL = this.maxDel;
        Globals.NO_HASHING = this.noHashing;
    }

    /**
     * Flats multiple fasta file and splits it files with 100 sequences.
     */
    private void flatAndSave() {
        Map<String, String> far = FastaParser.parseHaplotypeFile(input);
        StringBuilder sb = new StringBuilder();
        int x = 0;
        int i = 0;
        for (Map.Entry<String, String> entry : far.entrySet()) {
            sb.append(entry.getValue()).append("\n").append(entry.getKey()).append("\n");
            if (i++ % 1000 == 0) {
                Utils.saveFile(output + "flat-" + x++ + ".fasta", sb.toString());
                sb.setLength(0);
            }
        }
        if (sb.length() > 0) {
            Utils.saveFile(output + "flat-" + x + ".fasta", sb.toString());
        }
    }

    public static byte[] splitReadIntoByteArray(String read) {
        byte[] Rs = new byte[read.length()];
        char[] readSplit = read.toCharArray();
        int length = readSplit.length;
        for (int i = 0; i < length; i++) {
            switch ((short) readSplit[i]) {
                case 65:
                    Rs[i] = 0;
                    break;
                case 67:
                    Rs[i] = 1;
                    break;
                case 71:
                    Rs[i] = 2;
                    break;
                case 84:
                    Rs[i] = 3;
                    break;
                case 45:
                    Rs[i] = 4;
                    break;
                default:
                    break;
            }
        }
        return Rs;
    }

    public static String shortenSmall(double value) {
        String s;
        if (value < 1e-16) {
            s = "0      ";
        } else if (value == 1.0) {
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

    private boolean compute(Genome[] genomes) throws CmdLineException, IllegalStateException {
        if (this.freq) {
            double[][] allels = new double[genomes[0].getSequence().length()][5];
            for (Genome g : genomes) {
                byte[] a = splitReadIntoByteArray(g.getSequence());
                for (int j = 0; j < g.getSequence().length(); j++) {
                    allels[j][a[j]] += 1d / genomes.length;
                }
            }
            StringBuilder sb = new StringBuilder();
            sb.append("\tA\tC\tG\tT\t-\n");
            for (int j = 0; j < allels.length; j++) {
                sb.append(j);
                double sum = 0d;
                for (int v = 0; v < 5; v++) {
                    sum += allels[j][v] > 1e-5 ? allels[j][v] : 0;;
                }
                for (int v = 0; v < 5; v++) {
                    sb.append("\t").append(shortenSmall(allels[j][v] / sum));
                }
                sb.append("\n");
            }
            Utils.saveFile("Genome_allels.txt", sb.toString());
            return true;
        }
        for (Genome g : genomes) {
            g.split();
        }
        if (!new File(this.input).exists()) {
            throw new CmdLineException("Input file does not exist");
        }
        File readsSam = new File(this.output + "reads.sam");
        if (readsSam.exists()) {
            readsSam.delete();
        }
        if (this.inputReverse != null) {
            if (!new File(this.inputReverse).exists()) {
                throw new CmdLineException("Input reverse file does not exist");
            }
            new ProcessingIlluminaPaired(this.input, this.inputReverse);
        } else if (Utils.isFastaGlobalMatePairFormat(this.input)) {
            new ProcessingIlluminaSingle(this.input);
        } else if (Utils.isFastaFormat(this.input)) {
            new ProcessingFastaSingle(this.input);
        } else {
            new ProcessingSFFSingle(SFFParsing.parse(this.input));
        }
        return false;
    }
}
