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

import ch.ethz.bsse.indelfixer.minimal.processing.ProcessingSFFSingle;
import ch.ethz.bsse.indelfixer.minimal.processing.ProcessingFastaSingle;
import ch.ethz.bsse.indelfixer.minimal.processing.ProcessingIlluminaPaired;
import ch.ethz.bsse.indelfixer.sffParser.SFFParsing;
import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.utils.FastaParser;
import ch.ethz.bsse.indelfixer.utils.Utils;
import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.logging.Handler;
import java.util.logging.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
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
    private int overlap = 10;
    @Option(name = "-k")
    private int kmerLength = 21;
    @Option(name = "-fill")
    private boolean fill = false;
    @Option(name = "-t")
    private int threshold = 30;
    @Option(name = "-step")
    private int step = 1000;
    @Option(name = "-r")
    private String regions;
    @Option(name = "--flat")
    private boolean flat;
    @Option(name = "--count")
    private boolean count;
    @Option(name = "-l")
    private int minlength = 0;

    public static void main(String[] args) throws IOException {
        new Start().doMain(args);
    }

    public void doMain(String[] args) {
        CmdLineParser parser = new CmdLineParser(this);

        parser.setUsageWidth(80);
        try {
            parser.parseArgument(args);
            this.checkOutput();

            if (this.input == null && this.genome == null) {
                throw new CmdLineException("");
            }
            if (this.flat) {
                flatAndSave();
            } else {
                this.setGlobals();
                Genome[] genomes = parseGenome(this.genome);
                if (this.regions != null) {
                    Globals.RS = this.splitRegion();
                    this.cutGenomes(genomes);
                }

                if (this.inputReverse != null) {
                    new ProcessingIlluminaPaired(this.input, this.inputReverse);
                } else if (Utils.isFastaFormat(this.input)) {
                    new ProcessingFastaSingle(this.input);
                } else {
                    new ProcessingSFFSingle(SFFParsing.parse(this.input));
                }
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("java -jar InDelFixer.jar options...\n");
            System.err.println(" ------------------------");
            System.err.println(" === GENERAL options ===");
            System.err.println("  -o PATH\t\t: Path to the output directory (default: current directory)");
            System.err.println("  -i PATH\t\t: Path to the NGS input file (FASTA, FASTQ or SFF format) [REQUIRED]");
            System.err.println("  -ir PATH\t\t: Path to the second paired end file (FASTQ) [REQUIRED if first file is also fastq]");
            System.err.println("  -g PATH\t\t: Path to the reference genomes file (FASTA format) [REQUIRED]");
            System.err.println("  -r interval\t\t: Region on the reference genome (i.e. 342-944)");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" Example: java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta -r 342-944");
            System.err.println(" Example: java -jar InDelFixer.jar -i lib12_R1.fastq -ir lib12_R2.fastq -g referenceGenomes.fasta -r 342-944");
        }
    }

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

    private int[][] splitRegion() {
        String[] r = regions.split(",");
        int[][] rs = new int[r.length][2];
        int i = 0;
        for (String s : r) {
            String[] ss = s.split("-");
            rs[i][0] = Integer.parseInt(ss[0]) - 1;
            rs[i++][1] = Integer.parseInt(ss[1]);
        }
        return rs;
    }

    private void cutGenomes(Genome[] genomes) {
        int[][] rs = Globals.RS;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rs.length; i++) {
            for (Genome g : genomes) {
                sb.append(g.getHeader()).append("\n");
                sb.append(g.getSequence().substring(rs[i][0] - 1, rs[i][1] - 1)).append("\n");
            }
            Utils.saveFile(Globals.output + "ref_" + (rs[i][0]) + "-" + (rs[i][1] - 1) + ".fasta", sb.toString());
        }
    }

    private Genome[] parseGenome(String genomePath) {
        Map<String, String> haps = FastaParser.parseHaplotypeFile(genomePath);
        Genome[] gs = new Genome[haps.size()];
        int i = 0;
        for (Map.Entry<String, String> hap : haps.entrySet()) {
            gs[i++] = new Genome(hap);
        }
        Globals.GENOMES = gs;
        Globals.GENOME_COUNT = gs.length;
        Globals.GENOME_SEQUENCES = haps.keySet().toArray(new String[gs.length]);
        Globals.GENOME_LENGTH = Globals.GENOME_SEQUENCES[0].length();
        return gs;
    }

    private void setGlobals() {
        Globals.MIN_LENGTH = minlength;
        Globals.FILL = this.fill;
        Globals.STEPSIZE = this.step;
        Globals.THRESHOLD = this.threshold;
        Globals.KMER_OVERLAP = this.overlap;
        Globals.KMER_LENGTH = this.kmerLength;
    }

    {
        Logger rootLogger = Logger.getLogger("");
        Handler[] handlers = rootLogger.getHandlers();
        if (handlers.length > 0) {
            rootLogger.removeHandler(handlers[0]);
        }
    }

    private void flatAndSave() {
        String[] far = FastaParser.parseFarFile(input);
        StringBuilder sb = new StringBuilder();
        int x = 0;
        for (int i = 0; i < far.length; i++) {
            sb.append(">READ-").append(i).append("\n" + far[i] + "\n");
            if (i % 100 == 0) {
                Utils.saveFile(output + "flat-" + x++ + ".fasta", sb.toString());
                sb.setLength(0);
            }
        }
        if (sb.length() > 0) {
            Utils.saveFile(output + "flat-" + x + ".fasta", sb.toString());
        }
    }
}
