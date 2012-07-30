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
package ch.ethz.bsse.quasialign;

import ch.ethz.bsse.quasialign.sffParser.SFFParsing;
import ch.ethz.bsse.quasialign.stored.Globals;
import ch.ethz.bsse.quasialign.stored.Read;
import ch.ethz.bsse.quasialign.utils.Cutter;
import ch.ethz.bsse.quasialign.utils.FastaParser;
import ch.ethz.bsse.quasialign.utils.StatusUpdate;
import ch.ethz.bsse.quasialign.utils.Utils;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Startup {

    @Option(name = "-i")
    private String input;
    @Option(name = "-g")
    private String genome;
//    @Option(name = "-verbose")
//    private boolean verbose;
    @Option(name = "-o", usage = "Path to the output directory (default: current directory)", metaVar = "PATH")
    private String output;
    @Option(name = "-v")
    private int overlap = 1;
    @Option(name = "-save")
    private boolean save = false;
    @Option(name = "-regions")
    private String regions;
    @Option(name = "-heuristic")
    private boolean heuristic = false;
    @Option(name = "--fill")
    private boolean fill = false;

    public static void main(String[] args) throws IOException {
        new Startup().doMain(args);
    }

    public void doMain(String[] args) {
        CmdLineParser parser = new CmdLineParser(this);

        parser.setUsageWidth(80);
        try {
            parser.parseArgument(args);

            if (output == null) {
                this.output = System.getProperty("user.dir") + File.separator;
            } else {
                Globals.output = this.output;
            }
            if (!new File(this.output).exists()) {
                new File(this.output).mkdirs();

            }

            if (this.input == null && this.genome == null) {
                throw new CmdLineException("");
            }
            Globals.KMER_OVERLAP = this.overlap;
            Globals.FILL = this.fill;
            Globals.SAVE = this.save;
            Globals.HEURISTIC = this.heuristic;
            Workflow workflow;
            if (Utils.isFastaFormat(this.input)) {
                workflow = new Workflow(this.genome, FastaParser.parseFarFile(input));
            } else {
                workflow = new Workflow(this.genome, SFFParsing.parse(this.input));
            }

            if (regions == null) {
                StringBuilder sb = new StringBuilder();
                int i = 0;
                for (Read read : Globals.INFO_HOLDER.getReads()) {
                    if (read.isAligned()) {
                        sb.append(">READ").append(i++).append("_").append(read.getBegin()).append("-").append(read.getEnd()).append("\n");
                        sb.append(read.getAlignedRead()).append("\n");
                    }
                }
                Utils.saveFile(Globals.output + "region.fasta", sb.toString());
            } else {
                String[] r = regions.split(",");
                int[][] rs = new int[r.length][2];
                int i = 0;
                for (String s : r) {
                    String[] ss = s.split("-");
                    rs[i][0] = Integer.parseInt(ss[0]) - 1;
                    rs[i++][1] = Integer.parseInt(ss[1]);
                }
                saveReads(rs);
//                Map<String, List<Read>>[] windows = new HashMap[sbs.length];
//                for (int j = 0; j < sbs.length; j++) {
//                    windows[j] = new HashMap<>();
//                }
//                for (int j = 0; j < rs.length; j++) {
//                    for (int m = rs[j][0] + 300;; m += 100) {
//                        int k = m - 300;
//                        if (m + 100 > rs[j][1]) {
//                            m = rs[j][1];
//                        }
//                        windows[j].put(k + "-" + m, new ArrayList());
//                        if (m == rs[j][1]) {
//                            break;
//                        }
//                    }
//                }
//                for (Read read : Globals.INFO_HOLDER.getReads()) {
//                    if (read.isAligned()) {
//                        for (int j = 0; j < rs.length; j++) {
//                            for (int m = rs[j][0] + 300;; m += 100) {
//                                int k = m - 300;
//                                if (m + 100 > rs[j][1]) {
//                                    m = rs[j][1];
//                                }
//                                final int B = read.getBegin();
//                                Read readNew = null;
//                                final int E = read.getEnd();
//                                if (B < k && E > m) {
//                                    readNew = new Read(k, m, read.getAlignedRead().substring(k - B, m - B), Arrays.copyOfRange(read.getQuality(), k - B, m - B));
//                                } else if (B < k && E > k && E < m) {
//                                    readNew = new Read(k, E, read.getAlignedRead().substring(k - B, E - B), Arrays.copyOfRange(read.getQuality(), k - B, E - B));
//                                } else if (B > k && B < m && E > m) {
//                                    readNew = new Read(B, m, read.getAlignedRead().substring(0, m - B), Arrays.copyOfRange(read.getQuality(), 0, m - B));
//                                }
//                                if (readNew != null) {
//                                    windows[j].get(k + "-" + m).add(readNew);
//                                }
//                                if (m == rs[j][1]) {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//                for (int j = 0; j < sbs.length; j++) {
//                    new File(Globals.output + "window_" + (rs[j][0] + 1) + "-" + rs[j][1]).mkdirs();
//                    for (String w : windows[j].keySet()) {
//                        StringBuilder sbw = new StringBuilder();
//                        for (Read rw : windows[j].get(w)) {
//                            sbw.append(">READ").append(i++).append("_").append(rw.getBegin()).append("-").append(rw.getEnd()).append("\n");
//                            sbw.append(rw.getAlignedRead()).append("\n");
//                        }
//                        Utils.saveFile(Globals.output + "window_" + (rs[j][0] + 1) + "-" + rs[j][1] + File.separator + "w_" + w + ".fasta", sbw.toString());
//
//                    }
//                }
                StatusUpdate.println("Saving coverage");
                workflow.saveCoverage();
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("java -jar QuasiAlign.jar options...\n");
            System.err.println(" ------------------------");
            System.err.println(" === GENERAL options ===");
            System.err.println("  -o PATH\t\t: Path to the output directory (default: current directory)");
            System.err.println("  -i PATH\t\t: Path to the NGS input file (FASTA or SFF format) [REQUIRED]");
            System.err.println("  -g PATH\t\t: Path to the reference genome file (FASTA format) [REQUIRED]");
            System.err.println("  -r intervals\t\t: Regions on the reference genome (i.e. 342-944,1239-3889,9953-12293)");
//            System.err.println("  -verbose\t\t: Print debug information");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" Example: java -jar QuasiAlign.jar -i libCase102.sff -g referenceGenome.fasta -r 342-944,1239-3889");
        }


//        StatusUpdate.println("Init");

//        new Workflow(args[0]+"/Dropbox/QuasiAlign/src/main/resources/HIV1.fas", args[0]+"/Dropbox/QuasiAlign/src/main/resources/haplotypes/SRR069887.far");
    }

    private void saveReads(int[][] rs) {
        List<Read>[] maps = new ArrayList[rs.length];
        StringBuilder[] sbs = new StringBuilder[rs.length];
        for (int j = 0; j < sbs.length; j++) {
            sbs[j] = new StringBuilder();
            maps[j] = new ArrayList<>();
        }
//        StringBuilder[] shorah = new StringBuilder[rs.length];
//        for (int j = 0; j < shorah.length; j++) {
//            shorah[j] = new StringBuilder();
//        }
        int cuts = 0;
        int i = 0;
        for (Read read : Globals.INFO_HOLDER.getReads()) {
            for (int j = 0; j < rs.length; j++) {
                if (read.getEnd() > rs[j][0] && read.getBegin() < rs[j][1] && read.isAligned()) {
                    if (read.getEnd() > rs[j][0] && read.getBegin() < rs[j][0]) {
                        String s = read.getAlignedRead().substring(rs[j][0] - read.getBegin());
                        if (s.length() < 150) {
                            continue;
                        }
//                                read.setQuality(Arrays.copyOfRange(read.getQuality(), rs[j][0] - read.getBegin(), read.getQuality().length));
                        read.setAlignedRead(s);
                        read.setBegin(rs[j][0]);
                        cuts++;
                    }
                    if (read.getEnd() > rs[j][1] && read.getBegin() < rs[j][1]) {
                        String s = read.getAlignedRead().substring(0, read.getEnd() - rs[j][1]);
                        if (s.length() < 150) {
                            continue;
                        }
//                                read.setQuality(Arrays.copyOfRange(read.getQuality(), 0, read.getEnd() - rs[j][1]));
                        read.setAlignedRead(s);
                        read.setEnd(rs[j][1]);
                        cuts++;
                    }
                    maps[j].add(read);
                    sbs[j].append(">READ").append(i++).append("_").append(read.getBegin()).append("-").append(read.getEnd()).append("\n");
                    sbs[j].append(read.getAlignedRead()).append("\n");
//                    shorah[j].append(">READ").append(i++).append("_").append(read.getBegin()).append("-").append(read.getEnd()).append("\n");
//                    for (int k = rs[j][0]; k < read.getBegin(); k++) {
//                        shorah[j].append("-");
//                    }
//                    shorah[j].append(read.getAlignedRead()).append("\n");
//                    for (int k = read.getEnd(); k < rs[j][1]; k++) {
//                        shorah[j].append("-");
//                    }
//                            sbs[j].append("+\n");
//                            for (int k = 0; k < read.getQuality().length; k++) {
//                                sbs[j].append((char) k);
//                            }
//                            sbs[j].append("\n");
                    break;
                }
            }
        }
        StatusUpdate.println("Cuts:\t\t\t" + cuts);
        for (int j = 0; j < sbs.length; j++) {
            Utils.saveFile(Globals.output + "region_" + (rs[j][0] + 1) + "-" + rs[j][1] + ".fasta", sbs[j].toString());
//            Utils.saveFile(Globals.output + "shorah_region_" + (rs[j][0] + 1) + "-" + rs[j][1] + ".fasta", shorah[j].toString());
            Cutter.cut(this.genome, Globals.output, rs[j][0], rs[j][1]);
        }
    }
}
