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
package ch.ethz.bsse.indelfixer.minimal.processing.parallel;

import ch.ethz.bsse.indelfixer.stored.Genome;
import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.Kmer;
import ch.ethz.bsse.indelfixer.stored.Read;
import ch.ethz.bsse.indelfixer.stored.SequenceEntry;
import ch.ethz.bsse.indelfixer.utils.StatusUpdate;
import ch.ethz.bsse.indelfixer.utils.Utils;
import jaligner.Alignment;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.util.SequenceParser;
import jaligner.util.SequenceParserException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class FutureSequence implements Callable<Pair<String, Map<Integer, Map<Integer, Integer>>>> {

    private SequenceEntry watsonTriple;
    private Genome[] genome;
    private Matrix matrix;
    private int number;
    private Map<Integer, Map<Integer, Integer>> substitutionsForward = new HashMap<>();

    public FutureSequence(SequenceEntry watson, int number) {
        this.watsonTriple = watson;
        this.genome = Globals.GENOMES;
        this.matrix = Globals.MATRIX;
        this.number = number;
        this.initSubs();
    }

    @Override
    public Pair<String, Map<Integer, Map<Integer, Integer>>> call() {
        StringBuilder sb = new StringBuilder();
        if (this.watsonTriple != null) {
            Read watsonRead = map(createRead(watsonTriple, false));
            Read watsonRevRead = map(createRead(watsonTriple, true));
            try {
                Read watson = align(watsonRead.getMaximumHits() > watsonRevRead.getMaximumHits() ? watsonRead : watsonRevRead, this.substitutionsForward);
                if (watson != null) {
                    sb.append(toString(watson));
                    StatusUpdate.processReads();
                    if (watson.getAlignedRead().length() > Globals.MIN_LENGTH_ALIGNED) {
                        return Pair.with(sb.toString(), substitutionsForward);
                    }
                }
            } catch (Exception e) {
                System.err.println(e);
                Utils.error();
                System.exit(0);
            }
        }
        return null;
    }

    private void initSubs() {
        for (int v = 0; v < 6; v++) {
            substitutionsForward.put(v, new HashMap<Integer, Integer>());
            for (int b = 0; b < 6; b++) {
                substitutionsForward.get(v).put(b, 0);
            }
        }
    }

    private Read createRead(SequenceEntry entry, boolean reverse) {
        Read r;
        if (reverse) {
            r = new Read(Utils.reverseComplement(entry.sequence.toUpperCase()), 0, true);
        } else {
            r = new Read(entry.sequence.toUpperCase(), 0, false);
        }
        if (entry.tag != null) {
            r.setDescription(entry.tag);
        }
        if (entry.pairedNumber != 0) {
            r.setMatePair(entry.pairedNumber);
        }
        if (entry.quality != null) {
            r.setQuality(entry.quality);
        }
        return r;
    }

    private Read map(Read r) {
        List<Kmer> kmers = r.getKmers();
        for (int x = 0; x < this.genome.length; x++) {
            for (Kmer kmer : kmers) {
                if (this.genome[x].getKmerMap().containsKey(kmer.getSequence())) {
                    for (int pos : this.genome[x].getKmerMap().get(kmer.getSequence())) {
                        r.addHit(pos, pos + Globals.KMER_LENGTH, x);
                    }
                }
            }
        }

        r.findMaxHitRegion();
        int y = r.getBestFittingGenome();
        if (y != -1) {
            if (r.getRegion(y)[0] - r.getRead().length() < 0) {
                r.setBegin(0);
            } else {
                r.setBegin(r.getRegion(y)[0] - r.getRead().length());
            }
            if (r.getRegion(y)[1] + r.getRead().length() > genome[y].getSequence().length()) {
                r.setEnd(genome[y].getSequence().length());
            } else {
                r.setEnd(r.getRegion(y)[1] + r.getRead().length());
            }
        }
        return r;
    }

    private Read align(Read r, Map<Integer, Map<Integer, Integer>> sub) throws SequenceParserException {
        
        Alignment align;
        if (r.getBestFittingGenome() == -1 || r.getEnd() < 0) {
            return null;
        } else {
            if (r.getBegin() <= 0) {
                align = SmithWatermanGotoh.align(
                        SequenceParser.parse(Globals.GENOME_SEQUENCES[r.getBestFittingGenome()]),
                        SequenceParser.parse(r.getRead()),
                        matrix, Globals.GOP, Globals.GEX);
                r.setBegin(r.getBegin() + 1);
            } else {
                int readEnd = r.getEnd() >= Globals.GENOME_SEQUENCES[r.getBestFittingGenome()].length() ? Globals.GENOME_SEQUENCES[r.getBestFittingGenome()].length() : r.getEnd();
                align = SmithWatermanGotoh.align(
                        SequenceParser.parse(Globals.GENOME_SEQUENCES[r.getBestFittingGenome()].substring(r.getBegin() - 1, readEnd)),
                        SequenceParser.parse(r.getRead()),
                        matrix, Globals.GOP, Globals.GEX);
            }
        }
        StringBuilder sb = new StringBuilder();
        r.setBegin(align.getStart1() + r.getBegin());
        char[] m = align.getMarkupLine();
        char[] c = align.getSequence2();
        char[] g = align.getSequence1();
        double L = m.length;
        int ins = 0;
        int del = 0;
        int subs = 0;
        char[] cigar = new char[m.length];
        for (int j = 0; j < L; j++) {
            char currentConsensus = '*';
            try {
                if (m[j] == '|') {
                    currentConsensus = c[j];
                    cigar[j] = 'M';
                } else if (m[j] == ' ') {
                    if (isGAP(c[j]) && isGAP(g[j])) {
                        cigar[j] = 'D';
                        currentConsensus = '-';
                    } else if (isGAP(c[j])) {
                        if (Globals.ADJUST) {
                            cigar[j] = 'M';
                            currentConsensus = g[j];
                        } else {
                            del++;
                            currentConsensus = '-';
                            cigar[j] = 'D';
                        }
                    } else if (isGAP(g[j])) {
                        if (Globals.ADJUST) {
                            ins++;
                        } else {
                            cigar[j] = 'I';
                            currentConsensus = c[j];
                        }
                    }
                } else if (m[j] == '.') {
                    subs++;
                    if (isGAP(g[j])) {
                        if (!isGAP(c[j])) {
                            if (Globals.ADJUST) {
                                ins++;
                                currentConsensus = g[j];
                            } else {
                                currentConsensus = c[j];
                            }
                            if (currentConsensus == '-') {
                                cigar[j] = 'D';
                            } else {
                                cigar[j] = 'M';
                            }
                        } else {
                            System.out.println("Contact program creator");
                            System.exit(0);
                        }
                    } else if (isGAP(c[j])) {
                        if (Globals.ADJUST) {
                            del++;
                            cigar[j] = 'M';
                            currentConsensus = g[j];
                        } else {
                            cigar[j] = 'D';
                            currentConsensus = c[j];
                        }
                    } else {
                        cigar[j] = 'M';
                        currentConsensus = c[j];
                    }
                }

                if (c[j] == 'N') {
                    currentConsensus = g[j];
                    if (isGAP(g[j])) {
                        cigar[j] = 'I';
                    } else {
                        cigar[j] = 'M';
                    }

                }
                if (cigar[j] != 0) {
                    if (cigar[j] == 'M') {
                        sb.append(currentConsensus);
                    }
                    sub.get(convert(currentConsensus)).put(convert(g[j]), sub.get(convert(currentConsensus)).get(convert(g[j])) + 1);
//                    if (cigar[j] == 'M' || cigar[j] == 'I') {
//                        sb.append(currentConsensus);
//                    }
//                    sub.get(convert(currentConsensus)).put(convert(g[j]), sub.get(convert(currentConsensus)).get(convert(g[j])) + 1);
                }
            } catch (Exception e) {
                System.err.println(e);
            }
        }
        if (((del / L) > Globals.MAX_DEL) || ((ins / L) > Globals.MAX_INS) || ((subs / L) > Globals.MAX_SUB)) {
            return null;
        }

//        StringBuilder cigarSB = new StringBuilder();
//        char prev = 'x';
//        for (char h : cigar) {
//            if (h != 0) {
//                if (prev == 'x') {
//                    prev = h;
//                } else {
//                    if (h == prev) {
//                    } else {
//                        if (prev == 'I' && h == 'D') {
//                            return null;
//                        } else if (prev == 'D' && h == 'I') {
//                            return null;
//                        }
//                        prev = h;
//                    }
//                }
//            }
//        }

        int length = 0;
        for (char h : cigar) {
            if (h != 0) {
                if (h == 'M') {
                    length++;
                }
            }
        }
        if (sb.toString().length() != length) {
            System.out.println("kick out");
            return null;
        }

        r.setCigars(cigar);
        r.setAlignedRead(sb.toString());
        r.setEnd(r.getBegin() + sb.length());
        return r;
    }

    private boolean isGAP(char c) {
        return c == '-' || c == 'N';
    }

    private int convert(char c) {
        switch (c) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            case '-':
                return 4;
            case 'N':
                return 5;
            default:
                return 6;
        }
    }

    private StringBuilder toString(Read read) {
        StringBuilder sb = new StringBuilder();
        sb.append("R");
        if (read.getDescription() != null) {
            sb.append(read.getDescription());
        } else {
            sb.append(this.number);
        }
        sb.append("\t0\t").append(Globals.GENOMES[read.getBestFittingGenome()].getHeader());
        sb.append("\t").append(read.getBegin());
        sb.append("\t").append("255");
        sb.append("\t").append(read.getCigars());
        sb.append("\t").append("*");
        sb.append("\t").append("0");
        sb.append("\t").append("0");
        sb.append("\t").append(read.getAlignedRead());
        sb.append("\t").append("*");
//        for (int i = 0; i < read.getAlignedRead().length(); i++) {
//            sb.append("I");
//        }
        sb.append("\n");
//        sb.append(">READ").append(this.number).append("_").append(read.getBegin()).append("-").append(read.getEnd());
//        if (read.getDescription() != null) {
//            sb.append("|").append(read.getDescription()).append("/").append(read.getMatePair());
//        }
//        sb.append("\n");
//        sb.append(read.getAlignedRead()).append("\n");
        return sb;
    }
}
