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
package ch.ethz.bsse.indelfixer.utils;

import ch.ethz.bsse.indelfixer.stored.Globals;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class FastaParser {

    public static String[] parseFarFile(String location) {
        List<String> readList = new LinkedList<>();

        boolean cut = false;
        boolean first = true;
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            readList.add(sb.toString());
                            sb.setLength(0);
                        }
                    } else {
                        if (first) {
                            boolean upperCase = false;
                            boolean lowerCase = false;
                            for (char c : strLine.toCharArray()) {
                                if (!lowerCase) {
                                    if (Character.isLowerCase(c)) {
                                        lowerCase = true;
                                        continue;
                                    }
                                }
                                if (!upperCase) {
                                    if (Character.isUpperCase(c)) {
                                        upperCase = true;
                                        continue;
                                    }
                                }
                                if (upperCase && lowerCase) {
                                    cut = true;
                                    break;
                                }
                            }
                            first = false;
                        }
                        if (cut) {
                            for (char c : strLine.toCharArray()) {
                                if (Character.isUpperCase(c)) {
                                    sb.append(c);
                                }
                            }
                        } else {
                            sb.append(strLine);
                        }
                    }
                }
                readList.add(sb.toString());
            }
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return readList.toArray(new String[readList.size()]);
    }

    public static List<Triplet<String, String, Integer>> parseFastq(String location) {
        List<Triplet<String, String, Integer>> readList = new LinkedList<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (!strLine.startsWith("@")) {
                        continue;
                    } else {
                        String tag = strLine.split(" ")[0];
                        int pairedNumber = Integer.parseInt(strLine.split(" ")[1].split(":")[0]);
                        String seq = br.readLine();
                        br.readLine();
                        char[] quality = br.readLine().toCharArray();
                        int[] p = new int[quality.length];
                        int begin = -1;
                        boolean started = false;
                        int end = -1;
                        for (int i = 0; i < quality.length; i++) {
                            p[i] = ((int) quality[i]) - 33;
                            if (!started && p[i] >= Globals.THRESHOLD) {
                                started = true;
                                begin = i;
                            }
                        }
                        for (int i = quality.length - 1; i >= 0; i--) {
                            if (p[i] >= Globals.THRESHOLD) {
                                end = i;
                                break;
                            }
                        }
                        if (begin == -1) {
                            continue;
                        }
                        seq = seq.substring(begin, end + 1);
                        readList.add(Triplet.with(seq, tag, pairedNumber));
                    }
                }
            } catch (Exception e) {
                System.err.println("Error Fastq: " + e.getMessage());
                System.err.println(e);
                System.exit(0);
            }
        } catch (Exception e) {
            System.err.println("Error Fastq: " + e.getMessage());
            System.err.println(e);
            System.exit(0);
        }

//        FastqReader fastqReader = new IlluminaFastqReader();
//        
//        try {
//            Iterable<Fastq> reads = fastqReader.read(new File(location));
//            for (Fastq f : reads) {
//                //@generator-0_0.25_899_1068_0:0:0_0:0:0_0/2
//                String description = f.getDescription();
////                String tag = description.split(":")[0];
////                int pairedNumber = Integer.parseInt(description.split("/")[1]);
//                //@HWI-M00276:5:000000000-A1RAA:1:1101:15245:1321 1:N:0:GGACTCCTTATCCTCT
//                String tag = description.split(" ")[0];
//                int pairedNumber = Integer.parseInt(tag.split(":")[0]);
//                readList.add(Triplet.with(f.getSequence(), tag, pairedNumber));
//            }
//        } catch (IOException ex) {
//            Logger.getLogger(FastaParser.class.getName()).log(Level.SEVERE, null, ex);
//        }
        return readList;
    }
    
    public static Map<String, String> parseHaplotypeFile(String location) {
        Map<String, String> hapMap = new ConcurrentHashMap<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            String head = null;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            hapMap.put(sb.toString(), head);
                            sb.setLength(0);
                        }
                        head = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
                hapMap.put(sb.toString(), head);

            }
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return hapMap;
    }
}
