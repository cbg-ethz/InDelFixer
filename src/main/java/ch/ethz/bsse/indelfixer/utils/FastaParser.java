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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava3.sequencing.io.fastq.Fastq;
import org.biojava3.sequencing.io.fastq.FastqReader;
import org.biojava3.sequencing.io.fastq.IlluminaFastqReader;
import org.javatuples.Pair;
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class FastaParser {

    public static String[] parseFarFile(String location) {
        List<String> readList = new LinkedList<>();
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
                        sb.append(strLine);
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
        FastqReader fastqReader = new IlluminaFastqReader();
        List<Triplet<String, String, Integer>> readList = new LinkedList<>();
        try {
            Iterable<Fastq> reads = fastqReader.read(new File(location));
            for (Fastq f : reads) {
                //@generator-0_0.25_899_1068_0:0:0_0:0:0_0/2
                String description = f.getDescription();
                String tag = description.split(":")[0];
                int pairedNumber = Integer.parseInt(description.split("/")[1]);
                readList.add(Triplet.with(f.getSequence(), tag, pairedNumber));
            }
        } catch (IOException ex) {
            Logger.getLogger(FastaParser.class.getName()).log(Level.SEVERE, null, ex);
        }
        return readList;
    }
}
