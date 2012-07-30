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
package ch.ethz.bsse.quasialign.utils;

import java.io.*;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Cutter {

    public static void cut(String input, String output, int begin, int end) {
        try {
            // Create file 
            FileWriter outstream = new FileWriter(output+File.separator+"ref_"+(begin+1)+"-"+end+".fasta");
            try (BufferedWriter out = new BufferedWriter(outstream)) {
                FileInputStream fstream = new FileInputStream(input);
                StringBuilder sb;
                String head = null;
                int i = 0;
                try (DataInputStream in = new DataInputStream(fstream)) {
                    BufferedReader br = new BufferedReader(new InputStreamReader(in));
                    String strLine;
                    sb = new StringBuilder();

                    while ((strLine = br.readLine()) != null) {
                        if (strLine.startsWith(">")) {
                            if (sb.length() > 0) {
                                out.write(head);
                                out.write("\n");
                                out.write(sb.substring(begin, end));
                                out.write("\n");
                                sb.setLength(0);
                            }
                            head = strLine;
                        } else {
                            sb.append(strLine);
                        }
                    }

                }
                out.write(head);
                out.write("\n");
                end = end > sb.length() ? sb.length() : end;
                out.write(sb.substring(begin - 1, end));
                out.write("\n");
            }
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
    }
}
