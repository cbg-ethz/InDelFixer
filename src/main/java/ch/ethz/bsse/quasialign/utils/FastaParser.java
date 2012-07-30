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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;

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
}
