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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

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
