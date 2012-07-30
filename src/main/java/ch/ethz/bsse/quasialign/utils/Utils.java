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

import ch.ethz.bsse.quasialign.stored.Globals;
import ch.ethz.bsse.quasialign.stored.InformationHolder;
import java.io.*;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Utils {

    public static void saveFile(String path, String sb) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(path);
            try (BufferedWriter out = new BufferedWriter(fstream)) {
                out.write(sb);
            }
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error save file: ");
            System.err.println(path);
            System.err.println(sb);
        }
    }

    public static void saveInformationHolder(InformationHolder ih) {
        try {
            String s = Globals.output + "support" + File.separator;
            if (!new File(s).exists()) {
                new File(s).mkdirs();
            }
            s += "snapshot";
            FileOutputStream fos = new FileOutputStream(s);
            try (ObjectOutputStream out = new ObjectOutputStream(fos)) {
                out.writeObject(ih);
            }
        } catch (IOException ex) {
            System.out.println("Snapshot saving problem\n" + ex.getMessage());
        }
    }

    public static String reverseComplement(String s) {
        StringBuilder sb = new StringBuilder();
        for (char c : s.toUpperCase().toCharArray()) {
            switch (c) {
                case 'A':
                    sb.append("T");
                    break;
                case 'C':
                    sb.append("G");
                    break;
                case 'G':
                    sb.append("C");
                    break;
                case 'T':
                    sb.append("A");
                    break;
                case '-':
                    sb.append("-");
                    break;
            }
        }
        return sb.reverse().toString();
    }

    public static boolean isFastaFormat(String path) {
        try {
            FileInputStream fstream = new FileInputStream(path);
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        } catch (Exception e) {
            System.err.println("Error identifying format of input file: " + e.getMessage());
        }
        return false;
    }
}
