/**
 * Copyright (c) 2011-2012 Armin Töpfer
 *
 * This file is part of InDelFixer.
 *
 * InDelFixer is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * InDelFixer is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * InDelFixer. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.indelfixer.sffParser;

import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.sffParser.SFFParser;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SFFParsing {

    public static void main(String args[]) {
        parse("/Users/XLR/Documents/lib100-11_E2+NS5B_TV1-98NRBrescia.sff");
    }

    public static String[] parse(String path) {
        Globals.N = 0;
        List<SFFRead> sffReads = new ArrayList<>();
        String[] split = path.split(";");
        for (String s : split) {
            SFFParser sff = null;
            try {
                sff = new SFFParser(path);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(SFFParsing.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(SFFParsing.class.getName()).log(Level.SEVERE, null, ex);
            }
            Globals.N += sff.getNumberOfReads();
            sffReads.addAll(sff.getReads());
        }

        String[] reads = new String[sffReads.size()];
        int i = 0;
        for (SFFRead r : sffReads) {
            reads[i++] = r.getRead();
        }
        sffReads = null;
        System.gc();
        System.gc();
        return reads;
    }
}
