/**
 * Copyright (c) 2011-2013 Armin Töpfer
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
package ch.ethz.bsse.indelfixer.minimal.processing;

import ch.ethz.bsse.indelfixer.stored.Globals;
import ch.ethz.bsse.indelfixer.stored.GridOutput;
import ch.ethz.bsse.indelfixer.stored.Read;
import ch.ethz.bsse.indelfixer.stored.SequenceEntry;
import ch.ethz.bsse.indelfixer.utils.Utils;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ProcessResults implements Callable<Void> {

    private ProcessingGeneral pg;
    private List<Future<List<Object>>> newList;
    private Map<Integer, Read> undone = new HashMap<>();

    public ProcessResults(ProcessingGeneral pg, List<Future<List<Object>>> newList) {
        this.pg = pg;
        this.newList = newList;
    }

    @Override
    public Void call() throws Exception {
        StringBuilder samSB = new StringBuilder();
        List<SequenceEntry> trash = new LinkedList<>();
        for (Future<List<Object>> future : newList) {
            Object o = future.get();
            if (o == null) {
                continue;
            } else if (o instanceof List) {
                List<Object> list = (List) o;
                for (Object o2 : list) {
                    if (o2 instanceof SequenceEntry) {
                        SequenceEntry result = (SequenceEntry) o2;
                        trash.add(result);
                    } else if (o2 instanceof Read) {
                        Read r = (Read) o2;
                        if (r.getPairedNumber() != -1) {
                            if (undone.containsKey(r.getPairedNumber())) {
                                samSB.append(r.toString(undone.get(r.getPairedNumber())));
                                samSB.append(undone.get(r.getPairedNumber()).toString(r));
                                undone.remove(r.getPairedNumber());
                            } else {
                                undone.put(r.getPairedNumber(), r);
                            }
                        } else {
                            samSB.append(r.toString(null));
                        }
                        if (Globals.REFINE || Globals.CONSENSUS) {
                            int[] x = Utils.reverse(r.getAlignedRead());
                            int i = 0;
                            int y = 0;
                            for (char c : r.getCigarsPure()) {
                                switch (c) {
                                    case 'X':
                                    case 'M':
                                        try {
                                            int a = x[y++];
                                            pg.alignment[r.getBegin() - 1 + i++][a]++;
                                        } catch (ArrayIndexOutOfBoundsException e) {
                                            System.err.println(r.getCigarsPure().length);
                                            System.err.println(y);
                                            System.err.println("w00t");
                                            System.exit(9);
                                        }
                                        break;
                                    case 'D':
                                        pg.alignment[r.getBegin() - 1 + i++][4]++;
                                        break;
                                    case 'I':
                                        y++;
                                        break;
                                    default:
                                        System.err.println("why");
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (Read r : undone.values()) {
            samSB.append(r.toString(null));
        }
        Utils.appendFile(Globals.output + "reads.sam", samSB.toString());

        StringBuilder sb = new StringBuilder();
        for (SequenceEntry s : trash) {
            sb.append(s.header).append("\n");
            sb.append(s.sequence).append("\n");
            sb.append("+").append("\n");
            sb.append(s.quality).append("\n");
        }
        Utils.appendFile(Globals.output + "trash.fastq", sb.toString());
        return null;
    }
}
