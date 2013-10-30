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
package ch.ethz.bsse.indelfixer.stored;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class GapCosts {

    public float open;
    public float extend;

    public GapCosts(float open, float extend) {
        this.open = open;
        this.extend = extend;
    }

    @Override
    public String toString() {
        return "XO:i:"+(int) open + "\tXE:i:" + (int) extend;
    }
}
