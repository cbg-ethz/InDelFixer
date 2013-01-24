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
package ch.ethz.bsse.indelfixer.sffParser;

import ch.ethz.bsse.indelfixer.stored.Globals;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SFFRead {

    private String read;
    private int[] quality;
    private double[] errorProbability;
    private String description;
    private int clipQualLeft;
    private int clipQualRight;
    private boolean clipped;
    
    /**
     *
     * @param read
     * @param quality
     * @param description
     * @param clipQualLeft
     * @param clipQualRight
     */
    public SFFRead(String read, int[] quality, String description, int clipQualLeft, int clipQualRight) {
        this.read = read;
        this.quality = quality;
        this.description = description;
        this.clipQualLeft = clipQualLeft;
        this.clipQualRight = clipQualRight;
    }

    /**
     *
     */
    public void clip() {
        if (clipQualLeft != 0 && clipQualRight != 0) {
            if (clipQualLeft < Globals.CUT) {
                clipQualLeft = Globals.CUT;
            }
            this.read = this.read.substring(this.clipQualLeft - 1, this.clipQualRight);
            int[] qualityTmp = new int[this.read.length()];
            this.errorProbability = new double[this.read.length()];
            for (int i = 0; i < read.length(); i++) {
                try {
                    qualityTmp[i] = this.quality[clipQualLeft - 1 + i];
                    this.errorProbability[i] = Math.pow(10d, -(this.quality[clipQualLeft + i - 1] / 10d));
                } catch (Exception e) {
                }
            }
            this.quality = qualityTmp;
            this.clipped = true;
        } else {
        }
    }

    /**
     *
     * @return
     */
    public boolean isClipped() {
        return clipped;
    }

    /**
     *
     * @return
     */
    public int getClipQualLeft() {
        return clipQualLeft;
    }

    /**
     *
     * @param clipQualLeft
     */
    public void setClipQualLeft(int clipQualLeft) {
        this.clipQualLeft = clipQualLeft;
    }

    /**
     *
     * @return
     */
    public int getClipQualRight() {
        return clipQualRight;
    }

    /**
     *
     * @param clipQualRight
     */
    public void setClipQualRight(int clipQualRight) {
        this.clipQualRight = clipQualRight;
    }

    /**
     *
     * @return
     */
    public String getDescription() {
        return description;
    }

    /**
     *
     * @param description
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     *
     * @return
     */
    public String getQuality() {
        StringBuilder sb = new StringBuilder(this.quality.length);
        for (int q : quality) {
            sb.append(String.valueOf(q));
        }
        return sb.toString();
    }

    /**
     *
     * @param quality
     */
    public void setQuality(int[] quality) {
        this.quality = quality;
    }

    /**
     *
     * @return
     */
    public String getRead() {
        return read;
    }

    /**
     *
     * @param read
     */
    public void setRead(String read) {
        this.read = read;
    }

    /**
     *
     * @return
     */
    public boolean isEmpty() {
        return read == null;
    }
}
