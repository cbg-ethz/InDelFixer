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
package ch.ethz.bsse.quasialign.sffParser;

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

    public SFFRead(String read, int[] quality, String description,  int clipQualLeft, int clipQualRight) {
        this.read = read;
        this.quality = quality;
        this.description = description;
        this.clipQualLeft = clipQualLeft;
        this.clipQualRight = clipQualRight;
    }
    
    public void clip() {
        if (clipQualLeft!=0 && clipQualRight != 0) {
            this.read = this.read.substring(this.clipQualLeft-1, this.clipQualRight);
            int[] qualityTmp = new int[this.read.length()];
            this.errorProbability = new double[this.read.length()];
            for (int i = 0; i < read.length(); i++) {
                try {
                qualityTmp[i] = this.quality[clipQualLeft-1+i];
                this.errorProbability[i] = Math.pow(10d,-(this.quality[clipQualLeft+i-1]/10d));
                } catch (Exception e) {
                    System.out.println(clipQualLeft-1+i +"\t"+this.quality.length);
                }
            }
            this.quality = qualityTmp;
            this.clipped = true;
        } else {
            System.out.println(this.clipQualLeft+"\t"+this.clipQualRight);
        }
    }

    public boolean isClipped() {
        return clipped;
    }

    public int getClipQualLeft() {
        return clipQualLeft;
    }

    public void setClipQualLeft(int clipQualLeft) {
        this.clipQualLeft = clipQualLeft;
    }

    public int getClipQualRight() {
        return clipQualRight;
    }

    public void setClipQualRight(int clipQualRight) {
        this.clipQualRight = clipQualRight;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public int[] getQuality() {
        return quality;
    }

    public void setQuality(int[] quality) {
        this.quality = quality;
    }

    public String getRead() {
        return read;
    }

    public void setRead(String read) {
        this.read = read;
    }

    public boolean isEmpty() {
        return read == null;
    }
}
