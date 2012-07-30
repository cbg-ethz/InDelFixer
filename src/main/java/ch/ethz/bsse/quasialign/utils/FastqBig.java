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

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class FastqBig {
    private String read;
    private String quality;
    private String description;
    private String identifier;

    public FastqBig(String read, String quality, String description, String identifier) {
        this.read = read;
        this.quality = quality;
        this.description = description;
        this.identifier = identifier;
    }

    public FastqBig() {
    }

    public String getRead() {
        return read;
    }

    public String getQuality() {
        return quality;
    }

    public String getDescription() {
        return description;
    }

    public String getIdentifier() {
        return identifier;
    }
    
    public boolean isEmpty() {
        return read == null;
    }

    public void setRead(String read) {
        this.read = read;
    }

    public void setQuality(String quality) {
        this.quality = quality;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
}
