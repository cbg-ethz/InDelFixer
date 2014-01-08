/*
 * $Id: MatrixLoader.java,v 1.3 2006/01/18 20:16:37 ahmed Exp $
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package jaligner.matrix;

import jaligner.ui.filechooser.NamedInputStream;
import jaligner.util.Commons;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Scoring matrices loader from a jar file or a file system.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class MatrixLoader {
	/**
	 * The starter character of a comment line.
	 */
	private static final char COMMENT_STARTER = '#';

	/**
	 * The path to the matrices within the package.
	 */
	private static final String MATRICES_HOME = "jaligner/matrix/matrices/";

	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(MatrixLoader.class
			.getName());

	/**
	 * Loads scoring matrix from Jar file or file system.
	 * 
	 * @param matrix
	 *            to load
	 * @return loaded matrix
	 * @throws MatrixLoaderException
	 * @see Matrix
	 */
	public static Matrix load(String matrix) throws MatrixLoaderException {
		InputStream is = null;

		if (new StringTokenizer(matrix, Commons.getFileSeparator())
				.countTokens() == 1) {
			// Matrix does not include the path
			// Load the matrix from matrices.jar
			is = MatrixLoader.class.getClassLoader().getResourceAsStream(
					MATRICES_HOME + matrix);
		} else {
			// Matrix includes the path information
			// Load the matrix from the file system
			try {
				is = new FileInputStream(matrix);
			} catch (Exception e) {
				String message = "Failed opening input stream: "
						+ e.getMessage();
				logger.log(Level.SEVERE, message, e);
				throw new MatrixLoaderException(message);
			}
		}

		return load(new NamedInputStream(matrix, is));
	}

	/**
	 * Loads scoring matrix from {@link InputStream}
	 * 
	 * @param nis
	 *            named input stream
	 * @return loaded matrix
	 * @throws MatrixLoaderException
	 * @see Matrix
	 * @see NamedInputStream
	 */
	public static Matrix load(NamedInputStream nis)
			throws MatrixLoaderException {
		logger.info("Loading scoring matrix...");
		char[] acids = new char[Matrix.SIZE];

		// Initialize the acids array to null values (ascii = 0)
		for (int i = 0; i < Matrix.SIZE; i++) {
			acids[i] = 0;
		}

		float[][] scores = new float[Matrix.SIZE][Matrix.SIZE];

		BufferedReader reader = new BufferedReader(new InputStreamReader(nis
				.getInputStream()));

		String line;

		try {
			// Skip the comment lines
			while ((line = reader.readLine()) != null
					&& line.trim().charAt(0) == COMMENT_STARTER)
				;
		} catch (Exception e) {
			String message = "Failed reading from input stream: "
					+ e.getMessage();
			logger.log(Level.SEVERE, message, e);
			throw new MatrixLoaderException(message);
		}

		// Read the headers line (the letters of the acids)
		StringTokenizer tokenizer;
		tokenizer = new StringTokenizer(line.trim());
		for (int j = 0; tokenizer.hasMoreTokens(); j++) {
			acids[j] = tokenizer.nextToken().charAt(0);
		}

		try {
			// Read the scores
			while ((line = reader.readLine()) != null) {
				tokenizer = new StringTokenizer(line.trim());
				char acid = tokenizer.nextToken().charAt(0);
				for (int i = 0; i < Matrix.SIZE; i++) {
					if (acids[i] != 0) {
						scores[acid][acids[i]] = Float.parseFloat(tokenizer
								.nextToken());
					}
				}
			}
		} catch (Exception e) {
			String message = "Failed reading from input stream: "
					+ e.getMessage();
			logger.log(Level.SEVERE, message, e);
			throw new MatrixLoaderException(message);
		}
		logger.info("Finished loading scoring matrix");
		return new Matrix(nis.getName(), scores);
	}
}