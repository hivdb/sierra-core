/*

    Copyright (C) 2017 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.mutations;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import edu.stanford.hivdb.viruses.Virus;

public class MutationFileReader {

	/**
	 * Reads lists of mutations, one line at a time, from a comma-delimited file
	 *   and returns a list of mutation objects.
	 * Each mutation consists of a gene, position, and one or more amino acids.
	 *   There is no consensus amino acid preceding the position
	 *
	 * @param fileInputStream
	 * @return List<Mutation>
	 */
	public static <VirusT extends Virus<VirusT>> List<MutationSet<VirusT>> readMutationLists(InputStream fileInputStream, VirusT virusIns) {
		List<MutationSet<VirusT>> fileMuts = new ArrayList<>();

		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(fileInputStream));
			String line;
			while ((line = br.readLine()) != null) {
				if (shouldSkip(line)) continue;
				MutationSet<VirusT> lineMuts = virusIns.newMutationSet(line.trim());
				if (!lineMuts.isEmpty()) {
					fileMuts.add(lineMuts);
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return fileMuts;
	}

	protected static boolean shouldSkip(String line) {
		return line.isEmpty() || line.startsWith("#");
	}
}
