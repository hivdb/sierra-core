/*

    Copyright (C) 2019-2020 Stanford HIVDB team

    This file is part of Sierra.

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sierra.  If not, see <https://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.mutations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

/**
 * An implementation of the Mutation interface to accept
 * multiple codon reads. In this class, a codon is strictly
 * restricted to non-ambiguous nucleotide code (ACGT).
 */
public class MultiCodonsMutation<VirusT extends Virus<VirusT>> extends AAMutation<VirusT> {

	private static final int DEFAULT_MAX_DISPLAY_AAS = 6;

	private final long totalCount;
	private final String compatTriplet;

	public static <VirusT extends Virus<VirusT>> MultiCodonsMutation<VirusT> initUnsequenced(Gene<VirusT> gene, int position) {
		return new MultiCodonsMutation<>(gene, position, Collections.emptySet(), 0, "NNN");
	}

	private static <VirusT extends Virus<VirusT>> Set<Character>
	getAACharSet(PositionCodonReads<VirusT> posCodonReads, long minReads) {
		Set<Character> aaCharSet = new TreeSet<>();

		for (CodonReads<VirusT> codonReads : posCodonReads.getCodonReads()) {
			char aa = codonReads.getAminoAcid();
			if (aa == 'X') {
				continue;
			}
			long count = codonReads.getReads();
			if (count < minReads) {
				// remove minor variants below min-prevalence
				continue;
			}
			aaCharSet.add(aa);
		}
		return aaCharSet;
	}

	private static <VirusT extends Virus<VirusT>> String getCompatTriplet(
		PositionCodonReads<VirusT> posCodonReads, long minReads
	) {
		List<String> cleanCodons = new ArrayList<>();
		for (CodonReads<VirusT> codonReads : posCodonReads.getCodonReads()) {
			// Tolerant spaces and dashes
			String codon = codonReads.getCodon().replaceAll("[ -]", "");
			long count = codonReads.getReads();
			if (count <= minReads) {
				// remove minor variants below min-prevalence
				continue;
			}
			if (!codon.matches("^[ACGT]*$")) {
				// do not allow ambiguous codes
				continue;
			}
			int codonLen = codon.length();
			if (codonLen < 3 || codonLen > 5) {
				// skip indels
				continue;
			}
			cleanCodons.add(codon.substring(0, 3));
		}
		return CodonUtils.getMergedCodon(cleanCodons);
	}

	public static <VirusT extends Virus<VirusT>> MultiCodonsMutation<VirusT> fromPositionCodonReads(
		PositionCodonReads<VirusT> posCodonReads, double minPrevalence
	) {
		Gene<VirusT> gene = posCodonReads.getGene();
		int position = (int) posCodonReads.getPosition();
		long totalCount = posCodonReads.getTotalReads();
		long minReads = Math.round(totalCount * minPrevalence + 0.5);
		Set<Character> aaCharSet = getAACharSet(posCodonReads, minReads);
		char ref = gene.getRefChar(position);
		if (aaCharSet.isEmpty() ||
			(aaCharSet.size() == 1 && aaCharSet.contains(ref))
		) {
			return null;
		}
		String compatTriplet = getCompatTriplet(posCodonReads, minReads);
		return new MultiCodonsMutation<>(
			gene, position, aaCharSet, totalCount, compatTriplet);
	}

	private MultiCodonsMutation(
		Gene<VirusT> gene, int position,
		Set<Character> aaCharSet,
		long totalCount, String compatTriplet
	) {
		super(gene, position, aaCharSet, DEFAULT_MAX_DISPLAY_AAS);
		this.totalCount = totalCount;
		this.compatTriplet = compatTriplet;
	}

	/**
	 * Gets total read count of all codons (include codons of reference AA)
	 *
	 * @return a Long number
	 */
	public Long getTotalCount() { return totalCount; }

	@Override
	public boolean isUnsequenced() { return this.totalCount == 0; }

	@Override
	public String getTriplet() { return compatTriplet; }

	/**
	 * There's no way to tell about inserted NAs for multiple codons without
	 * an alignment tools. Therefore we simply returns an empty result.
	 */
	@Override
	public String getInsertedNAs() { return ""; }

	@Override
	public boolean hasBDHVN() {
		return getTriplet().matches(".*[BDHVN].*");
	}

}
