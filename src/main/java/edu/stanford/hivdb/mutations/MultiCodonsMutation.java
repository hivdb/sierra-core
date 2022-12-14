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

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

/**
 * An implementation of the Mutation interface to accept
 * multiple codon reads. In this class, a codon is strictly
 * restricted to non-ambiguous nucleotide code (ACGT).
 */
public class MultiCodonsMutation<VirusT extends Virus<VirusT>> extends AAMutation<VirusT> {

	public static class AAReads {
		private final Character aminoAcid;
		private final Long numReads;
		private final Double percent;
		
		public AAReads(char aminoAcid, long numReads, double percent) {
			this.aminoAcid = aminoAcid;
			this.numReads = numReads;
			this.percent = percent;
		}
		
		public Character getAminoAcid() {
			return aminoAcid;
		}
		
		public Long getNumReads() {
			return numReads;
		}
		
		public Double getPercent() {
			return percent;
		}
	}
	
	private static final int DEFAULT_MAX_DISPLAY_AAS = 6;

	private final long totalReads;
	private final List<AAReads> allAAReads;
	private final String compatTriplet;

	public static <VirusT extends Virus<VirusT>> MultiCodonsMutation<VirusT> initUnsequenced(Gene<VirusT> gene, int position) {
		return new MultiCodonsMutation<>(gene, position, Collections.emptyMap(), 0, "NNN");
	}

	private static <VirusT extends Virus<VirusT>> Map<Character, Long>
	getAACharReadsMap(List<CodonReads<VirusT>> codonReads) {
		Map<Character, Long> aaCharReadsMap = new TreeMap<>();

		for (CodonReads<VirusT> cdr : codonReads) {
			if (cdr.isAmbiguous() || cdr.isDelFrameshift()) {
				continue;
			}
			char aa = cdr.getAminoAcid();
			long count = cdr.getReads();
			aaCharReadsMap.put(aa, aaCharReadsMap.getOrDefault(aa, 0L) + count);
		}
		return aaCharReadsMap;
	}

	public static <VirusT extends Virus<VirusT>> MultiCodonsMutation<VirusT> fromPositionCodonReads(
		PositionCodonReads<VirusT> posCodonReads, double minPrevalence, long minCodonReads
	) {
		Gene<VirusT> gene = posCodonReads.getGene();
		int position = (int) posCodonReads.getPosition();
		List<CodonReads<VirusT>> codonReads = posCodonReads.getCodonReadsUsingThreshold(minPrevalence, minCodonReads);
		if (codonReads.isEmpty()) {
			return initUnsequenced(gene, position);
		}
		long totalCount = posCodonReads.getTotalReads();
		String compatTriplet = PositionCodonReads.getCodonConsensusWithoutIns(codonReads);
		Map<Character, Long> aaCharReadsMap = getAACharReadsMap(codonReads);
		char ref = gene.getRefChar(position);
		if (aaCharReadsMap.isEmpty() ||
			(aaCharReadsMap.size() == 1 && aaCharReadsMap.containsKey(ref))
		) {
			// don't create a mutation object when only reference AA is present
			return null;
		}
		return new MultiCodonsMutation<>(
			gene, position, aaCharReadsMap, totalCount, compatTriplet);
	}

	private MultiCodonsMutation(
		Gene<VirusT> gene, int position,
		Map<Character, Long> aaCharReadsMap,
		long totalReads, String compatTriplet
	) {
		super(gene, position, aaCharReadsMap.keySet(), DEFAULT_MAX_DISPLAY_AAS);
		this.allAAReads = (
			aaCharReadsMap.entrySet()
			.stream()
			.map(entry -> new AAReads(
				entry.getKey(),
				entry.getValue(),
				entry.getValue().doubleValue() / totalReads * 100
			))
			.collect(Collectors.toList())
		);
		this.totalReads = totalReads;
		this.compatTriplet = compatTriplet;
	}

	/**
	 * Gets total read count of all codons (include codons of reference AA)
	 *
	 * @return a Long number
	 */
	public Long getTotalReads() { return totalReads; }
	
	/**
	 * Gets list of read count / prevalence for each AA (include reference AA)
	 * 
	 * @return a List of AAReads objects
	 */
	public List<AAReads> getAllAAReads() { return allAAReads; }

	@Override
	public boolean isUnsequenced() { return this.totalReads == 0; }

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
