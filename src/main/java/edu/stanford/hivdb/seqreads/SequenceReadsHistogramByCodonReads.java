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
package edu.stanford.hivdb.seqreads;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import org.apache.commons.lang3.tuple.Pair;

import edu.stanford.hivdb.mutations.CodonReads;
import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogram.AggregationOption;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceReadsHistogramByCodonReads<VirusT extends Virus<VirusT>> {

	public static interface WithSequenceReadsHistogramByCodonReads<VirusT extends Virus<VirusT>> {
		
		public SequenceReadsHistogramByCodonReads<VirusT> getHistogramByCodonReads(
			final Long[] codonReadsCutoffs,
			final AggregationOption aggregatesBy);
		
	}
	
	public static class HistogramByCodonReadsBin {
		final public Long cutoff;
		final public Integer count;
		// final public String type;
		
		public HistogramByCodonReadsBin(long cutoff, int count) {
			this.cutoff = cutoff;
			this.count = count;
		}
		
		public Long getCutoff() {
			return cutoff;
		}
		
		public Integer getCount() {
			return count;
		}
	}
	
	final private List<GeneSequenceReads<VirusT>> allGeneSequenceReads;
	final private Integer numBins;
	final private Long[] codonReadsCutoffs;
	final private AggregationOption aggregatesBy;
	
	public SequenceReadsHistogramByCodonReads(
		List<GeneSequenceReads<VirusT>> allGeneSequenceReads,
		Long[] codonReadsCutoffs,
		AggregationOption aggregatesBy
	) {
		this.allGeneSequenceReads = allGeneSequenceReads;
		this.numBins = codonReadsCutoffs.length;
		this.codonReadsCutoffs = codonReadsCutoffs;
		this.aggregatesBy = aggregatesBy;
	}

	private List<HistogramByCodonReadsBin> getSites(
		Function<CodonReads<VirusT>, Boolean> filter
	) {
		/* initialize binsCount */
		int[] binsCount = new int[numBins];

		/* iterate all codon reads; remove filtered ones; save result in binsCount */
		for (GeneSequenceReads<VirusT> geneSeqReads : allGeneSequenceReads) {
			for (PositionCodonReads<VirusT> pcr : geneSeqReads.getAllPositionCodonReads()) {
				Set<Pair<Integer, String>> aggregatedSites = new HashSet<>();
				for (CodonReads<VirusT> cr : pcr.getCodonReads()) {
					long reads = cr.getReads();
					if (cr.isReference()) {
						continue;
					}
					if (filter.apply(cr)) {
						String key;
						switch (aggregatesBy) {
							case AminoAcid:
								key = cr.getAminoAcid().toString();
								break;
							case Position:
								key = "" + pcr.getPosition();
								break;
							default:
								key = cr.getCodon();
						}
						for (int idx = 0; idx < numBins; idx ++) {
							if (reads >= codonReadsCutoffs[idx]) {
								aggregatedSites.add(Pair.of(idx, key));
							}
						}
					}
				}
				for (Pair<Integer, String> idxKey : aggregatedSites) {
					int idx = idxKey.getLeft();
					binsCount[idx] ++;
				}
			}
		}
		
		List<HistogramByCodonReadsBin> result = new ArrayList<>();
		for (int idx = 0; idx < numBins; idx ++) {
			result.add(new HistogramByCodonReadsBin(
				codonReadsCutoffs[idx],
				binsCount[idx]));
		}
		return result;
	}
	
	public List<HistogramByCodonReadsBin> getUsualSites() {
		switch (aggregatesBy) {
			case Position:
				return getSites(cr -> !cr.isUnusual());
			case AminoAcid:
				return getSites(cr -> !cr.isUnusual());
			default:  // case Codon:
				return getSites(cr -> !cr.isUnusualByCodon());
		}
	}
	
	public List<HistogramByCodonReadsBin> getUnusualSites() {
		switch (aggregatesBy) {
			case Position:
			case AminoAcid:
				return getSites(cr -> cr.isUnusual());
			default:  // case Codon:
				return getSites(cr -> cr.isUnusualByCodon());
		}
	}
	
	public List<HistogramByCodonReadsBin> getUnusualApobecSites() {
		switch (aggregatesBy) {
			case Position:
			case AminoAcid:
				return getSites(cr -> cr.isUnusual() && cr.isApobecMutation());
			default:  // case Codon:
				return getSites(cr -> cr.isUnusualByCodon() && cr.isApobecMutation());
		}
	}

	public List<HistogramByCodonReadsBin> getUnusualNonApobecSites() {
		switch (aggregatesBy) {
			case Position:
			case AminoAcid:
				return getSites(cr -> cr.isUnusual() && !cr.isApobecMutation());
			default:  // case Codon:
				return getSites(cr -> cr.isUnusualByCodon() && !cr.isApobecMutation());
		}
	}
	
	public List<HistogramByCodonReadsBin> getApobecSites() {
		return getSites(cr -> cr.isApobecMutation());
	}
	
	public List<HistogramByCodonReadsBin> getApobecDrmSites() {
		return getSites(cr -> cr.isApobecDRM());
	}

	public List<HistogramByCodonReadsBin> getStopCodonSites() {
		return getSites(cr -> cr.hasStop());
	}
	
	public List<HistogramByCodonReadsBin> getDrmSites() {
		return getSites(cr -> cr.isDRM());
	}
	
	public Integer getNumPositions() {
		return (allGeneSequenceReads
				.stream()
				.mapToInt(gsr -> gsr.getNumPositions())
				.sum());
	}
		
}
