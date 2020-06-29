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

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.tuple.Pair;

import edu.stanford.hivdb.mutations.CodonReads;
import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceReadsHistogram<VirusT extends Virus<VirusT>> {
	
	private static double MIN_LOG10_LOWER_LIMIT = -100d;

	public static enum AggregationOption {
		Codon, AminoAcid, Position;
	}
	
	public static interface WithSequenceReadsHistogram<VirusT extends Virus<VirusT>> {
		
		public SequenceReadsHistogram<VirusT> getHistogram(
			final Double pcntLowerLimit,
			final Double pcntUpperLimit,
			final Integer numBins,
			final Boolean cumulative,
			final AggregationOption aggregatesBy);

		public SequenceReadsHistogram<VirusT> getHistogram(
			final Double pcntLowerLimit,
			final Double pcntUpperLimit,
			final Double[] binTicks,
			final Boolean cumulative,
			final AggregationOption aggregatesBy);
		
	}
	
	public static class HistogramBin {
		final public Double log10Percent;
		final public Double binWidth;
		final public Integer count;
		// final public String type;
		
		public HistogramBin(double log10Percent, double binWidth, int count/*, String type*/) {
			this.log10Percent = log10Percent;
			this.binWidth = binWidth;
			this.count = count;
			// this.type = type;
		}
		
		public Double getPercentStart() {
			return Math.pow(10, log10Percent) * 100;
		}
		
		public Double getPercentStop() {
			return Math.pow(10, log10Percent + binWidth) * 100;
		}
	}
	
	final private List<GeneSequenceReads<VirusT>> allGeneSequenceReads;
	final private Double log10PcntLowerLimit;
	final private Double log10PcntUpperLimit;
	final private Integer numBins;
	final private Double[] binSteps;
	final private Boolean cumulative;
	final private AggregationOption aggregatesBy;
	
	public SequenceReadsHistogram(
		List<GeneSequenceReads<VirusT>> allGeneSequenceReads,
		double pcntLowerLimit, double pcntUpperLimit, int numBins,
		boolean cumulative, AggregationOption aggregatesBy
	) {
		this.allGeneSequenceReads = allGeneSequenceReads;
		this.numBins = numBins;
		this.cumulative = cumulative;
		this.aggregatesBy = aggregatesBy;
		log10PcntLowerLimit = Math.max(MIN_LOG10_LOWER_LIMIT, Math.log10(pcntLowerLimit));
		log10PcntUpperLimit = Math.log10(pcntUpperLimit);
		double binWidth = (log10PcntUpperLimit - log10PcntLowerLimit) / numBins;

		Double[] binSteps = new Double[numBins + 1];
		for (int idx = 0; idx < numBins + 1; idx ++) {
			binSteps[idx] = log10PcntLowerLimit + idx * binWidth;
		}
		this.binSteps = binSteps;
	}

	public SequenceReadsHistogram(
		List<GeneSequenceReads<VirusT>> allGeneSequenceReads,
		double pcntLowerLimit, double pcntUpperLimit, Double[] binTicks,
		boolean cumulative, AggregationOption aggregatesBy
	) {
		this.allGeneSequenceReads = allGeneSequenceReads;
		this.numBins = binTicks.length;
	
		this.cumulative = cumulative;
		this.aggregatesBy = aggregatesBy;
		log10PcntLowerLimit = Math.max(MIN_LOG10_LOWER_LIMIT, Math.log10(pcntLowerLimit));
		log10PcntUpperLimit = Math.log10(pcntUpperLimit);

		Double[] binSteps = new Double[numBins + 1];
		for (int idx = 0; idx < numBins + 1; idx ++) {
			if (idx == numBins) {
				binSteps[idx] = log10PcntUpperLimit;
			}
			else {
				binSteps[idx] = Math.log10(binTicks[idx]);
			}
		}
		this.binSteps = binSteps;
	}
	
	private boolean testBetween(double val, double lowerLimit, double upperLimit) {
		return val >= lowerLimit && (cumulative || val <= upperLimit);
	}

	private List<HistogramBin> getSites(
		Function<CodonReads<VirusT>, Boolean> filter
	) {
		/* initialize binsCount */
		int[] binsCount = new int[numBins];

		/* iterate all codon reads; remove filtered ones; save result in binsCount */
		for (GeneSequenceReads<VirusT> geneSeqReads : allGeneSequenceReads) {
			for (PositionCodonReads<VirusT> pcr : geneSeqReads.getAllPositionCodonReads()) {
				Set<Pair<Integer, String>> aggregatedSites = new HashSet<>();
				long total = pcr.getTotalReads();
				double log10Total = Math.log10(total);
				for (CodonReads<VirusT> cr : pcr.getCodonReads()) {
					long reads = cr.getReads();
					double log10Reads = Math.log10(reads);
					double log10Pcnt = log10Reads - log10Total;
					if (cr.isReference()) {
						continue;
					}
					if (!testBetween(log10Pcnt, log10PcntLowerLimit, log10PcntUpperLimit)) {
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
							if (testBetween(log10Pcnt, binSteps[idx], binSteps[idx + 1])) {
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
		
		List<HistogramBin> result = new ArrayList<>();
		for (int idx = 0; idx < numBins; idx ++) {
			result.add(new HistogramBin(
				binSteps[idx],
				binSteps[idx + 1] - binSteps[idx],
				binsCount[idx]));
		}
		return result;
	}
	
	public List<HistogramBin> getUsualSites() {
		switch (aggregatesBy) {
			case Position:
				return getSites(cr -> !cr.isUnusual());
			case AminoAcid:
				return getSites(cr -> !cr.isUnusual());
			default:  // case Codon:
				return getSites(cr -> !cr.isUnusualByCodon());
		}
	}
	
	public List<HistogramBin> getUsualSites(String treatment, String subtype) {
		switch (aggregatesBy) {
			case Position:
				return getSites(cr -> !cr.isUnusual(treatment, subtype));
			case AminoAcid:
				return getSites(cr -> !cr.isUnusual(treatment, subtype));
			default:  // case Codon:
				throw new NotImplementedException("Unsupport aggregate by codon");
				// return getSites(cr -> !cr.isUnusualByCodon());
		}
		
	}
	
	public List<HistogramBin> getUnusualSites() {
		switch (aggregatesBy) {
			case Position:
			case AminoAcid:
				return getSites(cr -> cr.isUnusual());
			default:  // case Codon:
				return getSites(cr -> cr.isUnusualByCodon());
		}
	}
	
	public List<HistogramBin> getUnusualApobecSites() {
		switch (aggregatesBy) {
			case Position:
			case AminoAcid:
				return getSites(cr -> cr.isUnusual() && cr.isApobecMutation());
			default:  // case Codon:
				return getSites(cr -> cr.isUnusualByCodon() && cr.isApobecMutation());
		}
	}

	public List<HistogramBin> getUnusualNonApobecSites() {
		switch (aggregatesBy) {
			case Position:
			case AminoAcid:
				return getSites(cr -> cr.isUnusual() && !cr.isApobecMutation());
			default:  // case Codon:
				return getSites(cr -> cr.isUnusualByCodon() && !cr.isApobecMutation());
		}
	}
	
	public List<HistogramBin> getApobecSites() {
		return getSites(cr -> cr.isApobecMutation());
	}
	
	public List<HistogramBin> getApobecDrmSites() {
		return getSites(cr -> cr.isApobecDRM());
	}

	public List<HistogramBin> getStopCodonSites() {
		return getSites(cr -> cr.hasStop());
	}
	
	public List<HistogramBin> getDrmSites() {
		return getSites(cr -> cr.isDRM());
	}
	
	public Integer getNumPositions() {
		return (allGeneSequenceReads
				.stream()
				.mapToInt(gsr -> gsr.getNumPositions())
				.sum());
	}
		
}
