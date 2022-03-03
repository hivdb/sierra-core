/*

	Copyright (C) 2017-2020 Stanford HIVDB team

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
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import com.google.common.collect.Streams;

import edu.stanford.hivdb.genotypes.BoundGenotype;
import edu.stanford.hivdb.genotypes.GenotypeResult;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.UntranslatedRegion;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.seqreads.CutoffCalculator.CutoffKeyPoint;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogram.AggregationOption;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogram.WithSequenceReadsHistogram;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogramByCodonReads.WithSequenceReadsHistogramByCodonReads;
import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceReads<VirusT extends Virus<VirusT>>
implements WithSequenceReadsHistogram<VirusT>, WithSequenceReadsHistogramByCodonReads<VirusT> {

	private final Strain<VirusT> strain;
	private final Map<Gene<VirusT>, GeneSequenceReads<VirusT>> allGeneSequenceReads;
	private final List<UntranslatedRegion> untranslatedRegions;
	private final List<OneCodonReadsCoverage<VirusT>> codonReadsCoverage;
	private final CutoffCalculator<VirusT> cutoffObj;
	private final Double proportionTrimmedPositions;
	private final String name;
	private GenotypeResult<VirusT> subtypeResult;
	private MutationSet<VirusT> mutations;
	private String concatenatedSeq;
	private transient DescriptiveStatistics readDepthStats;
	private transient DescriptiveStatistics readDepthStatsDRP;
	private transient List<ValidationResult> validationResults;
	private transient Map<String, UntranslatedRegion> utrLookup;
	
	public static <VirusT extends Virus<VirusT>> SequenceReads<VirusT> fromCodonReadsTable(
			final String name,
			final Strain<VirusT> strain, List<PositionCodonReads<VirusT>> allReads,
			final List<UntranslatedRegion> untranslatedRegions,
			final Double maxMixtureRate,
			final Double minPrevalence,
			final Long minCodonReads,
			final Long minPositionReads
		) {
		// TODO: dynamic cutoff
		CutoffCalculator<VirusT> cutoffObj = new CutoffCalculator<>(
			allReads,
			maxMixtureRate,
			minPrevalence,
			minCodonReads,
			minPositionReads
		);
		
		// sort by genePosition first
		allReads.sort((o1, o2) -> o1.getGenePosition().compareTo(o2.getGenePosition()));
		
		List<PositionCodonReads<VirusT>> filteredAllReads = allReads.stream()
			// remove all codons with their read depth < minPositionReads
			.filter(read -> read.getTotalReads() >= cutoffObj.getMinPositionReads())
			.collect(Collectors.toList());

		Map<Gene<VirusT>, GeneSequenceReads<VirusT>> geneSequences = filteredAllReads
			.stream()
			.collect(
				Collectors.groupingBy(
					PositionCodonReads::getGene,
					TreeMap::new,
					Collectors.collectingAndThen(
						Collectors.toList(),
						list -> new GeneSequenceReads<>(list, cutoffObj)
					)
				)
			);
		List<OneCodonReadsCoverage<VirusT>> codonReadsCoverage = allReads.stream()
			.map(read -> new OneCodonReadsCoverage<>(
				read.getGene(),
				read.getPosition(),
				read.getTotalReads(),
				read.getTotalReads() < cutoffObj.getMinPositionReads()
			))
			.collect(Collectors.toList());
		double proportionTrimmedPositions = 1. - (double) filteredAllReads.size() / (double) allReads.size();

		// TODO: add support for HIV2
		return new SequenceReads<>(
				name,
				strain,
				geneSequences,
				untranslatedRegions,
				cutoffObj,
				codonReadsCoverage,
				proportionTrimmedPositions);
	}

	protected SequenceReads(
			final String name,
			final Strain<VirusT> strain,
			final Map<Gene<VirusT>, GeneSequenceReads<VirusT>> allGeneSequenceReads,
			final List<UntranslatedRegion> untranslatedRegions,
			final CutoffCalculator<VirusT> cutoffSuggestion,
			final List<OneCodonReadsCoverage<VirusT>> codonReadsCoverage,
			final double proportionTrimmedPositions) {
		this.name = name;
		this.strain = strain;
		this.allGeneSequenceReads = Collections.unmodifiableMap(allGeneSequenceReads);
		this.untranslatedRegions = untranslatedRegions;
		this.cutoffObj = cutoffSuggestion;
		this.codonReadsCoverage = Collections.unmodifiableList(codonReadsCoverage);
		this.proportionTrimmedPositions = proportionTrimmedPositions;
	}
	
	public List<OneCodonReadsCoverage<VirusT>> getCodonReadsCoverage() {
		return codonReadsCoverage;
	}
	
	public Map<String, UntranslatedRegion> getUntranslatedRegionLookup() {
		if (utrLookup == null) {
			utrLookup = (
				untranslatedRegions
				.stream()
				.collect(Collectors.toMap(utr -> utr.getName(), utr -> utr))
			);
		}
		return utrLookup;
	}
	
	public String getAssembledUnambiguousConsensus(boolean skipIns) {
		return (
			strain
			.getSequenceReadsAssembler()
			.assemble(
				/* allGeneSeqReads = */allGeneSequenceReads,
				/* utrLookup = */getUntranslatedRegionLookup(),
				skipIns,
				/* includeAmbiguousNA = */false
			)
		);
	}
	
	public String getAssembledConsensus(boolean skipIns) {
		return (
			strain
			.getSequenceReadsAssembler()
			.assemble(
				/* allGeneSeqReads = */allGeneSequenceReads,
				/* utrLookup = */getUntranslatedRegionLookup(),
				skipIns,
				/* includeAmbiguousNA = */true
			)
		);
	}
	
	public String getAssembledUnambiguousConsensus() {
		return getAssembledUnambiguousConsensus(false).replaceAll("(^N+|N+$)",  "");
	}
	
	public String getAssembledConsensus() {
		return getAssembledConsensus(false).replaceAll("(^N+|N+$)", "");
	}
	
	public Integer getSize() {
		int totalSize = 0;
		for (Gene<VirusT> gene : allGeneSequenceReads.keySet()) {
			GeneSequenceReads<VirusT> gsr = allGeneSequenceReads.get(gene);
			totalSize += gsr.getSize();
		}
		return totalSize;
	}

	public Double getProportionTrimmedPositions() { return proportionTrimmedPositions; }
	
	public String getName() { return name; }
	
	public Strain<VirusT> getStrain() { return strain; }

	public boolean isEmpty() { return allGeneSequenceReads.isEmpty(); }

	public Double getMaxMixtureRate() { return cutoffObj.getMaxMixtureRate(); }
	
	public Double getMixtureRate() { return cutoffObj.getActualMixtureRate(); }

	public Double getMinPrevalence() { return cutoffObj.getMinPrevalence(); }

	public Double getActualMinPrevalence() { return cutoffObj.getActualMinPrevalence(); }
	
	public Long getMinCodonReads() { return cutoffObj.getMinCodonReads(); }

	public Long getMinPositionReads() { return cutoffObj.getMinPositionReads(); }
	
	public List<CutoffKeyPoint> getCutoffKeyPoints() { return cutoffObj.getCutoffKeyPoints(); }
	
	public List<ValidationResult> getValidationResults() {
		if (validationResults == null) {
			validationResults = (
				strain.getVirusInstance()
				.validateSequenceReads(this)
			);
		}
		return validationResults;
	}

	public DescriptiveStatistics getReadDepthStats() {
		if (readDepthStats == null) {
			Optional<DoubleStream> readDepthStream = allGeneSequenceReads.values().stream()
				.map(gsr -> gsr.getAllPositionCodonReads())
				.map(pcrs -> pcrs.stream().mapToDouble(pcr -> pcr.getTotalReads()))
				.reduce((a, b) -> Streams.concat(a, b));
			if (readDepthStream.isPresent()) {
				double[] readDepthArray = readDepthStream.get().toArray();
				readDepthStats = new DescriptiveStatistics(readDepthArray);
			}
			else {
				readDepthStats = new DescriptiveStatistics(new double[] {0, 0, 0});
			}
		}
		return readDepthStats;
	}
	
	public DescriptiveStatistics getReadDepthStatsDRP() {
		if (readDepthStatsDRP == null) {
			Optional<DoubleStream> readDepthStream = allGeneSequenceReads.values().stream()
				.map(gsr -> gsr.getAllPositionCodonReads())
				.map(pcrs -> (
					pcrs.stream()
					.filter(pcr -> pcr.getGenePosition().isDrugResistancePosition())
					.mapToDouble(pcr -> pcr.getTotalReads()))
				)
				.reduce((a, b) -> Streams.concat(a, b));
			if (readDepthStream.isPresent()) {
				double[] readDepthArray = readDepthStream.get().toArray();
				readDepthStatsDRP = new DescriptiveStatistics(readDepthArray);
			}
			else {
				readDepthStatsDRP = new DescriptiveStatistics(new double[] {0, 0, 0});
			}
		}
		return readDepthStatsDRP;
	}

	public List<GeneSequenceReads<VirusT>> getAllGeneSequenceReads() {
		return new ArrayList<>(allGeneSequenceReads.values());
	}

	public GeneSequenceReads<VirusT> getGeneSequenceReads(Gene<VirusT> gene) {
		return allGeneSequenceReads.get(gene);
	}

	public List<Gene<VirusT>> getAvailableGenes() {
		return new ArrayList<>(allGeneSequenceReads.keySet());
	}

	public String getConcatenatedSeq() {
		if (concatenatedSeq == null) {
			StringBuilder concatSeq = new StringBuilder();
			for (Gene<VirusT> gene : strain.getGenes()) {
				GeneSequenceReads<VirusT> geneSeq = allGeneSequenceReads.get(gene);
				if (geneSeq == null) {
					concatSeq.append(StringUtils.repeat("...", gene.getAASize()));
				} else {
					concatSeq.append(geneSeq.getAlignedNAs(true));
				}
			}
			concatenatedSeq = concatSeq.toString();
		}
		return concatenatedSeq;
	}
	
	@Override
	public SequenceReadsHistogram<VirusT> getHistogram(
		final Double pcntLowerLimit,
		final Double pcntUpperLimit,
		final Integer numBins,
		final Boolean cumulative,
		final AggregationOption aggregatesBy) {
		return new SequenceReadsHistogram<>(
			getAllGeneSequenceReads(),
			pcntLowerLimit, pcntUpperLimit,
			numBins, cumulative, aggregatesBy);
	}

	@Override
	public SequenceReadsHistogram<VirusT> getHistogram(
		final Double pcntLowerLimit,
		final Double pcntUpperLimit,
		final Double[] binTicks,
		final Boolean cumulative,
		final AggregationOption aggregatesBy) {
		return new SequenceReadsHistogram<>(
			getAllGeneSequenceReads(),
			pcntLowerLimit / 100, pcntUpperLimit / 100,
			binTicks, cumulative, aggregatesBy);
	}

	@Override
	public SequenceReadsHistogramByCodonReads<VirusT> getHistogramByCodonReads(
		final Long[] codonReadsCutoffs,
		final AggregationOption aggregatesBy
	) {
		return new SequenceReadsHistogramByCodonReads<>(
			getAllGeneSequenceReads(),
			codonReadsCutoffs,
			aggregatesBy
		);
	}


	public MutationSet<VirusT> getMutations(final double minPrevalence, final long minCodonReads) {
		if (!isEmpty()) {
			return allGeneSequenceReads.values().stream()
				.map(gs -> gs.getMutations(minPrevalence, minCodonReads))
				.reduce((m1, m2) -> m1.mergesWith(m2))
				.get();
		} else {
			return new MutationSet<>();
		}
	}

	public MutationSet<VirusT> getMutations() {
		if (mutations == null) {
			mutations = getMutations(
				this.cutoffObj.getActualMinPrevalence(),
				this.cutoffObj.getMinCodonReads()
			);
		}
		return mutations;
	}

	public Long getMutationCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced());
	}

	public Long getUnusualMutationCount() {
		return getMutations().countIf(mut -> mut.isUnusual() && !mut.isUnsequenced());
	}

	public GenotypeResult<VirusT> getSubtypeResult() {
		if (!isEmpty() && subtypeResult == null) {
			subtypeResult = strain
				.getVirusInstance()
				.getGenotyper()
				.compareAll(
					getAssembledUnambiguousConsensus(true),
					strain.getAbsoluteFirstNA()
				);
		}
		return subtypeResult;
	}

	public BoundGenotype<VirusT> getBestMatchingSubtype() {
		if (isEmpty()) {
			return null;
		}
		return getSubtypeResult().getBestMatch();
	}

	public String getSubtypeText() {
		if (isEmpty()) {
			return "NA";
		}
		return getBestMatchingSubtype().getDisplay();
	}

}
