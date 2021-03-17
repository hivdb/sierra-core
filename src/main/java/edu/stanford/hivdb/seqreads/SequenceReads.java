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
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogram.AggregationOption;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogram.WithSequenceReadsHistogram;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogramByCodonReads.WithSequenceReadsHistogramByCodonReads;
import edu.stanford.hivdb.sequences.SeqUtils;
import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceReads<VirusT extends Virus<VirusT>> implements WithSequenceReadsHistogram<VirusT>, WithSequenceReadsHistogramByCodonReads<VirusT> {

	private final static int HXB2_PR_FIRST_NA = 2253;
	private final static double MIN_PREVALENCE_FOR_SUBTYPING = 0.2;
	private final Strain<VirusT> strain;
	private final Map<Gene<VirusT>, GeneSequenceReads<VirusT>> allGeneSequenceReads;
	private final List<OneCodonReadsCoverage<VirusT>> codonReadsCoverage;
	private final CutoffSuggestion<VirusT> cutoffSuggestion;
	private final Double proportionTrimmedPositions;
	private final String name;
	private GenotypeResult<VirusT> subtypeResult;
	private MutationSet<VirusT> mutations;
	private String concatenatedSeq;
	private Double mixturePcnt;
	private Double minPrevalence;
	private Long minCodonReads;
	private Long minPositionReads;
	private transient DescriptiveStatistics readDepthStats;
	private transient DescriptiveStatistics readDepthStatsDRP;
	private transient List<ValidationResult> validationResults;
	
	public static <VirusT extends Virus<VirusT>> SequenceReads<VirusT> fromCodonReadsTable(
			String name, Strain<VirusT> strain, List<PositionCodonReads<VirusT>> allReads,
			Double minPrevalence, Long minCodonReads, Long minPositionReads) {
		// TODO: dynamic cutoff
		CutoffSuggestion<VirusT> cutoffSuggestion = new CutoffSuggestion<>(allReads);
		double finalMinPrevalence = minPrevalence >= 0 ? minPrevalence : 0; // TODO: cutoffSuggestion.getStricterLimit()
		long finalMinCodonReads = minCodonReads > 0 ? minCodonReads : 1L;
		long finalMinPositionReads = minPositionReads > 0 ? minPositionReads : 1L;
		
		// sort by genePosition first
		allReads.sort((o1, o2) -> o1.getGenePositon().compareTo(o2.getGenePositon()));
		
		List<PositionCodonReads<VirusT>> filteredAllReads = allReads.stream()
			// remove all codons with their read depth < minPositionReads
			.filter(read -> read.getTotalReads() >= finalMinPositionReads)
			.collect(Collectors.toList());

		Map<Gene<VirusT>, GeneSequenceReads<VirusT>> geneSequences = filteredAllReads
			.stream()
			.collect(
				Collectors.groupingBy(
					PositionCodonReads::getGene,
					TreeMap::new,
					Collectors.collectingAndThen(
						Collectors.toList(),
						list -> new GeneSequenceReads<>(list, finalMinPrevalence, finalMinCodonReads)
					)
				)
			);
		List<OneCodonReadsCoverage<VirusT>> codonReadsCoverage = allReads.stream()
			.map(read -> new OneCodonReadsCoverage<>(
				read.getGene(),
				read.getPosition(),
				read.getTotalReads(),
				read.getTotalReads() < finalMinPositionReads
			))
			.collect(Collectors.toList());
		double proportionTrimmedPositions = 1. - (double) filteredAllReads.size() / (double) allReads.size();

		// TODO: add support for HIV2
		return new SequenceReads<>(
				name,
				strain,
				geneSequences,
				finalMinPrevalence,
				finalMinCodonReads,
				finalMinPositionReads,
				cutoffSuggestion,
				codonReadsCoverage,
				proportionTrimmedPositions);
	}

	protected SequenceReads(
			final String name,
			final Strain<VirusT> strain,
			final Map<Gene<VirusT>, GeneSequenceReads<VirusT>> allGeneSequenceReads,
			final double minPrevalence,
			final long minCodonReads,
			final long minPositionReads,
			final CutoffSuggestion<VirusT> cutoffSuggestion,
			final List<OneCodonReadsCoverage<VirusT>> codonReadsCoverage,
			final double proportionTrimmedPositions) {
		this.name = name;
		this.strain = strain;
		this.allGeneSequenceReads = Collections.unmodifiableMap(allGeneSequenceReads);
		this.minPrevalence = minPrevalence;
		this.minCodonReads = minCodonReads;
		this.minPositionReads = minPositionReads;
		this.cutoffSuggestion = cutoffSuggestion;
		this.codonReadsCoverage = Collections.unmodifiableList(codonReadsCoverage);
		this.proportionTrimmedPositions = proportionTrimmedPositions;
	}
	
	public List<OneCodonReadsCoverage<VirusT>> getCodonReadsCoverage() {
		return codonReadsCoverage;
	}
	
	public Integer getSize() {
		int totalSize = 0;
		for (Gene<VirusT> gene : allGeneSequenceReads.keySet()) {
			GeneSequenceReads<VirusT> gsr = allGeneSequenceReads.get(gene);
			totalSize += gsr.getSize();
		}
		return totalSize;
	}

	public Double getCutoffSuggestionLooserLimit() { return cutoffSuggestion.getLooserLimit(); }

	public Double getCutoffSuggestionStricterLimit() { return cutoffSuggestion.getStricterLimit(); }

	public Double getProportionTrimmedPositions() { return proportionTrimmedPositions; }
	
	public String getName() { return name; }
	
	public Strain<VirusT> getStrain() { return strain; }

	public boolean isEmpty() { return allGeneSequenceReads.isEmpty(); }

	public double getMinPrevalence() { return minPrevalence; }
	
	public double getMinCodonReads() { return minCodonReads; }

	public long getMinPositionReads() { return minPositionReads; }
	
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
					.filter(pcr -> pcr.getGenePositon().isDrugResistancePosition())
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
	
	public String getConcatenatedSeqForSubtyping() {
		StringBuilder concatSeq = new StringBuilder();
		for (Gene<VirusT> gene : strain.getGenes()) {
			GeneSequenceReads<VirusT> geneSeq = allGeneSequenceReads.get(gene);
			if (geneSeq == null) {
				concatSeq.append(StringUtils.repeat("...", gene.getAASize()));
			} else {
				concatSeq.append(
					geneSeq.getAlignedNAs(MIN_PREVALENCE_FOR_SUBTYPING, 0, true));
			}
		}
		return concatSeq.toString();
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
		if (!isEmpty() && mutations == null) {
			mutations = getMutations(this.minPrevalence, this.minCodonReads);
		}
		return mutations;
	}

	public GenotypeResult<VirusT> getSubtypeResult() {
		if (!isEmpty() && subtypeResult == null) {
			subtypeResult = strain.getVirusInstance().getGenotyper().compareAll(
				getConcatenatedSeqForSubtyping(), HXB2_PR_FIRST_NA);
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

	public double getMixturePcnt() {
		if (mixturePcnt == null) {
			StringBuilder concatSeq = new StringBuilder();
			for (GeneSequenceReads<VirusT> geneSeqReads : allGeneSequenceReads.values()) {
				concatSeq.append(geneSeqReads.getAlignedNAs(false));
			}
			mixturePcnt = SeqUtils.mixturePcnt(concatSeq.toString());
		}
		return mixturePcnt;
	}

}
