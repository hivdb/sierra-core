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
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.google.common.collect.Lists;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MultiCodonsMutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogram.AggregationOption;
import edu.stanford.hivdb.seqreads.SequenceReadsHistogram.WithSequenceReadsHistogram;
import edu.stanford.hivdb.sequences.GeneRegions;
import edu.stanford.hivdb.utilities.CodonUtils;

public class GeneSequenceReads<VirusT extends Virus<VirusT>> implements WithSequenceReadsHistogram<VirusT>, WithGene<VirusT> {

	private final Gene<VirusT> gene;
	private final int firstAA;
	private final int lastAA;
	private final List<PositionCodonReads<VirusT>> posCodonReads;
	private final CutoffCalculator<VirusT> cutoffObj;
	private MutationSet<VirusT> mutations;
	private transient DescriptiveStatistics readDepthStats;
	private transient GeneRegions<VirusT> unseqRegions;

	public GeneSequenceReads(
		final Gene<VirusT> gene,
		final List<PositionCodonReads<VirusT>> posCodonReads,
		final CutoffCalculator<VirusT> cutoffObj
	) {
		this.gene = gene;
		this.cutoffObj = cutoffObj;
		this.firstAA = Math.max(1, (int) posCodonReads.get(0).getPosition());
		this.lastAA = Math.min(gene.getAASize(), (int) posCodonReads.get(posCodonReads.size() - 1).getPosition());
		this.posCodonReads = Collections.unmodifiableList(
			posCodonReads
			.stream()
			// remove illegal positions
			.filter(pcr -> {
				long pos = pcr.getPosition();
				return pos >= this.firstAA && pos <= this.lastAA;
			})
			.collect(Collectors.toList())
		);
	}

	/** initializes GeneSequence without specify gene
	 *
	 * Warning: This constructor is only intended to use internally
	 *
	 * @param posCodonReads
	 * @param minPrevalence
	 */
	protected GeneSequenceReads(
		final List<PositionCodonReads<VirusT>> posCodonReads,
		final CutoffCalculator<VirusT> cutoffObj
	) {
		this(posCodonReads.get(0).getGene(), posCodonReads, cutoffObj);
	}

	@Override
	public Gene<VirusT> getGene() { return gene; }

	public int getFirstAA() { return firstAA; }
	public int getLastAA() { return lastAA; }
	public int getSize() { return lastAA - firstAA; }
	public int getNumPositions() { return posCodonReads.size(); }
	public List<PositionCodonReads<VirusT>> getAllPositionCodonReads() { return posCodonReads; }

	public MutationSet<VirusT> getMutations(final double minPrevalence, final long minCodonReads) {
		if (
			minPrevalence != cutoffObj.getActualMinPrevalence()||
			minCodonReads != cutoffObj.getMinCodonReads() ||
			mutations == null
		) {
			List<Mutation<VirusT>> myMutations = new ArrayList<>();
			long prevPos = firstAA - 1;
			for (PositionCodonReads<VirusT> pcr : posCodonReads) {
				long curPos = pcr.getPosition();
				for (Long pos = prevPos + 1; pos < curPos; pos ++) {
					// add unsequenced regions
					myMutations.add(MultiCodonsMutation.initUnsequenced(
						gene, pos.intValue()
					));
				}
				prevPos = curPos;
				Mutation<VirusT> mut = MultiCodonsMutation
					.fromPositionCodonReads(pcr, minPrevalence, minCodonReads);
				if (mut != null) {
					myMutations.add(mut);
				}
			}
			if (
				minPrevalence == cutoffObj.getActualMinPrevalence() &&
				minCodonReads == cutoffObj.getMinCodonReads()
			) {
				mutations = new MutationSet<>(myMutations);
			}
			else {
				return new MutationSet<>(myMutations);
			}
		}
		return mutations;
	}

	public Double getMedianReadDepth() {
		Median median = new Median();
		double[] ReadDepths = (
			posCodonReads.stream()
			.mapToDouble(read -> read.getTotalReads())
			.toArray()
		);
		double medianReadDepth = -1;
		if (ReadDepths.length > 0) {
			medianReadDepth = median.evaluate(ReadDepths);
		}
		return medianReadDepth;
	}

	@Override
	public SequenceReadsHistogram<VirusT> getHistogram(
		final Double pcntLowerLimit,
		final Double pcntUpperLimit,
		final Integer numBins,
		final Boolean cumulative,
		final AggregationOption aggregatesBy) {
		return new SequenceReadsHistogram<VirusT>(
			Lists.newArrayList(this),
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
		return new SequenceReadsHistogram<VirusT>(
			Lists.newArrayList(this),
			pcntLowerLimit / 100, pcntUpperLimit / 100,
			binTicks, cumulative, aggregatesBy);
	}

	public DescriptiveStatistics getReadDepthStats() {
		if (readDepthStats == null) {
			double[] readDepthArray = posCodonReads.stream()
				.mapToDouble(pcr -> pcr.getTotalReads())
				.toArray();
			
			if (readDepthArray.length > 2) {
				readDepthStats = new DescriptiveStatistics(readDepthArray);
			}
			else {
				readDepthStats = new DescriptiveStatistics(new double[] {0, 0, 0});
			}
		}
		return readDepthStats;
	}

	public MutationSet<VirusT> getMutations() {
		return getMutations(
			this.cutoffObj.getActualMinPrevalence(),
			this.cutoffObj.getMinCodonReads()
		);
	}
	
	public Long getMutationCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced());
	}

	public Long getUnusualMutationCount() {
		return getMutations().countIf(mut -> mut.isUnusual() && !mut.isUnsequenced());
	}

	/** Returns consensus sequence aligned to reference.
	 *  All insertions are removed from the result.
	 *
	 * @param autoComplete specify <tt>true</tt> to prepend and/or append
	 * wildcard "." to incomplete sequence
	 * @return the aligned consensus sequence
	 */
	public String getAlignedNAs(boolean autoComplete) {
		return getAlignedNAs(
			this.cutoffObj.getActualMinPrevalence(),
			this.cutoffObj.getMinCodonReads(),
			autoComplete
		);
	}

	/** Returns consensus sequence aligned to reference.
	 *  All insertions are removed from the result. 1a/1b frameshifts are preserved.
	 *
	 * @param threshold specify the minimal prevalence requirement for
	 * creating codon consensus
	 * @param autoComplete specify <tt>true</tt> to prepend and/or append
	 * wildcard "." to incomplete sequence
	 * @return the aligned consensus sequence
	 */
	public String getAlignedNAs(double pcntThreshold, long countThreshold, boolean autoComplete) {
		StringBuilder seq = new StringBuilder();
		if (autoComplete) {
			seq.append(StringUtils.repeat("...", firstAA - 1));
		}
		long prevPos = firstAA - 1;
		for (PositionCodonReads<VirusT> pcr : posCodonReads) {
			long curPos = pcr.getPosition();
			if (curPos - prevPos > 1) {
				seq.append(StringUtils.repeat("...", (int) (curPos - prevPos - 1)));
			}
			prevPos = curPos;
			seq.append(pcr.getCodonConsensusWithoutIns(pcntThreshold, countThreshold));
		}
		if (autoComplete) {
			seq.append(StringUtils.repeat("...", gene.getAASize() - lastAA));
		}
		return seq.toString();
	}
	
	public CutoffCalculator<VirusT> getCutoffObj() { return cutoffObj; }
	
	/** Returns consensus sequence aligned to subtype B reference without
	 *  initial and trailing "." for incomplete sequence. All insertions are
	 *  removed from the result. The result is equivalent to the result of
	 *  <tt>getAlignedSequence(false)</tt>.
	 *
	 * @return the aligned consensus NA sequence
	 */
	public String getAlignedNAs() {
		return getAlignedNAs(false);
	}

	/** Returns consensus sequence aligned to subtype B reference in amino
	 *  acid form. Unsequenced region(s) are ignored. All insertions are
	 *  removed from the result.
	 *
	 * @return the aligned consensus AA sequence
	 */
	public String getAlignedAAs() {
		return CodonUtils.simpleTranslate(
			this.getAlignedNAs(false), firstAA, gene.getRefSequence());
	}

	public GeneRegions<VirusT> getUnsequencedRegions() {
		if (unseqRegions == null) {
			unseqRegions = GeneRegions.newUnsequencedRegions(gene, firstAA, lastAA, getMutations());
		}
		return unseqRegions;
	}

}
