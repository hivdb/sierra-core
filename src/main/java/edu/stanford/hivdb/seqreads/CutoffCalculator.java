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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Pair;

import edu.stanford.hivdb.mutations.CodonReads;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.viruses.Virus;

public class CutoffCalculator<VirusT extends Virus<VirusT>> {
	
	public static final Double minKeyPointMixtureRateFoldDiff = 1.1;
	public static final Double minKeyPointMinPrevalenceRateFoldDiff = 1.1;
	
	public static class CutoffKeyPoint {
		private final Double mixtureRate;
		private final Double minPrevalence;
		private final Boolean aboveMixtureRateThreshold;
		private final Boolean belowMinPrevalenceThreshold;

		public CutoffKeyPoint(
			double mixtureRate,
			double minPrevalence,
			boolean aboveMixtureRateThreshold,
			boolean belowMinPrevalenceThreshold
		) {
			this.mixtureRate = mixtureRate;
			this.minPrevalence = minPrevalence;
			this.aboveMixtureRateThreshold = aboveMixtureRateThreshold;
			this.belowMinPrevalenceThreshold = belowMinPrevalenceThreshold;
		}
		
		public Double getMixtureRate() {
			return mixtureRate;
		}
		
		public Double getMinPrevalence() {
			return minPrevalence;
		}
		
		public Boolean isAboveMixtureRateThreshold() {
			return aboveMixtureRateThreshold;
		}
		
		public Boolean isBelowMinPrevalenceThreshold() {
			return belowMinPrevalenceThreshold;
		}
		
	}

	private static final int EMPTY = 0; 
	private static final int MIXTURE = -1;
	private final Double maxMixtureRate;
	private final Double minPrevalence;
	private final Long minCodonReads;
	private final Long minPositionReads;
	private final CutoffKeyPoint selectedCutoff;
	private final List<CutoffKeyPoint> cutoffKeyPoints;
	
	protected static List<Integer> mergeCodon(List<Integer> codonA, List<Integer> codonB) {
		if (codonA == null) {
			return codonB;
		}
		if (codonB == null) {
			return codonA;
		}
		List<Integer> newCodon = new ArrayList<>(List.of(EMPTY, EMPTY, EMPTY));
		int codonASize = codonA.size();
		int codonBSize = codonB.size();
		for (int i = 0; i < 3; i ++) {
			int naA = codonASize > i ? codonA.get(i) : EMPTY;
			int naB = codonBSize > i ? codonB.get(i) : EMPTY;
			if (naA == EMPTY) {
				newCodon.set(i, naB);
			}
			else if (naB == EMPTY || naA == naB) {
				newCodon.set(i, naA);
			}
			else {
				newCodon.set(i, MIXTURE);
			}
		}
		return newCodon;
	}
	
	protected static
	<VirusT extends Virus<VirusT>> List<CutoffKeyPoint>
	calcCutoffKeyPoints(
		List<PositionCodonReads<VirusT>> allReads,
		Double maxMixtureRate,
		Double minPrevalence,
		Long minCodonReads,
		Long minPositionReads
	) {
		List<Pair<GenePosition<VirusT>, CodonReads<VirusT>>> sortedCodonReads = allReads
			.stream()
			.filter(pcr -> pcr.getTotalReads() >= minPositionReads)
			.flatMap(
				pcr -> (
					pcr.getCodonReads()
					.stream()
					.filter(cdr -> cdr.getReads() >= minCodonReads)
					.map(cdr -> Pair.of(pcr.getGenePosition(), cdr))
				)
			)
			.sorted(
				(r1, r2) -> (
					r2.getRight().getProportion()
					.compareTo(
						r1.getRight().getProportion()
					)
				)
			)
			.collect(Collectors.toList());

		Map<GenePosition<VirusT>, List<Integer>> codonLookup = new HashMap<>();
		double numMixtures = 0;
		double numNAs = 0;
		double prevMixtureRate = 0.;
		double prevProportion = 1.;
		List<CutoffKeyPoint> cutoffKeyPoints = new ArrayList<>();
		CutoffKeyPoint prevKeyPoint = new CutoffKeyPoint(
			prevMixtureRate,
			prevProportion,
			prevMixtureRate > maxMixtureRate,
			prevProportion < minPrevalence
		);
		cutoffKeyPoints.add(prevKeyPoint);
		for (Pair<GenePosition<VirusT>, CodonReads<VirusT>> pcdr : sortedCodonReads) {
			GenePosition<VirusT> genePos = pcdr.getLeft();
			CodonReads<VirusT> cdr = pcdr.getRight();
			List<Integer> codon = (
				(cdr.getCodon() + "---").substring(0, 3).chars()
				.mapToObj(i -> i).collect(Collectors.toList())
			);
			List<Integer> prevCodon = codonLookup.getOrDefault(genePos, null);
			if (prevCodon == null) {
				numNAs += 3;
			}
			List<Integer> mergedCodon = mergeCodon(prevCodon, codon);
			codonLookup.put(genePos, mergedCodon);
			numMixtures += (
				- (prevCodon == null ? 0 : prevCodon.stream().filter(na -> na == MIXTURE).count())
				+ mergedCodon.stream().filter(na -> na == MIXTURE).count()
			);
			double curMixtureRate = numMixtures / numNAs;
			double curProportion = cdr.getProportion();
			if (
				(
					prevMixtureRate / prevKeyPoint.mixtureRate >= minKeyPointMixtureRateFoldDiff &&
					prevKeyPoint.minPrevalence / prevProportion >= minKeyPointMinPrevalenceRateFoldDiff
				) ||
				(
					curMixtureRate > maxMixtureRate ^  // xor
					curProportion < minPrevalence
				)
			) {
				prevKeyPoint = new CutoffKeyPoint(
					prevMixtureRate,
					prevProportion,
					prevMixtureRate > maxMixtureRate,
					prevProportion < minPrevalence
				);
				cutoffKeyPoints.add(prevKeyPoint);
			}
			prevMixtureRate = curMixtureRate;
			prevProportion = curProportion;
		}
		cutoffKeyPoints.add(new CutoffKeyPoint(
			prevMixtureRate,
			prevProportion,
			prevMixtureRate > maxMixtureRate,
			prevProportion < minPrevalence
		));
		return Collections.unmodifiableList(cutoffKeyPoints);
	}
	
	public CutoffCalculator(
		List<PositionCodonReads<VirusT>> allReads,
		Double maxMixtureRate,
		Double minPrevalence,
		Long minCodonReads,
		Long minPositionReads
	) {
		if (maxMixtureRate == null || maxMixtureRate < 0 || maxMixtureRate > 1) {
			maxMixtureRate = 1.;
		}
		if (minPrevalence == null || minPrevalence < 0 || minPrevalence > 1) {
			minPrevalence = 0.;
		}
		if (minCodonReads == null || minCodonReads < 1) {
			minCodonReads = 1L;
		}
		if (minPositionReads == null || minPositionReads < 1) {
			minPositionReads = 1L;
		}
		this.maxMixtureRate = maxMixtureRate;
		this.minPrevalence = minPrevalence;
		this.minCodonReads = minCodonReads;
		this.minPositionReads = minPositionReads;
		cutoffKeyPoints = calcCutoffKeyPoints(
			allReads,
			maxMixtureRate,
			minPrevalence,
			minCodonReads,
			minPositionReads
		);
		selectedCutoff = cutoffKeyPoints
			.stream()
			.filter(ckp -> (
				!ckp.isAboveMixtureRateThreshold() &&
				!ckp.isBelowMinPrevalenceThreshold()
			))
			.reduce((first, second) -> second)
			.get();
	}

	public Double getMaxMixtureRate() { return maxMixtureRate; }
	public Double getMinPrevalence() { return minPrevalence; }
	public Long getMinCodonReads() { return minCodonReads; }
	public Long getMinPositionReads() { return minPositionReads; }
	
	public Double getActualMixtureRate() { return selectedCutoff.getMixtureRate(); }
	public Double getActualMinPrevalence() { return selectedCutoff.getMinPrevalence(); }
	
	public List<CutoffKeyPoint> getCutoffKeyPoints() { return cutoffKeyPoints; }
	
}