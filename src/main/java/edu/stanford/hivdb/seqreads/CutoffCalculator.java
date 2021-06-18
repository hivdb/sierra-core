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

	private static final int EMPTY = 0; 
	private static final int MIXTURE = -1;
	private final Double actualMixturePcnt;
	private final Double actualMinPrevalence;
	private final Double maxMixturePcnt;
	private final Double minPrevalence;
	private final Long minCodonReads;
	private final Long minPositionReads;
	
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
	
	
	public CutoffCalculator(
		List<PositionCodonReads<VirusT>> allReads,
		Double maxMixturePcnt,
		Double minPrevalence,
		Long minCodonReads,
		Long minPositionReads
	) {
		if (maxMixturePcnt == null || maxMixturePcnt < 0 || maxMixturePcnt > 1) {
			maxMixturePcnt = 1.;
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
		this.maxMixturePcnt = maxMixturePcnt;
		this.minPrevalence = minPrevalence;
		this.minCodonReads = minCodonReads;
		this.minPositionReads = minPositionReads;
		
		List<Pair<GenePosition<VirusT>, CodonReads<VirusT>>> sortedCodonReads = allReads
			.stream()
			.filter(pcr -> pcr.getTotalReads() >= this.minPositionReads)
			.flatMap(
				pcr -> (
					pcr.getCodonReads()
					.stream()
					.filter(cdr -> cdr.getReads() >= this.minCodonReads)
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
		double prevProportion = 1.;
		double prevMixturePcnt = 0.;
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
			/* System.out.println("prevCodon: " + prevCodon + " " +
				(prevCodon == null ? 0 : prevCodon.stream().filter(na -> na == MIXTURE).count()) +
				"  mergedCodon: " + mergedCodon + " " +
				mergedCodon.stream().filter(na -> na == MIXTURE).count() +
				"  mixture: " +				
				numMixtures +
				"  total: " + numNAs); */
			double curMixturePcnt = numMixtures / numNAs;
			double curProportion = cdr.getProportion();
			if (curMixturePcnt > maxMixturePcnt || curProportion < minPrevalence) {
				// System.out.println(codonLookup);
				break;
			}
			else {
				prevMixturePcnt = curMixturePcnt;
				prevProportion = curProportion;
			}
		}
		this.actualMixturePcnt = prevMixturePcnt;
		this.actualMinPrevalence = prevProportion;
	}
	
	public Double getMaxMixturePcnt() { return maxMixturePcnt; }
	public Double getMinPrevalence() { return minPrevalence; }
	public Long getMinCodonReads() { return minCodonReads; }
	public Long getMinPositionReads() { return minPositionReads; }
	
	public Double getActualMixturePcnt() { return actualMixturePcnt; }
	public Double getActualMinPrevalence() { return actualMinPrevalence; }
	
}