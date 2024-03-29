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
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.base.Strings;

import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;


public class PositionCodonReads<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {

	private final Gene<VirusT> gene;
	private final int position;
	private final long totalReads;
	private final Map<String, Long> allCodonReads;

	public PositionCodonReads(
		final Gene<VirusT> gene,
		final int position,
		final long totalReads,
		final Map<String, Long> allCodonReads) {
		this.gene = gene;
		this.position = position;
		this.totalReads = totalReads;
		this.allCodonReads = Collections.unmodifiableMap(
			allCodonReads.entrySet().stream()
			.sorted(Comparator.comparingLong(Map.Entry<String, Long>::getValue).reversed())
			.collect(Collectors.toMap(
				e -> e.getKey(),
				e -> e.getValue(),
				(e1, e2) -> e1,
				LinkedHashMap::new))
		);
	}

	@Override
	public Strain<VirusT> getStrain() { return gene.getStrain(); }
	
	@Override
	public Gene<VirusT> getGene() { return gene; }
	
	@Override
	public String getAbstractGene() { return gene.getAbstractGene(); }

	public long getPosition() { return position; }
	
	public GenePosition<VirusT> getGenePosition() { return new GenePosition<VirusT>(gene, position); }
	public long getTotalReads() { return totalReads; }
	public List<CodonReads<VirusT>> getCodonReads() {
		return getCodonReads(false, 2., -1.);
	}
	public List<CodonReads<VirusT>> getCodonReads(
		boolean mutationOnly,
		double maxProportion,
		double minProportion
	) {
		return allCodonReads.entrySet().stream()
			.map(e -> new CodonReads<>(gene, position, e.getKey(), e.getValue(), totalReads))
			.filter(cr -> !cr.isAmbiguous())  // Illumina, ONT nor Ion Torrent would report ambiguous NAs
			.filter(cr -> mutationOnly ? !cr.isReference() : true)
			.filter(cr -> {
				double prop = cr.getProportion();
				return prop >= minProportion && prop <= maxProportion;
			})
			.collect(Collectors.toList());
	}

	protected long calcMinReads(double minPrevalence, long minCodonReads) {
		return Math.max(
			((Double) Math.ceil(totalReads * minPrevalence)).longValue(),
			minCodonReads
		);
	}

	public List<CodonReads<VirusT>> getCodonReadsUsingThreshold(
		double minPrevalence, long minCodonReads
	) {
		long minReads = calcMinReads(minPrevalence, minCodonReads);
		return (
			getCodonReads().stream()
			.filter(cdr -> cdr.getReads() >= minReads)
			.collect(Collectors.toList())
		);
		
	}
	
	public Map<String, Double> getCodonWithPrevalence(double minPrevalence, long minCodonReads) {
		return (
			getCodonReadsUsingThreshold(minPrevalence, minCodonReads)
			.stream()
			.collect(Collectors.toMap(
				cdr -> cdr.getCodon(),
				cdr -> cdr.getProportion(),
				(e1, e2) -> e1,
				LinkedHashMap::new))
		);
	}
	
	public static <VirusT extends Virus<VirusT>> String getCodonConsensusWithoutIns(List<CodonReads<VirusT>> codonReads) {
		return (getCodonConsensusWithIns(codonReads) + "---").substring(0, 3);
	}
	
	public static <VirusT extends Virus<VirusT>> String getCodonConsensusWithIns(List<CodonReads<VirusT>> codonReads) {
		List<String> codons = (
			codonReads
			.stream()
			.sorted((cdr1, cdr2) -> {
				int cmp = cdr2.getReads().compareTo(cdr1.getReads());
				if (cmp != 0) { return cmp; }
				return cdr1.getCodon().compareTo(cdr2.getCodon());
			})
			.map(cdr -> cdr.getCodon())
			.collect(Collectors.toList()));
		if (codons.isEmpty()) {
			// do not return null
			return "NNN";
		}
		String codon = CodonUtils.getMergedCodon(codons);
		if (codon.length() < 3) {
			codon = codon + Strings.repeat("-", 3 - codon.length());
		}
		return codon;
	}
	
	public static <VirusT extends Virus<VirusT>> String getMajorCodonWithoutIns(List<CodonReads<VirusT>> codonReads) {
		return (getMajorCodonWithIns(codonReads) + "---").substring(0, 3);
	}
	
	public static <VirusT extends Virus<VirusT>> String getMajorCodonWithIns(List<CodonReads<VirusT>> codonReads) {
		if (codonReads.isEmpty()) {
			// do not raise error
			return "NNN";
		}
		else {
			String codon = (
				codonReads
				.stream()
				.sorted((cdr1, cdr2) -> {
					int cmp = cdr2.getReads().compareTo(cdr1.getReads());
					if (cmp != 0) { return cmp; }
					return cdr1.getCodon().compareTo(cdr2.getCodon());
				})
				.map(cdr -> cdr.getCodon())
				.findFirst()
				.get()
			);
			if (codon.length() < 3) {
				codon = codon + Strings.repeat("-", 3 - codon.length());
			}
			return codon;
		}
	}
	
	public String getCodonConsensusWithoutIns(double minPrevalence, long minCodonReads) {
		return getCodonConsensusWithoutIns(
			getCodonReadsUsingThreshold(minPrevalence, minCodonReads)
		);
	}

	public String getCodonConsensusWithIns(double minPrevalence, long minCodonReads) {
		return getCodonConsensusWithIns(
			getCodonReadsUsingThreshold(minPrevalence, minCodonReads)
		);
	}
	
	public String getMajorCodonWithoutIns(long minCodonReads) {
		return getMajorCodonWithoutIns(getCodonReadsUsingThreshold(0., minCodonReads));
	}

	public String getMajorCodonWithIns(long minCodonReads) {
		return getMajorCodonWithIns(getCodonReadsUsingThreshold(0., minCodonReads));
	}
	
	/**
	 * Get Extended map
	 * 
	 * @param mutationOnly	Only mutations
	 * @param maxProp		Maximum proportion
	 * @param minProp		Minimum proportion
	 * @return 				Map&lt;String, Object&gt;
	 */
	public Map<String, Object> extMap(
		boolean mutationOnly, double maxProp, double minProp
	) {
		Map<String, Object> result = new LinkedHashMap<>();
		result.put("position", getPosition());
		result.put("totalReads", getTotalReads());
		result.put(
			"codonReads",
			getCodonReads(mutationOnly, maxProp, minProp)
			.stream().map(cr -> cr.extMap())
			.collect(Collectors.toList()));
		return result;
	}
}
