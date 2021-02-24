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

import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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
		this.allCodonReads = allCodonReads.entrySet().stream()
			.sorted(Comparator.comparingLong(Map.Entry<String, Long>::getValue).reversed())
			.collect(Collectors.toMap(
				e -> e.getKey(),
				e -> e.getValue(),
				(e1, e2) -> e1,
				LinkedHashMap::new));
	}

	@Override
	public Strain<VirusT> getStrain() { return gene.getStrain(); }
	
	@Override
	public Gene<VirusT> getGene() { return gene; }
	
	@Override
	public String getAbstractGene() { return gene.getAbstractGene(); }

	public long getPosition() { return position; }
	
	public GenePosition<VirusT> getGenePositon() { return new GenePosition<VirusT>(gene, position); }
	public long getTotalReads() { return totalReads; }
	public List<CodonReads<VirusT>> getCodonReads() {
		return getCodonReads(false, 1., .0);
	}
	public List<CodonReads<VirusT>> getCodonReads(
		boolean mutationOnly,
		double maxProportion,
		double minProportion
	) {
		return allCodonReads.entrySet().stream()
			.map(e -> new CodonReads<>(gene, position, e.getKey(), e.getValue(), totalReads))
			.filter(cr -> cr.getAminoAcid() != 'X')
			.filter(cr -> mutationOnly ? !cr.isReference() : true)
			.filter(cr -> {
				double prop = cr.getProportion();
				return prop > minProportion && prop < maxProportion;
			})
			.collect(Collectors.toList());
	}
	
	public Map<String, Double> getCodonWithPrevalence(double minPrevalence, long minCodonCount) {
		long minReads = Math.max(Math.round(totalReads * minPrevalence + 0.5), minCodonCount);
		return allCodonReads.entrySet().stream()
			.filter(e -> e.getValue() > minReads)
			.collect(Collectors.toMap(
				e -> e.getKey(),
				e -> Double.valueOf(e.getValue() / totalReads),
				(e1, e2) -> e1,
				LinkedHashMap::new));
	}

	public String getCodonConsensus(double minPrevalence, long minCodonCount) {
		List<String> codons = (
			getCodonWithPrevalence(minPrevalence, minCodonCount).keySet().stream()
			.map(cd -> cd.substring(0, 3))
			.collect(Collectors.toList()));
		if (codons.isEmpty()) {
			// do not return null
			return "NNN";
		}
		return CodonUtils.getMergedCodon(codons);
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
