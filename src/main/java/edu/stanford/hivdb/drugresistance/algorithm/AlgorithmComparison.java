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

package edu.stanford.hivdb.drugresistance.algorithm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.stanford.hivdb.drugresistance.GeneDR;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.MutationSet;

public class AlgorithmComparison<VirusT extends Virus<VirusT>> {

	protected static <VirusT extends Virus<VirusT>> List<GeneDR<VirusT>> calcGeneDRFromAlgorithms(
		Gene<VirusT> gene,
		MutationSet<VirusT> mutations,
		Collection<DrugResistanceAlgorithm<VirusT>> algorithms
	) {
		return algorithms
			.stream()
			.map(alg -> new GeneDR<>(gene, mutations, alg))
			.collect(Collectors.toList());
	}

	private final Collection<DrugResistanceAlgorithm<VirusT>> algorithms;
	private final List<ASIDrugSusc<VirusT>> comparisonResults = new ArrayList<>();
	private final transient Map<Gene<VirusT>, List<GeneDR<VirusT>>> geneDRMap;

	public List<ASIDrugSusc<VirusT>> getComparisonResults() { return comparisonResults; }

	public AlgorithmComparison (
		MutationSet<VirusT> allMutations,
		Collection<DrugResistanceAlgorithm<VirusT>> algorithms
	) {
		this.algorithms = algorithms;
		this.geneDRMap = new LinkedHashMap<>();
		Map<Gene<VirusT>, MutationSet<VirusT>> mutationsByGroup = allMutations.groupByGene();
		for (Gene<VirusT> gene : mutationsByGroup.keySet()) {
			final MutationSet<VirusT> mutations = mutationsByGroup.get(gene);
			geneDRMap.put(gene, calcGeneDRFromAlgorithms(gene, mutations, algorithms));
			compareResults(gene);
		}
	}

	private void compareResults(Gene<VirusT> gene) {
		for (GeneDR<VirusT> geneDR : this.geneDRMap.get(gene)) {
			comparisonResults.addAll(geneDR.getDrugSuscs());
		}
	}

	public List<GeneDR<VirusT>> getGeneDR(Gene<VirusT> gene) {
		return this.geneDRMap.get(gene);
	}
	
	public Collection<DrugResistanceAlgorithm<VirusT>> getAlgorithms() {
		return algorithms;
	}

}