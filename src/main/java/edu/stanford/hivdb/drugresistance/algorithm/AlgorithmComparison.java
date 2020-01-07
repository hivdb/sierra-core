/*

    Copyright (C) 2017 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.drugresistance.algorithm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.drugs.DrugResistanceAlgorithm;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.MutationSet;

public class AlgorithmComparison<VirusT extends Virus<VirusT>> {

	public static class ComparableDrugScore<VirusT extends Virus<VirusT>> {
		private Drug<VirusT> drug;
		private String algorithm;
		private SIREnum SIR;
		private String interpretation;
		private String explanation;

		public ComparableDrugScore(
				Drug<VirusT> drug, String alg, SIREnum SIR,
				String interpretation, String explanation) {
			this.drug = drug;
			this.algorithm = alg;
			this.SIR = SIR;
			this.interpretation = interpretation;
			this.explanation = explanation;
		}
		
		public Drug<VirusT> getDrug() { return drug; }
		public String getAlgorithm() { return algorithm; }
		public SIREnum getSIR() { return SIR; }
		public String getInterpretation() { return interpretation; }
		public String getExplanation() { return explanation; }

		@Override
		public String toString() {
			return String.format("%s (%s): %s", drug, algorithm, SIR);
		}

		@Override
			public boolean equals(Object o) {
				if (o == this) { return true; }
				if (o == null) { return false; }
				if (!(o instanceof ComparableDrugScore)) { return false;}
				ComparableDrugScore<?> ds = (ComparableDrugScore<?>) o;

				// isDeletion and isInsertion is related to aas
				return new EqualsBuilder()
					.append(drug, ds.drug)
					.append(algorithm, ds.algorithm)
					.append(SIR, ds.SIR)
					.append(interpretation, ds.interpretation)
					.append(explanation, ds.explanation)
					.isEquals();
			}

			@Override
			public int hashCode() {
				return new HashCodeBuilder(24757, 43569)
					.append(drug)
					.append(algorithm)
					.append(SIR)
					.append(interpretation)
					.append(explanation)
					.toHashCode();
			}
	}

	protected static <VirusT extends Virus<VirusT>> List<AsiResult<VirusT>> calcAsiListFromAlgorithms(
		Gene<VirusT> gene,
		MutationSet<VirusT> mutations,
		Collection<DrugResistanceAlgorithm<VirusT>> algorithms
	) {
		return algorithms
			.stream()
			.map(alg -> new AsiResult<>(gene, mutations, alg))
			.collect(Collectors.toList());
	}

	private final Collection<DrugResistanceAlgorithm<VirusT>> algorithms;
	private final List<ComparableDrugScore<VirusT>> comparisonResults = new ArrayList<>();
	private final transient Map<Gene<VirusT>, List<AsiResult<VirusT>>> asiListMap;

	public List<ComparableDrugScore<VirusT>> getComparisonResults() { return comparisonResults; }

	public AlgorithmComparison (
		MutationSet<VirusT> allMutations,
		Collection<DrugResistanceAlgorithm<VirusT>> algorithms
	) {
		this.algorithms = algorithms;
		this.asiListMap = new LinkedHashMap<>();
		Map<Gene<VirusT>, MutationSet<VirusT>> mutationsByGroup = allMutations.groupByGene();
		for (Gene<VirusT> gene : mutationsByGroup.keySet()) {
			final MutationSet<VirusT> mutations = mutationsByGroup.get(gene);
			asiListMap.put(gene, calcAsiListFromAlgorithms(gene, mutations, algorithms));
			compareResults(gene);
		}
	}

	private void compareResults (Gene<VirusT> gene) {

		for (AsiResult<VirusT> asiObj : this.asiListMap.get(gene)) {
			String algName = asiObj.getAlgorithmName();

			for (DrugClass<VirusT> drugClass : gene.getDrugClasses()) {

				for (Drug<VirusT> drug : drugClass.getDrugs()) {
					AsiDrugComparableResult result = asiObj.getDrugComparableResult(drug);
					comparisonResults.add(new ComparableDrugScore<>(
						drug, algName, result.SIR,
						result.Interpretation, result.Explanation));
				}
			}
		}
	}

	public List<AsiResult<VirusT>> getAsiList(Gene<VirusT> gene) {
		return this.asiListMap.get(gene);
	}
	
	public Collection<DrugResistanceAlgorithm<VirusT>> getAlgorithms() {
		return algorithms;
	}

}
