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
package edu.stanford.hivdb.drugresistance.algorithm;

import java.io.Serializable;
import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.builder.CompareToBuilder;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.Precision;

import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.viruses.Virus;

public class ASIDrugSusc<VirusT extends Virus<VirusT>> implements Comparable<ASIDrugSusc<VirusT>>, Serializable {

	private static final long serialVersionUID = 6192833587557227063L;
	private final Drug<VirusT> drug;
	private final DrugResistanceAlgorithm<VirusT> algorithm;
	private final Double score;
	private final Integer level;
	private final String levelText;
	private final SIREnum sir;
	private final Map<MutationSet<VirusT>, Double> partialScores;
	private final String statement;
	private final boolean triggered;
	
	public ASIDrugSusc(
		Drug<VirusT> drug,
		DrugResistanceAlgorithm<VirusT> algorithm,
		Double score,
		Integer level,
		String levelText,
		SIREnum sir,
		Map<MutationSet<VirusT>, Double> partialScores,
		String statement,
		boolean triggered
	) {
		this.drug = drug;
		this.algorithm = algorithm;
		this.score = score == null ? null : Precision.round(score, 2);
		this.level = level;
		this.levelText = levelText;
		this.sir = sir;
		this.partialScores = Collections.unmodifiableMap(partialScores);
		this.statement = statement;
		this.triggered = triggered;
	}
	
	public boolean isTriggered() {
		return triggered;
	}
	
	public Drug<VirusT> getDrug() {
		return drug;
	}
	
	public String getAlgorithm() {
		return algorithm.getName();
	}
	
	public DrugResistanceAlgorithm<VirusT> getAlgorithmObj() {
		return algorithm;
	}

	public Boolean drugIs(String name) {
		return getDrug().getName().equals(name);
	}
	
	public Boolean drugIs(Drug<VirusT> drug) {
		return getDrug() == drug;
	}
	
	public DrugClass<VirusT> getDrugClass() {
		return drug.getDrugClass();
	}
	
	public Boolean drugClassIs(String name) {
		return getDrugClass().getName().equals(name);
	}
	
	public Boolean drugClassIs(DrugClass<VirusT> drugClass) {
		return getDrugClass() == drugClass;
	}
	
	public Double getScore() {
		return score;
	}

	public Integer getLevel() {
		return level;
	}

	public String getLevelText() {
		return levelText;
	}
	
	public SIREnum getSIR() {
		return sir;
	}

	public String getInterpretation() {
		return levelText;
	}
	
	public String getExplanation() {
		String explanation = "";

		if (triggered) {
			String textPartialScores = (
				getParialScorePairs()
				.stream()
				.map(pair -> String.format(
					"%s (%.1f)",
					pair.getLeft().join(" + "),
					pair.getRight()
				))
				.collect(Collectors.joining(", "))
			);
			if (textPartialScores.isEmpty()) {
				// rule was triggered but not partial score description
				// describe the statement and result instead
				explanation = String.format(
					"%s (%s)",
					statement.replace("+", " + ").replace(",", ", "),
					levelText
				);
			}
			else {
				explanation = String.format(
					"Total score: %.1f\n%s",
					score,
					textPartialScores
				);
			}
		}
		else {
			explanation = "No rules were triggered";
		}
		return explanation;
	}

	public Map<MutationSet<VirusT>, Double> getPartialScores() {
		return partialScores;
	}
	
	public Double getPartialScore(MutationSet<VirusT> muts) {
		return partialScores.getOrDefault(muts, 0.);
	}
	
	public Set<Pair<MutationSet<VirusT>, Double>> getParialScorePairs() {
		return (
			partialScores
			.entrySet()
			.stream()
			.map(e -> Pair.of(e.getKey(), e.getValue()))
			.collect(Collectors.toCollection(TreeSet::new))
		);
	}
	
	public String getStatement() {
		return statement;
	}

	public Boolean hasSingleMutPartialScore() {
		return (
			partialScores
			.keySet()
			.stream()
			.anyMatch(muts -> muts.size() == 1)
		);
	}
	
	public Map<Mutation<VirusT>, Double> getSingleMutPartialScores() {
		return (
			partialScores
			.entrySet()
			.stream()
			.filter(e -> e.getKey().size() == 1)
			.collect(Collectors.toMap(
				e -> e.getKey().first(),
				e -> e.getValue(),
				(a, b) -> a,
				TreeMap::new
			))
		);
	}
	
	public Boolean hasMultiMutsPartialScore() {
		return (
			partialScores
			.keySet()
			.stream()
			.anyMatch(muts -> muts.size() > 1)
		);
	}
	
	public Map<MutationSet<VirusT>, Double> getMultiMutsPartialScores() {
		return (
			partialScores
			.entrySet()
			.stream()
			.filter(e -> e.getKey().size() > 1)
			.collect(Collectors.toMap(
				e -> e.getKey(),
				e -> e.getValue(),
				(a, b) -> a,
				TreeMap::new
			))
		);
	}


	@Override
	public String toString() {
		return String.format("%s (%s): %s", drug, algorithm.getName(), sir);
	}
	
	@Override
	public int compareTo(ASIDrugSusc<VirusT> other) {
		return new CompareToBuilder()
			.append(algorithm, other.algorithm)
			.append(drug, other.drug)
			.append(sir, other.sir)
			.append(-score, -other.score) // sort by desc(score)
			.append(-level, -other.level) // probably not necessary
			.append(levelText, other.levelText)
			.append(statement, other.statement)
			.toComparison();
	}

	@Override
	public boolean equals(final Object obj) {
		if (obj == this) { return true; }
		if (obj == null) { return false; }
		if (!(obj instanceof ASIDrugSusc<?>)) { return false; }
		final ASIDrugSusc<?> other = (ASIDrugSusc<?>) obj;
		
		return new EqualsBuilder()
			.append(algorithm, other.algorithm)
			.append(drug, other.drug)
			.append(sir, other.sir)
			.append(score, other.score)
			.append(level, other.level)
			.append(levelText, other.levelText)
			.append(statement, other.statement)
			.isEquals();
    }

	@Override
	public int hashCode() {
		return new HashCodeBuilder(347973, 875361)
			.append(algorithm)
			.append(drug)
			.append(sir)
			.append(score)
			.append(level)
			.append(levelText)
			.append(statement)
			.toHashCode();
	}

}
