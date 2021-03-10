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

package edu.stanford.hivdb.viruses;

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.StrainModifier;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.utilities.AssertUtils;
import edu.stanford.hivdb.utilities.Json;

public class Gene<VirusT extends Virus<VirusT>> implements Comparable<Gene<VirusT>> {

	private final String name;
	private final Integer ordinal;
	private final String strain;
	private final String abstractGene;
	private final String refSequence;
	private final Map<String, String> strainModifiers;
	private final List<String> mutationTypes;
	private final Integer nucaminoMinNumOfAA;

	private transient VirusT virusInstance;
	private transient Map<String, StrainModifier> strainModifierMap;
	private transient Set<MutationType<VirusT>> mutationTypeObjs;
	private transient Set<DrugClass<VirusT>> drugClasses;

	public static <VirusT extends Virus<VirusT>> Map<String, Gene<VirusT>> loadJson(String raw, VirusT virusIns) {
		Map<String, Gene<VirusT>> genes = new LinkedHashMap<>();
		List<Gene<VirusT>> geneList = Json.loads(
			raw, new TypeToken<List<Gene<VirusT>>>(){});
		for (Gene<VirusT> gene : geneList) {
			gene.virusInstance = virusIns;
			genes.put(gene.getName(), gene);
			genes.put(gene.getName().toUpperCase(), gene);
		}
		return Collections.unmodifiableMap(genes);
	}

	public VirusT getVirusInstance() {
		return virusInstance;
	}

	private Gene(
		String name, Integer ordinal, String strain, String abstractGene,
		String refSequence, List<String> mutationTypes,
		Map<String, String> strainModifiers,
		Integer nucaminoMinNumOfAA) {
		this.name = name;
		this.ordinal = ordinal;
		this.strain = strain;
		this.abstractGene = abstractGene;
		this.refSequence = refSequence;
		this.mutationTypes = mutationTypes;
		this.strainModifiers = strainModifiers;
		this.nucaminoMinNumOfAA = nucaminoMinNumOfAA;
	}
	
	public String name() {
		return name;
	}
	
	public String getName() {
		return name;
	}
	
	/**
	 * Get corresponding gene in main strain.
	 * 
	 * For example, for HIV2 virus, the main strain is HIV2A.
	 * Therefore for HIV2BRT, the main strain is HIV2ART. 
	 * 
	 * @return gene
	 */
	public Gene<VirusT> getMainStrainGene() {
		Strain<VirusT> mainStrain = virusInstance.getMainStrain();
		return getTargetStrainGene(mainStrain);
	}
	
	public Gene<VirusT> getTargetStrainGene(String targetStrainName) {
		Strain<VirusT> targetStrain = virusInstance.getStrain(targetStrainName);
		return getTargetStrainGene(targetStrain);
	}

	public Gene<VirusT> getTargetStrainGene(Strain<VirusT> targetStrain) {
		if (getStrain() == targetStrain) {
			return this;
		}
		else {
			return targetStrain.getGene(getAbstractGene());
		}
	}
	
	@Deprecated
	public int getLength() {
		// deprecated, use getAASize instead
		return refSequence.length();
	}
	
	public int getAASize() {
		return refSequence.length();
	}
	
	public int getNASize() {
		return refSequence.length() * 3;
	}
	
	public int getNucaminoMinNumOfAA() {
		return nucaminoMinNumOfAA;
	}

	public Character getRefChar(int pos) {
		try {
			return refSequence.charAt(pos - 1);
		}
		catch (StringIndexOutOfBoundsException exc) {
			throw new RuntimeException(String.format(
				"Position out of %s range: %d",
				name, pos
			));
		}
	}

	public String getRefSequence(int pos, int length) {
		return refSequence.substring(pos - 1, pos - 1 + length);
	}

	public String getRefSequence() {
		return refSequence;
	}

	public Strain<VirusT> getStrain() {
		return virusInstance.getStrain(strain);
	}

	public String getAbstractGene() {
		return abstractGene;
	}
	
	/**
	 * Get StrainModifier where targetStrain=mainStrain
	 */
	public StrainModifier getMainStrainModifier() {
		return getTargetStrainModifier(virusInstance.getMainStrain());
	}

	public StrainModifier getTargetStrainModifier(Gene<?> targetGene) {
		return getTargetStrainModifier(targetGene.getStrain().getName());
	}
	
	public StrainModifier getTargetStrainModifier(Strain<?> targetStrain) {
		return getTargetStrainModifier(targetStrain.getName());
	}
	
	public StrainModifier getTargetStrainModifier(String targetStrainText) {
		if (strainModifierMap == null) {
			Map<String, StrainModifier> strainModifierMap = (
				strainModifiers
				.entrySet()
				.stream()
				.collect(Collectors.toMap(
					e -> e.getKey(),
					e -> new StrainModifier(e.getKey(), e.getValue())
				))
			);
			strainModifierMap.put(
				this.strain,
				new StrainModifier(
					this.strain,
					String.format("%dM", this.getAASize())
				)
			);
			this.strainModifierMap = Collections.unmodifiableMap(strainModifierMap);
		}
		return AssertUtils.notNull(
			strainModifierMap.get(targetStrainText),
			"Strain modifier for target strain %s is not defined", targetStrainText
		);
	}

	public Set<DrugClass<VirusT>> getDrugClasses() {
		if (drugClasses == null) {
			Strain<VirusT> strain = getStrain();
			Set<DrugClass<VirusT>> drugClasses = (
				virusInstance
				.getDrugClasses()
				.stream()
				.filter(dc -> (
					dc.getStrains().contains(strain) &&
					dc.getAbstractGene().equals(abstractGene)
				))
				.collect(Collectors.toCollection(LinkedHashSet::new))
			);
			this.drugClasses = Collections.unmodifiableSet(drugClasses);
		}
		return drugClasses;
	}
	
	public Set<MutationType<VirusT>> getMutationTypes() {
		if (mutationTypeObjs == null) {
			Set<MutationType<VirusT>> mutationTypeObjs = (
				mutationTypes.stream()
				.map(mtype -> virusInstance.getMutationType(mtype))
				.collect(Collectors.toCollection(LinkedHashSet::new))
			);
			this.mutationTypeObjs = Collections.unmodifiableSet(mutationTypeObjs);
		}
		return mutationTypeObjs;
	}

	public Collection<GenePosition<VirusT>> getGenePositionsBetween(int startPos, int endPos) {
		Set<GenePosition<VirusT>> genePositions = new LinkedHashSet<>();
		startPos = Math.max(startPos, 1);
		endPos = Math.min(endPos, getLength());
		for (int pos = startPos; pos <= endPos; pos ++) {
			genePositions.add(new GenePosition<>(this, pos));
		}
		return genePositions;
	}

	@Override
	public String toString() {
		return name();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// genes are singletons
		return false;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(46548647, 769543361)
			.append(name)
			.toHashCode();
	}

	@Override
	public int compareTo(Gene<VirusT> o) {
		return ordinal.compareTo(o.ordinal);
	}

}
