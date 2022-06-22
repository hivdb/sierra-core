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

package edu.stanford.hivdb.drugs;

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class DrugClass<VirusT extends Virus<VirusT>> implements Comparable<DrugClass<VirusT>> {

	private final String name;
	private final Integer ordinal;
	private final String fullName;
	private final String abstractGene;
	private final List<String> strains;
	private final List<String> synonyms;
	private final List<String> mutationTypes;

	private transient VirusT virusInstance;
	private transient Set<Strain<VirusT>> strainObjs;
	private transient Set<MutationType<VirusT>> mutationTypeObjects;
	private transient Collection<Drug<VirusT>> allDrugs;

	public static <VirusT extends Virus<VirusT>> Map<String, DrugClass<VirusT>> loadJson(String raw, VirusT virusIns) {
		Map<String, DrugClass<VirusT>> drugClasses = new LinkedHashMap<>();
		List<DrugClass<VirusT>> dcList = Json.loads(
			raw, new TypeToken<List<DrugClass<VirusT>>>(){});
		for (DrugClass<VirusT> dc : dcList) {
			dc.virusInstance = virusIns;
			drugClasses.put(dc.getName(), dc);
			drugClasses.put(dc.getFullName(), dc);
			for (String synonym : dc.getSynonyms()) {
				drugClasses.put(synonym, dc);
			}
		}
		return Collections.unmodifiableMap(drugClasses);

	}

	private void initStrains() {
		Set<Strain<VirusT>> strainObjs = (
			strains
			.stream()
			.map(strainText -> virusInstance.getStrain(strainText))
			.collect(Collectors.toCollection(LinkedHashSet::new))
		);
		this.strainObjs = Collections.unmodifiableSet(strainObjs);
	}

	private DrugClass(
		String name, Integer ordinal, String fullName,
		String abstractGene, List<String> strains,
		List<String> synonyms, List<String> mutationTypes
	) {
		this.name = name;
		this.ordinal = ordinal;
		this.fullName = fullName;
		this.abstractGene = abstractGene;
		this.strains = strains;
		this.synonyms = synonyms;
		this.mutationTypes = mutationTypes;
	}

	public List<String> getSynonyms() {
		return synonyms;
	}

	public Collection<Drug<VirusT>> getDrugs() {
		if (allDrugs == null) {
			allDrugs = Collections.unmodifiableCollection(virusInstance.getDrugs(this));
		}
		return allDrugs;
	}

	public String getName() {
		return name;
	}

	public String name() {
		return name;
	}

	public String getFullName() {
		return fullName;
	}

	public String getAbstractGene() {
		return abstractGene;
	}

	public Set<Strain<VirusT>> getStrains() {
		if (strainObjs == null) {
			initStrains();
		}
		return strainObjs;
	}

	public boolean supportStrain(Strain<VirusT> strain) {
		if (strainObjs == null) {
			initStrains();
		}
		return strainObjs.contains(strain);
	}

	public Set<MutationType<VirusT>> getMutationTypes() {
		if (mutationTypeObjects == null) {
			Set<MutationType<VirusT>> mutationTypeObjects = (
				mutationTypes.stream()
				.map(mt -> virusInstance.getMutationType(mt))
				.collect(Collectors.toCollection(TreeSet::new))
			);
			this.mutationTypeObjects = Collections.unmodifiableSet(mutationTypeObjects);
		}
		return mutationTypeObjects;
	}

	public MutationSet<VirusT> getDrugResistMutations() {
		return virusInstance.getDrugResistMutations(this);
	}

	public MutationSet<VirusT> getSurveilDrugResistMutations() {
		return virusInstance.getSurveilDrugResistMutations(this);
	}

	public MutationSet<VirusT> getRxSelectedMutations() {
		return virusInstance.getRxSelectedMutations(this);
	}

	public boolean hasDrugResistMutations() {
		MutationSet<VirusT> drms = virusInstance.getDrugResistMutations(this);
		return drms != null && !drms.isEmpty();
	}

	public boolean hasSurveilDrugResistMutations() {
		MutationSet<VirusT> sdrms = virusInstance.getSurveilDrugResistMutations(this);
		return sdrms != null && !sdrms.isEmpty();
	}

	public boolean hasRxSelectedMutations() {
		MutationSet<VirusT> tsms = virusInstance.getRxSelectedMutations(this);
		return tsms != null && !tsms.isEmpty();
	}

	@Override
	public String toString() {
		return name();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// drug classes are singletons
		return false;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(577251, 467488937)
			.append(name)
			.toHashCode();
	}

	@Override
	public int compareTo(DrugClass<VirusT> o) {
		return ordinal.compareTo(o.ordinal);
	}

}
