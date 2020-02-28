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
import edu.stanford.hivdb.utilities.Json;

public class Strain<VirusT extends Virus<VirusT>> implements Comparable<Strain<VirusT>> {
	
	private final String name;
	private final Integer ordinal;
	private final String displayText;
	private final String nucaminoProfile;
	private final String nucaminoGene;
	private final Integer nucaminoGeneOffset;
	private final Integer absoluteFirstNA;

	private transient VirusT virusInstance;
	private transient Set<DrugClass<VirusT>> drugClasses;

	public static <VirusT extends Virus<VirusT>> Map<String, Strain<VirusT>> loadJson(String raw, VirusT virusIns) {
		Map<String, Strain<VirusT>> strains = new LinkedHashMap<>();
		List<Strain<VirusT>> strainList = Json.loads(
			raw, new TypeToken<List<Strain<VirusT>>>(){});
		for (Strain<VirusT> strain : strainList) {
			strain.virusInstance = virusIns;
			strains.put(strain.getName(), strain);
		}
		return Collections.unmodifiableMap(strains);
	}
	
	private Strain(
		String name, Integer ordinal, String displayText,
		String nucaminoProfile, String nucaminoGene,
		Integer nucaminoGeneOffset, Integer offsetNA
	) {
		this.name = name;
		this.ordinal = ordinal;
		this.displayText = displayText;
		this.nucaminoProfile = nucaminoProfile;
		this.nucaminoGene = nucaminoGene;
		this.absoluteFirstNA = offsetNA;
		this.nucaminoGeneOffset = nucaminoGeneOffset;
	}
	
	public VirusT getVirusInstance() {
		return virusInstance;
	}
	
	public String name() {
		return name;
	}
	
	public String getName() {
		return name;
	}

	public String getDisplayText() {
		return displayText;
	}

	public String getNucaminoProfile() {
		return nucaminoProfile;
	}
	
	public String getNucaminoGene() {
		return nucaminoGene;
	}
	
	public Integer getNucaminoGeneOffset() {
		return nucaminoGeneOffset;
	}
	
	public Integer getAbsoluteFirstNA() {
		return absoluteFirstNA;
	}
	
	public Gene<VirusT> getGene(String name) {
		return virusInstance.getGene(this.name + name);
	}

	public Collection<Gene<VirusT>> getGenes() {
		return virusInstance.getGenes(this);
	}
	
	public Set<DrugClass<VirusT>> getDrugClasses() {
		if (drugClasses == null) {
			Set<DrugClass<VirusT>> drugClasses = (
				virusInstance
				.getDrugClasses()
				.stream()
				.filter(dc -> dc.getStrains().contains(this))
				.collect(Collectors.toCollection(LinkedHashSet::new))
			);
			this.drugClasses = Collections.unmodifiableSet(drugClasses);
		}
		return drugClasses;
	}
	
	@Override
	public String toString() {
		return name();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// strains are singletons
		return false;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(643277, 541754393)
			.append(name)
			.toHashCode();
	}

	@Override
	public int compareTo(Strain<VirusT> o) {
		return ordinal.compareTo(o.ordinal);
	}
	

}
