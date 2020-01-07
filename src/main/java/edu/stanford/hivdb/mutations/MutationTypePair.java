/*

    Copyright (C) 2019 Stanford HIVDB team

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

package edu.stanford.hivdb.mutations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.utilities.AAUtils;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class MutationTypePair<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	private final String strain;
	private final String gene;
	private final String drugClass;
	private final Integer position;
	private final String aas;
	private final String mutationType;
	private final Boolean isUnusual;
	private transient Mutation<VirusT> mutObj;

	private transient String name;
	private transient VirusT virusInstance;
	
	public static <VirusT extends Virus<VirusT>> List<MutationTypePair<VirusT>> loadJson(String raw, VirusT virusIns) {
		List<MutationTypePair<VirusT>> mutTypePairs = new ArrayList<>();
		List<MutationTypePair<VirusT>> mutTypePairList = Json.loads(
			raw, new TypeToken<List<MutationTypePair<VirusT>>>(){});
		for (MutationTypePair<VirusT> mutTypePair : mutTypePairList) {
			mutTypePair.virusInstance = virusIns;
			mutTypePairs.add(mutTypePair);
		}
		return Collections.unmodifiableList(mutTypePairs);
	}
	
	private MutationTypePair(
		final String strain,
		final String gene,
		final String drugClass,
		final int position, final String aas,
		final String mutationType,
		final Boolean isUnusual
	) {
		this.strain = strain;
		this.gene = gene;
		this.drugClass = drugClass;
		this.position = position;
		this.aas = aas;
		this.mutationType = mutationType;
		this.isUnusual = isUnusual;
	}
	
	@Override
	public Strain<VirusT> getStrain() {
		return virusInstance.getStrain(strain);
	}

	@Override
	public Gene<VirusT> getGene() {
		return virusInstance.getGene(strain + gene);
	}
	
	@Override
	public String getAbstractGene() {
		return gene;
	}
	
	public GenePosition<VirusT> getGenePosition() {
		return new GenePosition<>(getGene(), position);
	}
	
	public DrugClass<VirusT> getDrugClass() {
		return virusInstance.getDrugClass(drugClass);
	}
	public Integer getPosition() { return position; }
	public String getTriggeredAAs() { return aas; }
	public MutationType<VirusT> getMutationType() {
		return virusInstance.getMutationType(mutationType);
	}
	public boolean isUnusual() { return isUnusual; }
	public Mutation<VirusT> getMutObj() {
		if (mutObj == null) {
			mutObj = new AAMutation<>(getGene(), position, aas.toCharArray(), 0xff);
		}
		return mutObj;
	}

	public boolean isMutationMatched(Mutation<VirusT> targetMut) {
		return getMutObj().containsSharedAA(targetMut);
	}
	
	public String getName() {
		if (name == null) {
			StringBuilder typeStr = new StringBuilder();
			if (drugClass.equals(mutationType)) {
				typeStr.append(mutationType);
			}
			else {
				typeStr.append(drugClass);
				typeStr.append(mutationType);
			}
			name = String.format(
				"%s%s_POS%d%s_%s", strain, gene, position,
				AAUtils.toASIFormat(aas), typeStr.toString()
			);
		}
		return name;
	}
	
	public String name() {
		return getName();
	}
	
	@Override
	public String toString() {
		return getName();
	}

}