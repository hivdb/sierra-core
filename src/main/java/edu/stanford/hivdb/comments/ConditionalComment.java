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
package edu.stanford.hivdb.comments;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.utilities.AAUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class ConditionalComment<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	final private Strain<VirusT> strain;
	final private String commentName;
	final private DrugClass<VirusT> drugClass;
	final private ConditionType conditionType;
	final private Map<String, ?> conditionValue;
	final private String comment;
	
	final private transient VirusT virusInstance;

	protected ConditionalComment(
			Strain<VirusT> strain, String commentName,
			DrugClass<VirusT> drugClass, ConditionType conditionType,
			Map<String, ?> conditionValue, String comment) {
		this.strain = strain;
		this.commentName = commentName;
		this.drugClass = drugClass;
		this.conditionType = conditionType;
		this.conditionValue = conditionValue;
		this.comment = comment;
		this.virusInstance = strain.getVirusInstance();
	}
	
	public Gene<VirusT> getMutationGene() {
		if (conditionType != ConditionType.MUTATION) {
			return null;
		}
		return strain.getGene((String) conditionValue.get("gene"));
	}

	public Integer getMutationPosition() {
		if (conditionType != ConditionType.MUTATION) {
			return null;
		}
		return ((Double) conditionValue.get("pos")).intValue();
	}

	public String getMutationAAs() {
		if (conditionType != ConditionType.MUTATION) {
			return null;
		}
		return AAUtils.normalizeAAs((String) conditionValue.get("aas"));
	}
	
	public GenePosition<VirusT> getMutationGenePosition() {
		if (conditionType != ConditionType.MUTATION) {
			return null;
		}
		return new GenePosition<>(getMutationGene(), getMutationPosition());
	}

	public Map<Drug<VirusT>, List<Integer>> getDrugLevels() {
		Map<Drug<VirusT>, List<Integer>> drugLevels = new LinkedHashMap<>();
		if (conditionType != ConditionType.DRUGLEVEL) {
			return drugLevels;
		}
		if (conditionValue.containsKey("and")) {
			for (Object dlevel : ((List<?>) conditionValue.get("and"))) {
				Map<?, ?> dlevelMap = ((Map<?, ?>) dlevel);
				Drug<VirusT> drug = virusInstance.getDrug((String) dlevelMap.get("drug"));
				List<Integer> levels =
					((List<?>) dlevelMap.get("levels"))
					.stream().map(i -> ((Double) i).intValue())
					.collect(Collectors.toList());
				drugLevels.put(drug, levels);
			}
		}
		else {
			Drug<VirusT> drug = virusInstance.getDrug((String) conditionValue.get("drug"));
			List<Integer> levels =
				((List<?>) conditionValue.get("levels"))
				.stream().map(i -> ((Double) i).intValue())
				.collect(Collectors.toList());
			drugLevels.put(drug, levels);
		}
		return drugLevels;
	}

	public String getDrugLevelsText() {
		StringBuilder text = new StringBuilder();
		Map<Drug<VirusT>, List<Integer>> drugLevels = getDrugLevels();
		String delimiter = "";
		for (Map.Entry<Drug<VirusT>, List<Integer>> e : drugLevels.entrySet()) {
			Drug<VirusT> drug = e.getKey();
			List<Integer> levels = e.getValue();
			text.append(String.format(
				"%s%s: %s", delimiter, drug,
				levels.stream().map(l -> l.toString())
				.collect(Collectors.joining(", "))
			));
			delimiter = "; ";
		}
		return text.toString();
	}

	public String getName() { return commentName; }
	public String getText() { return comment; }
	public DrugClass<VirusT> getDrugClass() { return drugClass; }
	public ConditionType getConditionType() { return conditionType; }
	
	@Override
	public Gene<VirusT> getGene() {
		return strain.getGene(drugClass.getAbstractGene());
	}
}