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

package edu.stanford.hivdb.comments;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.fstrf.stanfordAsiInterpreter.resistance.definition.CommentDefinition;
import org.fstrf.stanfordAsiInterpreter.resistance.definition.Definition;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedDrugLevelCondition;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedResultCommentRule;

import com.google.common.primitives.Chars;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.utilities.Json;

public class ConditionalComments<VirusT extends Virus<VirusT>> {

	private static final String WILDCARD_REGEX = "\\$listMutsIn\\{.+?\\}";

	private final Map<String, ConditionalComment<VirusT>> condCommentMap;
	private final VirusT virusInstance;
	private final transient Map<GenePosition<VirusT>, List<ConditionalComment<VirusT>>> condCommentsByGenePos;
	
	@SuppressWarnings("unchecked")
	public ConditionalComments(String jsonRaw, VirusT virusIns) {
		List<Map<String, Object>> condCommentMaps = Json.loads(jsonRaw, new TypeToken<List<Map<String, Object>>>() {});
		Map<String, ConditionalComment<VirusT>> condComments = new LinkedHashMap<>();
		for (Map<String, Object> condCmtMap: condCommentMaps) {
			String commentName = (String) condCmtMap.get("commentName");
			condComments.put(commentName, new ConditionalComment<VirusT>(
				virusIns.getStrain((String) condCmtMap.get("strain")),
				commentName,
				virusIns.getDrugClass((String) condCmtMap.get("drugClass")),
				ConditionType.valueOf((String) condCmtMap.get("conditionType")),
				(Map<String, ?>) condCmtMap.get("conditionValue"),
				(String) condCmtMap.get("comment")
			));
		}
		virusInstance = virusIns;
		condCommentMap = Collections.unmodifiableMap(condComments);
		condCommentsByGenePos = Collections.unmodifiableMap(
			condComments.values()
			.stream()
			.filter(cmt -> cmt.getConditionType() == ConditionType.MUTATION)
			.collect(Collectors.groupingBy(ConditionalComment::getMutationGenePosition))
		);
	}
	
	private static <VirusT extends Virus<VirusT>> BoundComment<VirusT> findMutationComment(
			Mutation<VirusT> mutation, ConditionalComment<VirusT> cc) {

		List<Character> aaChars = Chars.asList(
			cc.getMutationAAs().toCharArray()
		);
		Mutation<VirusT> resultMut = mutation.intersectsWith(aaChars);
		if (resultMut == null) {
			return null;
		}
		List<String> highlight = new ArrayList<>();
		highlight.add(resultMut.getHumanFormat());

		return new BoundComment<>(
			cc.getStrain(), cc.getName(), cc.getDrugClass(),
			CommentType.fromMutType(resultMut.getPrimaryType()),
			cc.getText().replaceAll(WILDCARD_REGEX, resultMut.getHumanFormat()),
			highlight,
			resultMut
		);
	}
	
	public List<BoundComment<VirusT>> getComments(Mutation<VirusT> mutation) {
		GenePosition<VirusT> genePos = mutation.getGenePosition();
		List<BoundComment<VirusT>> comments = new ArrayList<>();
		List<ConditionalComment<VirusT>> candidates = condCommentsByGenePos.get(genePos);
		if (candidates == null) {
			return comments;
		}
		for (ConditionalComment<VirusT> cc : candidates) {
			BoundComment<VirusT> comment = findMutationComment(mutation, cc);
			if (comment != null) {
				comments.add(comment);
			}
		}
		return comments;
	}

	public List<BoundComment<VirusT>> fromAsiMutationComments(
		Collection<?> defs, MutationSet<VirusT> muts
	) {
		List<BoundComment<VirusT>> results = new ArrayList<>();
		for (Object def : defs) {
			CommentDefinition cmtDef = (CommentDefinition) def;
			String commentName = cmtDef.getId();
			ConditionalComment<VirusT> condComment = condCommentMap.get(commentName);
			Mutation<VirusT> mut;
			if (condComment.getConditionType() == ConditionType.MUTATION) {
				mut = (
					muts
					.get(condComment.getMutationGene(), condComment.getMutationPosition())
					.intersectsWith(
						condComment.getMutationAAs().chars()
						.mapToObj(e -> (char) e)
						.collect(Collectors.toList())
					)
				);
			}
			else {
				throw new RuntimeException(
					String.format("Invalid comment name: %s", commentName)
				);
			}
			
			List<String> highlight = new ArrayList<>();
			highlight.add(mut.getHumanFormat());
			results.add(new BoundComment<>(
				condComment.getStrain(), condComment.getName(),
				condComment.getDrugClass(), CommentType.fromMutType(mut.getPrimaryType()),
				cmtDef.getText().replaceAll(WILDCARD_REGEX, mut.getHumanFormat()),
				highlight,
				mut
			));
		}
		results.sort((a, b) -> a.getBoundMutation().compareTo(b.getBoundMutation())); 
		return results;
	}
	
	public List<BoundComment<VirusT>> fromAsiDrugLevelComments(
		Collection<EvaluatedResultCommentRule> resultComments
	) {
		List<BoundComment<VirusT>> results = new ArrayList<>();
		for (EvaluatedResultCommentRule resultComment : resultComments) {
			if (!resultComment.getResult()) {
				continue;
			}
			List<Drug<VirusT>> drugs = new ArrayList<>();
			for (
				EvaluatedDrugLevelCondition cond :
				resultComment.getEvaluatedDrugLevelConditions()
			) {
				String drugName = cond.getDrug();
				drugs.add(virusInstance.getDrug(drugName));
			}
			for (Definition def : resultComment.getDefinitions()) {
				CommentDefinition cmtDef = (CommentDefinition) def;
				String commentName = cmtDef.getId();
				ConditionalComment<VirusT> condComment = condCommentMap.get(commentName);
				Set<String> highlights = drugs.stream()
					.map(d -> d.getDisplayAbbr())
					.collect(Collectors.toSet());
				highlights.addAll(
					drugs.stream()
					.map(d -> d.name())
					.collect(Collectors.toSet())
				);
				results.add(new BoundComment<>(
					condComment.getStrain(), condComment.getName(),
					condComment.getDrugClass(), CommentType.Dosage,
					cmtDef.getText(), highlights, null
				));
			}
		}
		return results;
	}
}
