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

package edu.stanford.hivdb.comments;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedResultCommentRule;

import com.google.common.primitives.Chars;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugresistance.algorithm.ASIResultHandler;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.utilities.Json;

public class ConditionalComments<VirusT extends Virus<VirusT>> {

	private final Map<String, ConditionalComment<VirusT>> condCommentMap;
	private final VirusT virusInstance;
	private final transient Map<GenePosition<VirusT>, List<ConditionalComment<VirusT>>> condCommentsByGenePos;

	private static final String WILDCARD_REGEX = "\\$listMutsIn\\{.+?\\}";
	
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
		Map<GenePosition<VirusT>, List<ConditionalComment<VirusT>>> selfCondCommentsByGenePos = new HashMap<>();
		for (ConditionalComment<VirusT> cmt : condComments.values()) {
			if (cmt.getConditionType() != ConditionType.MUTATION) {
				continue;
			}
			for (GenePosition<VirusT> genePos : cmt.getMutationLookup().keySet()) {
				selfCondCommentsByGenePos.putIfAbsent(genePos, new ArrayList<>());
				selfCondCommentsByGenePos.get(genePos).add(cmt);
			}
		}
		condCommentsByGenePos = Collections.unmodifiableMap(selfCondCommentsByGenePos);
	}
	
	public static String getCommentWildcardRegex() {
		return WILDCARD_REGEX;
	}
	
	public ConditionalComment<VirusT> get(String name) { return condCommentMap.get(name); }
	
	private static <VirusT extends Virus<VirusT>> BoundComment<VirusT> findMutationComment(
			Mutation<VirusT> mutation, ConditionalComment<VirusT> cc) {

		List<Character> aaChars = Chars.asList(
			cc.getMutationAAs(mutation.getGenePosition()).toCharArray()
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

	@Deprecated
	public List<BoundComment<VirusT>> fromAsiMutationComments(
		Collection<?> defs, MutationSet<VirusT> muts
	) {
		return ASIResultHandler.extractMutationComments(virusInstance, defs, muts);
	}
	
	@Deprecated
	public List<BoundComment<VirusT>> fromAsiDrugLevelComments(
		Collection<EvaluatedResultCommentRule> resultComments
	) {
		return ASIResultHandler.extractDrugLevelComments(virusInstance, resultComments);
	}
}
