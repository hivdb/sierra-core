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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.fstrf.stanfordAsiInterpreter.resistance.definition.CommentDefinition;
import org.fstrf.stanfordAsiInterpreter.resistance.definition.Definition;
import org.fstrf.stanfordAsiInterpreter.resistance.definition.LevelDefinition;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedCondition;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedDrug;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedDrugClass;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedDrugLevelCondition;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedGene;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedResultCommentRule;
import org.fstrf.stanfordAsiInterpreter.resistance.grammar.AsiGrammarAdapter.ScoredItem;
import org.fstrf.stanfordAsiInterpreter.resistance.grammar.MutationComparator;
import org.fstrf.stanfordAsiInterpreter.resistance.grammar.StringMutationComparator;

import edu.stanford.hivdb.comments.BoundComment;
import edu.stanford.hivdb.comments.CommentType;
import edu.stanford.hivdb.comments.ConditionType;
import edu.stanford.hivdb.comments.ConditionalComment;
import edu.stanford.hivdb.comments.ConditionalComments;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.StrainModifier;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

public class ASIResultHandler {
	
	/**
	 * Convert an ASI EvaluatedDrug to Drug<VirusT>
	 * 
	 * @param <T>
	 * @param virusIns
	 * @param evalDrug
	 * @return
	 */
	public final static <T extends Virus<T>> Drug<T> convertDrug(T virusIns, EvaluatedDrug evalDrug) {
		String drugName = evalDrug.getDrug().toString();
		return virusIns.getDrug(drugName);
	}
	
	public final static <T extends Virus<T>> ASIDrugSusc<T> extractDrugSusc(
		T virusIns,
		Gene<T> targetGene,
		EvaluatedDrug evalDrug,
		DrugResistanceAlgorithm<T> algorithm
	) {
		Gene<T> srcGene = targetGene.getMainStrainGene();
		StrainModifier strainModifier = srcGene.getTargetStrainModifier(targetGene);
		
		String drugName = evalDrug.getDrug().toString();
		Drug<T> drug = virusIns.getDrug(drugName);
		Map<String, MutationSet<T>> gpMutations = new TreeMap<>();
		Map<String, Double> gpPartialScores = new HashMap<>();

		if (drug == null) {
			// skip unknown drug
			return null;
		}

		LevelDefinition levelDef = evalDrug.getHighestLevelDefinition();
		Double highestTotalScore = Double.NEGATIVE_INFINITY;
		String ruleStatement = "";
		boolean isTriggered = false;

		for(Object condObj : evalDrug.getEvaluatedConditions()) {
			EvaluatedCondition evalCond = (EvaluatedCondition) condObj;
			Object evaluatedResult = evalCond.getEvaluator().getResult();
			String tmpRuleStatement = (
				evalCond
				.getRuleCondition().toString()
				.replaceAll("\\s+", " ")
			);
			if (evaluatedResult instanceof Double) {
				// score rules available
				Double totalScore = (Double) evaluatedResult;

				if (totalScore <= highestTotalScore) {
					continue;
				}
				Collection<?> scoredItemObjs = evalCond.getEvaluator().getScoredItems();
				if (scoredItemObjs.size() > 0) {
					isTriggered = true;
					highestTotalScore = totalScore;
					ruleStatement = tmpRuleStatement;
				}

				for (Object scoredItemObj : scoredItemObjs) {
					ScoredItem scoredItem = (ScoredItem) scoredItemObj;

					Set<?> muts = scoredItem.getMutations();
					MutationSet<T> mutations = MutationSet.parseString(
						srcGene, muts.stream()
						.map(m -> (String) m)
						.collect(Collectors.toSet()))
						.displayAmbiguities();
					// use modifyMutationSet to convert all ASI-compat mutations back
					mutations = strainModifier.modifyMutationSet(srcGene, targetGene, mutations);

					// aggregate scores by positions instead of mutations
					String gpKey = StringUtils.join(mutations.getPositions(), "+");
					MutationSet<T> prevMuts = gpMutations.getOrDefault(gpKey, new MutationSet<>());
					mutations = mutations.mergesWith(prevMuts);
     
					gpMutations.put(gpKey, mutations);
					gpPartialScores.put(
						gpKey,
						gpPartialScores.getOrDefault(gpKey, .0) + scoredItem.getScore()
					);
				}
				
			}
			else if (evaluatedResult instanceof Boolean) {
				// level rules available
				Boolean tmpTriggered = (Boolean) evaluatedResult;
				if (highestTotalScore > Double.NEGATIVE_INFINITY) {
					/**
					 * don't save level rules for a drug already with scores
					 *
					 * The rationale:
					 * For Rega HIV1, the NRTIs have rules. The NNRTIs, PIs, and INSTIs have scores.
					 * However, the PIs also have an extra rule "SELECT ATLEAST 0 FROM (1P)" leading
					 * to level 2.
					 */
					continue;
				}
				if (!tmpTriggered) {
					continue;
				}
				isTriggered = true;
				ruleStatement = tmpRuleStatement;
			}
		}

		Map<MutationSet<T>, Double> partialScores = new TreeMap<>();
		for (String gpKey : gpMutations.keySet()) {
			partialScores.put(gpMutations.get(gpKey), gpPartialScores.get(gpKey));
		}

		return new ASIDrugSusc<T>(
			drug,
			algorithm,
			highestTotalScore == Double.NEGATIVE_INFINITY ? 0 : highestTotalScore,
			levelDef == null ? 1 : levelDef.getOrder(),
			levelDef == null ? algorithm.getOriginalLevelText() : levelDef.getText(),
			levelDef == null ? algorithm.getOriginalLevelSIR() : SIREnum.valueOf(levelDef.getSir()),
			partialScores,
			ruleStatement,
			isTriggered
		);
	}

	public final static <T extends Virus<T>> EvaluatedGene evalutateGeneMutations(
		Gene<T> srcGene,
		MutationSet<T> mutations,
		DrugResistanceAlgorithm<T> algorithm
	) {
		StrainModifier strainModifier = srcGene.getMainStrainModifier();
		Gene<T> targetGene = srcGene.getMainStrainGene();

		// use modifyMutationSet to convert all mutations to ASI-compat mutations
		List<String> asiMutations = strainModifier.modifyMutationSet(
			srcGene, targetGene, mutations
		).toASIFormat();

		org.fstrf.stanfordAsiInterpreter.resistance.definition.Gene asiGene = algorithm.getASIGene(srcGene);

		MutationComparator mutationComparator = new StringMutationComparator(false);
		if (!mutationComparator.areMutationsValid(asiMutations)){
			throw new RuntimeException(
				String.format("Invalid list of mutations: %s",
				asiMutations.toString()
			));
		}

		try {
			return asiGene.evaluate(asiMutations, mutationComparator);
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	public final static <T extends Virus<T>> SortedSet<ASIDrugSusc<T>> extractDrugSuscs(
		Gene<T> targetGene, EvaluatedGene evaluatedGene, DrugResistanceAlgorithm<T> algorithm
	) {
		T virusIns = targetGene.getVirusInstance();
		SortedSet<ASIDrugSusc<T>> asiDrugSuscs = new TreeSet<>();
		for(Object drugClassObj : evaluatedGene.getEvaluatedDrugClasses()) {
			EvaluatedDrugClass evalDrugClass = (EvaluatedDrugClass) drugClassObj;

			for (Object drugObj : evalDrugClass.getEvaluatedDrugs()) {
				EvaluatedDrug evalDrug = (EvaluatedDrug) drugObj;
				ASIDrugSusc<T> drugSusc = extractDrugSusc(
					virusIns, targetGene, evalDrug, algorithm
				);
				if (drugSusc != null) {
					asiDrugSuscs.add(drugSusc);
				}
			}
			
		}
		return asiDrugSuscs;
	}

	public static <T extends Virus<T>> List<BoundComment<T>> extractMutationComments(
		T virusIns, Collection<?> defs, MutationSet<T> muts
	) {
		ConditionalComments<T> condComments = virusIns.getConditionalComments();
		
		List<BoundComment<T>> results = new ArrayList<>();
		for (Object def : defs) {
			CommentDefinition cmtDef = (CommentDefinition) def;
			String commentName = cmtDef.getId();
			ConditionalComment<T> condComment = condComments.get(commentName);
			Mutation<T> mut;
			if (condComment.getConditionType() != ConditionType.MUTATION) {
				throw new RuntimeException(
					String.format("Invalid comment name: %s", commentName)
				);
			}
			
			Mutation<T> matchedMut = muts.get(condComment.getMutationGenePosition());
			if (matchedMut != null) {
				mut = matchedMut.intersectsWith(condComment.getMutationAAs().chars()
						.mapToObj(e -> (char) e)
						.collect(Collectors.toList())
						);
				if (mut == null) {
					throw new IllegalArgumentException(String.format(
							"Mutation %s is not match with comment definition %s.",
							matchedMut.getASIFormat(), commentName));
				}
				List<String> highlight = new ArrayList<>();
				highlight.add(mut.getHumanFormat());
				results.add(new BoundComment<T>(
					condComment.getStrain(), condComment.getName(),
					condComment.getDrugClass(), CommentType.fromMutType(mut.getPrimaryType()),
					cmtDef.getText().replaceAll(
						ConditionalComments.getCommentWildcardRegex(),
						mut.getHumanFormat()
					),
					highlight,
					mut
				));
			}

		}
		results.sort((a, b) -> a.getBoundMutation().compareTo(b.getBoundMutation())); 
		return results;
	}
	
	public static <T extends Virus<T>> List<BoundComment<T>> extractDrugLevelComments(
		T virusIns, Collection<EvaluatedResultCommentRule> resultComments
	) {
		ConditionalComments<T> condComments = virusIns.getConditionalComments();
		List<BoundComment<T>> results = new ArrayList<>();
		for (EvaluatedResultCommentRule resultComment : resultComments) {
			if (!resultComment.getResult()) {
				continue;
			}
			List<Drug<T>> drugs = new ArrayList<>();
			for (
				EvaluatedDrugLevelCondition cond :
				resultComment.getEvaluatedDrugLevelConditions()
			) {
				String drugName = cond.getDrug();
				drugs.add(virusIns.getDrug(drugName));
			}
			for (Definition def : resultComment.getDefinitions()) {
				CommentDefinition cmtDef = (CommentDefinition) def;
				String commentName = cmtDef.getId();
				ConditionalComment<T> condComment = condComments.get(commentName);
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