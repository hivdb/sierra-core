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
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.fstrf.stanfordAsiInterpreter.resistance.definition.LevelDefinition;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedCondition;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedDrug;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedDrugClass;
import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedGene;
import org.fstrf.stanfordAsiInterpreter.resistance.grammar.AsiGrammarAdapter.ScoredItem;
import org.fstrf.stanfordAsiInterpreter.resistance.grammar.MutationComparator;
import org.fstrf.stanfordAsiInterpreter.resistance.grammar.StringMutationComparator;

import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.drugs.DrugResistanceAlgorithm;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationMapUtils;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationMapUtils.SortOrder;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class AsiResult<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {

	// Use canonical name of HIVDB Gene class since
	// FSTRFAsi code has its own Gene class we call asiGene
	protected final edu.stanford.hivdb.viruses.Gene<VirusT> gene;
	protected final List<String> mutations;

	private final DrugResistanceAlgorithm<VirusT> algorithm;

	// Objects created by the ASI code.
	protected org.fstrf.stanfordAsiInterpreter.resistance.definition.Gene asiGene;
	protected EvaluatedGene evaluatedGene;

	protected MutationSet<VirusT> triggeredMutations;
	protected Map<DrugClass<VirusT>, MutationSet<VirusT>>
		triggeredMutationsByDrugClass = new TreeMap<>();

	// drugLevel: 1 to 5 corresponding to drugLevel
	// drugLevelText: Susceptible, Potential Low Level Res., Low-Level Res, Intermediate Res, High-level Res.
	// drugLevelSir: S (Susceptible), I (Intermediate), R (Resistance)
	protected Map<Drug<VirusT>, Integer> drugLevel = new TreeMap<>();
	protected Map<Drug<VirusT>, String> drugLevelText = new TreeMap<>();
	protected Map<Drug<VirusT>, String> drugLevelSir = new TreeMap<>();

	protected Map<Drug<VirusT>, Double> totalDrugScores = new TreeMap<>();
	protected Map<Drug<VirusT>, Map<Mutation<VirusT>, Double>> drugMutScores = new TreeMap<>();
	protected Map<Drug<VirusT>, Map<MutationSet<VirusT>, Double>> drugComboMutScores = new TreeMap<>();

	// List of rules that triggered for each drug mapped to the level of the corresponding action
	protected Map<Drug<VirusT>, Map<String, String>> triggeredDrugRules = new TreeMap<>();

	/**
	 * Instantiated by a gene and list of mutations in a sequence.
	 * asiGene is an object created by the XMLAsiTransformer that contains all of the
	 *   specifications and rules in the algorithm for the submitted gene
	 * evaluatedGene is an object that contains all of the results obtained by applying
	 *   the asiGene rules to the submitted mutations
	 *
	 * @param submittedGene
	 * @param mutations
	 * @param algorithm
	 */
	public AsiResult (
		final edu.stanford.hivdb.viruses.Gene<VirusT> submittedGene,
		final MutationSet<VirusT> mutations,
		final DrugResistanceAlgorithm<VirusT> algorithm
	) {
		this.gene = submittedGene;
		this.mutations = mutations.toASIFormat();
		this.algorithm = algorithm;

		asiGene = algorithm.getASIGene(submittedGene);

		MutationComparator mutationComparator = new StringMutationComparator(false);
		if (!mutationComparator.areMutationsValid(this.mutations)){
			throw new RuntimeException("Invalid list of mutations: " + this.mutations.toString());
		}

		try {
			this.evaluatedGene = asiGene.evaluate(this.mutations, mutationComparator);
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}
		populateEvaluatedResults();
	}
	
	@Override
	public edu.stanford.hivdb.viruses.Gene<VirusT> getGene() {
		return this.gene;
	}

	/**
	 * @param drug
	 * @return
	 */
	public final int getDrugLevel(Drug<VirusT> drug) {
		if (!drugLevel.containsKey(drug)) {
			return 1;
		} else {
			return drugLevel.get(drug);
		}
	}

	public String getDrugLevelText(Drug<VirusT> drug) {
		if (!drugLevelText.containsKey(drug)) {
			return algorithm.getOriginalLevelText();
		} else {
			return drugLevelText.get(drug);
		}
	}

	public final String getDrugLevelSir(Drug<VirusT> drug) {
		if (!drugLevelSir.containsKey(drug)) {
			return algorithm.getOriginalLevelSIR();
		} else {
			return drugLevelSir.get(drug);
		}
	}

	public final Double getTotalScore(Drug<VirusT> drug) {
		if (!totalDrugScores.containsKey(drug)) {
			return .0;
		} else {
			return totalDrugScores.get(drug);
		}
	}

	private final <T> Map<DrugClass<VirusT>, Map<Drug<VirusT>, T>> groupingByDrugClass(Map<Drug<VirusT>, T> input) {
		return input
			.entrySet()
			.stream()
			.collect(Collectors.groupingBy(
				e -> e.getKey().getDrugClass(),
				() -> new TreeMap<>(),
				Collectors.toMap(
					e -> e.getKey(),
					e -> e.getValue(),
					(v1, v2) -> v1,
					() -> new TreeMap<>()
				)
			));
	}

	private final <T> Map<Drug<VirusT>, T> filterByDrugClass(Map<Drug<VirusT>, T> input, DrugClass<VirusT> drugClass) {
		return input
			.entrySet()
			.stream()
			.filter(e -> e.getKey().getDrugClass() == drugClass)
			.collect(Collectors.toMap(
				e -> e.getKey(),
				e -> e.getValue()
			));
	}

	/**
	 * Data structure:
	 *   DrugClass => Drug => totalScore (obtained from adding up the individual and combination scores
	 *   for each mutation in a sequence
	 *
	 * @return Map: DrugClass => Drug => totalScore
	 */
	public final Map<DrugClass<VirusT>, Map<Drug<VirusT>, Double>> getDrugClassTotalDrugScores() {
		return groupingByDrugClass(totalDrugScores);
	}

	public final Map<Drug<VirusT>, Double> getDrugClassTotalDrugScores(DrugClass<VirusT> drugClass) {
		return filterByDrugClass(totalDrugScores, drugClass);
	}

	public final MutationSet<VirusT> getTriggeredMutations() {
		if (triggeredMutations == null) {
			MutationSet<VirusT> allMuts = new MutationSet<>();
			for (Map<Mutation<VirusT>, Double> mutScores : drugMutScores.values()) {
				allMuts = allMuts.mergesWith(mutScores.keySet());
			}
			for (Map<MutationSet<VirusT>, Double> comboMutScores : drugComboMutScores.values()) {
				for (MutationSet<VirusT> muts : comboMutScores.keySet()) {
					allMuts = allMuts.mergesWith(muts);
				}
			}
			allMuts.displayAmbiguities();
			triggeredMutations = allMuts;
		}
		return triggeredMutations;
	}

	public final MutationSet<VirusT> getTriggeredMutations(DrugClass<VirusT> drugClass) {
		if (!triggeredMutationsByDrugClass.containsKey(drugClass)) {
			MutationSet<VirusT> allMuts = new MutationSet<>();
			for (Drug<VirusT> drug : drugClass.getDrugs()) {
				allMuts = allMuts.mergesWith(drugMutScores
					.getOrDefault(drug, new HashMap<>()).keySet());
				for (MutationSet<VirusT> muts : drugComboMutScores
					 .getOrDefault(drug, new HashMap<>()).keySet()) {
					allMuts = allMuts.mergesWith(muts);
				}
			}
			allMuts = allMuts.displayAmbiguities();
			triggeredMutationsByDrugClass.put(drugClass, allMuts);
		}
		return triggeredMutationsByDrugClass.get(drugClass);
	}

	/**
	 * This data map only has entries for Drugs with scored Mutations
	 * @return
	 */
	public final Map<DrugClass<VirusT>, Map<Drug<VirusT>, Map<Mutation<VirusT>, Double>>> getDrugClassDrugMutScores() {
		return groupingByDrugClass(drugMutScores);
	}

	/**
	 * This data map only has entries for Drugs with scored Mutation Combinations
	 * @return
	 */
	public final Map<DrugClass<VirusT>, Map<Drug<VirusT>, Map<MutationSet<VirusT>, Double>>> getDrugClassDrugComboMutScores() {
		return groupingByDrugClass(drugComboMutScores);
	}

	public final Map<Drug<VirusT>, Map<String, String>> getTriggeredDrugRules() {
		return triggeredDrugRules;
	}

	public final Map<Drug<VirusT>, Map<Mutation<VirusT>, Double>> getDrugMutScores() {
		return drugMutScores;
	}

	public final Map<Drug<VirusT>, Map<MutationSet<VirusT>, Double>> getDrugComboMutScores() {
		return drugComboMutScores;
	}

	public AsiDrugComparableResult getDrugComparableResult(Drug<VirusT> drug) {
		String sir = getDrugLevelSir(drug);
		String level = getDrugLevelText(drug);
		String explanation = getDrugExplanation(drug);
		return new AsiDrugComparableResult(sir, level, explanation);
	}

	protected String getDrugExplanation(Drug<VirusT> drug) {
		String explanation = "";
		Boolean hasMuts = false;
		String individualMuts = "";
		String comboMuts = "";
		Double totalScore = totalDrugScores.get(drug);

		Map<Mutation<VirusT>, Double> mutScores =
			drugMutScores.getOrDefault(drug, new HashMap<>());

		Map<MutationSet<VirusT>, Double> comboMutScores =
			drugComboMutScores.getOrDefault(drug, new HashMap<>());

		if (!mutScores.isEmpty()) {
			Map<Mutation<VirusT>, Double> mutScoresSortedByScore =
					MutationMapUtils.sortByComparator(mutScores, SortOrder.DESC);
			individualMuts = MutationMapUtils.printMutScoresAsInts(mutScoresSortedByScore);
			hasMuts = true;
		}
		if (!comboMutScores.isEmpty()) {
			Map<MutationSet<VirusT>, Double> comboMutsSortedByScore =
					MutationMapUtils.sortByComparator(comboMutScores, SortOrder.DESC);
			comboMuts = MutationMapUtils.printMutSetScoresAsInts(comboMutsSortedByScore);
			hasMuts = true;
		}

		if (hasMuts) {
			explanation = String.format(
				"Total score: %d\n%s %s",
				totalScore.intValue(), individualMuts, comboMuts);
		}
		else if (triggeredDrugRules.containsKey(drug)) {
			Map<String, String> rules = triggeredDrugRules.get(drug);
			List<String> explanationList = new ArrayList<>();
			for (String rule : rules.keySet()) {
				explanationList.add(String.format("%s (%s)",
					rule.replace("+", " + ").replace(",", ", "), rules.get(rule)));
			}
			explanation = String.join("\n", explanationList);
		} else {
			explanation = "No rules were triggered";
		}
		return explanation;
	}

	private final DrugClass<VirusT> convertDrugClass(EvaluatedDrugClass evalDrugClass) {
		String drugClassName = evalDrugClass.getDrugClass().getClassName();
		return gene.getVirusInstance().getDrugClass(drugClassName);
	}

	private final Drug<VirusT> convertDrug(EvaluatedDrug evalDrug) {
		String drugName = evalDrug.getDrug().toString();
		return gene.getVirusInstance().getDrug(drugName);
	}

	protected void scoreHandler(
			DrugClass<VirusT> drugClass, Drug<VirusT> drug, double score, MutationSet<VirusT> mutations) {
		// Populate map for individual mutation scores
		if (mutations.size() == 1) {
			Mutation<VirusT> scoredMut = mutations.first();
			drugMutScores.putIfAbsent(drug, new LinkedHashMap<Mutation<VirusT>, Double>());
			double drugMutScore = drugMutScores.get(drug).getOrDefault(scoredMut, -99.0);
			drugMutScores.get(drug).put(scoredMut, Math.max(score, drugMutScore));

		// Populate map for combination mutation scores
		} else {
			drugComboMutScores.putIfAbsent(drug, new LinkedHashMap<>());
			double drugComboMutScore = drugComboMutScores.get(drug).getOrDefault(mutations, -99.0);
			drugComboMutScores.get(drug).put(mutations, Math.max(score, drugComboMutScore));
		}
	}

	protected void susceptibilityHandler(
			DrugClass<VirusT> drugClass, Drug<VirusT> drug, String condition, String susceptibility) {
		triggeredDrugRules.putIfAbsent(drug, new LinkedHashMap<String, String>());
		triggeredDrugRules.get(drug).putIfAbsent(condition, susceptibility);
	}
	
	public EvaluatedGene getEvaluatedGene() {
		return evaluatedGene;
	}

	// drugClasstotalDrugScores: DrugClass => Drug => totalScore
	// drugClassDrugMutScores: DrugClass => Drug => Mutation=> score
	// drugClassDrugComboMutScores DrugClass => Drug => Mutation combination (List<Mutation>) => score
	protected void populateEvaluatedResults() {
		for(Object drugClassObj : evaluatedGene.getEvaluatedDrugClasses()) {
			EvaluatedDrugClass evalDrugClass = (EvaluatedDrugClass) drugClassObj;
			DrugClass<VirusT> drugClass = convertDrugClass(evalDrugClass);

			for (Object drugObj : evalDrugClass.getEvaluatedDrugs()) {
				EvaluatedDrug evalDrug = (EvaluatedDrug) drugObj;
				Drug<VirusT> drug = convertDrug(evalDrug);
				if (drug == null) {
					// skip unknown drug
					continue;
				}

				LevelDefinition levelDef = evalDrug.getHighestLevelDefinition();

				int levelDefOrder = levelDef == null ? 1 : levelDef.getOrder();
				String levelDefText = levelDef == null ? algorithm.getOriginalLevelText() : levelDef.getText();
				String levelDefSir = levelDef == null ? algorithm.getOriginalLevelSIR() : levelDef.getSir();

				drugLevel.put(drug, levelDefOrder);
				drugLevelText.put(drug, levelDefText);
				drugLevelSir.put(drug, levelDefSir);

				for(Object condObj : evalDrug.getEvaluatedConditions()) {
					EvaluatedCondition evalCond = (EvaluatedCondition) condObj;
					Object o = evalCond.getEvaluator().getResult();
					if (o instanceof Double) {
						totalDrugScores.put(drug, (Double) o);
						for (Object scoredItemObj : evalCond.getEvaluator().getScoredItems()) {
							ScoredItem scoredItem = (ScoredItem) scoredItemObj;
							Set<?> muts = scoredItem.getMutations();
							MutationSet<VirusT> mutations = MutationSet.parseString(
								gene, muts.stream()
								.map(m -> (String) m)
								.collect(Collectors.toSet()))
								.displayAmbiguities();
							this.scoreHandler(drugClass, drug, scoredItem.getScore(), mutations);
						}
					}
					else if (o instanceof Boolean) {
						String ruleCondition = evalCond
							.getRuleCondition().toString()
							.replaceAll("\\s+", " ");
						Boolean triggered = (Boolean) o;
						if (!triggered) {
							continue;
						}
						for (Object defObj : evalCond.getDefinitions()) {
							String ruleSusceptibilityText = "";
							if (defObj instanceof LevelDefinition) {
								LevelDefinition ruleLevelDefinition = (LevelDefinition) defObj;
								ruleSusceptibilityText = ruleLevelDefinition.getText();
							}
							this.susceptibilityHandler(
								drugClass, drug, ruleCondition, ruleSusceptibilityText);
						}
					}
				}
			}
		}
	}

	public String getAlgorithmName() {
		return algorithm.getName();
	}
}
