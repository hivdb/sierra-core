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

package edu.stanford.hivdb.drugresistance;

import java.util.Collection;
import java.util.Collections;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.stream.Collectors;

import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedGene;

import edu.stanford.hivdb.comments.BoundComment;
import edu.stanford.hivdb.comments.CommentType;
import edu.stanford.hivdb.drugresistance.algorithm.ASIDrugSusc;
import edu.stanford.hivdb.drugresistance.algorithm.ASIResultHandler;
import edu.stanford.hivdb.drugresistance.algorithm.DrugResistanceAlgorithm;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.seqreads.GeneSequenceReads;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.utilities.AssertUtils;
import edu.stanford.hivdb.utilities.MySetUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.Mutation;

public class GeneDR<VirusT extends Virus<VirusT>> {

	// Required data structures used to instantiate the class
	private final Gene<VirusT> gene;
	private final MutationSet<VirusT> mutations;
	private final DrugResistanceAlgorithm<VirusT> algorithm;

	// Optional data structures used to instantiate the class
	private final transient EvaluatedGene evalGene;
	private final transient SortedSet<ASIDrugSusc<VirusT>> drugSuscs;
	private transient SortedMap<MutationType<VirusT>, MutationSet<VirusT>> mutTypes;
	private transient SortedMap<CommentType, List<BoundComment<VirusT>>> commentsByTypes;

	public static <VirusT extends Virus<VirusT>> SortedMap<Gene<VirusT>, GeneDR<VirusT>> newFromAlignedGeneSeqs(
		List<AlignedGeneSeq<VirusT>> alignedGeneSeqs, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		return alignedGeneSeqs
			.parallelStream()
			.collect(Collectors.toMap(
				geneSeq -> geneSeq.getGene(),
				geneSeq -> new GeneDR<>(geneSeq.getGene(), geneSeq.getMutations(), algorithm),
				(g1, g2) -> g1,
				TreeMap::new
			));
	}

	public static <VirusT extends Virus<VirusT>> SortedMap<Gene<VirusT>, GeneDR<VirusT>> newFromGeneSequenceReads(
		List<GeneSequenceReads<VirusT>> allGeneSeqReads, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		return allGeneSeqReads
			.parallelStream()
			.collect(Collectors.toMap(
				geneReads -> geneReads.getGene(),
				geneReads -> new GeneDR<>(geneReads.getGene(), geneReads.getMutations(), algorithm),
				(r1, r2) -> r1,
				TreeMap::new
			));
	}

	public static <VirusT extends Virus<VirusT>> SortedMap<MutationSet<VirusT>, GeneDR<VirusT>> newFromMutationSets(
		Gene<VirusT> gene, Set<MutationSet<VirusT>> allMuts, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		return allMuts
			.parallelStream()
			.collect(Collectors.toMap(
				muts -> muts,
				muts -> new GeneDR<>(gene, muts, algorithm),
				(m1, m2) -> m1,
				TreeMap::new
			));
	}

	public GeneDR(Gene<VirusT> gene, AlignedGeneSeq<VirusT> seq, DrugResistanceAlgorithm<VirusT> algorithm) {
		this(gene, seq.getMutations(), algorithm);
	}

	public GeneDR(Gene<VirusT> gene, MutationSet<VirusT> mutations, DrugResistanceAlgorithm<VirusT> algorithm) {
		mutations = mutations.filterBy(mut -> !mut.isUnsequenced());
		this.gene = gene;
		this.mutations = mutations;
		this.algorithm = algorithm;
		evalGene = ASIResultHandler.evalutateGeneMutations(gene, mutations, algorithm);
		drugSuscs = ASIResultHandler.extractDrugSuscs(gene, evalGene, algorithm);
		
	}

	public final Gene<VirusT> getGene() { return gene; }

	public DrugResistanceAlgorithm<VirusT> getAlgorithm() {
		return algorithm;
	}
	
	@Deprecated
	public DrugResistanceAlgorithm<VirusT> getVersion() {
		// use getAlgorithm() instead
		return algorithm;
	}

	public List<BoundComment<VirusT>> getComments() {
		VirusT virusIns = gene.getVirusInstance();
		List<BoundComment<VirusT>> comments = (
			ASIResultHandler.extractMutationComments(
				virusIns,
				(Collection<?>) evalGene.getGeneCommentDefinitions(),
				getMutations()
			)
		);
		comments.addAll(
			ASIResultHandler.extractDrugLevelComments(
				virusIns,
				evalGene.getEvaluatedResultCommentRules()
			)
		);
		return comments;
	}

	public List<BoundComment<VirusT>> getComments(CommentType commentType) {
		return groupCommentsByTypes().getOrDefault(commentType, Collections.emptyList());
	}

	@Deprecated
	public Map<Drug<VirusT>, Double> getTotalDrugScores() {
		return (
			getDrugSuscs()
			.stream()
			.collect(Collectors.toMap(
				ds -> ds.getDrug(),
				ds -> ds.getScore()
			))
		);
	}
	
	@Deprecated
	public Map<Drug<VirusT>, Double> getTotalDrugScores(DrugClass<VirusT> drugClass) {
		return (
			getDrugSuscs(ds -> ds.drugClassIs(drugClass))
			.stream()
			.collect(Collectors.toMap(
				ds -> ds.getDrug(),
				ds -> ds.getScore()
			))
		);
	}

	@Deprecated
	public Double getTotalDrugScore(Drug<VirusT> drug) {
		return getDrugSusc(drug).getScore();
	}

	@Deprecated
	public Integer getDrugLevel(Drug<VirusT> drug) {
		return getDrugSusc(drug).getLevel();
	}

	@Deprecated
	public String getDrugLevelText(Drug<VirusT> drug) {
		return getDrugSusc(drug).getLevelText();
	}

	@Deprecated
	public String getDrugLevelSIR(Drug<VirusT> drug) {
		return getDrugSusc(drug).getSIR().toString();
	}
	
	public SortedMap<MutationType<VirusT>, MutationSet<VirusT>> groupMutationsByTypes() {
		if (mutTypes == null) {
			mutTypes = mutations.groupByMutType(gene);
		}
		return mutTypes;
	}

	public MutationSet<VirusT> getMutations() { return mutations; }

	public MutationSet<VirusT> getMutations(MutationType<VirusT> mutType) {
		return groupMutationsByTypes().getOrDefault(mutType, new MutationSet<>());
	}
	
	public SortedMap<CommentType, List<BoundComment<VirusT>>> groupCommentsByTypes() {
		if (commentsByTypes == null) {
			commentsByTypes = getComments()
				.stream()
				.collect(Collectors.groupingBy(
					cmt -> cmt.getType(),
					TreeMap::new,
					Collectors.toList()
				));
		}
		return commentsByTypes;
	}

	public final SortedSet<MutationSet<VirusT>> getScoredMutations() {
		return (
			drugSuscs
			.parallelStream()
			.flatMap(ds -> ds.getPartialScores().keySet().stream())
			.collect(Collectors.toCollection(TreeSet::new))
		);
	}

	public final SortedSet<MutationSet<VirusT>> getScoredMutations(Predicate<? super ASIDrugSusc<VirusT>> predicate) {
		return (
			drugSuscs
			.parallelStream()
			.filter(predicate)
			.flatMap(ds -> ds.getPartialScores().keySet().stream())
			.collect(Collectors.toCollection(TreeSet::new))
		);
	}
	
	public final SortedSet<ASIDrugSusc<VirusT>> getDrugSuscs() {
		return new TreeSet<>(drugSuscs);
	}
	
	public final SortedSet<ASIDrugSusc<VirusT>> getDrugSuscs(Predicate<? super ASIDrugSusc<VirusT>> predicate) {
		return MySetUtils.filter(drugSuscs, predicate);
	}
	
	public final ASIDrugSusc<VirusT> getDrugSusc(Drug<VirusT> drug) {
		AssertUtils.isTrue(
			drug.getDrugClass().getAbstractGene().equals(gene.getAbstractGene()),
			"The input drug %s is for gene %s, but this GeneDR object is for %s",
			drug.getName(), drug.getDrugClass().getAbstractGene(), gene.getAbstractGene()
		);
		SortedSet<ASIDrugSusc<VirusT>> drugSusc = MySetUtils.filter(drugSuscs, ds -> ds.drugIs(drug));
		AssertUtils.isTrue(
			drugSusc.size() == 1,
			RuntimeException.class,
			"Expect one DrugSusc data for %s but received %d",
			drug, drugSusc.size()
		);
		return drugSusc.first();
	}

	@Deprecated
	public final Map<Mutation<VirusT>, Double> getMutScores(Drug<VirusT> drug) {
		return getDrugSusc(drug).getSingleMutPartialScores();
	}

	@Deprecated
	public Map<Drug<VirusT>, Map<MutationSet<VirusT>, Double>> getDrugComboMutScores(DrugClass<VirusT> drugClass) {
		Map<Drug<VirusT>, Map<MutationSet<VirusT>, Double>> multiMutsPartialScores = new TreeMap<>();
		SortedSet<ASIDrugSusc<VirusT>> partialDrugSuscs = MySetUtils.filter(drugSuscs, ds -> ds.drugClassIs(drugClass));
		for (ASIDrugSusc<VirusT> drugSusc : partialDrugSuscs) {
			multiMutsPartialScores.put(drugSusc.getDrug(), drugSusc.getMultiMutsPartialScores());
		}
		return multiMutsPartialScores;
	}
	
	@Deprecated
	public final Map<MutationSet<VirusT>, Double> getComboMutScores(Drug<VirusT> drug) {
		Map<MutationSet<VirusT>, Double> multiMutsPartialScores = new TreeMap<>();
		ASIDrugSusc<VirusT> drugSusc = getDrugSusc(drug);
		multiMutsPartialScores.putAll(drugSusc.getMultiMutsPartialScores());
		return multiMutsPartialScores;
	}
	
	@Deprecated
	public final boolean hasScoredMuts(DrugClass<VirusT> drugClass) {
		return hasScoredIndividualMuts(drugClass) || hasScoredComboMuts(drugClass);
	}

	@Deprecated
	public final boolean hasScoredMuts(Drug<VirusT> drug) {
		return hasScoredIndividualMuts(drug) || hasScoredComboMuts(drug);
	}
	
	@Deprecated
	public final boolean hasScoredIndividualMuts(DrugClass<VirusT> drugClass) {
		return (
			MySetUtils.anyMatch(
				drugSuscs,
				ds -> (
					ds.drugClassIs(drugClass) &&
					ds.hasSingleMutPartialScore()
				)
			)
		);
	}

	@Deprecated
	public final boolean hasScoredIndividualMuts(Drug<VirusT> drug) {
		return (
			MySetUtils.anyMatch(
				drugSuscs,
				ds -> (
					ds.drugIs(drug) &&
					ds.hasSingleMutPartialScore()
				)
			)
		);
	}

	@Deprecated
	public final boolean hasScoredComboMuts(DrugClass<VirusT> drugClass) {
		return (
			MySetUtils.anyMatch(
				drugSuscs,
				ds -> (
					ds.drugClassIs(drugClass) &&
					ds.hasMultiMutsPartialScore()
				)
			)
		);
	}

	@Deprecated
	public final boolean hasScoredComboMuts(Drug<VirusT> drug) {
		return (
			MySetUtils.anyMatch(
				drugSuscs,
				ds -> (
					ds.drugIs(drug) &&
					ds.hasMultiMutsPartialScore()
				)
			)
		);
	}

}
