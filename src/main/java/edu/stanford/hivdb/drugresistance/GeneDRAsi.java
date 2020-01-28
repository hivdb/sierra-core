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

package edu.stanford.hivdb.drugresistance;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.fstrf.stanfordAsiInterpreter.resistance.evaluate.EvaluatedGene;

import edu.stanford.hivdb.comments.BoundComment;
import edu.stanford.hivdb.drugresistance.algorithm.AsiResult;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.drugs.DrugResistanceAlgorithm;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.seqreads.GeneSequenceReads;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

/**
 * Initialize with Gene, mutTypes, mutCommentsAsi, drugClassTotalDrugScoresAsi,
 * drugClassDrugMutScoresAsi, and drugClassDrugComboMutScoresAsi
 *
 * drugClassDrugMutScoresAsi, and drugClassDrugComboMutScoresAsi differ from those that instantiate GeneDRHivdb
 * as they have keys for each drugClass whether or not there are DRMs associated with that class
 *
 * drugClassTotalDrugScoresAsi is supplied by AsiHivdb but not by any of the classes that directly
 * query HIVDB_Scores
 *
 * drugClassDrugMutScoresAsi, and drugClassDrugComboMutScoresAsi need to be manipulated to four new maps
 * to provide all data required for reports. These new maps should be identical to the four maps created in
 * GeneDRHivdb
 */
public class GeneDRAsi<VirusT extends Virus<VirusT>> extends GeneDR<VirusT> implements WithGene<VirusT> {

	protected final AsiResult<VirusT> asiObject;
	private DrugResistanceAlgorithm<VirusT> algorithm;
	
	public static <VirusT extends Virus<VirusT>> Map<Gene<VirusT>, GeneDR<VirusT>> getResistanceByGeneFromAlignedGeneSeqs(
		List<AlignedGeneSeq<VirusT>> alignedGeneSeqs, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		Map<Gene<VirusT>, GeneDR<VirusT>> resistanceForSequence = new TreeMap<>();
		for (AlignedGeneSeq<VirusT> geneSeq : alignedGeneSeqs) {
			Gene<VirusT> gene = geneSeq.getGene();
			final GeneDR<VirusT> geneDR = new GeneDRAsi<>(gene, geneSeq.getMutations(), algorithm);
			resistanceForSequence.put(gene, geneDR);
		}
		return resistanceForSequence;
	}

	public static <VirusT extends Virus<VirusT>> Map<Gene<VirusT>, GeneDR<VirusT>> getResistanceByGeneFromReads(
		List<GeneSequenceReads<VirusT>> allGeneSeqReads, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		Map<Gene<VirusT>, GeneDR<VirusT>> resistanceForSequence = new TreeMap<>();
		for (GeneSequenceReads<VirusT> geneSeqReads : allGeneSeqReads) {
			Gene<VirusT> gene = geneSeqReads.getGene();
			final GeneDR<VirusT> geneDR = new GeneDRAsi<>(gene, geneSeqReads.getMutations(), algorithm);
			resistanceForSequence.put(gene, geneDR);
		}
		return resistanceForSequence;
	}

	public GeneDRAsi(Gene<VirusT> gene, AlignedGeneSeq<VirusT> seq, DrugResistanceAlgorithm<VirusT> algorithm) {
		this(gene, seq.getMutations(), algorithm);
	}

	public GeneDRAsi(Gene<VirusT> gene, MutationSet<VirusT> mutations, DrugResistanceAlgorithm<VirusT> algorithm) {
		super(gene, mutations);
		this.algorithm = algorithm;
		this.asiObject = new AsiResult<>(gene, mutations, algorithm);

		drugClassDrugMutScores = asiObject.getDrugClassDrugMutScores();
		drugClassDrugComboMutScores = asiObject.getDrugClassDrugComboMutScores();
		postConstructor();
	}

	public static <VirusT extends Virus<VirusT>> Map<MutationSet<VirusT>, GeneDR<VirusT>> parallelConstructor(
		Gene<VirusT> gene, Set<MutationSet<VirusT>> allMuts, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		return allMuts
			.parallelStream()
			.collect(Collectors.toMap(
				muts -> muts,
				muts -> new GeneDRAsi<>(gene, muts, algorithm)
			));
	}
	
	public AsiResult<VirusT> getAsiObject() {
		return asiObject;
	}
	
	public DrugResistanceAlgorithm<VirusT> getAlgorithm() {
		return algorithm;
	}
	
	public DrugResistanceAlgorithm<VirusT> getVersion() {
		return algorithm;
	}

	@Override
	public List<BoundComment<VirusT>> getAllComments() {
		EvaluatedGene evaluatedGene = asiObject.getEvaluatedGene();
		List<BoundComment<VirusT>> comments = (
			gene.getVirusInstance().getConditionalComments().fromAsiMutationComments(
				(Collection<?>) evaluatedGene.getGeneCommentDefinitions(),
				getMutations()
			)
		);
		comments.addAll(
			gene.getVirusInstance().getConditionalComments().fromAsiDrugLevelComments(
				evaluatedGene.getEvaluatedResultCommentRules()
			)
		);
		return comments;
	}

	@Override
	public Map<Drug<VirusT>, Double> getDrugClassTotalDrugScores(DrugClass<VirusT> drugClass) {
		return asiObject.getDrugClassTotalDrugScores(drugClass);
	}

	@Override
	public Double getTotalDrugScore(Drug<VirusT> drug) {
		return asiObject.getTotalScore(drug);
	}

	@Override
	public Integer getDrugLevel(Drug<VirusT> drug) {
		return asiObject.getDrugLevel(drug);
	}

	@Override
	public String getDrugLevelText(Drug<VirusT> drug) {
		return asiObject.getDrugLevelText(drug);
	}

	@Override
	public String getDrugLevelSIR(Drug<VirusT> drug) {
		return asiObject.getDrugLevelSir(drug);
	}
}
