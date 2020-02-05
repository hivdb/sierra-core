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

package edu.stanford.hivdb.reports;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import edu.stanford.hivdb.drugresistance.GeneDR;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.drugs.DrugResistanceAlgorithm;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.utilities.TSV;


/**
 * Print the sequence name, list of genes, "major" and "minor" mutations for each gene,
 *   resistance scores and levels for selected drugs
 * Currrently the header and list of drugs for a class (not the complete list) are constants
 *
 */
public class ResistanceSummaryTSV<VirusT extends Virus<VirusT>> {

	private static final Map<String, ResistanceSummaryTSV<? extends Virus<?>>> singletons = new HashMap<>();

	@SuppressWarnings("unchecked")
	public static <VirusT extends Virus<VirusT>> ResistanceSummaryTSV<VirusT> getInstance(VirusT virusIns) {
		String name = virusIns.getName();
		if (!singletons.containsKey(name)) {
			singletons.put(name, new ResistanceSummaryTSV<>(virusIns));
		}
		return (ResistanceSummaryTSV<VirusT>) singletons.get(name);
	}

	private final List<String> headerFields;
	private final VirusT virusIns;

	private ResistanceSummaryTSV(VirusT virusIns) {
		this.virusIns = virusIns;

		List<String> hFields = new ArrayList<>();
		hFields.add("Sequence Name");
		hFields.add("Strain");
		hFields.add("Genes");
		for (DrugClass<VirusT> drugClass : virusIns.getDrugClasses()) {
			for (MutationType<VirusT> mtype : drugClass.getMutationTypes()) {
				String mtypeText = mtype.getName();
				if (mtypeText.equals("Other")) {
					continue;
				}
				hFields.add(mtype.getFullName(drugClass));
			}
			for (Drug<VirusT> drug : drugClass.getDrugs()) {
				hFields.add(drug.getDisplayAbbr() + " Score");
				hFields.add(drug.getDisplayAbbr() + " Level");
			}
		}
		hFields.add("Algorithm Name");
		hFields.add("Algorithm Version");
		hFields.add("Algorithm Date");
		headerFields = Collections.unmodifiableList(hFields);
	}

	public String makeReport(
		List<AlignedSequence<VirusT>> alignedSeqs, List<Map<Gene<VirusT>,
		GeneDR<VirusT>>> allResistanceResults, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		List<List<String>> sequenceRows = new ArrayList<>();
		Map<String, Map<String, String>> tabularResults = new TreeMap<>();

		int numSeqs = alignedSeqs.size();
		for (int i=0; i < numSeqs; i ++) {
			AlignedSequence<VirusT> alignedSeq = alignedSeqs.get(i);
			Map<Gene<VirusT>, GeneDR<VirusT>> resistanceResults = allResistanceResults.get(i);

			String seqName = alignedSeq.getInputSequence().getHeader();
			tabularResults.put(seqName, new TreeMap<String, String>());

			List<Gene<VirusT>> geneList = Lists.newArrayList(resistanceResults.keySet());
			geneList.sort(Gene::compareTo);
			String genes = (
				geneList.stream()
				.map(gene -> gene.getAbstractGene())
				.collect(Collectors.joining(","))
			);
			Strain<VirusT> strain = alignedSeq.getStrain();

			List<String> sequenceRecord = new ArrayList<>();
			sequenceRecord.add(seqName);
			sequenceRecord.add(strain.getName());
			sequenceRecord.add(genes);
			for (String absGene : virusIns.getAbstractGenes()) {
				Gene<VirusT> gene = strain.getGene(absGene);
				GeneDR<VirusT> result = resistanceResults.getOrDefault(gene, null);
				for (DrugClass<VirusT> drugClass : gene.getDrugClasses()) {
					sequenceRecord.addAll(getScoredMutations(drugClass, result));
					sequenceRecord.addAll(getScores(drugClass, result));
				}
			}

			sequenceRecord.add(algorithm.getFamily());
			sequenceRecord.add(algorithm.getVersion());
			sequenceRecord.add(algorithm.getPublishDate());
			sequenceRows.add(sequenceRecord);

			for (int j=0; j<headerFields.size(); j++) {
				String field = headerFields.get(j);
				String dataItem = sequenceRecord.get(j);
				tabularResults.get(seqName).put(field, dataItem);
			}
		}
		return TSV.dumps(headerFields, sequenceRows);
	}

	private List<String> getScores(DrugClass<VirusT> drugClass, GeneDR<VirusT> geneDR) {
		List<String> resistanceScoresAndLevels = new ArrayList<>();

		for (Drug<VirusT> drug : drugClass.getDrugs()) {
			if (geneDR == null) {
				resistanceScoresAndLevels.add("NA");
				resistanceScoresAndLevels.add("NA");
			}
			else {
				int score = geneDR.getTotalDrugScore(drug).intValue();
				int level = geneDR.getDrugLevel(drug);
				resistanceScoresAndLevels.add(Integer.toString(score));
				resistanceScoresAndLevels.add(Integer.toString(level));
			}
		}
		return resistanceScoresAndLevels;
	}


	private List<String> getScoredMutations(DrugClass<VirusT> drugClass, GeneDR<VirusT> geneDR) {
		List<String> scoredMutations = new ArrayList<>();
		Map<MutationType<VirusT>, MutationSet<VirusT>> mutationsByTypes;
		if (geneDR == null) {
			mutationsByTypes = Collections.emptyMap();
		}
		else {
			mutationsByTypes = geneDR.groupMutationsByTypes();
		}

		for (MutationType<VirusT> mtype : drugClass.getMutationTypes()) {
			if (mtype.toString().equals("Other")) {
				continue;
			}
			String mutationsText = "NA";
			if (mutationsByTypes.containsKey(mtype)) {
				MutationSet<VirusT> muts = mutationsByTypes.get(mtype);
				if (muts.size() > 0) {
					mutationsText = muts.join();
				}
				else {
					mutationsText = "None";
				}
			}
			scoredMutations.add(mutationsText);
		}
		return scoredMutations;
	}


}
