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

package edu.stanford.hivdb.reports;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.stream.Collectors;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import edu.stanford.hivdb.drugresistance.GeneDR;
import edu.stanford.hivdb.drugresistance.algorithm.DrugResistanceAlgorithm;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
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
		MutationType<VirusT> otherMutType = virusIns.getOtherMutationType();

		List<String> hFields = new ArrayList<>();
		hFields.add("Sequence Name");
		hFields.add("Strain");
		hFields.add("Genes");
		for (DrugClass<VirusT> drugClass : virusIns.getDrugClasses()) {
			for (MutationType<VirusT> mtype : drugClass.getMutationTypes()) {
				if (mtype == otherMutType) {
					continue;
				}
				hFields.add(mtype.getFullName(drugClass));
			}
			for (Drug<VirusT> drug : drugClass.getDrugs()) {
				hFields.add(String.format("%s Score", drug.getDisplayAbbr()));
				hFields.add(String.format("%s Level", drug.getDisplayAbbr()));
			}
		}
		hFields.add("Algorithm Name");
		hFields.add("Algorithm Version");
		hFields.add("Algorithm Date");
		headerFields = Collections.unmodifiableList(hFields);
	}

	protected List<String> getHeaderFields() {
		return headerFields;
	}

	protected List<Map<String, String>> getReportRows(
		List<AlignedSequence<VirusT>> alignedSeqs, List<Map<Gene<VirusT>,
		GeneDR<VirusT>>> allResistanceResults, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		List<Map<String, String>> sequenceRows = new ArrayList<>();

		int numSeqs = alignedSeqs.size();
		for (int i=0; i < numSeqs; i ++) {
			AlignedSequence<VirusT> alignedSeq = alignedSeqs.get(i);
			Map<Gene<VirusT>, GeneDR<VirusT>> resistanceResults = allResistanceResults.get(i);

			String seqName = alignedSeq.getInputSequence().getHeader();

			List<Gene<VirusT>> geneList = Lists.newArrayList(resistanceResults.keySet());
			geneList.sort(Gene::compareTo);
			String genes = (
				geneList.stream()
				.map(gene -> gene.getAbstractGene())
				.collect(Collectors.joining(","))
			);
			Strain<VirusT> strain = alignedSeq.getStrain();

			Map<String, String> sequenceRecord = new HashMap<>();
			sequenceRecord.put("Sequence Name", seqName);
			sequenceRecord.put("Strain", strain.getName());
			sequenceRecord.put("Genes", genes);
			for (String absGene : virusIns.getAbstractGenes()) {
				Gene<VirusT> gene = strain.getGene(absGene);
				GeneDR<VirusT> result = resistanceResults.getOrDefault(gene, null);
				sequenceRecord.putAll(getScoreDetails(result));
			}

			sequenceRecord.put("Algorithm Name", algorithm.getFamily());
			sequenceRecord.put("Algorithm Version", algorithm.getVersion());
			sequenceRecord.put("Algorithm Date", algorithm.getPublishDate());
			sequenceRows.add(sequenceRecord);

		}
		return sequenceRows;
	}

	public String getReport(
		List<AlignedSequence<VirusT>> alignedSeqs, List<Map<Gene<VirusT>,
		GeneDR<VirusT>>> allResistanceResults, DrugResistanceAlgorithm<VirusT> algorithm
	) {
		return TSV.dumpMaps(
			headerFields,
			getReportRows(alignedSeqs, allResistanceResults, algorithm),
			"NA");
	}

	private Map<String, String> getScoreDetails(GeneDR<VirusT> geneDR) {
		if (geneDR == null) {
			return Collections.emptyMap();
		}
		MutationType<VirusT> otherMutType = virusIns.getOtherMutationType();
		Map<String, String> results = new HashMap<>();

		for (DrugClass<VirusT> drugClass : geneDR.getGene().getDrugClasses()) {
			for (Drug<VirusT> drug : drugClass.getDrugs()) {
				int score = geneDR.getDrugSusc(drug).getScore().intValue();
				int level = geneDR.getDrugSusc(drug).getLevel();
				results.put(
					String.format("%s Score", drug.getDisplayAbbr()),
					Integer.toString(score));
				results.put(
					String.format("%s Level", drug.getDisplayAbbr()),
					Integer.toString(level));
			}
			for (MutationType<VirusT> mtype : drugClass.getMutationTypes()) {
				if (mtype == otherMutType) {
					continue;
				}
				results.put(
					mtype.getFullName(drugClass),
					geneDR.getMutations(mtype).join()
				);
			}
		}
		return results;
	}

}
