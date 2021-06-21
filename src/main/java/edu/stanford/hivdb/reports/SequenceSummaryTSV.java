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
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.stanford.hivdb.utilities.TSV;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.utilities.NumberFormats;


public class SequenceSummaryTSV<VirusT extends Virus<VirusT>> {
	
	private static final Map<String, SequenceSummaryTSV<? extends Virus<?>>> singletons = new HashMap<>();

	@SuppressWarnings("unchecked")
	public static <VirusT extends Virus<VirusT>> SequenceSummaryTSV<VirusT> getInstance(VirusT virusIns) {
		String name = virusIns.getName();
		if (!singletons.containsKey(name)) {
			singletons.put(name, new SequenceSummaryTSV<>(virusIns));
		}
		return (SequenceSummaryTSV<VirusT>) singletons.get(name);
	}

	private final List<String> headerFields;
	
	private SequenceSummaryTSV(VirusT virusIns) {
		List<String> hFields = new ArrayList<>();
		hFields.add("Sequence Name");
		hFields.add("Genes");
		for (String geneText : virusIns.getAbstractGenes()) {
			hFields.add(String.format("%s Start", geneText));
			hFields.add(String.format("%s End", geneText));
		}
		hFields.add("Subtype (%)");
		hFields.add("Pcnt Mix");
		Strain<VirusT> mainStrain = virusIns.getMainStrain();
		MutationType<VirusT> otherMutType = virusIns.getOtherMutationType();
		for (Gene<VirusT> gene : mainStrain.getGenes()) {
			String otherMutTypeField = otherMutType.getName(); 
			for (DrugClass<VirusT> drugClass : gene.getDrugClasses()) {
				for (MutationType<VirusT> mutType : drugClass.getMutationTypes()) {
					if (mutType == otherMutType) {
						otherMutTypeField = mutType.getFullName(drugClass);
						continue;
					}
					hFields.add(mutType.getFullName(drugClass));
				}
			}
			hFields.add(otherMutTypeField);
		}
		for (DrugClass<VirusT> drugClass : virusIns.getDrugClasses()) {
			hFields.add(String.format("%s SDRMs", drugClass));
		}
		for (DrugClass<VirusT> drugClass : virusIns.getDrugClasses()) {
			hFields.add(String.format("%s TSMs", drugClass));
		}
		hFields.add("Num Frame Shifts");
		hFields.add("Frame Shifts");
		hFields.add("Num Insertions");
		hFields.add("Insertions");
		hFields.add("Num Deletions");
		hFields.add("Deletions");
		hFields.add("Num Stop Codons");
		hFields.add("StopCodons");
		hFields.add("Num BDHVN");
		hFields.add("BDHVN");
		hFields.add("Num Apobec Mutations");
		hFields.add("Apobec Mutations");
		hFields.add("Num Unusual Mutations");
		hFields.add("UnusualMutations");
		this.headerFields = Collections.unmodifiableList(hFields);
	}
	
	protected List<String> getHeaderFields() {
		return headerFields;
	}
	
	protected List<Map<String, String>> getReportRows(List<AlignedSequence<VirusT>> overallResults) {

		List<Map<String, String>> sequenceRows = new ArrayList<>();

		for (AlignedSequence<VirusT> alignedSeq : overallResults) {
			Map<String, String> sequenceRecord = new HashMap<>();

			List<Gene<VirusT>> geneList = alignedSeq.getAvailableGenes();
			MutationSet<VirusT> seqMutations = alignedSeq.getMutations();
			String seqName = alignedSeq.getInputSequence().getHeader();

			String genes = (
				geneList
				.stream()
				.map(g -> g.getAbstractGene())
				.collect(Collectors.joining(","))
			);

			// sequenceName
			sequenceRecord.put("Sequence Name", seqName);

			// Genes
			sequenceRecord.put("Genes", genes);

			// PRStart, PREnd, RTStart, RTEnd, INStart, INEnd
			sequenceRecord.putAll(determineGeneBoundaries(alignedSeq));

			// Subtype (%)
			sequenceRecord.put("Subtype (%)", alignedSeq.getGenotypeText());

			// PcntMix
			sequenceRecord.put(
				"Pcnt Mix",
				NumberFormats.prettyDecimalAsString(alignedSeq.getMixtureRate() * 100));

			// %(DrugClass)s %(MutationType)s
			sequenceRecord.putAll(determineMutLists(alignedSeq));

		  	// PI SDRMs, NRTI SDRMs, NNRTI SDRMs, INSTI SDRMs
			sequenceRecord.putAll(determineSdrms(alignedSeq));

		  	// PI-TSMs, NRTI-TSMs, NNRTI-TSMs, INSTI-TSMs
			sequenceRecord.putAll(determineNonDrmTsms(alignedSeq));

		  	// NumFS, FrameShifts
			sequenceRecord.putAll(determineFrameShiftText(alignedSeq));
			// NumIns, Insertions
			sequenceRecord.putAll(determineSeqInsertions(seqMutations));
			// NumDel, Deletions
			sequenceRecord.putAll(determineSeqDeletions(seqMutations));
			// NumStops, StopCodons
			sequenceRecord.putAll(determineSeqStopCodons(seqMutations));
			// NumBDHVN, BDHVN
			sequenceRecord.putAll(determineSeqBDHVN(seqMutations));
			// NumApobec, ApobecMuts
			sequenceRecord.putAll(determineApobecFields(seqMutations));
			// NumUnusual, UnusualMuts
			sequenceRecord.putAll(determineSeqUnusualMuts(seqMutations));
			sequenceRows.add(sequenceRecord);
		}

		return sequenceRows;
	}

	public String getReport(List<AlignedSequence<VirusT>> overallResults) {
		return TSV.dumpMaps(headerFields, getReportRows(overallResults), "NA");
	}

	private Map<String, String> mutationListToTabularResult(MutationSet<VirusT> mutations, String numField, String mutsField) {
		Map<String, String> fields = new HashMap<>();
		String text = "None";
		if (mutations.size() > 0) {
			text = mutations.join();
		}
		fields.put(numField, "" + mutations.size());
		fields.put(mutsField, text);
		return fields;
	}

	private Map<String, String> determineApobecFields(MutationSet<VirusT> mutations) {
		MutationSet<VirusT> apobecMuts = mutations.getApobecMutations();
		return mutationListToTabularResult(apobecMuts, "Num Apobec Mutations", "Apobec Mutations");
	}

	private Map<String, String> determineSeqUnusualMuts(MutationSet<VirusT> mutations) {
		MutationSet<VirusT> unusualMuts = mutations.getUnusualMutations();
		return mutationListToTabularResult(unusualMuts, "Num Unusual Mutations", "UnusualMutations");
	}

	// TODO: What if bdhvn does not affect the amino acid. Check how this is handled
	private Map<String, String> determineSeqBDHVN(MutationSet<VirusT> mutations) {
		MutationSet<VirusT> bdhvnMuts = mutations.getAmbiguousCodons();
		return mutationListToTabularResult(bdhvnMuts, "Num BDHVN", "BDHVN");
	}

	private Map<String, String> determineSeqStopCodons(MutationSet<VirusT> mutations) {
		MutationSet<VirusT> stopCodons = mutations.getStopCodons();
		return mutationListToTabularResult(stopCodons, "Num Stop Codons", "StopCodons");
	}

	private Map<String, String> determineSeqDeletions(MutationSet<VirusT> mutations) {
		MutationSet<VirusT> deletions = mutations.getDeletions();
		return mutationListToTabularResult(deletions, "Num Deletions", "Deletions");
	}

	private Map<String, String> determineSeqInsertions(MutationSet<VirusT> mutations) {
		MutationSet<VirusT> insertions = mutations.getInsertions();
		return mutationListToTabularResult(insertions, "Num Insertions", "Insertions");
	}


	private Map<String, String> determineFrameShiftText(AlignedSequence<VirusT> alignedSeq) {
		List<FrameShift<VirusT>> frameShifts = alignedSeq.getFrameShifts();
		Map<String, String> frameShiftFields = new HashMap<>();
		String frameShiftsString = FrameShift.joinFrameShifts(frameShifts);
		frameShiftFields.put(
			"Num Frame Shifts",
			Integer.toString(frameShifts.size()));
		frameShiftFields.put("Frame Shifts", frameShiftsString);
		return frameShiftFields;
	}

	private Map<String, String> determineSdrms(AlignedSequence<VirusT> alignedSeq) {
		return (
			alignedSeq
			.getMutations()
			.filterAndGroupBy(Mutation::isSDRM, Mutation::getSDRMDrugClass)
			.entrySet()
			.stream()
			.collect(Collectors.toMap(
				e -> String.format("%s SDRMs", e.getKey()),
				e -> e.getValue().join()
			))
		);
	}

	public Map<String, String> determineMutLists(AlignedSequence<VirusT> alignedSeq) {
		Map<String, String> mutListStrings = new HashMap<>();
		for (AlignedGeneSeq<VirusT> geneSeq : alignedSeq.getAlignedGeneSequences()) {
			for (DrugClass<VirusT> drugClass : geneSeq.getGene().getDrugClasses()) {
				for (MutationType<VirusT> mutType : drugClass.getMutationTypes()) {
					MutationSet<VirusT> mutTypeMutations = geneSeq.getMutationsByMutType(mutType);
					mutListStrings.put(mutType.getFullName(drugClass), mutTypeMutations.join());
				}
			}
		}
		return mutListStrings;
	}

	private Map<String, String> determineNonDrmTsms(AlignedSequence<VirusT> alignedSeq) {
		return (
			alignedSeq
			.getMutations()
			.filterAndGroupBy(Mutation::isTSM, Mutation::getTSMDrugClass)
			.entrySet()
			.stream()
			.collect(Collectors.toMap(
				e -> String.format("%s TSMs", e.getKey()),
				e -> e.getValue().join()
			))
		);
	}

	private Map<String, String> determineGeneBoundaries(AlignedSequence<VirusT> alignedSeq) {
		Map<String, String> geneBoundaries = new HashMap<>();
		for (AlignedGeneSeq<VirusT> geneSeq : alignedSeq.getAlignedGeneSequences()) {
			geneBoundaries.put(
				String.format("%s Start", geneSeq.getAbstractGene()),
				"" + geneSeq.getFirstAA()
			);
			geneBoundaries.put(
				String.format("%s End", geneSeq.getAbstractGene()),
				"" + geneSeq.getLastAA()
			);
		}
		return geneBoundaries;
	}

}
