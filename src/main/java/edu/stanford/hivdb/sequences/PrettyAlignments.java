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

package edu.stanford.hivdb.sequences;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.StrainModifier;
import edu.stanford.hivdb.utilities.TSV;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class PrettyAlignments<VirusT extends Virus<VirusT>> {
	private transient Gene<VirusT> targetGene;
	private transient List<AlignedGeneSeq<VirusT>> alignedGeneSeqs = new ArrayList<>();
	private Map<String, Map<Integer, String>> sequenceAllPosAAs = new LinkedHashMap<>();

	/**
	 * This receives a list of alignedGeneSeqs for a single gene
	 * 
	 * TODO: by providing targetGene, this class need to apply StrainModifier to each of the pretty alignment results 
	 * 
	 * @param targetGene		Target gene
	 * @param alignedGeneSeqs	Aligned gene sequences
	 */
	public PrettyAlignments (Gene<VirusT> targetGene, List<AlignedGeneSeq<VirusT>> alignedGeneSeqs) {
		this.targetGene = targetGene;
		this.alignedGeneSeqs = alignedGeneSeqs;
		createSequencePosAAsMap();
	}

	public int getFirstAA() { return 1; }
	public int getLastAA() { return targetGene.getAASize(); }
	public Map<String, Map<Integer, String>> getSequenceAllPosAAs() {
		return sequenceAllPosAAs;
	}

	@Override
	public String toString() {
		List<String> headers = new ArrayList<>();
		List<List<String>> rows = new ArrayList<>();

		if (sequenceAllPosAAs.isEmpty()) { return ""; }

		headers.add("Sequence Names");
		int firstAA = getFirstAA();
		int lastAA = getLastAA();

		List<String> consRow = new ArrayList<>();
		String cons = targetGene.getRefSequence();
		consRow.add("Consensus");

		for (int pos=firstAA; pos <= lastAA; pos ++) {
			headers.add("" + pos);
			consRow.add("" + cons.charAt(pos - 1));
		}

		rows.add(consRow);
		rows.addAll(
			sequenceAllPosAAs
			.entrySet()
			.stream()
			.map(e -> {
				String seqName = e.getKey();
				Map<Integer, String> posAAs = e.getValue();

				List<String> row = new ArrayList<>();
				row.add(seqName);
				for (int pos=firstAA; pos <= lastAA; pos ++) {
					row.add(posAAs.get(pos));
				}
				return row;
			})
			.collect(Collectors.toList())
		);
		return TSV.dumps(headers, rows);
	}

	private void createSequencePosAAsMap() {
		int firstAA = getFirstAA();
		int lastAA = getLastAA();
		Strain<VirusT> targetStrain = targetGene.getStrain();
		for (AlignedGeneSeq<VirusT> alignedGeneSeq : alignedGeneSeqs) {
			String seqName = alignedGeneSeq.getSequence().getHeader();
			int seqFirstPos = alignedGeneSeq.getFirstAA();
			int seqLastPos = alignedGeneSeq.getLastAA();
			Gene<VirusT> gene = alignedGeneSeq.getGene();
			StrainModifier strainModifier = gene.getTargetStrainModifier(targetStrain);
			MutationSet<VirusT> mutations = strainModifier.modifyMutationSet(
				gene, targetGene, alignedGeneSeq.getMutations());
			sequenceAllPosAAs.put(seqName, new HashMap<Integer, String>());

			for (int pos=firstAA; pos<=lastAA; pos++) {
				String aas;
				if (pos < seqFirstPos || pos > seqLastPos) {
					aas = ".";
				} else {
					Mutation<VirusT> mut = mutations.get(targetGene, pos);
					if (mut == null) {
						aas = "-";
					} else {
						aas = mut.getDisplayAAs().replace("-", "del");
					}
				}
				sequenceAllPosAAs.get(seqName).put(pos, aas);
			}
		}
	}


}
