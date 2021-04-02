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

import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.utilities.TSV;
import edu.stanford.hivdb.viruses.Virus;

public class TabularAlgorithmsComparison<VirusT extends Virus<VirusT>> {

	private List<String> headers;
	private List<AlignedSequence<VirusT>> alignedSeqs;
	private List<DrugResistanceAlgorithm<VirusT>> algorithms;

	public TabularAlgorithmsComparison(
		List<AlignedSequence<VirusT>> alignedSeqs,
		List<DrugResistanceAlgorithm<VirusT>> algorithms
	) {
		this.alignedSeqs = alignedSeqs;
		this.algorithms = algorithms;
		headers = new ArrayList<>();
		headers.add("sequenceName");
		headers.add("gene");
		headers.add("drugName");
		for (DrugResistanceAlgorithm<VirusT> alg : algorithms) {
			headers.add(alg.getName() + ".drugLevel");
		}
	}

	@Override
	public String toString() {
		List<Map<String, String>> allRows = new ArrayList<>();
		for (AlignedSequence<VirusT> alignedSeq : alignedSeqs) {
			Map<Drug<VirusT>, Map<String, String>> rows = new LinkedHashMap<>();
			MutationSet<VirusT> allMuts = alignedSeq.getMutations();
			AlgorithmComparison<VirusT> algCmp = new AlgorithmComparison<>(allMuts, algorithms);
			List<ASIDrugSusc<VirusT>> cmpResults = algCmp.getComparisonResults();
			// an assumption was made here that algorithms are in the same order as the input
			for (ASIDrugSusc<VirusT> result : cmpResults) {
				Drug<VirusT> drug = result.getDrug();
				if (!rows.containsKey(drug)) {
					Map<String, String> row = new HashMap<>();
					row.put("sequenceName", alignedSeq.getInputSequence().getHeader());
					row.put("gene", drug.getDrugClass().getAbstractGene());
					row.put("drugName", drug.getDisplayAbbr());
					rows.put(drug, row);
				}
				Map<String, String> row = rows.get(drug);
				SIREnum sir = result.getSIR();
				DrugResistanceAlgorithm<VirusT> alg = result.getAlgorithmObj();
				row.put(alg.getName() + ".drugLevel", sir == null ? "-" : sir.toString());
			}
			allRows.addAll(rows.values());
		}
		return TSV.dumpMaps(headers, allRows, "-");
	}

}
