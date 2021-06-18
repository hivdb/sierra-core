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
package edu.stanford.hivdb.genotypes;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.viruses.Virus;

public class Genotyper<VirusT extends Virus<VirusT>> {

	private final Integer treeFirstNA;
	private final Integer treeLastNA;
	private final int[][][] referenceMismatchTree;
	private final VirusT virusInstance;
	private transient Map<Integer, Set<String>> sdrmCodonMap;
	// private final int codonNAOffset = 2253;

	private static <VirusT extends Virus<VirusT>> int[][][] buildReferenceMismatchTree(
		List<GenotypeReference<VirusT>> references, int firstNA, int lastNA
	) {
		int seqLen = lastNA - firstNA + 1;
		int numRefs = references.size();
		// int[NAPosition][NA] rootNode = [...mismatchedRefs]
		int[][][] rootNode = new int[seqLen][4][numRefs + 1];
		int[][] rootNodePointer = new int[seqLen][4];
		for (int i = 0; i < seqLen; i ++) {
			// initialize rootNode; -1 indicates the end of matched refIdx
			for (int naIdx = 0; naIdx < 4; naIdx ++) {
				rootNode[i][naIdx][0] = -1;
			}
		}
		for (int refIdx = 0; refIdx < numRefs; refIdx ++) {
			GenotypeReference<VirusT> ref = references.get(refIdx);
			String sequence = ref.getSequence().toUpperCase();
			// build tree for each position
			for (int i = 0; i < seqLen; i ++) {
				char na = sequence.charAt(i);	
				// should search for mismatched refs but not matched
				int[] inverseNAIndice = getInverseNAIndice(na);
				for (int naIdx : inverseNAIndice) {
					int[] naNode = rootNode[i][naIdx];
					// add the newly found refIdx and shift end by 1
					naNode[rootNodePointer[i][naIdx] ++] = refIdx;
					naNode[rootNodePointer[i][naIdx]] = -1;
				}
			}
		}
		return rootNode;
	}

	private static int[] getNAIndice(char na) {
		List<Integer> naIndice = new ArrayList<>();
		String unNA = CodonUtils.expandAmbiguityNA(na);
		for (char sNA : unNA.toCharArray()) {
			switch(sNA) {
				case 'A':
					naIndice.add(0);
					break;
				case 'C':
					naIndice.add(1);
					break;
				case 'G':
					naIndice.add(2);
					break;
				case 'T':
					naIndice.add(3);
					break;
			}
		}
		return Ints.toArray(naIndice);
	}
	
	private static int[] getInverseNAIndice(char na) {
		int[] naIndice = getNAIndice(na);
		List<Integer> inverseNAIndice = new ArrayList<>();
		for (int i = 0; i < 4; i ++) {
			boolean found = false;
			for (int naIdx : naIndice) {
				if (i == naIdx) {
					found = true;
					break;
				}
			}
			if (!found) {
				inverseNAIndice.add(i);
			}
		}
		return Ints.toArray(inverseNAIndice);
	}
	
	private static void appendCodonDiscordance(
		int codonStartNAPos, String curCodon,
		Map<Integer, List<Integer>> discordanceListPerRef,
		Map<Integer, List<Integer>> curCodonDiscordancePerRef,
		Map<Integer, Set<String>> ignoredCodons) {
		Set<String> codons = ignoredCodons.get(codonStartNAPos);
		if (codons == null || !codons.contains(curCodon)) {
			// keep the result if the current codon is not a SDRM
			for (Map.Entry<Integer, List<Integer>> entry :
					curCodonDiscordancePerRef.entrySet()) {
				int refIdx = entry.getKey();
				if (!discordanceListPerRef.containsKey(refIdx)) {
					discordanceListPerRef.put(refIdx, new ArrayList<>());
				}
				List<Integer> discordanceList = discordanceListPerRef.get(refIdx);
				discordanceList.addAll(entry.getValue());
			}
		}
	}

	public Genotyper(VirusT virusIns) {
		virusInstance = virusIns;
		
		List<GenotypeReference<VirusT>> references = virusIns.getGenotypeReferences(); 
		
		Integer firstNA = null, lastNA = null;
		
		if (references.size() == 0) {
			treeFirstNA = 0;
			treeLastNA = 0;
			referenceMismatchTree = new int[0][0][0];
			return;
		}
		
		for (GenotypeReference<VirusT> ref : references) {
			if (firstNA == null) {
				firstNA = ref.getFirstNA();
			}
			if (lastNA == null) {
				lastNA = ref.getLastNA();
			}
			if (!firstNA.equals(ref.getFirstNA()) || !lastNA.equals(ref.getLastNA())) {
				throw new IllegalArgumentException(String.format(
					"Reference %s has a different NA boundary (%d - %d) " +
					"like other references (%d - %d).",
					ref.getAccession(),
					ref.getFirstNA(), ref.getLastNA(),
					firstNA, lastNA
				));
			}
		}
		treeFirstNA = firstNA;
		treeLastNA = lastNA;
		
		referenceMismatchTree = buildReferenceMismatchTree(references, firstNA, lastNA);
		
	}
	
	private Map<Integer, Set<String>> getSDRMCodonMap() {
		if (this.sdrmCodonMap == null) {
			Map<DrugClass<VirusT>, MutationSet<VirusT>> sdrmsMap = virusInstance.getSurveilDrugResistMutations();
			Map<Integer, Set<String>> sdrmCodonMap = new HashMap<>();
			for (MutationSet<VirusT> sdrms : sdrmsMap.values()) {
				for (Mutation<VirusT> mut : sdrms) {
					int absPos = mut.getGenePosition().getPositionInStrain();
					for (char aa : mut.getAAChars()) {
						if (aa == '-' || aa == '_') {
							continue;
						}
						sdrmCodonMap.put(
							absPos * 3 - 3,
							Collections.unmodifiableSet(
								Sets.newHashSet(CodonUtils.translateAA(aa))
							)
						);
					}
				}
			}
			this.sdrmCodonMap = Collections.unmodifiableMap(sdrmCodonMap);
		}
		return this.sdrmCodonMap;
	}
	
	/** compare given sequence with a set of references (provided by mismatchTree)
	 *
	 * @param sequence a string of DNA sequence
	 * @param seqFirstNA starting position of the given sequence
	 * @param seqLastNA ending position of the given sequence
	 * 
	 * @return A list of discordance for each given reference
	 * 
	 */
	protected Map<Integer, List<Integer>> compareWithSearchTree(
		String sequence, int seqFirstNA, int seqLastNA
	) {
		int maxFirstNA = Math.max(seqFirstNA, treeFirstNA);
		int minLastNA = Math.min(seqLastNA, treeLastNA);
		int treeOffset = maxFirstNA - treeFirstNA;
		int seqOffset = maxFirstNA - seqFirstNA;
		int compareLength = minLastNA - maxFirstNA + 1;
		Map<Integer, Set<String>> ignoredCodons = getSDRMCodonMap();
		Map<Integer, List<Integer>> discordanceListPerRef = new HashMap<>();
		Map<Integer, List<Integer>> curCodonDiscordancePerRef = new HashMap<>();
		StringBuffer curCodon = new StringBuffer();
		for (int i = 0; i < compareLength; i ++) {
			if ((maxFirstNA + i) % 3 == 0) {
				// to check if the current position is the beginning of a codon
				appendCodonDiscordance(
					/* codonStartNAPos */ maxFirstNA + i - 3,
					curCodon.toString(),
					discordanceListPerRef, curCodonDiscordancePerRef,
					ignoredCodons
				);
				curCodon.setLength(0);
				curCodonDiscordancePerRef.clear();
			}
			int[][] treeNAs = referenceMismatchTree[treeOffset + i];
			char seqNA = sequence.charAt(seqOffset + i);
			if (seqNA == '.') {
				// no need for further processing if seqNA == '.'
				curCodon.append(seqNA);
				continue;
			}
			int[] naIndice = getNAIndice(seqNA);
			Map<Integer, Integer> mismatchRefs = new HashMap<>();
			for (int naIdx : naIndice) {
				for (int mismatchRef : treeNAs[naIdx]) {
					if (mismatchRef < 0) {
						// the end
						break;
					}
					int mismatchCount = mismatchRefs.getOrDefault(mismatchRef, 0) + 1;
					mismatchRefs.put(mismatchRef, mismatchCount);
				}
			}
			for (Map.Entry<Integer, Integer> e : mismatchRefs.entrySet()) {
				if (e.getValue() < naIndice.length) {
					// only get counted as discordance when no unambiguous NA was matched
					continue;
				}
				int mismatchRef = e.getKey();
				if (!curCodonDiscordancePerRef.containsKey(mismatchRef)) {
					curCodonDiscordancePerRef.put(mismatchRef, new ArrayList<>());
				}
				curCodonDiscordancePerRef.get(mismatchRef).add(maxFirstNA + i);
			}
			curCodon.append(seqNA);
		}
		appendCodonDiscordance(
			/* codonStartNAPos */ minLastNA - 3,
			curCodon.toString(),
			discordanceListPerRef, curCodonDiscordancePerRef,
			ignoredCodons
		);
		return discordanceListPerRef;
	}
	
	public GenotypeResult<VirusT> compareAll(String sequence, int firstNA) {
		int lastNA = firstNA + sequence.length() - 1;
		return compareAll(sequence, firstNA, lastNA);
	}

	protected GenotypeResult<VirusT> compareAll(String sequence, int firstNA, int lastNA) {
		List<GenotypeReference<VirusT>> references = virusInstance.getGenotypeReferences();
		Map<Integer, List<Integer>> discordanceListPerRef = compareWithSearchTree(
			sequence, firstNA, lastNA);
		int numRefs = references.size();
		List<BoundGenotype<VirusT>> results = new ArrayList<>();
		for (int refIdx = 0; refIdx < numRefs; refIdx ++) {
			GenotypeReference<VirusT> ref = references.get(refIdx);
			results.add(ref.getBoundGenotype(
				sequence, firstNA, lastNA,
				discordanceListPerRef.getOrDefault(refIdx, new ArrayList<>())
			));
		}
		return new GenotypeResult<>(results);
	}

}
