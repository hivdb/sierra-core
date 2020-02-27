/*

    Copyright (C) 2019 Stanford HIVDB team

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

package edu.stanford.hivdb.mutations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import edu.stanford.hivdb.sequences.AlignedSite;
import edu.stanford.hivdb.utilities.AssertUtils;
import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

public class StrainModifier {
	
	public static final Character WILDCARD = '.';
	public static final Pattern TRIM_WILDCARD_PATTERN = Pattern.compile("(^\\.+|\\.+$)");
	
	public static enum CIGARFlag {
		M("Match"),
		I("Insertion"),
		D("Deletion");
		// S, H, =, X, and N are not supported 

		private String fullName; 
		
		CIGARFlag(String fullName) {
			this.fullName = fullName;
		}
		
		public String getFullName() {
			return fullName;
		}
	}
	
	public static class PosModifier {
		private Integer targetPos;
		private Integer size;
		private CIGARFlag flag;
		
		public PosModifier(Integer targetPos, Integer size, CIGARFlag flag) {
			this.targetPos = targetPos;
			this.size = size;
			this.flag = flag;
		}
		
		public Integer getTargetPos() {
			return targetPos;
		}
		
		public Integer getSize() {
			return size;
		}
		
		public CIGARFlag getFlag() {
			return flag;
		}
	}

	private static final Pattern CIGAR_PATTERN = Pattern.compile("(\\d+)([MID])");
	
	private String targetStrain;
	private final String cigar;

	private transient List<Pair<Integer, CIGARFlag>> cigarList;
	private transient Map<Integer, PosModifier> posModifierMap;
	
	public StrainModifier(
		String targetStrain, String cigar
	) {
		this.targetStrain = targetStrain;
		this.cigar = cigar;
	}
	
	public String getTargetStrain() {
		return targetStrain;
	}
	
	public List<Pair<Integer, CIGARFlag>> getCIGAR() {
		if (this.cigarList == null) {
			List<Pair<Integer, CIGARFlag>> cigarList = new ArrayList<>();
			Matcher m = CIGAR_PATTERN.matcher(cigar);
			while (m.find()) {
				Integer offset = Integer.parseInt(m.group(1));
				CIGARFlag flag = CIGARFlag.valueOf(m.group(2));
				cigarList.add(Pair.of(offset, flag));
			}
			this.cigarList = Collections.unmodifiableList(cigarList);
		}
		return this.cigarList;
	}
	
	public Map<Integer, PosModifier> getPosModifiers() {
		if (this.posModifierMap == null) {
			List<Pair<Integer, CIGARFlag>> cigarList = getCIGAR();
			int srcPos = 1;
			int tgtPos = 1;
			Map<Integer, PosModifier> posMods = new LinkedHashMap<>();
			for (Pair<Integer, CIGARFlag> cigarMod : cigarList) {
				int srcPosEnd;
				int size = cigarMod.getLeft();
				CIGARFlag flag = cigarMod.getRight();
				switch (flag) {
					case M:
						srcPosEnd = srcPos + size;
						while (srcPos < srcPosEnd) {
							posMods.put(srcPos, new PosModifier(tgtPos, 1, flag));
							srcPos ++;
							tgtPos ++;
						}
						break;
					case I:
						// use previous M to store insertion info
						posMods.put(srcPos - 1, new PosModifier(tgtPos - 1, size, flag));
						tgtPos += size;
						break;
					case D:
						srcPosEnd = srcPos + size;
						while (srcPos < srcPosEnd) {
							posMods.put(srcPos, new PosModifier(tgtPos - 1, 1, flag));
							srcPos ++;
						}
						break;
				}
			}
			this.posModifierMap = Collections.unmodifiableMap(posMods);
		}
		return this.posModifierMap;
	}

	// TODO: do we really need modify alignedSites?
	public List<AlignedSite> modifyAlignedSites(List<AlignedSite> srcAlignedSites) {
		Map<Integer, PosModifier> posModifiers = getPosModifiers();
		List<AlignedSite> targetAlignedSites = new ArrayList<>();
		for (AlignedSite site : srcAlignedSites) {
			int srcPos = site.getPosAA();
			PosModifier posMod = posModifiers.get(srcPos);
			int tgtPos = posMod.getTargetPos();
			int posNA = site.getPosNA();
			int lenNA = site.getLengthNA();
			switch(posMod.getFlag()) {
				case M:
					targetAlignedSites.add(new AlignedSite(tgtPos, posNA, lenNA));
					break;
				case I:
					targetAlignedSites.add(new AlignedSite(tgtPos, posNA, lenNA));
					int insSize = posMod.getSize();
					for (int i = 0; i < insSize; i ++) {
						targetAlignedSites.add(new AlignedSite(tgtPos + i + 1, posNA, 0));
					}
					break;
				case D:
					int lastTgtSiteIdx = targetAlignedSites.size() - 1;
					AlignedSite lastTgtSite = targetAlignedSites.get(lastTgtSiteIdx);
					AssertUtils.isTrue(
						lastTgtSite.getPosAA() == tgtPos,
						"Argument alignedSites is not continuous around position %d",
						tgtPos
					);
					lastTgtSite = new AlignedSite(
						lastTgtSite.getPosAA(),
						lastTgtSite.getPosNA(),
						lastTgtSite.getLengthNA() + 3);
					targetAlignedSites.set(lastTgtSiteIdx, lastTgtSite);
					break;
			}
		}
		return targetAlignedSites;
	}
	
	private <VirusS extends Virus<VirusS>, VirusT extends Virus<VirusT>>
	List<Mutation<VirusT>> modifyMutation(Mutation<VirusS> srcMutation, Gene<VirusT> targetGene, PosModifier posMod) {
		int tgtPos = posMod.getTargetPos();
		List<Mutation<VirusT>> targetMutations = new ArrayList<>();

		if (srcMutation instanceof CodonMutation) {
			String aas = srcMutation.getAAs();
			String triplet = srcMutation.getTriplet();
			String insertedNAs = srcMutation.getInsertedNAs();
			int maxDisplayAAs = ((CodonMutation<VirusS>) srcMutation).getMaxDisplayAAs();
			switch(posMod.getFlag()) {
				case M:
					targetMutations.add(new CodonMutation<>(
						targetGene, tgtPos, aas, triplet, insertedNAs, maxDisplayAAs));
					break;
				case I:
					targetMutations.add(new CodonMutation<>(
						targetGene, tgtPos, aas, triplet, insertedNAs, maxDisplayAAs));
					int insSize = posMod.getSize();
					for (int i = 0; i < insSize; i ++) {
						targetMutations.add(new CodonMutation<>(
						targetGene, tgtPos + i + 1, "-", "---", "", maxDisplayAAs));
					}
					break;
				case D:
					int lastTgtMutIdx = targetMutations.size() - 1;
					Mutation<VirusT> lastTgtMut = targetMutations.get(lastTgtMutIdx);
					if (lastTgtMut.getPosition() == tgtPos) {
						insertedNAs = lastTgtMut.getInsertedNAs() + triplet + insertedNAs;
						aas = lastTgtMut.getAAs() + '_' + CodonUtils.simpleTranslate(insertedNAs);
						targetMutations.set(lastTgtMutIdx, new CodonMutation<>(
							targetGene, tgtPos, aas, lastTgtMut.getTriplet(), insertedNAs, maxDisplayAAs
						));
					}
					break;
			}

		}
		else if (srcMutation instanceof MultiCodonsMutation) {
			throw new NotImplementedException("Modifying MultiCodonsMutation is not yet supported");
		}
		else { // default to AAMutation
			Set<Character> aas = srcMutation.getAAChars();
			switch(posMod.getFlag()) {
				case M:
					targetMutations.add(new AAMutation<>(targetGene, tgtPos, aas));
					break;
				case I:
					targetMutations.add(new AAMutation<>(targetGene, tgtPos, aas));
					int insSize = posMod.getSize();
					for (int i = 0; i < insSize; i ++) {
						targetMutations.add(new AAMutation<>(targetGene, tgtPos + i + 1, '-'));
					}
					break;
				case D:
					int lastTgtMutIdx = targetMutations.size() - 1;
					Mutation<VirusT> lastTgtMut = targetMutations.get(lastTgtMutIdx);
					if (lastTgtMut.getPosition() == tgtPos) {
						targetMutations.set(lastTgtMutIdx, new AAMutation<>(targetGene, tgtPos, '_'));
					}
					break;
			}
		}
		return targetMutations;
	}

	public <VirusS extends Virus<VirusS>, VirusT extends Virus<VirusT>>
	MutationSet<VirusT> modifyMutationSet(Gene<VirusS> srcGene, Gene<VirusT> targetGene, MutationSet<VirusS> srcMutations) {
		AssertUtils.notNull(srcGene, "Argument srcGene should not be null");
		AssertUtils.notNull(targetGene, "Argument targetGene should not be null");
		AssertUtils.notNull(srcMutations, "Argument srcMutations should not be null");
		AssertUtils.isTrue(
			srcGene.getAbstractGene().equals(targetGene.getAbstractGene()),
			"Arguments srcGene is '%s', however targetGene is '%s'",
			srcGene.getAbstractGene(), targetGene.getAbstractGene());
		Map<Integer, PosModifier> posModifiers = getPosModifiers();
		List<Mutation<VirusT>> targetMutations = new ArrayList<>();
		for (Mutation<VirusS> srcMut : srcMutations) {
			if (srcMut.getGene() != srcGene) {
				continue;
			}
			int srcPos = srcMut.getPosition();
			PosModifier posMod = posModifiers.get(srcPos);
			targetMutations.addAll(modifyMutation(srcMut, targetGene, posMod));
		}
		return new MutationSet<>(targetMutations);
	}

	/**
	 * Re-align amino acid sequence with target strain
	 *
	 * This function allows to re-align the sequences from one strain
	 * into another strain's numbering system using predefined CIGAR
	 * string associated with given targetStrain.
	 *
	 * @param srcGene
	 * @param aaseq
	 * @param firstAA
	 * @param lastAA
	 * @return String of adjusted AA alignment in full gene length
	 */
	public String modifyAASeq(
		Gene<?> srcGene, String aaseq, int firstAA, int lastAA
	) {
		int numPrefixAAs = firstAA - 1;
		int numSuffixAAs = srcGene.getAASize() - aaseq.length() - numPrefixAAs;
		aaseq =
			StringUtils.repeat(WILDCARD, numPrefixAAs) +
			aaseq +
			StringUtils.repeat(WILDCARD, numSuffixAAs);
		StringBuilder targetAASeq = new StringBuilder();
		for (Pair<Integer, CIGARFlag> mod : getCIGAR()) {
			Integer offset = mod.getLeft();
			CIGARFlag flag = mod.getRight();
			switch (flag) {
				case M:
					targetAASeq.append(aaseq.substring(0, offset));
					aaseq = aaseq.substring(offset);
					break;
				case I:
					targetAASeq.append(StringUtils.repeat(WILDCARD, offset));
					break;
				case D:
					aaseq = aaseq.substring(offset);
					break;
			}
		}
		return targetAASeq.toString();
	}

	/**
	 * Re-align nucleotide sequence with target strain
	 *
	 * This function allows to re-align the sequences from one strain
	 * into another strain's numbering system using predefined CIGAR
	 * string associated with given targetStrain.
	 *
	 * @param srcGene
	 * @param naseq
	 * @param firstAA
	 * @param lastAA
	 * @return String of adjusted AA alignment in full gene length
	 */
	public String modifyNASeq(
		Gene<?> srcGene, String naseq, int firstAA, int lastAA
	) {
		int numPrefixNAs = (firstAA - 1) * 3;
		int numSuffixNAs = srcGene.getNASize() - naseq.length() - numPrefixNAs;
		naseq =
			StringUtils.repeat(WILDCARD, numPrefixNAs) +
			naseq +
			StringUtils.repeat(WILDCARD, numSuffixNAs);
		StringBuilder targetNASeq = new StringBuilder();
		for (Pair<Integer, CIGARFlag> mod : getCIGAR()) {
			Integer offset = mod.getLeft();
			CIGARFlag flag = mod.getRight();
			switch (flag) {
				case M:
					targetNASeq.append(naseq.substring(0, offset * 3));
					naseq = naseq.substring(offset);
					break;
				case I:
					targetNASeq.append(StringUtils.repeat(WILDCARD, offset * 3));
					break;
				case D:
					naseq = naseq.substring(offset * 3);
					break;
			}
		}
		return targetNASeq.toString();
	}

}
