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
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Pair;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.StrainModifier;
import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.utilities.PrettyPairwise;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;


/**
 * Result object of data from {@link edu.stanford.hivdb.sequences.PostAlignAligner}.
 *
 */
public class AlignedGeneSeq<VirusT extends Virus<VirusT>> implements WithGene<VirusT>, WithSequenceStat<VirusT> {

	// Variables assigned by CLapAlign
	private final Gene<VirusT> gene;
	private final Sequence sequence;
	private final int firstAA;
	private final int lastAA;
	private final int firstNA;
	private final int lastNA;
	private Long aaSize;
	private Double matchPcnt;
	private transient List<AlignedSite> alignedSites;
	private transient PrettyPairwise<VirusT> prettyPairwise;
	private transient GeneRegions<VirusT> unseqRegions;

	// Variables assigned by NucAminoAligner and changed by the methods in this class
	private String alignedNAs;

	// This information is generated by processSequence.
	private MutationSet<VirusT> mutations;
	private MutationSet<VirusT> sequencedMutations;
	private List<FrameShift<VirusT>> frameShifts = new ArrayList<>();
	private transient MutationSet<VirusT> unusualMutations;
	private transient MutationSet<VirusT> sdrms;
	private transient Map<DrugClass<VirusT>, MutationSet<VirusT>> nonDrmTsms;

	// This may help for debugging and may be preferable to serialize than mutations
	@SuppressWarnings("unused")
	private String mutationListString;

	private transient Map<MutationType<VirusT>, MutationSet<VirusT>> mutationsGroupingByMutType;
	
	/**
	 *
	 * @param sequence			Aligned gene sequence
	 * @param gene				Gene related to this sequenc
	 * @param firstAA			First amino acid position
	 * @param lastAA			Last amino acid position
	 * @param firstNA			First nucleic acid position
	 * @param lastNA			Last nucleic acid position
	 * @param alignedSites		Aligned sequence sites
	 * @param mutations			Mutations
	 * @param frameShifts		List&lt;FrameShift&gt;
	 * @param sequenceReversed	If sequence should be reversed
	 */
	public AlignedGeneSeq(
			Sequence sequence, Gene<VirusT> gene,
			final int firstAA, final int lastAA,
			final int firstNA, final int lastNA,
			List<AlignedSite> alignedSites,
			Collection<Mutation<VirusT>> mutations,
			List<FrameShift<VirusT>> frameShifts) {
		this.sequence = sequence;
		this.gene = gene;
		this.firstAA = firstAA;
		this.lastAA = lastAA;
		this.firstNA = firstNA;
		this.lastNA = lastNA;
		this.aaSize = -1L;
		this.matchPcnt = -1.;

		alignedSites = alignedSites.stream()
			.filter(m -> {
				int posAA = m.getPosAA();
				return posAA >= firstAA && posAA <= lastAA;
			})
			.collect(Collectors.toList());
		mutations = mutations.stream()
			.filter(m -> {
				int posAA = m.getPosition();
				return posAA >= firstAA && posAA <= lastAA;
			})
			.collect(Collectors.toList());
		frameShifts = frameShifts.stream()
			.filter(fs -> {
				int posAA = fs.getPosition();
				return posAA >= firstAA && posAA <= lastAA;
			})
			.collect(Collectors.toList());


		this.alignedSites = Collections.unmodifiableList(alignedSites);
		this.mutations = new MutationSet<>(mutations);
		this.frameShifts = Collections.unmodifiableList(frameShifts);
		mutationListString = getMutationListString();
	}


	@Override
	public Gene<VirusT> getGene() { return gene; }

	@Override
	public Strain<VirusT> getStrain() {return gene.getStrain(); }

	@Override
	public String getAbstractGene() { return gene.getAbstractGene(); }

	public Sequence getSequence() {	return sequence; }
	public List<AlignedSite> getAlignedSites() { return alignedSites; }

	/**
	 * Retrieve amino acids size of alignment.
	 * @return integer
	 */
	public long getSize() {
		if (aaSize == -1L) {
			aaSize = (long) lastAA - firstAA + 1;
			GeneRegions<VirusT> unseqRegions = getUnsequencedRegions();
			aaSize -= unseqRegions.size();
		}
		return aaSize;
	}
	
	protected Map<Integer, String> getCodonLookup() {
		String naSeq = sequence.getSequence();

		Map<Integer, String> allCodons = new LinkedHashMap<>();

		for (AlignedSite site : alignedSites) {
			List<Integer> posNAs = site.getPosNAs();
			StringBuilder codon = new StringBuilder();
			for (Integer posNA : posNAs) {
				if (posNA == null) {
					codon.append('-');
				}
				else {
					codon.append(naSeq.charAt(posNA - 1));
				}
			}
			allCodons.put(site.getPosAA(), codon.toString());
		}
		
		return allCodons;
	}
	
	public List<String> getAlignedCodons(boolean skipIns) {
		int aaSize = getGene().getAASize();
		List<String> codons = new ArrayList<>();
		Map<Integer, String> codonLookup = getCodonLookup();
		for (int pos = 1; pos <= aaSize; pos ++) {
			String posCodon = codonLookup.getOrDefault(pos, "NNN");
			if (skipIns) {
				// skip insertions by only using first 3 chars
				codons.add(posCodon.substring(0, 3));
			}
			else {
				codons.add(posCodon.replace("-", ""));
			}
		}
		return codons;
	}
	
	public String getAlignedNAs() {
		if (this.alignedNAs == null) {
			List<String> codons = getAlignedCodons(true);
			codons = codons.subList(firstAA - 1, lastAA);
			this.alignedNAs = String.join("", codons);
		}
		return this.alignedNAs;
	} // Need

	public String getAlignedNAsNoTrim() {
		if (this.alignedNAs == null) {
			List<String> codons = getAlignedCodons(true);
			this.alignedNAs = String.join("", codons);
		}
		return this.alignedNAs;
	} // Need
	
	public String getAlignedAAs() {
		return CodonUtils.simpleTranslate(
			this.getAlignedNAs(), firstAA, gene.getRefSequence());
	}

	public String getAdjustedAlignedNAs() {
		String targetStrain = gene.getVirusInstance().getMainStrain().getName();
		StrainModifier strainModifier = gene.getTargetStrainModifier(targetStrain);
		return strainModifier.modifyNASeq(gene, getAlignedNAs(), firstAA, lastAA);
	}

	public String getAdjustedAlignedNAs(String targetStrain) {
		StrainModifier strainModifier = gene.getTargetStrainModifier(targetStrain);
		return strainModifier.modifyNASeq(gene, getAlignedNAs(), firstAA, lastAA);
	}
	
	public String getAdjustedAlignedNAsNoTrim(String targetStrain) {
		StrainModifier strainModifier = gene.getTargetStrainModifier(targetStrain);
		return strainModifier.modifyNASeq(gene, getAlignedNAsNoTrim(), 1, gene.getAASize());
	}

	public String getAdjustedAlignedAAs() {
		String targetStrain = gene.getVirusInstance().getMainStrain().getName();
		StrainModifier strainModifier = gene.getTargetStrainModifier(targetStrain);
		return strainModifier.modifyAASeq(gene, getAlignedAAs(), firstAA, lastAA);
	}

	public String getAdjustedAlignedAAs(String targetStrain) {
		StrainModifier strainModifier = gene.getTargetStrainModifier(targetStrain);
		return strainModifier.modifyAASeq(gene, getAlignedAAs(), firstAA, lastAA);
	}

	public int getFirstNA() { return firstNA; }
	public int getLastNA() { return lastNA; }
	public int getFirstAA() { return firstAA;}  // Need
	public int getLastAA() { return lastAA; } // Need

	protected int getNumDiscordantNAs() {
		int numDiscordantNAs = 0;
			GeneRegions<VirusT> unseqRegions = getUnsequencedRegions();
		for (Mutation<VirusT> mut : mutations) {
			if (mut.isUnsequenced(unseqRegions)) {
				// NNN doesn't count
				continue;
			}
			else if (mut.isDeletion()) {
				numDiscordantNAs += 3;
			} else {
				numDiscordantNAs += 3; //CodonTranslation.getMinimalNAChanges(mut.getTriplet(), mut.getConsensus());
			}
		}
		return numDiscordantNAs;

	}

	public Double getMatchPcnt() {
		if (matchPcnt == -1) {
			GeneRegions<VirusT> unseqRegions = getUnsequencedRegions();
			int numNAs = (lastAA - firstAA + 1) * 3;
			for (Mutation<VirusT> mut : mutations) {
				if (mut.isUnsequenced(unseqRegions)) {
					// unsequenced region doesn't count
					numNAs -= 3;
				}
				else if (mut.isDeletion()) {
					numNAs += 3;
				}
			}
			matchPcnt = 100 - 100 * Double.valueOf(getNumDiscordantNAs()) / Double.valueOf(numNAs);
		}
		return matchPcnt;
	}

	public MutationSet<VirusT> getMutations() { return mutations; }  // Need
	
	@Override
	public MutationSet<VirusT> getSequencedMutations() {
		if (sequencedMutations == null) {
			GeneRegions<VirusT> unseqRegions = getUnsequencedRegions();
			sequencedMutations = mutations.filterByNoSplit(mut -> !mut.isUnsequenced(unseqRegions));
		}
		return sequencedMutations;
	}  // Need

	@Override
	public List<FrameShift<VirusT>> getFrameShifts() { return frameShifts; } // Need

	public PrettyPairwise<VirusT> getPrettyPairwise() {
		if (prettyPairwise == null) {
			prettyPairwise = new PrettyPairwise<>(
				gene, getAlignedNAs(), firstAA,
				getMutations(), Collections.emptyList());
		}
		return prettyPairwise;
	}

	public String getMutationListString() { return getSequencedMutations().join(); } // Need

	public MutationSet<VirusT> getInsertions() { return getSequencedMutations().getInsertions(); }
	public MutationSet<VirusT> getDeletions() { return getSequencedMutations().getDeletions(); }
	public MutationSet<VirusT> getStopCodons() { return getSequencedMutations().getStopCodons(); }
	public MutationSet<VirusT> getHighlyAmbiguousCodons() {
		return getSequencedMutations().getAmbiguousCodons();
	}

	public Map<MutationType<VirusT>, MutationSet<VirusT>> groupMutationsByMutType() {
		if (mutationsGroupingByMutType == null) {
			mutationsGroupingByMutType = Collections.unmodifiableMap(
				getSequencedMutations().groupByMutType(gene));
		}
		return mutationsGroupingByMutType;
	}

	public MutationSet<VirusT> getMutationsByMutType(MutationType<VirusT> mutType) {
		return groupMutationsByMutType().getOrDefault(mutType, new MutationSet<>());
	}

	public MutationSet<VirusT> getUnusualMutations() {
		if (unusualMutations == null) {
			unusualMutations = getSequencedMutations().getUnusualMutations();
		}
		return unusualMutations;
	}
	
	public MutationSet<VirusT> getUnusualMutationsAtDrugResistancePositions() {
		return getUnusualMutations().getAtDRPMutations();
	}

	public MutationSet<VirusT> getSdrms() {
		if (sdrms == null) {
			sdrms = getSequencedMutations().getSDRMs();
		}
		return sdrms;
	}

	public Map<DrugClass<VirusT>, MutationSet<VirusT>> getNonDrmTsms() {
		if (nonDrmTsms == null) {
			MutationSet<VirusT> mutations = getSequencedMutations();
			nonDrmTsms = (
				gene.getVirusInstance()
				.getRxSelectedMutations()
				.entrySet()
				.stream()
				.filter(e -> (
					e.getKey().getAbstractGene()
					.equals(getAbstractGene())))
				.map(e -> Pair.of(
					e.getKey(),
					e.getValue().intersectsWith(mutations)
				))
				.collect(Collectors.toMap(
					e -> e.getKey(),
					e -> e.getValue()
				))
			);
		}
		return nonDrmTsms;
	}
	
	public GeneRegions<VirusT> getUnsequencedRegions() {
		if (unseqRegions == null) {
			unseqRegions = GeneRegions.newUnsequencedRegions(gene, firstAA, lastAA, mutations);
		}
		return unseqRegions;
	}
	
	public GeneRegions<VirusT> getUnsequencedRegions(Gene<VirusT> gene) {
		if (gene == this.gene) {
			return getUnsequencedRegions();
		}
		return null;
	}

}
