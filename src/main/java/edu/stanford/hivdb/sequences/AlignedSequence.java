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
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.Collections;

import org.apache.commons.lang3.StringUtils;

import edu.stanford.hivdb.genotypes.BoundGenotype;
import edu.stanford.hivdb.genotypes.GenotypeResult;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.StrainModifier;
import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class AlignedSequence<VirusT extends Virus<VirusT>> {

	private final Strain<VirusT> strain;
	private final Sequence inputSequence;

	private List<ValidationResult> validationResults;
	private Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> alignedGeneSequenceMap;
	private Map<Gene<VirusT>, String> discardedGenes;
	private String concatenatedSequence;
	private MutationSet<VirusT> mutations;
	private MutationSet<VirusT> sdrms;
	private GenotypeResult<VirusT> genotypeResult;
	private Double mixturePcnt;
	private transient List<FrameShift<VirusT>> frameShifts;
	private final Boolean isReverseComplement;
	private final Boolean isEmpty;
	private transient Integer numMatchedNAs;
	private transient VirusT virusInstance;

	// private static final Logger LOGGER = LogManager.getLogger();

	public AlignedSequence(
			final Strain<VirusT> strain,
			final Sequence unalignedSequence,
			final Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> alignedGeneSequenceMap,
			final Map<Gene<VirusT>, String> discardedGenes,
			final boolean sequenceReversed) {
		inputSequence = unalignedSequence;
		this.strain = strain;
		this.virusInstance = strain.getVirusInstance();
		this.alignedGeneSequenceMap = alignedGeneSequenceMap;
		this.discardedGenes = discardedGenes;
		isReverseComplement = sequenceReversed;
		isEmpty = alignedGeneSequenceMap.isEmpty();
		numMatchedNAs = null;
	}

	public boolean isEmpty() {
		return isEmpty;
	}
	
	public Strain<VirusT> getStrain() {
		return strain;
	}

	public boolean isReverseComplement() {
		return isReverseComplement;
	}

	public Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> getAlignedGeneSequenceMap() {
		return alignedGeneSequenceMap;
	}

	/**
	 * Return an instance of AlignedGeneSeq for giving gene or
	 * {@code null} if the instance doesn't exist.
	 *
	 * @param gene the gene to be fetched
	 * @return corresponding instance of AlignedGeneSeq
	 */
	public AlignedGeneSeq<VirusT> getAlignedGeneSequence(Gene<VirusT> gene) {
		return alignedGeneSequenceMap.get(gene);
	}
	
	public AlignedGeneSeq<VirusT> getAlignedGeneSequence(String abstractGene) {
		Gene<VirusT> gene = strain.getGene(abstractGene);
		return alignedGeneSequenceMap.get(gene);
	}

	public List<AlignedGeneSeq<VirusT>> getAlignedGeneSequences() {
		return Collections.unmodifiableList(
			new ArrayList<>(alignedGeneSequenceMap.values()));
	}

	public Sequence getInputSequence() {
		return inputSequence;
	}

	public List<Gene<VirusT>> getAvailableGenes() {
		return new ArrayList<>(alignedGeneSequenceMap.keySet());
	}

	public List<ValidationResult> getValidationResults() {
		if (validationResults == null) {
			validationResults = (
				strain.getVirusInstance()
				.validateSequence(this)
			);
		}
		return validationResults;
	}

	public Map<Gene<VirusT>, String> getDiscardedGenes() {
		return discardedGenes;
	}

	public String getConcatenatedSeq() {
		if (concatenatedSequence == null) {
			concatenatedSequence = concatAlignments(true, strain);
		}
		return concatenatedSequence;
	}
	
	@SuppressWarnings("unchecked")
	private <T extends Virus<T>> String concatAlignments(boolean trimResult, Strain<T> targetStrain) {
		String geneSeqNAs;
		AlignedGeneSeq<VirusT> geneSeq;
		StringBuilder concatSeq = new StringBuilder();
		if (targetStrain == null) {
			targetStrain = (Strain<T>) strain;
		}
		String targetStrainText = targetStrain.getName();

		for (Gene<VirusT> srcGene : strain.getGenes()) {
			Gene<T> targetGene = targetStrain.getGene(srcGene.getAbstractGene());
			geneSeq = alignedGeneSequenceMap.get(srcGene);
			if (geneSeq == null) {
				concatSeq.append(StringUtils.repeat(StrainModifier.WILDCARD, targetGene.getNASize()));
				continue;
			}
			StrainModifier strainModifier = srcGene.getTargetStrainModifier(targetStrainText);
			geneSeqNAs = geneSeq.getAlignedNAs();
			geneSeqNAs = strainModifier.modifyNASeq(
				srcGene, geneSeqNAs, geneSeq.getFirstAA(), geneSeq.getLastAA());
			concatSeq.append(geneSeqNAs);
		}
		if (trimResult) {
			Matcher matcher = StrainModifier.TRIM_WILDCARD_PATTERN.matcher(concatSeq.toString());
			return matcher.replaceAll("");
		} else {
			return concatSeq.toString();
		}
	}

	public MutationSet<VirusT> getMutations() {
		if (mutations == null) {
			mutations = new MutationSet<>();
			alignedGeneSequenceMap.values()
				.stream()
				.forEach(geneSeq -> {
					mutations = mutations.mergesWith(geneSeq.getMutations());
				});
		}
		return mutations;
	}

	public MutationSet<VirusT> getSdrms() {
		if (sdrms == null) {
			sdrms = getMutations().getSDRMs();
		}
		return sdrms;
	}

	public GenotypeResult<VirusT> getGenotypeResult() {
		if (!isEmpty && genotypeResult == null) {
			Strain<VirusT> targetStrain = virusInstance.getMainStrain();
			// Tweak the alignment first
			String compatConcatSeq = concatAlignments(false, targetStrain);
			int absFirstNA = targetStrain.getAbsoluteFirstNA();
			genotypeResult = (
				virusInstance
				.getGenotyper()
				.compareAll(compatConcatSeq, absFirstNA)
			);
		}
		return genotypeResult;
	}
	
	public GenotypeResult<VirusT> getSubtypeResult() {
		return getGenotypeResult();
	}

	public BoundGenotype<VirusT> getBestMatchingGenotype() {
		if (isEmpty) {
			return null;
		}
		return getGenotypeResult().getBestMatch();
	}
	
	public BoundGenotype<VirusT> getBestMatchingSubtype() {
		return getBestMatchingGenotype();
	}

	public String getGenotypeText() {
		if (isEmpty) {
			return "NA";
		}
		return getBestMatchingGenotype().getDisplay();
	}
	
	public String getSubtypeText() {
		return getGenotypeText();
	}

	public double getMixturePcnt() {
		if (mixturePcnt == null) {
			mixturePcnt = SeqUtils.mixturePcnt(getConcatenatedSeq());
		}
		return mixturePcnt;
	}

	public List<FrameShift<VirusT>> getFrameShifts() {
		if (frameShifts == null) {
			frameShifts = new ArrayList<>();
			for (AlignedGeneSeq<VirusT> seq: alignedGeneSequenceMap.values()) {
				frameShifts.addAll(seq.getFrameShifts());
			}
		}
		return frameShifts;
	}
	
	public int getNumMatchedNAs() {
		if (numMatchedNAs == null) {
			numMatchedNAs = 0;
			for (Gene<VirusT> gene : alignedGeneSequenceMap.keySet()) {
				numMatchedNAs += gene.getNASize();
				AlignedGeneSeq<VirusT> geneSeq = alignedGeneSequenceMap.get(gene);
				numMatchedNAs -=
					geneSeq.getFirstAA() * 3 - 3 + // left missing NAs
					geneSeq.getNumDiscordantNAs() +
					gene.getNASize() - geneSeq.getLastAA() * 3; // right missing NAs
			}
		}
		return numMatchedNAs;
	}
}
