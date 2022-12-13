/*

    Copyright (C) 2022 Stanford HIVDB team

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

package edu.stanford.hivdb.viruses;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.Lists;

import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.sequences.SequenceValidator;
import edu.stanford.hivdb.sequences.GeneRegions;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.utilities.MyStringUtils;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;

public class DefaultSequenceValidator<VirusT extends Virus<VirusT>> implements SequenceValidator<VirusT> {

	protected DefaultSequenceValidator() {}

	@Override
	public List<ValidationResult> validate(AlignedSequence<VirusT> alignedSequence, Collection<String> includeGenes) {
		List<ValidationResult> results = new ArrayList<>();
		results.addAll(validateNotEmpty(alignedSequence, includeGenes));
		if (results.size() > 0) {
			return results;
		}
		results.addAll(validateReverseComplement(alignedSequence));
		results.addAll(validateNoMissingPositions(alignedSequence, includeGenes));
		results.addAll(validateLongGap(alignedSequence, includeGenes));
		results.addAll(validateNAs(alignedSequence));
		results.addAll(validateGaps(alignedSequence, includeGenes));
		results.addAll(validateNoStopCodons(alignedSequence, includeGenes));
		return results;
	}

	protected static <VirusT extends Virus<VirusT>> List<ValidationResult> validateNotEmpty(
		AlignedSequence<VirusT> alignedSequence,
		Collection<String> includeGenes
	) {
		VirusT virusIns = alignedSequence.getStrain().getVirusInstance();
		boolean isNotEmpty = alignedSequence.getAvailableGenes().stream()
			.anyMatch(gene -> includeGenes.contains(gene.getAbstractGene()));
		if (!isNotEmpty) {
			List<String> includeGeneDisplays = includeGenes.stream()
				.map(geneName -> virusIns.getGeneDisplay(geneName))
				.collect(Collectors.toList());
			return Lists.newArrayList(
				DefaultValidationMessage.NoGeneFound.format(
					MyStringUtils.andListFormat(includeGeneDisplays)
				)
			);
		}
		return Collections.emptyList();
	}

	protected static <VirusT extends Virus<VirusT>>List<ValidationResult> validateReverseComplement(AlignedSequence<?> alignedSequence) {
		if (alignedSequence.isReverseComplement()) {
			return Lists.newArrayList(DefaultValidationMessage.FASTAReverseComplement.format());
		}
		return Collections.emptyList();
	}

	protected static <VirusT extends Virus<VirusT>> List<ValidationResult> validateNoMissingPositions(
		final Set<GenePosition<VirusT>> needGenePositions,
		final Set<GenePosition<VirusT>> needDRGenePositions,
		final Set<GenePosition<VirusT>> availableGenePositions
	) {
		List<ValidationResult> results = new ArrayList<>();

		List<GenePosition<VirusT>> missingPositions = needGenePositions.stream()
				.filter(gp -> !availableGenePositions.contains(gp))
				.collect(Collectors.toList());
		long numMissingPositions = missingPositions.size();

		List<GenePosition<VirusT>> missingDRPs = needDRGenePositions.stream()
				.filter(gp -> !availableGenePositions.contains(gp))
				.collect(Collectors.toList());
		long numMissingDRPs = missingDRPs.size();
		
		String textMissingPositions = StringUtils.join(
			GeneRegions.newListOfGeneRegions(missingPositions),
			"; "
		);

		String textMissingDRPs = StringUtils.join(
			GeneRegions.newListOfGeneRegions(missingDRPs),
			"; "
		);
		
		if (numMissingDRPs > 1) {
			results.add(DefaultValidationMessage.MultiplePositionsMissingWithMultipleDRPs.formatWithLevel(
				numMissingDRPs > 5 ? ValidationLevel.SEVERE_WARNING : ValidationLevel.WARNING,
				numMissingPositions,
				textMissingPositions,
				numMissingDRPs,
				textMissingDRPs
			));
		}
		else if (numMissingDRPs > 0 && numMissingPositions > 1) {
			results.add(DefaultValidationMessage.MultiplePositionsMissingWithSingleDRP.formatWithLevel(
				ValidationLevel.WARNING,
				numMissingPositions,
				textMissingPositions,
				textMissingDRPs
			));
		}
		else if (numMissingDRPs > 0) {
			results.add(DefaultValidationMessage.SingleDRPMissing.formatWithLevel(
				ValidationLevel.NOTE,
				textMissingDRPs
			));
		}
		else if (numMissingPositions > 1) {
			results.add(DefaultValidationMessage.MultiplePositionsMissingWithoutDRP.formatWithLevel(
				ValidationLevel.WARNING,
				numMissingPositions,
				textMissingPositions
			));
		}
		else if (numMissingPositions > 0) {
			results.add(DefaultValidationMessage.SinglePositionMissingWithoutDRP.formatWithLevel(
				ValidationLevel.NOTE,
				textMissingPositions
			));
		}
		return results;
	}

	protected static <VirusT extends Virus<VirusT>> List<ValidationResult> validateNoMissingPositions(
		AlignedSequence<VirusT> alignedSequence,
		Collection<String> includeGenes
	) {
		List<AlignedGeneSeq<VirusT>> geneSeqs = alignedSequence.getAlignedGeneSequences(includeGenes);
		if (geneSeqs.isEmpty()) {
			return Collections.emptyList();
		}
		AlignedGeneSeq<VirusT> geneSeq = geneSeqs.get(0);
		GenePosition<VirusT> leftMost = new GenePosition<>(geneSeq.getGene(), 1);
		geneSeq = geneSeqs.get(geneSeqs.size() - 1);
		GenePosition<VirusT> rightMost = new GenePosition<>(geneSeq.getGene(), geneSeq.getGene().getAASize());
		Strain<VirusT> strain = alignedSequence.getStrain();
		Map<Gene<VirusT>, GeneRegions<VirusT>> unseqRegions = includeGenes.stream()
			.map(absGene -> strain.getGene(absGene))
			.collect(Collectors.toMap(
				gene -> gene,
				gene -> {
					AlignedGeneSeq<VirusT> gs = alignedSequence.getAlignedGeneSequence(gene);
					return gs == null ? (
						GeneRegions.newGeneRegions(gene, 1, gene.getAASize())
					) : gs.getUnsequencedRegions();
				}
			));
			
		Set<GenePosition<VirusT>> needGenePositions = GenePosition
			.getGenePositionsBetween(leftMost, rightMost, includeGenes);
		
		// For DRPs, the leftMost must be the begining of the first gene and the rightMost must be the ending of the last gene
		Set<GenePosition<VirusT>> needDRGenePositions = GenePosition
			.getDRGenePositionsBetween(leftMost, rightMost, includeGenes);

		Set<GenePosition<VirusT>> availableGenePositions = needGenePositions.stream()
			.filter(gpos -> {
				Gene<VirusT> gene = gpos.getGene();
				GeneRegions<VirusT> geneUnseqRegions = unseqRegions.get(gene);
				if (geneUnseqRegions == null) {
					return true;
				}
				if (geneUnseqRegions.contains(gpos.getPosition())) {
					return false;
				}
				return true;
			})	
			.collect(Collectors.toSet());
		return validateNoMissingPositions(
			needGenePositions,
			needDRGenePositions,
			availableGenePositions
		);
	}

	protected static <VirusT extends Virus<VirusT>> List<ValidationResult> validateLongGap(
		AlignedSequence<VirusT> alignedSequence,
		Collection<String> includeGenes
	) {
		int gapLenThreshold = 20;
		int totalIndels = 0;
		List<ValidationResult> result = new ArrayList<>();
		for (Mutation<VirusT> mut : alignedSequence.getMutations()) {
			if (!includeGenes.contains(mut.getAbstractGene())) {
				continue;
			}
			if (totalIndels > gapLenThreshold) {
				result.add(DefaultValidationMessage.FASTAGapTooLong.format());
				break;
			}
			if (mut.getInsertedNAs().length() > gapLenThreshold * 3) {
				result.add(DefaultValidationMessage.FASTAGapTooLong.format());
				break;
			}
			if (mut.isDeletion()) {
				totalIndels ++;
			}
			else if (mut.isInsertion()) {
				totalIndels += Math.round(mut.getInsertedNAs().length() / 3);
			}
			else {
				totalIndels = 0;
			}
		}
		return result;
	}

	protected static <VirusT extends Virus<VirusT>> List<ValidationResult> validateNAs(AlignedSequence<VirusT> alignedSequence) {
		List<String> invalids =
			alignedSequence.getInputSequence().removedInvalidChars()
			.stream().map(c -> "" + c)
			.collect(Collectors.toList());
		List<ValidationResult> result = new ArrayList<>();
		if (!invalids.isEmpty()) {
			result.add(DefaultValidationMessage.FASTAInvalidNAsRemoved.format(
				Json.dumps(String.join("", invalids))
			));
		}
		return result;
	}

	protected static <VirusT extends Virus<VirusT>> List<ValidationResult> validateNoStopCodons(
		AlignedSequence<VirusT> alignedSequence,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		MutationSet<VirusT> stopCodons = (
			alignedSequence.getMutations()
			.getStopCodons()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()))
		);
		for (Map.Entry<Gene<VirusT>, MutationSet<VirusT>> entry : stopCodons.groupByGene().entrySet()) {
			String geneDisplay = entry.getKey().getDisplay();
			MutationSet<VirusT> geneStopCodons = entry.getValue();
			int numGeneStopCodons = geneStopCodons.size();
			String geneStopText = geneStopCodons.join(", ", Mutation::getHumanFormat);
			if (numGeneStopCodons > 1) {
				results.add(DefaultValidationMessage.MultipleStopCodons.formatWithLevel(
					ValidationLevel.SEVERE_WARNING,
					numGeneStopCodons,
					geneDisplay,
					geneStopText
				));
			} else if (numGeneStopCodons > 0) {
				results.add(DefaultValidationMessage.SingleStopCodon.formatWithLevel(
					ValidationLevel.NOTE,
					geneDisplay,
					geneStopText
				));
			}
		}
		
		return results;
	}

	private static <VirusT extends Virus<VirusT>>List<ValidationResult> validateGaps(
		AlignedSequence<VirusT> alignedSequence,
		Collection<String> includeGenes
	) {
		Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> alignedGeneSeqs = alignedSequence.getAlignedGeneSequenceMap();
		List<Gene<VirusT>> seqGenes = alignedSequence.getAvailableGenes();
		List<ValidationResult> results = new ArrayList<>();

		for (Gene<VirusT> gene : seqGenes) {
			String geneName = gene.getAbstractGene();
			if (!includeGenes.contains(geneName)) {
				continue;
			}
			String geneDisplay = gene.getDisplay();
			AlignedGeneSeq<VirusT> alignedGeneSeq = alignedGeneSeqs.get(gene);
			List<FrameShift<VirusT>> frameShifts = alignedGeneSeq.getFrameShifts();
			MutationSet<VirusT> insertions = alignedGeneSeq.getInsertions();
			MutationSet<VirusT> deletions = alignedGeneSeq.getDeletions();
			MutationSet<VirusT> unusualInsertions = insertions.getUnusualMutations();
			MutationSet<VirusT> unusualDeletions = deletions.getUnusualMutations();
			MutationSet<VirusT> unusualIndels = unusualInsertions.mergesWith(unusualDeletions);
			int numTotal = frameShifts.size() + unusualInsertions.size() + unusualDeletions.size();
			String frameShiftListText = FrameShift.joinFrameShifts(frameShifts);
			String unusualIndelsListText = unusualIndels.join(", ");

			if (numTotal > 1) {
				if (frameShifts.size() > 0 && unusualIndels.size() > 0) {
					results.add(DefaultValidationMessage.MultipleUnusualIndelsAndFrameshifts.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						geneDisplay,
						numTotal,
						unusualIndelsListText,
						frameShiftListText
					));
				} else if (frameShifts.size() > 0) {
					results.add(DefaultValidationMessage.MultipleFrameShifts.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						geneDisplay,
						numTotal,
						frameShiftListText
					));
				} else {
					results.add(DefaultValidationMessage.MultipleUnusualIndels.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						geneDisplay,
						numTotal,
						unusualIndelsListText
					));
				}

			} else if (numTotal >0 ) {
				if (frameShifts.size() > 0) {
					results.add(DefaultValidationMessage.SingleFrameshift.formatWithLevel(
						ValidationLevel.WARNING,
						geneDisplay,
						frameShiftListText
					));
				} else {
					results.add(DefaultValidationMessage.SingleUnusualIndel.formatWithLevel(
						ValidationLevel.WARNING,
						geneDisplay,
						unusualIndelsListText
					));
				}

			}
		}
		return results;
	}

}
