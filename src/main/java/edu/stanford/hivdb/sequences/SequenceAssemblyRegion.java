package edu.stanford.hivdb.sequences;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import edu.stanford.hivdb.viruses.AssemblyRegion;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceAssemblyRegion<VirusT extends Virus<VirusT>> extends AssemblyRegion<
	VirusT,
	AlignedGeneSeq<VirusT>,
	SequenceAssemblyRegion<VirusT>
> {

	public SequenceAssemblyRegion(
		final String name,
		final AssemblyRegionType type,
		final Long refStart,
		final Long refEnd,
		final Gene<VirusT> gene,
		final List<Long> trim
	) {
		super(name, type, refStart, refEnd, gene, trim);
	}

	@Override
	public String assemble(AlignedGeneSeq<VirusT> geneSeq, Strain<VirusT> targetStrain, boolean skipIns) {
		List<String> codons;
		if (geneSeq == null) {
			int aaSize = getGene().getAASize();
			codons = (
				Stream.generate(() -> "NNN")
				.limit(aaSize)
				.collect(Collectors.toList())
			);
		}
		else {
			String alignedNAs;
			if (targetStrain != geneSeq.getStrain()) {
				alignedNAs = geneSeq.getAdjustedAlignedNAsNoTrim(targetStrain.getName());
			}
			else {
				alignedNAs = geneSeq.getAlignedNAsNoTrim();
			}
			codons = Lists.newArrayList(Splitter.fixedLength(3).split(alignedNAs));
		}
		codons = this.trimCodonList(codons);
		return String.join("", codons);
	}

}
