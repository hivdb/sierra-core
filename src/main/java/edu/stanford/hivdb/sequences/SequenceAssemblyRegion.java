package edu.stanford.hivdb.sequences;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import edu.stanford.hivdb.viruses.AssemblyRegion;
import edu.stanford.hivdb.viruses.Gene;
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
	public String assemble(AlignedGeneSeq<VirusT> geneSeq, boolean skipIns) {
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
			codons = geneSeq.getAlignedCodons(skipIns);
		}
		codons = this.trimCodonList(codons);
		return String.join("", codons);
	}

}
