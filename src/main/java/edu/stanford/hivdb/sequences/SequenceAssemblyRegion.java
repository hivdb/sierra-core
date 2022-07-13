package edu.stanford.hivdb.sequences;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
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
		int aaSize = getGene().getAASize();
		List<String> seq;
		if (geneSeq == null) {
			seq = (
				Stream.generate(() -> "NNN")
				.limit(aaSize)
				.collect(Collectors.toCollection(ArrayList::new))
			);
		}
		else {
			seq = new ArrayList<>();
			// TODO: merge this part with geneSeq.getAlignedNAs()
			Map<Integer, String> codonLookup = geneSeq.getCodonLookup();
			for (int pos = 1; pos <= aaSize; pos ++) {
				String posCodon = codonLookup.getOrDefault(pos, "NNN");
				if (skipIns) {
					seq.add(posCodon.substring(0, 3));
				}
				else {
					seq.add(posCodon.replace("-", ""));
				}
			}
		}
		seq = this.trimCodonList(seq);
		return String.join("", seq);
	}

}
