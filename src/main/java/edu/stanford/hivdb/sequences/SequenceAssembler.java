package edu.stanford.hivdb.sequences;

import java.util.List;
import java.util.Map;

import edu.stanford.hivdb.viruses.Assembler;
import edu.stanford.hivdb.viruses.AssemblyRegion.AssemblyRegionType;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.UntranslatedRegion;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceAssembler<VirusT extends Virus<VirusT>> extends Assembler<
	VirusT,
	AlignedGeneSeq<VirusT>,
	SequenceAssemblyRegion<VirusT>,
	SequenceAssembler<VirusT>
> {

	@SuppressWarnings("unchecked")
	public static <VirusT extends Virus<VirusT>>
	Map<Strain<VirusT>, SequenceAssembler<VirusT>> loadJson(
		String raw,
		VirusT virusIns
	) {
		return Assembler.loadJson(SequenceAssembler.class, SequenceAssemblyRegion.class, raw, virusIns);
	}

	public SequenceAssembler(List<SequenceAssemblyRegion<VirusT>> regions) {
		super(regions);
	}
	
	public String assemble(
		final Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> alignedGeneSeqs,
		final Map<String, UntranslatedRegion> utrLookup,
		final boolean strictAlign
	) {
		StringBuilder cons = new StringBuilder();
		for (SequenceAssemblyRegion<VirusT> abr : getRegions()) {
			if (abr.getType() == AssemblyRegionType.UNTRANS_REGION) {
				/* System.out.printf(
					"name: %s; length: %d; actual: %d\n",
					abr.getName(),
					abr.getRefEnd() - abr.getRefStart() + 1,
					abr.assemble(utrLookup.get(abr.getName()), skipIns).length()
				); */
				cons.append(
					abr.assemble(utrLookup.get(abr.getName()), strictAlign)
				);
			}
			else { // type == GENE
				/* System.out.printf(
					"name: %s; length: %d; actual: %d\n",
					abr.getName(),
					abr.getGene().getNASize(),
					abr.assemble(
						alignedGeneSeqs.get(abr.getGene()),
						skipIns
					).length()
				); */
				cons.append(
					abr.assemble(
						alignedGeneSeqs.get(abr.getGene()),
						strictAlign
					)
				);
			}
		}
		if (strictAlign) {
			return cons.toString();
		}
		else {
			return cons.toString().replace("-", "").replaceAll("^N+|N+$", "");
		}
	}
	
}