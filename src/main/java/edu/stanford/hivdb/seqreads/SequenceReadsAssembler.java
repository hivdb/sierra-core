package edu.stanford.hivdb.seqreads;

import java.util.List;
import java.util.Map;

import edu.stanford.hivdb.viruses.Assembler;
import edu.stanford.hivdb.viruses.AssemblyRegion.AssemblyRegionType;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.UntranslatedRegion;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceReadsAssembler<VirusT extends Virus<VirusT>> extends Assembler<
	VirusT,
	GeneSequenceReads<VirusT>,
	SequenceReadsAssemblyRegion<VirusT>,
	SequenceReadsAssembler<VirusT>
> {
	
	@SuppressWarnings("unchecked")
	public static <VirusT extends Virus<VirusT>>
	Map<Strain<VirusT>, SequenceReadsAssembler<VirusT>> loadJson(
		String raw,
		VirusT virusIns
	) {
		return Assembler.loadJson(SequenceReadsAssembler.class, SequenceReadsAssemblyRegion.class, raw, virusIns);
	}

	public SequenceReadsAssembler(List<SequenceReadsAssemblyRegion<VirusT>> regions) {
		super(regions);
	}
	
	public String assemble(
		final Map<Gene<VirusT>, GeneSequenceReads<VirusT>> allGeneSequenceReads,
		final Map<String, UntranslatedRegion> utrLookup,
		final boolean skipIns,
		final boolean includeAmbiguousNA
	) {
		StringBuilder cons = new StringBuilder();
		for (SequenceReadsAssemblyRegion<VirusT> abr : getRegions()) {
			if (abr.getType() == AssemblyRegionType.UNTRANS_REGION) {
				cons.append(
					abr.assemble(utrLookup.get(abr.getName()), skipIns)
				);
			}
			else { // type == GENE
				cons.append(
					abr.assemble(allGeneSequenceReads.get(abr.getGene()), skipIns, includeAmbiguousNA)
				);
			}
		}
		if (skipIns) {
			return cons.toString();
		}
		else {
			return cons.toString().replace("-", "");
		}
	}

}