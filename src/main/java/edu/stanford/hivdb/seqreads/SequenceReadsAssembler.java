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

	public SequenceReadsAssembler(Strain<VirusT> strain, List<SequenceReadsAssemblyRegion<VirusT>> regions) {
		super(strain, regions);
	}
	
	public String assemble(
		final Map<Gene<VirusT>, GeneSequenceReads<VirusT>> allGeneSequenceReads,
		final Map<String, UntranslatedRegion> utrLookup,
		final boolean strictAlign,
		final boolean includeAmbiguousNA
	) {
		StringBuilder cons = new StringBuilder();
		for (SequenceReadsAssemblyRegion<VirusT> abr : getRegions()) {
			if (abr.getType() == AssemblyRegionType.UNTRANS_REGION) {
				cons.append(
					abr.assemble(utrLookup.get(abr.getName()), strictAlign)
				);
			}
			else { // type == GENE
				cons.append(
					abr.assemble(allGeneSequenceReads.get(abr.getGene()), getStrain(), strictAlign, includeAmbiguousNA)
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