package edu.stanford.hivdb.seqreads;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.viruses.AssemblyRegion;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

public class SequenceReadsAssemblyRegion<VirusT extends Virus<VirusT>> extends AssemblyRegion<
	VirusT,
	GeneSequenceReads<VirusT>,
	SequenceReadsAssemblyRegion<VirusT>
> {

	public SequenceReadsAssemblyRegion(
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
	public String assemble(GeneSequenceReads<VirusT> geneSeqReads, boolean skipIns) {
		int aaSize = getGene().getAASize();
		List<String> seq;
		if (geneSeqReads == null) {
			seq = (
				Stream.generate(() -> "NNN")
				.limit(aaSize)
				.collect(Collectors.toCollection(ArrayList::new))
			);
		}
		else {
			int firstAA = geneSeqReads.getFirstAA();
			int lastAA = geneSeqReads.getLastAA();
			List<PositionCodonReads<VirusT>> posCodonReads = geneSeqReads.getAllPositionCodonReads();
			CutoffCalculator<VirusT> cutoffObj = geneSeqReads.getCutoffObj();
			double pcntThreshold = cutoffObj.getActualMinPrevalence();
			long countThreshold = cutoffObj.getMinCodonReads();
			seq = new ArrayList<>();
			for (int pos = 1; pos < firstAA; pos ++) {
				seq.add("NNN");
			}
			long prevPos = firstAA - 1;
			for (PositionCodonReads<VirusT> pcr : posCodonReads) {
				long curPos = pcr.getPosition();
				for (long pos = prevPos + 1; pos < curPos; pos ++) {
					seq.add("NNN");
				}
				prevPos = curPos;
				if (skipIns) {
					seq.add(pcr.getCodonConsensusWithoutIns(pcntThreshold, countThreshold));
				}
				else {
					seq.add(pcr.getCodonConsensusWithIns(pcntThreshold, countThreshold));
				}
			}
			for (int pos = lastAA + 1; pos < aaSize; pos ++) {
				seq.add("NNN");
			}
		}
		seq = this.trimCodonList(seq);
		return String.join("", seq);
	}
	
}