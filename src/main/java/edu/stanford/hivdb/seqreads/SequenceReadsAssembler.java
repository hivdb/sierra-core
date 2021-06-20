package edu.stanford.hivdb.seqreads;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class SequenceReadsAssembler<VirusT extends Virus<VirusT>> {

	public static enum AssemblyRegionType {
		UNTRANS_REGION, GENE
	}
	
	public static class AssemblyRegion<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	
		private final String name;
		private final AssemblyRegionType type;
		private final Long refStart;
		private final Long refEnd;
		private final String geneName;
		private final List<Long> trim;
	
		private transient VirusT virusIns;
		private transient Gene<VirusT> gene;

		private AssemblyRegion(
			final String name,
			final AssemblyRegionType type,
			final Long refStart,
			final Long refEnd,
			final String geneName,
			final List<Long> trim
		) {
			this.name = name;
			this.type = type;
			this.refStart = refStart;
			this.refEnd = refEnd;
			this.geneName = geneName;
			this.trim = trim;
		}
		
		public void setVirusInstance(VirusT virusIns) {
			this.virusIns = virusIns;
		}
		
		@Override
		public Gene<VirusT> getGene() {
			if (gene == null) {
				gene = virusIns.getGene(geneName); 
			}
			return gene;
		}
		
		public String getName() { return name; }
		public AssemblyRegionType getType() { return type; }
		public Long getRefStart() { return refStart; }
		public Long getRefEnd() { return refEnd; }

		public List<String> trimCodonList(List<String> codonList) { 
			List<Long> trim = this.trim == null ? Lists.newArrayList() : new ArrayList<>(this.trim);
			trim.sort(Long::compare);
			for (long napos0 : Lists.reverse(trim)) {
				napos0 --;
				int aapos0 = Math.toIntExact(napos0 / 3);
				int bp0 = Math.toIntExact(napos0 % 3);
				String codon = codonList.get(aapos0);
				codon = codon.substring(0, bp0) + codon.substring(bp0 + 1);
				codonList.set(aapos0, codon);
			}
			return codonList;
		}
		
		public String assemble(GeneSequenceReads<VirusT> geneSeqReads) {
			int aaSize = gene.getAASize();
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
					seq.add(pcr.getCodonConsensusWithIns(pcntThreshold, countThreshold));
				}
				for (int pos = lastAA + 1; pos < aaSize; pos ++) {
					seq.add("NNN");
				}
			}
			seq = this.trimCodonList(seq);
			return String.join("", seq);
		}
		
		public String assemble(UntranslatedRegion utr) {
			if (utr == null) {
				return Strings.repeat("N", Math.toIntExact(refEnd - refStart + 1));
			}
			int left = Math.toIntExact(refStart - utr.getRefStart());
			int right = Math.toIntExact(utr.getRefEnd() - refEnd);
			int trimLeft = left < 0 ? - left : 0;
			int trimRight = right < 0 ? - right : 0;
			int addLeft = left > 0 ? left : 0;
			int addRight = right > 0 ? right : 0;
			String consensus = utr.getConsensus();
			return (
				Strings.repeat("N", addLeft) +
				consensus.substring(trimLeft, consensus.length() - trimRight) +
				Strings.repeat("N", addRight)
			);
		}
	}
		
	public static <VirusT extends Virus<VirusT>>
	Map<Strain<VirusT>, SequenceReadsAssembler<VirusT>> loadJson(String raw, VirusT virusIns) {
		return (
			Json.loads(raw, new TypeToken<Map<String, List<AssemblyRegion<VirusT>>>>(){})
			.entrySet()
			.stream()
			.map(e -> {
				e.getValue()
				.stream()
				.forEach(r -> r.setVirusInstance(virusIns));
				return e;
			})
			.collect(Collectors.toMap(
				e -> virusIns.getStrain(e.getKey()),
				e -> new SequenceReadsAssembler<>(e.getValue())
			))
		);
	}
	
	private final List<AssemblyRegion<VirusT>> regions;
	
	private SequenceReadsAssembler(List<AssemblyRegion<VirusT>> regions) {
		this.regions = Collections.unmodifiableList(regions);
	}
	
	public List<AssemblyRegion<VirusT>> getRegions() { return regions; }
	
	public String assemble(
		final Map<Gene<VirusT>, GeneSequenceReads<VirusT>> allGeneSequenceReads,
		final Map<String, UntranslatedRegion> utrLookup
	) {
		StringBuilder cons = new StringBuilder();
		for (AssemblyRegion<VirusT> abr : regions) {
			if (abr.getType() == AssemblyRegionType.UNTRANS_REGION) {
				cons.append(
					abr.assemble(utrLookup.get(abr.getName()))
				);
			}
			else { // type == GENE
				cons.append(
					abr.assemble(allGeneSequenceReads.get(abr.getGene()))
				);
			}
		}
		return cons.toString().replace("-", "");
	}
	
}