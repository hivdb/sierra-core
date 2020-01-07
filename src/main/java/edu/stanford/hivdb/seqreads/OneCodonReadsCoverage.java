package edu.stanford.hivdb.seqreads;

import java.util.LinkedHashMap;
import java.util.Map;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.GenePosition;

public class OneCodonReadsCoverage<VirusT extends Virus<VirusT>> {
	private final Gene<VirusT> gene;
	private final long position;
	private final long totalReads;
	private final boolean isTrimmed;
	
	public OneCodonReadsCoverage(
		final Gene<VirusT> gene,
		final long position,
		final long totalReads,
		final boolean isTrimmed
	) {
		this.gene = gene;
		this.position = position;
		this.totalReads = totalReads;
		this.isTrimmed = isTrimmed;
	}
	
	public Gene<VirusT> getGene() {
		return gene;
	}
	
	public Long getPosition() {
		return position;
	}
	
	public GenePosition<VirusT> getGenePosition() {
		return new GenePosition<VirusT>(gene, (int) position);
	}
	
	public Long getTotalReads() {
		return totalReads;
	}
	
	public Boolean isTrimmed() {
		return isTrimmed;
	}
	
	public Integer getPositionInStrain() {
		return getGenePosition().getPositionInStrain();
	}

	public Map<String, Object> extMap() {
		Map<String, Object> result = new LinkedHashMap<>();
		Map<String, Object> geneMap = new LinkedHashMap<>();
		geneMap.put("name", gene.getName());
		result.put("gene", geneMap);
		result.put("position", position);
		result.put("totalReads", totalReads);
		result.put("isTrimmed", isTrimmed);
		return result;
	}
	
}