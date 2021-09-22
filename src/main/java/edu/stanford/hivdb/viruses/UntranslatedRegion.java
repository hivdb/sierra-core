package edu.stanford.hivdb.viruses;

public class UntranslatedRegion {

	private final String name;
	private final Long refStart;
	private final Long refEnd;
	private final String consensus;
	
	public UntranslatedRegion(
		final String name,
		final Long refStart,
		final Long refEnd,
		final String consensus
	) {
		this.name = name;
		this.refStart = refStart;
		this.refEnd = refEnd;
		this.consensus = consensus;
	}
	
	public String getName() { return name; }
	public Long getRefStart() { return refStart; }
	public Long getRefEnd() { return refEnd; }
	public String getConsensus() { return consensus; }
}
