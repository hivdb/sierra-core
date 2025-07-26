package edu.stanford.hivdb.motif;

public class MotifMatch {
	public final long firstAA;
	public final long lastAA;
	public final String matched;
	public final String pattern;
	public MotifMatch(long firstAA, long lastAA, String matched, String pattern) {
		this.firstAA = firstAA;
		this.lastAA = lastAA;
		this.matched = matched;
		this.pattern = pattern;
	}
	
	@Override
	public String toString() {
		return String.format("MotifMatch<%d-%d:%s>", firstAA, lastAA, matched);
	}
}
