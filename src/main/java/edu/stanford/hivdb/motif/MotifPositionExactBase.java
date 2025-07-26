package edu.stanford.hivdb.motif;

public final class MotifPositionExactBase implements MotifPosition {
	private final char base;
	
	protected MotifPositionExactBase(char base) {
		this.base = base;
	}
	
	@Override
	public String matches(char[] aas) {
		for (char aa : aas) {
			if (Character.toUpperCase(aa) == base) {
				return "" + base;
			}
		}
		return "";
	}
	
	@Override
	public String toString() {
		return String.valueOf(base);
	}
}
