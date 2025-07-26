package edu.stanford.hivdb.motif;

public class MotifPositionWildcard implements MotifPosition {

	@Override
	public String matches(char[] aas) {
		return new String(aas);
	}
	
	@Override
	public String toString() {
		return "X";
	}
}
