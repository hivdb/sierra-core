package edu.stanford.hivdb.motif;

import java.util.Set;
import java.util.stream.Collectors;


public final class MotifPositionAmbiguousBase implements MotifPosition {
	
	private final Set<Character> bases;
	
	protected MotifPositionAmbiguousBase(Set<Character> bases) {
		this.bases = bases;
	}
	
	@Override
	public String matches(char[] aas) {
		StringBuilder matches = new StringBuilder();
		for (char aa : aas) {
			aa = Character.toUpperCase(aa);
			if (bases.contains(aa)) {
				matches.append(aa);
			}
		}
		return matches.toString();
	}
	
	@Override
	public String toString() {
		String joined = bases.stream().map(String::valueOf).collect(Collectors.joining());
		if (bases.size() > 1) {
			return "[" + joined + "]";
		}
		else {
			return joined;
		}
	}

}
