package edu.stanford.hivdb.motif;

import java.util.Set;
import java.util.stream.Collectors;

public final class MotifPositionNotBase implements MotifPosition {
	
	private final Set<Character> excluded;
	
	protected MotifPositionNotBase(Set<Character> excluded) {
		this.excluded = excluded;
	}
	
	@Override
	public String matches(char[] aas) {
		StringBuilder matches = new StringBuilder();
		for (char aa : aas) {
			aa = Character.toUpperCase(aa);
			if (!excluded.contains(aa)) {
				// any AA matches will yield true result
				matches.append(aa);
			}
		}
		return matches.toString();
	}
	
	@Override
	public String toString() {
		String joined = excluded.stream().sorted().map(String::valueOf).collect(Collectors.joining());
		if (excluded.size() > 1) {
			return "~[" + joined + "]";
		}
		else {
			return "~" + joined;
		}
	}

}
