package edu.stanford.hivdb.motif;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class MotifUtils {
	
	private static Set<Character> parseAmbiguous(String motif, int offset) {
		int closeIdx = motif.indexOf(']', offset);
		if (closeIdx == -1) throw new IllegalArgumentException("Unclosed [ at position " + offset);
		Set<Character> bases = new LinkedHashSet<>();
		for (int j = offset + 1; j < closeIdx; j ++) {
			bases.add(motif.charAt(j));
		}
		return bases;
	}

	public static List<MotifPosition> parseMotif(String motif) {
		List<MotifPosition> positions = new ArrayList<>();
		int i = 0;
		int stringLen = motif.length();
		
		while (i < stringLen) {
			char c = motif.charAt(i);
			
			if (c == '-') {
				// Optional delimiter, skipped
				i ++;
			}
			else if (c == 'X' || c == '*' || c == '.') {
				// wildcard for everything
				positions.add(new MotifPositionWildcard());
				i ++;
			}
			else if (c == '[') {
				// Positive set, e.g. [ST]
				Set<Character> bases = parseAmbiguous(motif, i);
				positions.add(new MotifPositionAmbiguousBase(bases));
				i += bases.size() + 2;
			}
			else if (c == '~') {
				// Negative set, e.g. ~[ST], ~P
				i ++;
				if (i >= stringLen) throw new IllegalArgumentException("Dangling ~ at end of motif");
				c = motif.charAt(i);
				if (c == '[') {
					Set<Character> excluded = parseAmbiguous(motif, i);
					positions.add(new MotifPositionNotBase(excluded));
					i += excluded.size() + 2;
				}
				else {
					positions.add(new MotifPositionNotBase(Set.of(c)));
					i ++;
				}
			}
			else {
				// Exact single base
				positions.add(new MotifPositionExactBase(c));
				i ++;
			}
		}
		return positions;
	}
	
	public static List<MotifMatch> findMotifMatches(
		Map<? extends Number, String> aaLookup,
		String motifPattern
	) {
		List<MotifPosition> motif = parseMotif(motifPattern);
		List<MotifMatch> matches = new ArrayList<>();
		int motifSize = motif.size();
		
		// convert whatever key to Long
		Map<Long, String> longAALookup = aaLookup.entrySet().stream().collect(
			Collectors.toMap(
				e -> e.getKey().longValue(),
				Map.Entry::getValue
			)
		);		
		
		long minPos = longAALookup.keySet().stream().min(Long::compare).get();
		long maxPos = longAALookup.keySet().stream().max(Long::compare).get();
		
		for (long i = minPos; i <= maxPos - motifSize + 1; i ++) {
			boolean match = true;
			StringBuilder matched = new StringBuilder();
			for (int j = 0; j < motifSize; j ++) {
				String aas = longAALookup.getOrDefault(i + j, null);
				if (aas == null) {
					match = false;
					break;
				}
				String matchedAAs = motif.get(j).matches(aas.toCharArray());
				int matchedLen = matchedAAs.length();
				if (matchedLen == 0) {
					match = false;
					break;
				}
				if (j > 0) {
					matched.append('-');
				}
				if (matchedLen == 1) {
					matched.append(matchedAAs);
				}
				else {
					matched.append('[');
					matched.append(matchedAAs);
					matched.append(']');
				}
			}
			if (match) {
				matches.add(new MotifMatch(i, i + motifSize - 1, matched.toString(), motifPattern));
			}
		}
		return matches;
	}
	
}
