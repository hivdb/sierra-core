/*

    Copyright (C) 2017-2020 Stanford HIVDB team

    This file is part of Sierra.

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sierra.  If not, see <https://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.sequences;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public class AlignedSite {
	private final int posAA;
	private final List<Integer> posNAs;
	private final int lengthNA;
	
	private static List<Integer> inferPosNAs(int firstPosNA, int lengthNA) {
		List<Integer> posNAs = new ArrayList<>();
		for (int offset = 0; offset < lengthNA; offset ++) {
			posNAs.add(firstPosNA + offset);
		}
		for (int offset = posNAs.size(); offset < 3; offset ++) {
			posNAs.add(null);
		}
		return posNAs;
	}
	
	public AlignedSite(int posAA, int firstPosNA, int lengthNA) {
		this(posAA, inferPosNAs(firstPosNA, lengthNA), lengthNA);
	}

	public AlignedSite(int posAA, List<Integer> posNAs, int lengthNA) {
		this.posAA = posAA;
		this.posNAs = Collections.unmodifiableList(posNAs);
		this.lengthNA = lengthNA;
	}

	public int getPosAA() { return posAA; }

	public List<Integer> getPosNAs() { return posNAs; }

	public Optional<Integer> getFirstPosNA() {
		return (
			posNAs
			.stream()
			.filter(pos -> pos != null)
			.findFirst()
		);
	}

	public Optional<Integer> getLastPosNA() {
		List<Integer> nonNulls = (
			posNAs
			.stream()
			.filter(pos -> pos != null)
			.collect(Collectors.toList())
		);
		long count = nonNulls.size();
		return nonNulls.stream().skip(Math.max(0, count - 1)).findFirst();
	}

	public int getLengthNA() { return lengthNA; }
	
}
