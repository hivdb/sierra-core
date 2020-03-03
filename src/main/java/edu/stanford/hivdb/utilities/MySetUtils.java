/*

    Copyright (C) 2019-2020 Stanford HIVDB team

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
package edu.stanford.hivdb.utilities;

import java.util.LinkedHashSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class MySetUtils {
	
	public static <E> Boolean anyMatch(Set<E> unfiltered, Predicate<? super E> predicate) {
		return (
			unfiltered
			.parallelStream()
			.anyMatch(predicate)
		);
	}
	
	public static <E> SortedSet<E> filter(SortedSet<E> unfiltered, Predicate<? super E> predicate) {
		return (
			unfiltered
			.parallelStream()
			.filter(predicate)
			.collect(Collectors.toCollection(TreeSet::new))
		);
	}

	public static <E> Set<E> filter(Set<E> unfiltered, Predicate<? super E> predicate) {
		return (
			unfiltered
			.parallelStream()
			.filter(predicate)
			.collect(Collectors.toCollection(LinkedHashSet::new))
		);
	}

}
