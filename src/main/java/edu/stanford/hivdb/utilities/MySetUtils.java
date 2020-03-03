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
