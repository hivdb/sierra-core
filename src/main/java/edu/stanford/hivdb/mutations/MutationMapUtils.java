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

package edu.stanford.hivdb.mutations;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.stanford.hivdb.viruses.Virus;

public class MutationMapUtils {

	public enum SortOrder {ASC, DESC};

	/**
	 * @param <T>			Object type
	 * @param unsortedMap	unsorted Objects with doubles
	 * @param sortOrder		Sort order
	 * @return 				sorted map
	 */
	public static<T> Map<T, Double> sortByComparator(Map<T, Double> unsortedMap, SortOrder sortOrder) {
		List<Map.Entry<T, Double>> list = new LinkedList<Map.Entry<T, Double>>(unsortedMap.entrySet());

		// Sort list with comparator
		Collections.sort(list, new Comparator<Map.Entry<T, Double>>(){
			public int compare(Map.Entry<T, Double> o1, Map.Entry<T, Double> o2) {
				if (sortOrder.equals(SortOrder.DESC)) {
					return (o2.getValue()).compareTo(o1.getValue());
				} else return (o1.getValue()).compareTo(o2.getValue());
			}
		});

		// Convert sorted map back to a Map
		Map<T, Double> sortedMap = new LinkedHashMap<>();
		for (Iterator<Map.Entry<T, Double>> it = list.iterator(); it.hasNext();) {
			Map.Entry<T, Double> entry = it.next();
			sortedMap.put(entry.getKey(), entry.getValue());
		}

		return sortedMap;
	}

	/**
	 *
	 * @param <VirusT>	Virus subtype
	 * @param map		Mutations with scores
	 * @return 			a map in which the value is an integer rather than a double
	 */
	public static <VirusT extends Virus<VirusT>> Map<Mutation<VirusT>, Integer> convertMutScoresToInts(Map<Mutation<VirusT>, Double> map) {
		Map<Mutation<VirusT>, Integer> newMutScores = new HashMap<>();
		for (Mutation<VirusT> mut : map.keySet()) {
			int value = map.get(mut).intValue();
			newMutScores.put(mut, value);
		}
		return newMutScores;
	}

	/**
	 * Used only by getScoredMuts in Algorithm comparison
	 * 
	 * @param <VirusT> 	Virus subtype
	 * @param map		Mutations with scores
	 * @return			Formatted strint for print
	 */
	public static <VirusT extends Virus<VirusT>> String printMutScoresAsInts(Map<Mutation<VirusT>, Double> map) {
		StringBuffer output = new StringBuffer();
		for (Mutation<VirusT> mut : map.keySet()) {
			int value = map.get(mut).intValue();
			output.append(mut.getHumanFormat() + " (" + value + "), ");
		}
		if (output.length() > 0) {
			output.setLength(output.length() - 2);
		}
		return output.toString();
	}

	public static <VirusT extends Virus<VirusT>> String printMutSetScoresAsInts(Map<MutationSet<VirusT>, Double> comboMutsSortedByScore) {
		StringBuffer output = new StringBuffer();
		for (MutationSet<VirusT> mutList : comboMutsSortedByScore.keySet()) {
			int value = comboMutsSortedByScore.get(mutList).intValue();
			String mutListOutput = mutList.join(" + ");
			output.append(mutListOutput + " (" + value + "), ");
		}
		output.setLength(output.length() - 2);
		return output.toString();
	}

	public static <VirusT extends Virus<VirusT>> String printMutScoresAsDouble(Map<Mutation<VirusT>, Double> map) {
		StringBuffer output = new StringBuffer();
		for (Mutation<VirusT> mut : map.keySet()) {
			Double value = map.get(mut);
			output.append(mut.getHumanFormat() + " (" + value + "), ");
		}
		output.setLength(output.length() - 2);
		return output.toString();
	}

	public static <VirusT extends Virus<VirusT>> String printMutSetScoresAsDouble(Map<MutationSet<VirusT>, Double> comboMutsSortedByScore) {
		StringBuffer output = new StringBuffer();
		for (MutationSet<VirusT> mutSet : comboMutsSortedByScore.keySet()) {
			Double value = comboMutsSortedByScore.get(mutSet);
			String mutSetOutput = mutSet.join(" + ");
			output.append(mutSetOutput + " (" + value + "), ");
		}
		output.setLength(output.length() - 2);
		return output.toString();
	}
}
