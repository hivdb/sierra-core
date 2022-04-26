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

package edu.stanford.hivdb.utilities;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class MyStringUtils {

	/**
	 * Only for static access. DO NOT instantiate this class
	 */
	private MyStringUtils() {}

	/**
	 * Returns true if two strings share one or more characters
	 * @param aas1 	String waiting to be compare
	 * @param aas2 	String waiting to be compare
	 * @return		Has shared characters or not
	 */
	public static boolean hasSharedChar(String aas1, String aas2) {
		for (char aa1 : aas1.toCharArray()) {
			for (char aa2 : aas2.toCharArray()) {
				if (aa1 == aa2) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Returns alphabetically sorted string of input.
	 *
	 * Example:
	 *   "gfedcba" =&gt; "abcdefg"
	 *
	 * @param input	String waiting to be sorted
	 * @return		Alphabetically sorted string
	 */
	public static String sortAlphabetically(String input) {
		char[] chars = input.toCharArray();
		Arrays.sort(chars);
		return new String(chars);
	}
	
	/**
	 * Format a list in English.
	 * 
	 * This function is a simple implementation of JavaScript's Intl.ListFormat
	 * for English language. The rule is to join each item except the last one by
	 * commas, and join the last one with previous ones by text " and " or " or ".
	 * 
	 * For example, ["Apple", "Orange", "Banana"] will become
	 * "Apple, Orange, and Banana".
	 *  
	 * @param iter the list to be formatted
	 * @param isDisjunction true for the last joining text to be " or " and false to be " and " 
	 * @return the formatted result.
	 */
	public static String listFormat(Iterable<?> iter, boolean isDisjunction) {
		List<String> lst = Lists.newArrayList(iter)
			.stream()
			.map(Object::toString)
			.collect(Collectors.toList());
		String lastJoiner = isDisjunction ? " or " : " and ";
		int size = lst.size();
		switch (size) {
			case 0:
				return "";
			case 1:
				return lst.get(0);
			case 2:
				return String.join(lastJoiner, lst);
			default:
				return String.join(
					", ", lst.subList(0, size - 1)
				) + "," + lastJoiner + lst.get(size - 1);
		}
	}
	
	public static String orListFormat(Iterable<?> iter) {
		return listFormat(iter, true);
	}
	
	public static String andListFormat(Iterable<?> iter) {
		return listFormat(iter, false);
	}

}
