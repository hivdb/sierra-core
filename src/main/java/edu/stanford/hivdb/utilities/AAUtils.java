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

import java.util.Arrays;

public class AAUtils {

	private AAUtils() {}

	/**
	 * Normalize the input AAs.
	 *
	 * The code explains the normalization rules.
	 */
	public static String normalizeAAs(String aas) {
		if (aas == null) return null;

		aas = aas.replaceAll("^[dD]elet(e|ion)|d(el)?|~$", "-")
			     .replaceAll("^[iI]nsert(ion)?|i(ns)?$|#", "_")
			     .replaceAll("[.Z]", "*");

		if (aas.length() > 1 && !aas.contains("_")) {
			char[] aaChars = aas.toUpperCase().toCharArray();
			Arrays.sort(aaChars);
			return new String(aaChars);
		}

		return aas.toUpperCase();
	}

	public static String toHIVDBFormat(String input) {
		return input
			.replace("Insertion", "#")
			.replace("Deletion", "~")
			.replace("ins", "#")
			.replace("del", "~")
			.replace('i', '#')
			.replace('_', '#')
			.replace('d', '~')
			.replace('-', '~');
	}

	public static String toInternalFormat(String input) {
		return input
			.replace("Insertion", "_")
			.replace("Deletion", "-")
			.replace("ins", "_")
			.replace("del", "-")
			.replace('i', '_')
			.replace('#', '_')
			.replace('d', '-')
			.replace('~', '-');
	}

	public static String toASIFormat(String input) {
		return input
			.replace("Insertion", "i")
			.replace("Deletion", "d")
			.replace("ins", "i")
			.replace("del", "d")
			.replace('#', 'i')
			.replace('_', 'i')
			.replace('~', 'd')
			.replace('-', 'd');
	}

}


