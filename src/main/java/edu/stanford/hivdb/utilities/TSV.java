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

import java.util.Collection;
import java.util.List;
import java.util.Map;

public abstract class TSV {

	public static String dumpsHeader(String[] headers) {
		return String.join("\t", headers);
	}

	public static String dumpsHeader(List<String> headers) {
		return dumpsHeader(headers.toArray(new String[0]));
	}

	public static String dumpsBody(String[][] rows) {
		StringBuilder output = new StringBuilder();
		for (String[] row : rows) {
			output.append(String.join("\t", row));
			output.append('\n');
		}
		output.setLength(output.length() - 1);
		return output.toString();
	}

	public static String dumpsBody(List<List<String>> rows) {
		return dumpsBody(rows.stream().map(
			row -> row.toArray(new String[0])).toArray(String[][]::new));
	}

	public static String dumps(String[] headers, String[][] rows) {
		StringBuilder output = new StringBuilder();
		output.append(dumpsHeader(headers));
		output.append('\n');
		output.append(dumpsBody(rows));
		return output.toString();
	}

	public static String dumps(String[] headers, Collection<List<String>> rows) {
		return dumps(
			headers,
			rows.stream().map(
				row -> row.toArray(new String[0])).toArray(String[][]::new));
	}

	public static String dumps(List<String> headers, Collection<List<String>> rows) {
		return dumps(
			headers.toArray(new String[0]),
			rows.stream().map(
				row -> row.toArray(new String[0])).toArray(String[][]::new));
	}

	public static String dumpMaps(List<String> headers, Collection<Map<String, String>> rows, String naText) {
		return dumps(
			headers.toArray(new String[0]),
			rows.stream().map(
				row -> headers.stream().map(h -> row.getOrDefault(h, naText)).toArray(String[]::new)
			).toArray(String[][]::new)
		);
	}
}
