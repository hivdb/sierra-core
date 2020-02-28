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

public class NumberFormats {

	private NumberFormats() {}

	/**
	 * Formats the number of decimal places in a Double according to a set of rules
	 * that we use frequently:
	 * <ul>
	 * 	<li>&gt;=10: decimal = 0</li>
	 * 	<li>1.0 to 9.9: one decimal place</li>
	 * 	<li>&lt;1.0: one significant figure.</li>
	 * </ul>
	 * 
	 * @param number 	Real number
	 * @return 			formatted output
	 */
	public static final String prettyDecimalAsString(double number) {
		final String fmt;
		if (number >= 10.0) {
			fmt = "%.0f";
		} else if (number >= 1) {
			fmt = "%.1f";
		} else if (number >= 0.1) {
			fmt = "%.2f";
		} else if (number >= 0.001) {
			fmt = "%.3f";
		} else {
			fmt = "%.1g";
		}
		return String.format(fmt, number);
	}

	public static final Double prettyDecimal(double number) {
		return Double.parseDouble(prettyDecimalAsString(number));
	}
}
