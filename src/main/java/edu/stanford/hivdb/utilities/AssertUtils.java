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

public class AssertUtils {
	
	/**
	 * Assert if input value is null
	 * 
	 * @param <T>	Java Object Type
	 * @param value	Any Java Object
	 * @param msg	Assertion error message
	 * @param args	Arguments for formatting message string, using String.format
	 * @return		value is not null
	 * 
	 * @throws IllegalArgumentException value is null
	 */
	public static <T> T notNull(T value, String msg, Object... args) {
		if (value == null) {
			throw new IllegalArgumentException(String.format(msg, args));
		}
		return value;
	}
	
	public static void isTrue(Boolean value, String msg, Object... args) {
		if (!value) {
			throw new IllegalArgumentException(String.format(msg, args));
		}
	}
	
}
