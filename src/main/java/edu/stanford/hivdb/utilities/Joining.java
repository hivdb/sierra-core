/*

    Copyright (C) 2019 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.utilities;

import java.util.Collection;
import java.util.function.Function;
import java.util.stream.Collectors;

public interface Joining {

	public static <T extends Object> String join(
		final Collection<T> objects, final String delimiter,
		final Function<T, String> toStringFunc
	) {
		return (
			objects
			.stream()
			.map(toStringFunc)
			.collect(Collectors.joining(delimiter))
		);
	}

	public static <T extends Object> String join(
		final Collection<T> objects, final Character delimiter,
		final Function<T, String> toStringFunc
	) {
		return join(objects, "" + delimiter, toStringFunc);
	}

	public static <T extends Object> String join(
		final Collection<T> objects,
		final Function<T, String> toStringFunc
	) {
		return join(objects, ",", toStringFunc);
	}

	public static <T extends Object> String join(
		final Collection<T> objects
	) {
		return join(objects, ",", T::toString);
	}

	public static <T extends Object> String join(
		final Collection<T> objects, final String delimiter
	) {
		return join(objects, delimiter, T::toString);
	}

	public static <T extends Object> String join(
		final Collection<T> objects, final Character delimiter) {
		return join(objects, "" + delimiter, T::toString);
	}
}
