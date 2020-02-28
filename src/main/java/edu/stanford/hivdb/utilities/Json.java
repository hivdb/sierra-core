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

import java.io.Reader;
import java.lang.reflect.Type;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;

public class Json {

	private static final Gson gson;
	private static final Gson uglyGson;

	static {
		gson = new GsonBuilder()
			.serializeNulls().setPrettyPrinting().create();
		uglyGson = new GsonBuilder()
			.serializeNulls().create();
	}

	public static String dumps(Object object) {
		return gson.toJson(object);
	}
	
	public static String dumpsUgly(Object object) {
		return uglyGson.toJson(object);
	}

	public static <T> T loads(String json, Class<T> type) {
		return gson.fromJson(json, type);
	}

	public static <T> T loads(String json, Type type) {
		return gson.fromJson(json, type);
	}
	
	public static <T> T loads(String json, TypeToken<T> typeToken) {
		return gson.fromJson(json, typeToken.getType());
	}

	public static <T> T loads(Reader json, Class<T> type) {
		return gson.fromJson(json, type);
	}

	public static <T> T loads(Reader json, Type type) {
		return gson.fromJson(json, type);
	}
	
	public static <T> T loads(Reader json, TypeToken<T> typeToken) {
		return gson.fromJson(json, typeToken.getType());
	}
}
