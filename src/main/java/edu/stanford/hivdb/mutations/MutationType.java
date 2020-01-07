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

package edu.stanford.hivdb.mutations;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.builder.HashCodeBuilder;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Virus;

public class MutationType<VirusT extends Virus<VirusT>> implements Comparable<MutationType<VirusT>> {
	
	private final String name;
	private final Integer ordinal;
	
	public static <VirusT extends Virus<VirusT>> Map<String, MutationType<VirusT>> loadJson(String raw, VirusT virusIns) {
		Map<String, MutationType<VirusT>> mutTypes = new LinkedHashMap<>();
		List<MutationType<VirusT>> mutTypeList = Json.loads(
			raw, new TypeToken<List<MutationType<VirusT>>>(){});
		for (MutationType<VirusT> mutType : mutTypeList) {
			mutTypes.put(mutType.getName(), mutType);
		}
		if (!mutTypes.containsKey("Other")) {
			throw new IllegalStateException("Required but lack of \"Other\" mutation type.");
		}
		return Collections.unmodifiableMap(mutTypes);
	}

	
	private MutationType(String name, Integer ordinal) {
		this.name = name;
		this.ordinal = ordinal;
	}
	
	public String name() {
		return name;
	}
	
	public String getName() {
		return name;
	}
	
	public String getFullName(DrugClass<VirusT> drugClass) {
		if (name.equals(drugClass.getName())) {
			return name;
		}
		else {
			return String.format("%s %s", drugClass.getName(), name);
		}
	}
	
	@Override
	public String toString() {
		return name();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// genes are singletons
		return false;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(195300943, 5629821)
			.append(name)
			.toHashCode();
	}

	@Override
	public int compareTo(MutationType<VirusT> o) {
		return ordinal.compareTo(o.ordinal);
	}
}
