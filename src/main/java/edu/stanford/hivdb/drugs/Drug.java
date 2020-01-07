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

package edu.stanford.hivdb.drugs;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Virus;

public class Drug<VirusT extends Virus<VirusT>> implements Comparable<Drug<VirusT>> {

	private final String name;
	private final String fullName;
	private final String displayAbbr;
	private final String drugClass;

	private transient VirusT virusInstance;

	public static <VirusT extends Virus<VirusT>> Map<String, Drug<VirusT>> loadJson(String raw, VirusT virusIns) {
		Map<String, Drug<VirusT>> drugMap = new TreeMap<>();
		List<Drug<VirusT>> drugs = Json.loads(
			raw, new TypeToken<List<Drug<VirusT>>>() {});
		for (Drug<VirusT> drug : drugs) {
			drug.virusInstance = virusIns;
			drugMap.put(drug.getName(), drug);
			drugMap.put(drug.getFullName(), drug);
			drugMap.put(drug.getDisplayAbbr(), drug);
		}
		return Collections.unmodifiableMap(drugMap);
	}
	
	private Drug(
		String name, String fullName,
		String displayAbbr, String drugClass
	) {
		this.name = name;
		this.fullName = fullName;
		this.displayAbbr = displayAbbr;
		this.drugClass = drugClass;
	}

	public void setVirusInstance(VirusT virusInstance) {
		this.virusInstance = virusInstance;
	}
	
	private void checkVirusInstance() {
		if (virusInstance == null) {
			throw new ExceptionInInitializerError(
				"Object not properly initialzed: virusInstance is empty"
			);
		}
	}
	
	public DrugClass<VirusT> getDrugClass() {
		checkVirusInstance();
		return virusInstance.getDrugClass(drugClass);
		
	}
	
	public String name() {
		return name;
	}
	
	public String getName() {
		return name;
	}

	public String getFullName() {
		return fullName;
	}

	public String getDisplayAbbr() {
		return displayAbbr;
	}

	@Override
	public String toString() {
		return name();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// drugs are singletons
		return false;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(855342123, 6357125)
			.append(name)
			.toHashCode();
	}

	@Override
	public int compareTo(Drug<VirusT> o) {
		return name.compareTo(o.name);
	}

}
