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

package edu.stanford.hivdb.viruses;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.utilities.Json;

public class Strain<VirusT extends Virus<VirusT>> implements Comparable<Strain<VirusT>> {
	
	private final String name;
	private final Integer ordinal;
	private final String displayText;
	private final Integer absoluteFirstNA;
	private final Integer nucaminoStopCodonPenalty;
	private final Integer nucaminoGapOpeningPenalty;
	private final Integer nucaminoGapExtensionPenalty;
	private final Integer nucaminoIndelCodonOpeningBonus;
	private final Integer nucaminoIndelCodonExtensionBonus;
	private transient String nucaminoProfileString;
	private transient File nucaminoProfileFile;

	private transient VirusT virusInstance;
	private transient Set<DrugClass<VirusT>> drugClasses;

	public static <VirusT extends Virus<VirusT>> Map<String, Strain<VirusT>> loadJson(String raw, VirusT virusIns) {
		Map<String, Strain<VirusT>> strains = new LinkedHashMap<>();
		List<Strain<VirusT>> strainList = Json.loads(
			raw, new TypeToken<List<Strain<VirusT>>>(){});
		for (Strain<VirusT> strain : strainList) {
			strain.virusInstance = virusIns;
			strains.put(strain.getName(), strain);
		}
		return Collections.unmodifiableMap(strains);
	}
	
	private Strain(
		String name, Integer ordinal, String displayText,
		Integer offsetNA,
		Integer nucaminoStopCodonPenalty,
		Integer nucaminoGapOpeningPenalty,
		Integer nucaminoGapExtensionPenalty,
		Integer nucaminoIndelCodonOpeningBonus,
		Integer nucaminoIndelCodonExtensionBonus
	) {
		this.name = name;
		this.ordinal = ordinal;
		this.displayText = displayText;
		this.absoluteFirstNA = offsetNA;
		this.nucaminoStopCodonPenalty = nucaminoStopCodonPenalty;
		this.nucaminoGapOpeningPenalty = nucaminoGapOpeningPenalty;
		this.nucaminoGapExtensionPenalty = nucaminoGapExtensionPenalty;
		this.nucaminoIndelCodonOpeningBonus = nucaminoIndelCodonOpeningBonus;
		this.nucaminoIndelCodonExtensionBonus = nucaminoIndelCodonExtensionBonus;
	}
	
	public VirusT getVirusInstance() {
		return virusInstance;
	}
	
	public String name() {
		return name;
	}
	
	public String getName() {
		return name;
	}

	public String getDisplayText() {
		return displayText;
	}
	
	public Map<String, Gene<VirusT>> getNucaminoGeneMap() {
		Map<String, Gene<VirusT>> geneMap = new LinkedHashMap<>();
		for (Gene<VirusT> gene : getGenes()) {
			geneMap.put(gene.getAbstractGene().toUpperCase(), gene);
		}
		return geneMap;
	}
	
	public String makeNucaminoProfileString() {
		if (nucaminoProfileString == null) {
			Map<String, Object> profile = new LinkedHashMap<>();
			profile.put("StopCodonPenalty", nucaminoStopCodonPenalty);
			profile.put("GapOpeningPenalty", nucaminoGapOpeningPenalty);
			profile.put("GapExtensionPenalty", nucaminoGapExtensionPenalty);
			profile.put("IndelCodonOpeningBonus", nucaminoIndelCodonOpeningBonus);
			profile.put("IndelCodonExtensionBonus", nucaminoIndelCodonExtensionBonus);
			Map<String, String> refSeqs = new LinkedHashMap<>();
			for (Gene<VirusT> gene : getGenes()) {
				refSeqs.put(gene.getAbstractGene().toUpperCase(), gene.getRefSequence());
			}
			profile.put("ReferenceSequences", refSeqs);
			// TODO: support PositionalIndelScores
			nucaminoProfileString = Json.dumpsUgly(profile);;
		}
		return nucaminoProfileString;
	}
	
	public File makeNucaminoProfileFile() {
		if (nucaminoProfileFile == null || ! nucaminoProfileFile.exists()) {
			try {
				nucaminoProfileFile = File.createTempFile("nucamino-" + name, ".yaml");
				BufferedWriter bw = new BufferedWriter(new FileWriter(nucaminoProfileFile));
				bw.write(makeNucaminoProfileString());
				bw.close();
				nucaminoProfileFile.setReadOnly();
				nucaminoProfileFile.deleteOnExit();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		return nucaminoProfileFile;
	}

	public Integer getAbsoluteFirstNA() {
		return absoluteFirstNA;
	}
	
	public Gene<VirusT> getGene(String name) {
		return virusInstance.getGene(this.name + name);
	}

	public Collection<Gene<VirusT>> getGenes() {
		return virusInstance.getGenes(this);
	}
	
	public Set<DrugClass<VirusT>> getDrugClasses() {
		if (drugClasses == null) {
			Set<DrugClass<VirusT>> drugClasses = (
				virusInstance
				.getDrugClasses()
				.stream()
				.filter(dc -> dc.getStrains().contains(this))
				.collect(Collectors.toCollection(LinkedHashSet::new))
			);
			this.drugClasses = Collections.unmodifiableSet(drugClasses);
		}
		return drugClasses;
	}
	
	@Override
	public String toString() {
		return name();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// strains are singletons
		return false;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(643277, 541754393)
			.append(name)
			.toHashCode();
	}

	@Override
	public int compareTo(Strain<VirusT> o) {
		return ordinal.compareTo(o.ordinal);
	}
	

}
