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
package edu.stanford.hivdb.seqreads;

import java.util.LinkedHashMap;
import java.util.Map;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;
import edu.stanford.hivdb.mutations.GenePosition;

public class OneCodonReadsCoverage<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	private final Gene<VirusT> gene;
	private final long position;
	private final long totalReads;
	private final boolean isTrimmed;
	
	public OneCodonReadsCoverage(
		final Gene<VirusT> gene,
		final long position,
		final long totalReads,
		final boolean isTrimmed
	) {
		this.gene = gene;
		this.position = position;
		this.totalReads = totalReads;
		this.isTrimmed = isTrimmed;
	}
	
	public Gene<VirusT> getGene() {
		return gene;
	}
	
	public Long getPosition() {
		return position;
	}
	
	public GenePosition<VirusT> getGenePosition() {
		return new GenePosition<VirusT>(gene, (int) position);
	}
	
	public Long getTotalReads() {
		return totalReads;
	}
	
	public Boolean isTrimmed() {
		return isTrimmed;
	}
	
	public Integer getAbsoluteNAPosition() {
		return getGenePosition().getAbsoluteNAPosition();
	}

	public Map<String, Object> extMap() {
		Map<String, Object> result = new LinkedHashMap<>();
		Map<String, Object> geneMap = new LinkedHashMap<>();
		geneMap.put("name", gene.getAbstractGene());
		geneMap.put("nameWithStrain", gene.getName());
		result.put("gene", geneMap);
		result.put("position", position);
		result.put("totalReads", totalReads);
		result.put("isTrimmed", isTrimmed);
		return result;
	}
	
}