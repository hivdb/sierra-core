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

package edu.stanford.hivdb.mutations;

import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class CodonPercent<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	final private Gene<VirusT> gene;
	final private Integer position;
	final private String codon;
	final private Double percent;
	final private Integer count;
	final private Integer total;
	final transient private Character aa;

	protected CodonPercent(
			Gene<VirusT> gene, int position, String codon,
			double percent, int count, int total) {
		this.gene = gene;
		this.position = position;
		this.codon = codon;
		if (codon.equals("ins")) {
			this.aa = '_';
		}
		else if (codon.equals("del")) {
			this.aa = '-';
		}
		else {
			this.aa = CodonUtils.translateNATriplet(codon).charAt(0);
		}
		this.percent = percent;
		this.count = count;
		this.total = total;
	}

	@Override
	public Gene<VirusT> getGene() {
		return gene;
	}
	
	@Override
	public Strain<VirusT> getStrain() { return gene.getStrain(); }
	
	@Override
	public String getAbstractGene() { return gene.getAbstractGene(); }

	final public Integer getPosition() { return position; }
	final public String getCodon() { return codon; }
	final public Character getAA() { return aa; }
	final public Double getPercent() { return percent; }
	final public Integer getCount() { return count; }
	final public Integer getTotal() { return total; }

	public GenePosition<VirusT> getGenePosition() {
		return new GenePosition<>(gene, position);
	}
}
