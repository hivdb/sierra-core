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

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class AminoAcidPercent<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	final private Gene<VirusT> gene;
	final private Integer position;
	final private Character aa;
	final private Double percent;
	final private Integer count;
	final private Integer total;
	final private String reason;
	final private Boolean isUnusual;
	
	protected AminoAcidPercent (
		Gene<VirusT> gene, int position, char aa, double percent,
		int count, int total, String reason, boolean isUnusual
	) {
		this.gene = gene;
		this.position = position;
		this.aa = aa;
		this.percent = percent;
		this.count = count;
		this.total = total;
		this.reason = reason;
		this.isUnusual = isUnusual;
	}

	@Override
	public Strain<VirusT> getStrain() {
		return gene.getStrain();
	}

	@Override
	public Gene<VirusT> getGene() {
		return gene;
	}
	
	@Override
	public String getAbstractGene() {
		return gene.getAbstractGene();
	}
	
	public Integer getPosition() {
		return position;
	}
	
	public Character getAA() {
		return aa;
	}
	
	public Character getRefChar() {
		return gene.getRefChar(position);
	}
	
	public Mutation<VirusT> getMutation() {
		return new AAMutation<>(gene, position, aa);
	}
	
	public Double getPercent() {
		return percent;
	}
	
	public Integer getCount() {
		return count;
	}
	
	public Integer getTotal() {
		return total;
	}
	
	public String getReason() {
		return reason;
	}
	
	public Boolean isUnusual() {
		return isUnusual;
	}

	public GenePosition<VirusT> getGenePosition() {
		return new GenePosition<>(getGene(), position);
	}

}
