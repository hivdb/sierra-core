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

import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

/**
 * Helper class mainly used to build mutation search index.
 *
 * Many mutation-related indices use gene and position as their index
 * key. This class instantiates hashable and comparable objects using
 * value gene and pos.
 */
public class GenePosition<VirusT extends Virus<VirusT>> implements Comparable<GenePosition<VirusT>>, WithGene<VirusT> {

	private final Gene<VirusT> gene;
	private final Integer position;

	private transient Boolean isDrugResistancePosition;

	public static <VirusT extends Virus<VirusT>> Set<GenePosition<VirusT>> getGenePositionsBetween(
		final GenePosition<VirusT> start,
		final GenePosition<VirusT> end
	) {
		Gene<VirusT> startGene = start.getGene();
		Gene<VirusT> endGene = end.getGene();
		if (startGene.getStrain() != endGene.getStrain()) {
			throw new IllegalArgumentException(
				"Virus strain of `start` and `end` positions must be the same."
			);
		}
		Strain<VirusT> strain = startGene.getStrain();
		Set<GenePosition<VirusT>> genePositions = new LinkedHashSet<>();
		for (Gene<VirusT> gene : strain.getGenes()) {
			int startPos = 1;
			int endPos = gene.getAASize();
			if (gene.compareTo(startGene) == 0) {
				startPos = start.getPosition();
			}
			if (gene.compareTo(endGene) == 0) {
				endPos = end.getPosition();
			}
			if (gene.compareTo(startGene) >= 0 &&
				gene.compareTo(endGene) <= 0) {
				genePositions.addAll(
					gene.getGenePositionsBetween(startPos, endPos));
			}
		}
		return genePositions;
	}

	public static <VirusT extends Virus<VirusT>> Set<GenePosition<VirusT>> getGenePositionsBetween(
		final GenePosition<VirusT> start,
		final GenePosition<VirusT> end,
		final Collection<String> includeGenes
	) {
		return getGenePositionsBetween(start, end).stream()
			.filter(gp -> includeGenes.contains(gp.getAbstractGene()))
			.collect(Collectors.toCollection(LinkedHashSet::new));
	}

	public static <VirusT extends Virus<VirusT>> Set<GenePosition<VirusT>> getDRGenePositionsBetween(
		final GenePosition<VirusT> start,
		final GenePosition<VirusT> end
	) {
		return (
			getGenePositionsBetween(start, end)
			.stream()
			.filter(gp -> gp.isDrugResistancePosition())
			.collect(Collectors.toCollection(LinkedHashSet::new))
		);
	}

	public static <VirusT extends Virus<VirusT>> Set<GenePosition<VirusT>> getDRGenePositionsBetween(
		final GenePosition<VirusT> start,
		final GenePosition<VirusT> end,
		final Collection<String> includeGenes
	) {
		return (
			getGenePositionsBetween(start, end)
			.stream()
			.filter(gp -> gp.isDrugResistancePosition() && includeGenes.contains(gp.getAbstractGene()))
			.collect(Collectors.toCollection(LinkedHashSet::new))
		);
	}

	public GenePosition(final Gene<VirusT> gene, final int position) {
		this.gene = gene;
		this.position = position;
	}

	public GenePosition(final Gene<VirusT> gene, final Integer position) {
		this.gene = gene;
		this.position = position;
	}

	public Integer getPosition() {
		return position;
	}
	
	public Integer getPositionInStrain() {
		int absPos = 0;
		Strain<VirusT> strain = gene.getStrain();
		for (Gene<VirusT> prevGene : strain.getGenes()) {
			if (prevGene == gene) {
				absPos += position;
				break;
			}
			else {
				absPos += prevGene.getAASize();
			}
		}
		return absPos;
	}
	
	public boolean isDrugResistancePosition() {
		if (isDrugResistancePosition == null) {
			isDrugResistancePosition = (
				gene.getVirusInstance().getDrugResistMutations()
				.values().stream()
				.anyMatch(drms -> drms.get(this) != null)
			);
		}
		return isDrugResistancePosition;
	}

	@Override
	public Gene<VirusT> getGene() {
		return gene;
	}
	
	@Override
	public Strain<VirusT> getStrain() {
		return gene.getStrain();
	}
	
	@Override
	public String getAbstractGene() {
		return gene.getAbstractGene();
	}
	
	public List<MutationType<VirusT>> getMutationTypes() {
		return gene.getVirusInstance().getMutationTypePairs()
			.stream()
			.filter(mc -> mc.getGenePosition().equals(this))
			.map(mc -> mc.getMutationType())
			.collect(Collectors.toList());
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) { return false; }
		if (obj == this) { return true; }
		if (obj.getClass() != getClass()) { return false; }
		GenePosition<?> gp = (GenePosition<?>) obj;
		return new EqualsBuilder()
			.append(getGene(), gp.getGene())
			.append(getPosition(), gp.getPosition())
			.isEquals();
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(63261, 362788935)
			.append(gene).append(position).toHashCode();
	}

	@Override
	public int compareTo(GenePosition<VirusT> o) {
		if (o == null) throw new NullPointerException("Null is incomprable.");
		int cmp = getGene().compareTo(o.getGene());
		if (cmp == 0) {
			cmp = getPosition().compareTo(o.getPosition());
		}
		return cmp;
	}

	
	public String toStringWithAbstractGene() {
		return String.format("%s:%d", getAbstractGene(), getPosition());
	}
	
	@Override
	public String toString() {
		return String.format("%s:%d", getGene(), getPosition());
	}
}
