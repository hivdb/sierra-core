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

import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

import org.apache.commons.lang3.builder.EqualsBuilder;

public class FrameShift<VirusT extends Virus<VirusT>> implements Comparable<FrameShift<VirusT>>, WithGene<VirusT> {
	public static enum Type {INSERTION, DELETION};
	private Gene<VirusT> gene;
	private int position;
	private int size;
	private String nas;
	private Type type;

	public static <VirusT extends Virus<VirusT>> String joinFrameShifts(List<FrameShift<VirusT>> frameShifts) {
		StringBuilder output = new StringBuilder();
		if (frameShifts.size() == 0) {
			output.append("None");
		} else {
			for (FrameShift<VirusT> fs : frameShifts) {
				output.append(fs.getHumanFormat() + ", ");
			}
			output.setLength(output.length() - 2);
		}
		return output.toString();
	}

	private FrameShift(Gene<VirusT> gene, int position, int size, String nas, Type type) {
		this.gene = gene;
		this.position = position;
		this.nas = nas;
		this.type = type;
		this.size = Math.abs(size);
	}

	public int compareTo (FrameShift<VirusT> fs) {
		return new Integer(position).compareTo(new Integer(fs.position));
	}

	public boolean isInsertion() {
		return type == Type.INSERTION;
	}

	public boolean isDeletion() {
		return type == Type.DELETION;
	}


	public String getHumanFormat() {
		String output;
		if (getType() == Type.INSERTION) {
			output = String.format(
				"%s%dins%dbp_%s", gene.getName(), position, size, nas
			);
		} else {
			output = String.format(
				"%s%ddel%dbp", gene.getName(), position, size
			);
		}
		return output;
	}

	public static <VirusT extends Virus<VirusT>> FrameShift<VirusT> createDeletion(Gene<VirusT> gene, int aaPosition, int size) {
		return new FrameShift<>(gene, aaPosition, size, "", Type.DELETION);
	}

	public static <VirusT extends Virus<VirusT>> FrameShift<VirusT> createInsertion(Gene<VirusT> gene, int aaPosition, int size, String nas) {
		return new FrameShift<>(gene, aaPosition, size, nas, Type.INSERTION);
	}

	public static <VirusT extends Virus<VirusT>> FrameShift<VirusT> fromNucAminoFrameShift(Gene<VirusT> gene, int aaStart, Map<?, ?> fs) {
		return new FrameShift<>(
			gene,
			/* aaPosition */ ((Double) fs.get("Position")).intValue() - aaStart + 1,
			/* size */ ((Double) fs.get("GapLength")).intValue(),
			/* nas */ (String) fs.get("NucleicAcidsText"),
			/* type */ (Boolean) fs.get("IsInsertion") ? Type.INSERTION : Type.DELETION
		);
	}

	@Override
	public Strain<VirusT> getStrain() { return gene.getStrain(); }

	@Override
	public Gene<VirusT> getGene() {return gene; }

	@Override
	public String getAbstractGene() { return gene.getAbstractGene(); }

	public int getPosition() { return position; }
	public Type getType() { return type; }
	public int getSize() { return size; }
	public String getNAs() { return nas; }

	public String getText() {
		// same as `toString()`, only for GraphQL
		return getHumanFormat();
	}

	@Override
	public String toString() {
		return getHumanFormat();
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(6434, 15675) // any two random prime numbers
			.append(gene)
			.append(position)
			.append(type)
			.append(size)
			.append(nas)
			.toHashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof FrameShift)) {
			return false;
		}
		FrameShift<?> other = (FrameShift<?>) obj;
		return new EqualsBuilder()
			.append(gene, other.getGene())
			.append(position, other.getPosition())
			.append(type, other.getType())
			.append(size, other.getSize())
			.append(nas, other.getNAs())
			.isEquals();
	}

}
