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

import edu.stanford.hivdb.viruses.Virus;

/*
 * Class representing the frequency of mutations for treated and naive individuals
 * for one particular mutation and subtype.
 */
public class MutationPrevalence<VirusT extends Virus<VirusT>> {
	private final Mutation<VirusT> mutation;
	private final String subtype;
	private final Integer totalNaive;
	private final Integer totalTreated;
	private final Integer frequencyNaive;
	private final Integer frequencyTreated;
	private final Double percentageNaive;
	private final Double percentageTreated;

	private MutationPrevalence(
			Mutation<VirusT> mutation, String subtype,
			int totalNaive, int freqNaive, double pcntNaive,
			int totalTreated, int freqTreated, double pcntTreated) {
		this.mutation = mutation;
		this.subtype = subtype;
		this.totalNaive = totalNaive;
		this.totalTreated = totalTreated;
		this.frequencyNaive = freqNaive;
		this.frequencyTreated = freqTreated;
		this.percentageNaive = pcntNaive;
		this.percentageTreated = pcntTreated;
	}
	
	public MutationPrevalence(
		String subtype,
		AminoAcidPercent<VirusT> naiveAAPcnt,
		AminoAcidPercent<VirusT> artAAPcnt
	) {
		this(
			naiveAAPcnt.getMutation(), subtype,
			naiveAAPcnt.getTotal(), naiveAAPcnt.getCount(), naiveAAPcnt.getPercent(),
			artAAPcnt.getTotal(), artAAPcnt.getCount(), artAAPcnt.getPercent()
		);
	}
	
	public Mutation<VirusT> getMutation() {
		return mutation;
	}
	
	public String getSubtype() {
		return subtype;
	}
	
	public Integer getTotalNaive() {
		return totalNaive;
	}
	
	public Integer getTotalTreated() {
		return totalTreated;
	}
	
	public Integer getFrequencyNaive() {
		return frequencyNaive;
	}
	
	public Integer getFrequencyTreated() {
		return frequencyTreated;
	}
	
	public Double getPercentageNaive() {
		return percentageNaive;
	}
	
	public Double getPercentageTreated() {
		return percentageTreated;
	}

	public String getAA() {
		return mutation.getAAs();
	}

	@Override
	public String toString() {
		return String.format(
			"%s %s %d %d %d %d %f %f",
			mutation, subtype, totalNaive, totalTreated, frequencyNaive,
			frequencyTreated, percentageNaive, percentageTreated);
	}

	public boolean isRare() {
		return this.percentageNaive < 0.1 && this.percentageTreated < 0.1;
	}
}

