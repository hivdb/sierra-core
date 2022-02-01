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
package edu.stanford.hivdb.genotypes;

import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;

import edu.stanford.hivdb.viruses.Virus;

public class GenotypeResult<VirusT extends Virus<VirusT>> {

	private List<BoundGenotype<VirusT>> genotypes;

	public GenotypeResult(List<BoundGenotype<VirusT>> genotypes) {
		genotypes.sort((g1, g2) -> {
			int cmp = g1.getDistance().compareTo(g2.getDistance());
			if (cmp == 0) {
				List<Genotype<VirusT>> g1Parents = g1.getParentGenotypes();
				List<Genotype<VirusT>> g2Parents = g2.getParentGenotypes();
				// child genotype should be placed before parents
				if (g1Parents != null && g1Parents.contains(g2.getGenotype())) {
					cmp = -1;
				}
				else if (g2Parents != null && g2Parents.contains(g1.getGenotype())) {
					cmp = 1;
				}
			}
			return cmp;
		});
		this.genotypes = genotypes;
	}

	public BoundGenotype<VirusT> getFirstMatch() {
		if (genotypes.isEmpty()) {
			return null;
		}
		return genotypes.get(0);
	}

	public BoundGenotype<VirusT> getParentFallbackMatch(BoundGenotype<VirusT> firstBg) {
		try {
			return genotypes.stream()
				.filter(bg -> (
					bg != firstBg &&
					firstBg.getGenotype().hasParentGenotypes() &&
					firstBg.getGenotype().getParentGenotypes().contains(bg.getGenotype())
				))
				.findFirst().get();
		} catch (NoSuchElementException e) {
			return null;
		} catch (NullPointerException e) {
			throw new NullPointerException(
				String.format(
					"Genotyping reference %s doesn't have a valid genotype, " +
					"which triggers following NullPointerException: %s",
					firstBg.getReference().getAccession(),
					e
				)
			);
		}
	}

	public BoundGenotype<VirusT> getChildFallbackMatch(BoundGenotype<VirusT> firstBg) {
		try {
			return genotypes.stream()
				.filter(bg -> (
					bg != firstBg &&
					bg.getGenotype().hasParentGenotypes() &&
					bg.getGenotype().getParentGenotypes().contains(firstBg.getGenotype())
				))
				.findFirst().get();
		} catch (NoSuchElementException e) {
			return null;
		}
	}
	
	public List<BoundGenotype<VirusT>> getAllMatches() {
		return Collections.unmodifiableList(genotypes);
	}

	public BoundGenotype<VirusT> getBestMatch() {
		BoundGenotype<VirusT> first = getFirstMatch();
		BoundGenotype<VirusT> parentFb = getParentFallbackMatch(first);
		BoundGenotype<VirusT> childFb = getChildFallbackMatch(first);
		
		if (first == null) {
			return null;
		}
		
		return first.useOrFallbackTo(parentFb, childFb);
	}

}