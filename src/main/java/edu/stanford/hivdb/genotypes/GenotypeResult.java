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

import edu.stanford.hivdb.viruses.Virus;

public class GenotypeResult<VirusT extends Virus<VirusT>> {

	private List<BoundGenotype<VirusT>> genotypes;

	public GenotypeResult(List<BoundGenotype<VirusT>> genotypes) {
		genotypes.sort((g1, g2) -> g1.getDistance().compareTo(g2.getDistance()));
		this.genotypes = genotypes;
	}

	public BoundGenotype<VirusT> getFirstMatch() {
		return genotypes.get(0);
	}

	public BoundGenotype<VirusT> getFallbackMatch() {
		return genotypes.stream()
			.filter(bg -> !bg.getGenotype().hasParentGenotypes())
			.findFirst().get();
	}

	public List<BoundGenotype<VirusT>> getAllMatches() {
		return Collections.unmodifiableList(genotypes);
	}

	public BoundGenotype<VirusT> getBestMatch() {
		BoundGenotype<VirusT> first = getFirstMatch();
		BoundGenotype<VirusT> fallback = getFallbackMatch();

		return first.shouldFallbackTo(fallback) ? fallback : first;
	}

}
