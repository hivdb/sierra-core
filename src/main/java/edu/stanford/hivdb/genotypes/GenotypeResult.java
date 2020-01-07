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
