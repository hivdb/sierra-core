package edu.stanford.hivdb.sequences;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class UnsequencedRegions<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {

	public static class UnsequencedRegion {
		private final Long posStart;
		private final Long posEnd;

		protected UnsequencedRegion(long posStart, long posEnd) {
			this.posStart = posStart;
			this.posEnd = posEnd;
		}

		public Long getPosStart() { return posStart; }

		public Long getPosEnd() { return posEnd; }

		public Long getSize() { return posEnd - posStart + 1; }

		public Long size() { return getSize(); }

		@Override
		public String toString() {
			return String.format("%d-%d", posStart, posEnd);
		}


	}

	private static <T extends Virus<T>> List<UnsequencedRegion> findUnseqRegionsFromMutations(Gene<T> gene, MutationSet<T> mutations) {
		List<UnsequencedRegion> regions = new ArrayList<>();
		MutationSet<T> unseqs = mutations.getGeneMutations(gene).filterBy(Mutation::isUnsequenced);
		Long[] positions = (
			unseqs.getPositions().stream().map(gp -> gp.getPosition().longValue())
			.toArray(Long[]::new)
		);
		for (int i = 0; i < positions.length; i ++) {
			long posStart = positions[i];
			long posEnd = posStart;
			for (int j = i + 1; j < positions.length; j ++) {
				long newPosEnd = positions[j];
				if (newPosEnd - posEnd == 1) {
					posEnd = newPosEnd;
					i ++;
				}
				else {
					break;
				}
			}
			regions.add(new UnsequencedRegion(posStart, posEnd));
		}
		return regions;
	}

	private final Gene<VirusT> gene;
	private final List<UnsequencedRegion> regions;

	public UnsequencedRegions(Gene<VirusT> gene, long firstAA, long lastAA, MutationSet<VirusT> mutations) {
		List<UnsequencedRegion> regions = new ArrayList<>();
		if (firstAA > 1) {
			regions.add(new UnsequencedRegion(1, firstAA - 1));
		}
		regions.addAll(findUnseqRegionsFromMutations(gene, mutations));
		if (lastAA < gene.getAASize()) {
			regions.add(new UnsequencedRegion(lastAA + 1, gene.getAASize()));
		}
		this.gene = gene;
		this.regions = Collections.unmodifiableList(regions);
	}

	@Override
	public Gene<VirusT> getGene() { return gene; }

	public List<UnsequencedRegion> getRegions() { return regions; }

	public Long getSize() {
		return regions.stream().map(r -> r.getSize()).reduce(0L, (a, b) -> a + b);
	}

	public Long size() {
		return getSize();
	}

}
