package edu.stanford.hivdb.sequences;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.LongStream;

import com.google.common.collect.Streams;

import edu.stanford.hivdb.mutations.CodonMutation;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.utilities.Json;
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
		
		public boolean isUnsequenced(long position) {
			return posStart <= position && position <= posEnd;
		}

		@Override
		public String toString() {
			return String.format("%d-%d", posStart, posEnd);
		}


	}

	private static <T extends Virus<T>> List<UnsequencedRegion> findUnseqRegionsFromMutations(Gene<T> gene, long firstAA, long lastAA, MutationSet<T> mutations) {
		List<UnsequencedRegion> regions = new ArrayList<>();
		MutationSet<T> geneMuts = mutations.getGeneMutationsNoSplit(gene); 
		final MutationSet<T> definitiveUnseqs = geneMuts.filterByNoSplit(Mutation::isUnsequenced);
		final MutationSet<T> conditionalUnseqs = geneMuts.filterByNoSplit(mut -> {
			// mutation is considered "conditionally unsequenced"
			// if its "N" part is adjacent to a "definitively unsequenced" mutation
			if (mut.isUnsequenced()) {
				return false;
			}
			if (mut.isInsertion()) {
				return false;
			}
			if (mut instanceof CodonMutation) {
				String triplet = mut.getTriplet();
				int pos = mut.getPosition();
				if (triplet.endsWith("N") && (definitiveUnseqs.get(gene, pos + 1) != null || pos + 1 > lastAA)) {
					// right N adjacent to a "definitively unsequenced" mutation
					return true;
				}
				else if (triplet.startsWith("N") && (definitiveUnseqs.get(gene, pos - 1) != null || pos - 1 < firstAA)) {
					// left N adjacent to a "definitively unsequenced" mutation
					return true;
				}
			}
			return false;
		});
		MutationSet<T> unseqs = definitiveUnseqs.mergesWith(conditionalUnseqs);
	
		long[] positions = Streams.concat(
			LongStream.range(1, firstAA),
			unseqs.getPositions().stream().mapToLong(gp -> gp.getPosition().longValue()),
			LongStream.range(lastAA + 1, gene.getAASize() + 1)
		).toArray();
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
			if (posStart == firstAA) {
				posStart = 1;
			}
			if (posEnd == lastAA) {
				posEnd = gene.getAASize();
			}
			regions.add(new UnsequencedRegion(posStart, posEnd));
		}
		return regions;
	}

	private final Gene<VirusT> gene;
	private final List<UnsequencedRegion> regions;

	public UnsequencedRegions(Gene<VirusT> gene, long firstAA, long lastAA, MutationSet<VirusT> mutations) {
		this.gene = gene;
		this.regions = Collections.unmodifiableList(
			findUnseqRegionsFromMutations(gene, firstAA, lastAA, mutations)
		);
	}

	@Override
	public Gene<VirusT> getGene() { return gene; }

	public List<UnsequencedRegion> getRegions() { return regions; }

	public boolean isUnsequenced(Gene<VirusT> gene, long position) {
		if (gene != this.gene) {
			throw new RuntimeException("The input gene doesn't match these unsequenced regions.");
		}
		for (UnsequencedRegion region : regions) {
			if (region.isUnsequenced(position)) {
				return true;
			}
		}
		return false;
	}
	
	public Long getSize() {
		return regions.stream().map(r -> r.getSize()).reduce(0L, (a, b) -> a + b);
	}

	public Long size() {
		return getSize();
	}

}
