package edu.stanford.hivdb.sequences;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.Streams;

import edu.stanford.hivdb.mutations.CodonMutation;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class GeneRegions<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {

	public static class GeneRegion {
		private final Long posStart;
		private final Long posEnd;

		protected GeneRegion(long posStart, long posEnd) {
			this.posStart = posStart;
			this.posEnd = posEnd;
		}

		public Long getPosStart() { return posStart; }

		public Long getPosEnd() { return posEnd; }

		public Long getSize() { return posEnd - posStart + 1; }

		public Long size() { return getSize(); }
		
		public boolean contains(long position) {
			return posStart <= position && position <= posEnd;
		}

		@Override
		public String toString() {
			if (posStart.equals(posEnd)) {
				return String.format("%d", posStart);
			}
			return String.format("%d-%d", posStart, posEnd);
		}

	}
	protected static <T extends Virus<T>> GeneRegions<T> unsafeNewGeneRegions(Gene<T> gene, long[] positions) {
		List<GeneRegion> regions = new ArrayList<>();
		int size = positions.length;
		for (int i = 0; i < size; i ++) {
			long posStart = positions[i];
			long posEnd = posStart;
			for (int j = i + 1; j < size; j ++) {
				long newPosEnd = positions[j];
				if (newPosEnd - posEnd == 1) {
					posEnd = newPosEnd;
					i ++;
				}
				else {
					break;
				}
			}
			regions.add(new GeneRegion(posStart, posEnd));
		}
		return new GeneRegions<>(gene, regions);
	}
	
	
	public static <T extends Virus<T>> GeneRegions<T> newGeneRegions(Gene<T> gene, long posStart, long posEnd) {
		List<GeneRegion> regions = new ArrayList<>();
		regions.add(new GeneRegion(posStart, posEnd));
		return new GeneRegions<>(gene, regions);
	}

	public static <T extends Virus<T>> GeneRegions<T> newGeneRegions(Gene<T> gene, List<Number> positions) {
		return GeneRegions.unsafeNewGeneRegions(
			gene,
			positions.stream()
				.mapToLong(pos -> pos.longValue())
				.toArray()
		);
	}
	
	public static <T extends Virus<T>> List<GeneRegions<T>> newListOfGeneRegions(List<GenePosition<T>> genePos) {
		return genePos
			.stream()
			.collect(Collectors.groupingBy(
				GenePosition::getGene,
				TreeMap::new,
				Collectors.toCollection(TreeSet::new)
			))
			.entrySet()
			.stream()
			.map(e -> GeneRegions.unsafeNewGeneRegions(
				e.getKey(),
				e.getValue().stream()
					.mapToLong(gp -> gp.getPosition())
					.toArray()
			))
			.collect(Collectors.toList());
	}
	
	public static <T extends Virus<T>> GeneRegions<T> newUnsequencedRegions(Gene<T> gene, long firstAA, long lastAA, MutationSet<T> mutations) {
		MutationSet<T> geneMuts = mutations.getGeneMutationsNoSplit(gene); 
		final MutationSet<T> definitiveUnseqs = geneMuts.filterByNoSplit(Mutation::isUnsequenced);
		final List<Pair<Integer, Integer>> consecutiveDels = new ArrayList<>();
		for (Mutation<T> mut : mutations) {
			int pos = mut.getPosition();
			boolean breakFlag = false;
			for (int i = 0; i < consecutiveDels.size(); i ++) {
				Pair<Integer, Integer> range = consecutiveDels.get(i);
				if (pos == range.getRight()) {
					consecutiveDels.set(i, Pair.of(range.getLeft(), pos + 1));
					breakFlag = true;
					break;
				}
			}
			if (!breakFlag) {
				consecutiveDels.add(Pair.of(pos, pos + 1));
			}
		}
		final Set<Integer> unseqDels = consecutiveDels
			.stream()
			.filter(pair -> pair.getRight() - pair.getLeft() >= 20)
			.flatMapToInt(pair -> IntStream.range(pair.getLeft(), pair.getRight()))
			.mapToObj(i -> i)
			.collect(Collectors.toSet());
		
		final MutationSet<T> conditionalUnseqs = geneMuts.filterByNoSplit(mut -> {
			// mutation is considered "conditionally unsequenced"
			// if its "N" part is adjacent to a "definitively unsequenced" mutation
			if (mut.isUnsequenced()) {
				return false;
			}
			if (mut.isInsertion()) {
				return false;
			}
			if (mut.isDeletion() && unseqDels.contains(mut.getPosition())) {
				// long-stretch del are unsequenced regions
				return true;
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
		
		return unsafeNewGeneRegions(
			gene,
			Streams.concat(
				LongStream.range(1, firstAA),
				unseqs.getPositions().stream().mapToLong(gp -> gp.getPosition().longValue()),
				LongStream.range(lastAA + 1, gene.getAASize() + 1)
			).toArray()
		);
	}

	private final Gene<VirusT> gene;
	private final List<GeneRegion> regions;
	
	private GeneRegions(final Gene<VirusT> gene, final List<GeneRegion> regions) {
		this.gene = gene;
		this.regions = regions;
	}

	@Override
	public Gene<VirusT> getGene() { return gene; }

	public List<GeneRegion> getRegions() { return regions; }

	public boolean contains(long position) {
		for (GeneRegion region : regions) {
			if (region.contains(position)) {
				return true;
			}
		}
		return false;
	}

	public boolean contains(Gene<VirusT> gene, long position) {
		if (gene != this.gene) {
			throw new RuntimeException("The input gene doesn't match this GeneRegions.");
		}
		return contains(position);
	}

	public boolean contains(GenePosition<VirusT> genePos) {
		return contains(genePos.getGene(), genePos.getPosition());
	}
	
	/**
	 * Get the amino acid size of all regions.
	 * 
	 * @return int
	 */
	public Long getSize() {
		return regions.stream().map(r -> r.getSize()).reduce(0L, (a, b) -> a + b);
	}

	public Long size() {
		return getSize();
	}

	@Override
	public String toString() {
		return String.format("%s %s", getGeneDisplay(), StringUtils.join(regions, ", "));
	}

}