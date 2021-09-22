package edu.stanford.hivdb.viruses;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;

public abstract class AssemblyRegion<
	VirusT extends Virus<VirusT>,
	GeneSeqT extends WithGene<VirusT>,
	T extends AssemblyRegion<VirusT, GeneSeqT, T>
> implements WithGene<VirusT> {

	public static enum AssemblyRegionType {
		UNTRANS_REGION, GENE
	}

	private final String name;
	private final AssemblyRegionType type;
	private final Long refStart;
	private final Long refEnd;
	private final Gene<VirusT> gene;
	private final List<Long> trim;

	public static <
		VirusT extends Virus<VirusT>,
		GeneSeqT extends WithGene<VirusT>,
		T extends AssemblyRegion<VirusT, GeneSeqT, T>
	>
	T of(
		final Class<T> clazz,
		final String name,
		final AssemblyRegionType type,
		final Long refStart,
		final Long refEnd,
		final Gene<VirusT> gene,
		final List<Long> trim
	) {
		try {
			T instance = clazz
				.getDeclaredConstructor(
					String.class,
					AssemblyRegionType.class,
					Long.class,
					Long.class,
					Gene.class,
					List.class
				)
				.newInstance(name, type, refStart, refEnd, gene, trim);
			return instance;
		} catch (
			InstantiationException
			| IllegalAccessException
			| IllegalArgumentException
			| InvocationTargetException
			| NoSuchMethodException
			| SecurityException e
		) {
			throw new RuntimeException(e);
		}
	}
	
	public AssemblyRegion(
		final String name,
		final AssemblyRegionType type,
		final Long refStart,
		final Long refEnd,
		final Gene<VirusT> gene,
		final List<Long> trim
	) {
		this.name = name;
		this.type = type;
		this.refStart = refStart;
		this.refEnd = refEnd;
		this.gene = gene;
		this.trim = trim;
	}
	
	@Override
	public Gene<VirusT> getGene() { return gene; }
	
	public String getName() { return name; }
	public AssemblyRegionType getType() { return type; }
	public Long getRefStart() { return refStart; }
	public Long getRefEnd() { return refEnd; }

	public List<String> trimCodonList(List<String> codonList) { 
		List<Long> trim = this.trim == null ? Lists.newArrayList() : new ArrayList<>(this.trim);
		trim.sort(Long::compare);
		for (long napos0 : Lists.reverse(trim)) {
			napos0 --;
			int aapos0 = Math.toIntExact(napos0 / 3);
			int bp0 = Math.toIntExact(napos0 % 3);
			String codon = codonList.get(aapos0);
			codon = codon.substring(0, bp0) + codon.substring(bp0 + 1);
			codonList.set(aapos0, codon);
		}
		return codonList;
	}
	
	public String assemble(UntranslatedRegion utr, boolean skipIns) {
		if (type != AssemblyRegionType.UNTRANS_REGION) {
			throw new RuntimeException(
				"UntranslatedRegion object passed to Non-UTR AssemblyRegion " + name
			);
		}
		if (utr == null) {
			return Strings.repeat("N", Math.toIntExact(refEnd - refStart + 1));
		}
		int left = Math.toIntExact(refStart - utr.getRefStart());
		int right = Math.toIntExact(utr.getRefEnd() - refEnd);
		int trimLeft = left < 0 ? - left : 0;
		int trimRight = right < 0 ? - right : 0;
		int addLeft = left > 0 ? left : 0;
		int addRight = right > 0 ? right : 0;
		String consensus = utr.getConsensus();
		return (
			Strings.repeat("N", addLeft) +
			(skipIns ?
				consensus.substring(
					trimLeft,
					Math.toIntExact(trimLeft + utr.getRefEnd() - utr.getRefStart() + 1)
				) :
				consensus.substring(trimLeft, consensus.length() - trimRight)
			) +
			Strings.repeat("N", addRight)
		);
	}

	public abstract String assemble(GeneSeqT geneSeq, boolean skipIns);
}