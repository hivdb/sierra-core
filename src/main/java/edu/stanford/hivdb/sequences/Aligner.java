package edu.stanford.hivdb.sequences;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.stanford.hivdb.viruses.Virus;

public interface Aligner<VirusT extends Virus<VirusT>> {

	public final static Map<String, Aligner<?>> singletons = new HashMap<>(); 

	public static class MisAlignedException extends IllegalArgumentException {
		/**
		 *
		 */
		private static final long serialVersionUID = 46495128315347L;
		private final boolean suppressible;

		public MisAlignedException(String message, boolean suppressible) {
			super(message);
			this.suppressible = suppressible;
		}

		public boolean isSuppressible() { return suppressible; }
	}

	
	@SuppressWarnings("unchecked")
	public static <AlignerT extends Aligner<VirusT>, VirusT extends Virus<VirusT>> AlignerT getInstance(VirusT virusIns) {
		String className = virusIns.getClass().getName();
		if (!singletons.containsKey(className)) {
			AlignmentConfig<VirusT> alignConfig = virusIns.getAlignmentConfig();
			switch(alignConfig.getMethod()) {
			case NucAmino:
				singletons.put(className, new NucAminoAligner<>(virusIns));
				break;
			case PostAlign:
				singletons.put(className, new PostAlignAligner<>(virusIns));
				break;
			default:
				throw new NullPointerException("alignmentConfig.method can not be null");
			}
		}
		return (AlignerT) singletons.get(className);
	}

	/**
	 * Receives a sequence and aligns it to each gene by given method.
	 *
	 * @param sequence	Sequence waiting to be aligned.
	 * @return 			an AlignedSequence object
	 */
	public AlignedSequence<VirusT> align(Sequence sequence);

	/**
	 * Receives set of sequences and aligns them to each gene in parallel.
	 *
	 * @param sequences		Sequence list waiting to be aligned
	 * @return 				list of AlignedSequence objects
	 */
	public List<AlignedSequence<VirusT>> parallelAlign(Collection<Sequence> sequences);
	
}
