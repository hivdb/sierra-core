package edu.stanford.hivdb.sequences;

import java.util.List;

import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.viruses.Virus;

public interface WithSequenceStat<VirusT extends Virus<VirusT>> {
	
	public MutationSet<VirusT> getSequencedMutations();
	
	public List<FrameShift<VirusT>> getFrameShifts();
	
	default public Long getMutationCount() {
		return getSequencedMutations().count();
	}

	default public Long getUnusualMutationCount() {
		return getSequencedMutations().countIf(Mutation::isUnusual);
	}
	
	default public Long getFrameShiftCount() {
		return Long.valueOf(getFrameShifts().size());
	}
	
	default public Long getInsertionCount() {
		return getSequencedMutations().countIf(Mutation::isInsertion);
	}
	
	default public Long getDeletionCount() {
		return getSequencedMutations().countIf(Mutation::isDeletion);
	}

	default public Long getStopCodonCount() {
		return getSequencedMutations().countIf(Mutation::hasStop);
	}

	default public Long getAmbiguousMutationCount() {
		return getSequencedMutations().countIf(Mutation::isAmbiguous);
	}
	
	default public Long getApobecMutationCount() {
		return getSequencedMutations().countIf(Mutation::isApobecMutation);
	}

	default public Long getApobecDRMCount() {
		return getSequencedMutations().countIf(Mutation::isApobecDRM);
	}

}
