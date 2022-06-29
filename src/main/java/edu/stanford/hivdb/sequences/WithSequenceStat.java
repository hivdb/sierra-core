package edu.stanford.hivdb.sequences;

import java.util.List;

import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.viruses.Virus;

public interface WithSequenceStat<VirusT extends Virus<VirusT>> {
	
	public MutationSet<VirusT> getMutations();
	
	public List<FrameShift<VirusT>> getFrameShifts();

	default public Long getMutationCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced());
	}

	default public Long getUnusualMutationCount() {
		return getMutations().countIf(mut -> mut.isUnusual() && !mut.isUnsequenced());
	}
	
	default public Long getFrameShiftCount() {
		return Long.valueOf(getFrameShifts().size());
	}
	
	default public Long getInsertionCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced() && mut.isInsertion());
	}
	
	default public Long getDeletionCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced() && mut.isDeletion());
	}

	default public Long getStopCodonCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced() && mut.hasStop());
	}

	default public Long getAmbiguousMutationCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced() && mut.isAmbiguous());
	}
	
	default public Long getApobecMutationCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced() && mut.isApobecMutation());
	}

	default public Long getApobecDRMCount() {
		return getMutations().countIf(mut -> !mut.isUnsequenced() && mut.isApobecDRM());
	}

}
