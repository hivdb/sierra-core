package edu.stanford.hivdb.sequences;

import java.util.List;

import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Virus;

public interface SequenceValidator<VirusT extends Virus<VirusT>> {
	
	public List<ValidationResult> validate(AlignedSequence<VirusT> alignedSeq);

}
