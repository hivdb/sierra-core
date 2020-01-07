package edu.stanford.hivdb.seqreads;

import java.util.List;

import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Virus;

public interface SequenceReadsValidator<VirusT extends Virus<VirusT>> {
	
	public List<ValidationResult> validate(SequenceReads<VirusT> seqReads);

}
