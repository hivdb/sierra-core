package edu.stanford.hivdb.mutations;

import java.util.List;

import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Virus;

public interface MutationsValidator<VirusT extends Virus<VirusT>> {
	
	public List<ValidationResult> validate(MutationSet<VirusT> mutations);

}
