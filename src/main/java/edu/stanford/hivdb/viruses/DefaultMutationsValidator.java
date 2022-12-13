/*

    Copyright (C) 2022 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.viruses;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationsValidator;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;

public class DefaultMutationsValidator<VirusT extends Virus<VirusT>> implements MutationsValidator<VirusT> {

	@Override
	public List<ValidationResult> validate(MutationSet<VirusT> mutations, Collection<String> includeGenes) {
		List<ValidationResult> validationResults = new ArrayList<>();
		validationResults.addAll(validateNoStopCodons(mutations, includeGenes));
		return validationResults;
	}

	private static <VirusT extends Virus<VirusT>> List<ValidationResult> validateNoStopCodons(
		MutationSet<VirusT> mutations,
		Collection<String> includeGenes
	) {
		List<ValidationResult> validationResults = new ArrayList<>();
		MutationSet<VirusT> stopCodons = mutations
			.getStopCodons()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()));
		for (Map.Entry<Gene<VirusT>, MutationSet<VirusT>> entry : stopCodons.groupByGene().entrySet()) {
			String geneDisplay = entry.getKey().getDisplay();
			MutationSet<VirusT> geneStopCodons = entry.getValue();
			int numGeneStopCodons = geneStopCodons.size();
			String geneStopText = geneStopCodons.join(", ", Mutation::getHumanFormatWithAbstractGene);
			if (numGeneStopCodons > 1) {
				validationResults.add(DefaultValidationMessage.MultipleStopCodons.formatWithLevel(
					ValidationLevel.SEVERE_WARNING,
					numGeneStopCodons,
					geneDisplay,
					geneStopText
				));
			} else if (numGeneStopCodons > 0) {
				validationResults.add(DefaultValidationMessage.SingleStopCodon.formatWithLevel(
					ValidationLevel.WARNING,
					geneDisplay,
					geneStopText
				));
			}
		}
		
		return validationResults;
	}

}