/*

    Copyright (C) 2019-2020 Stanford HIVDB team

    This file is part of Sierra.

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sierra.  If not, see <https://www.gnu.org/licenses/>.
*/
package edu.stanford.hivdb.seqreads;

import java.util.Collection;
import java.util.List;

import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Virus;

public interface SequenceReadsValidator<VirusT extends Virus<VirusT>> {
	
	public List<ValidationResult> validate(SequenceReads<VirusT> seqReads, Collection<String> includeGenes);

}
