/*

    Copyright (C) 2017-2020 Stanford HIVDB team

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

package edu.stanford.hivdb.comments;

import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.viruses.Virus;

public enum CommentType {
	// TODO: change this enum to class; the enum elements should be decided by virus
	Major, Accessory, NRTI, NNRTI, Uncertain, Other, Dosage;

	public static CommentType fromMutType(MutationType<? extends Virus<?>> mutType) {
		return CommentType.valueOf(mutType.name());
	}
}
