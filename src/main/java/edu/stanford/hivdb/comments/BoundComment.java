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
package edu.stanford.hivdb.comments;

import java.util.Collection;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class BoundComment<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	final private Strain<VirusT> strain;
	final private DrugClass<VirusT> drugClass;
	final private String commentName;
	final private CommentType commentType;
	final private String comment;
	final private Collection<String> highlightText;
	final private Mutation<VirusT> mutation;

	public BoundComment(
			Strain<VirusT> strain, String commentName,
			DrugClass<VirusT> drugClass, CommentType commentType,
			String comment, Collection<String> highlightText,
			Mutation<VirusT> mutation) {
		this.strain = strain;
		this.drugClass = drugClass;
		this.commentName = commentName;
		this.commentType = commentType;
		this.comment = comment;
		this.highlightText = highlightText;
		this.mutation = mutation;
	}

	public String getName() { return commentName; }
	public CommentType getType() { return commentType; }
	public String getText() { return comment; }
	public Collection<String> getHighlightText() { return highlightText; }
	public Mutation<VirusT> getBoundMutation() { return mutation; }

	@Override
	public Gene<VirusT> getGene() { return strain.getGene(drugClass.getAbstractGene()); }
	public DrugClass<VirusT> drugClass() { return drugClass; }
}