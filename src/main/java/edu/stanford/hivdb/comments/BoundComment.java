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

	protected BoundComment(
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