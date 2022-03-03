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

package edu.stanford.hivdb.mutations;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import edu.stanford.hivdb.comments.BoundComment;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.sequences.UnsequencedRegions;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public interface Mutation<VirusT extends Virus<VirusT>> extends Comparable<Mutation<VirusT>>, Cloneable, WithGene<VirusT> {

	/**
	 * Exception for invalid mutation
	 */
	public static class InvalidMutationException extends RuntimeException {
		private static final long serialVersionUID = 4271016133470715497L;

		public InvalidMutationException(String message, Exception e) {
			super(message, e);
		}

		public InvalidMutationException(String message) {
			super(message);
		}
	}

	/**
	 * Compare two mutations by gene, position and aas.
	 * 
	 * @param <VirusT>
	 * @param mutA
	 * @param mutB
	 * @return int
	 */
	public static <VirusT extends Virus<VirusT>> int compare(Mutation<VirusT> mutA, Mutation<VirusT> mutB) {
		return mutA.compareTo(mutB);
	}

	// /**
	//  * Split a mixture into a set of mutations.
	//  *
	//  * @return A set of mutations
	//  */
	// public Set<Mutation> split();

	/**
	 * Merges with another Mutation object and returns the merged mutation.
	 *
	 * @param another Mutation
	 * 
	 * @return A new merged Mutation object
	 * 
	 * @throws IllegalArgumentException if another is null or the gene/position didn't match.
	 */
	@Deprecated
	public Mutation<VirusT> mergesWith(Mutation<VirusT> another);

	/**
	 * Merges with another Collection&lt;Character&gt; of
	 * AAs and returns the merged mutation.
	 *
	 * @param otherAAChars Amino acid characters
	 * 
	 * @return A new merged Mutation object
	 */
	public Mutation<VirusT> mergesWith(Collection<Character> otherAAChars);


	/**
	 * Subtracts by another Mutation object and returns the result mutation.
	 *
	 * @param another Mutation
	 * 
	 * @return A new Mutation object
	 * 
	 * @throws IllegalArgumentException if another is null or the gene/position didn't match.
	 */
	@Deprecated
	public Mutation<VirusT> subtractsBy(Mutation<VirusT> another);

	/**
	 * Subtracts by another Collection&lt;Character&gt; of
	 * AAs and returns the result mutation.
	 *
	 * @param otherAAChars Amino acid characters
	 * 
	 * @return A new Mutation object
	 */
	public Mutation<VirusT> subtractsBy(Collection<Character> otherAAChars);

	/**
	 * Intersects with another Mutation object and returns the result mutation.
	 *
	 * @param another Mutation
	 * 
	 * @return A new merged Mutation object
	 * 
	 * @throws IllegalArgumentException if another is null or the gene/position didn't match.
	 */
	@Deprecated
	public Mutation<VirusT> intersectsWith(Mutation<VirusT> another);

	/**
	 * Intersects with another Collection&lt;Character&gt; of
	 * AAs and returns the result mutation.
	 *
	 * @param otherAAChars Amino acid characters	
	 * 
	 * @return A new merged Mutation object
	 */
	public Mutation<VirusT> intersectsWith(Collection<Character> otherAAChars);

	/**
	 * Checks if the mutation is at a drug resistance position
	 *
	 * @return true if the mutation is at DRM position
	 */
	public boolean isAtDrugResistancePosition();

	/**
	 * Checks if the position is unsequenced.
	 *
	 * @return true if the position is unsequenced.
	 */
	public boolean isUnsequenced();

	/**
	 * Checks if the position is unsequenced by given unsequenced regions.
	 *
	 * @return true if the position is unsequenced.
	 */
	public boolean isUnsequenced(UnsequencedRegions<VirusT> unseqRegions);

	/**
	 * Checks if the mutation is a DRM
	 *
	 * @return true if the mutation is a DRM
	 */
	public boolean isDRM();

	/**
	 * Retrieves drug class of the mutation if it's a DRM
	 * 
	 * @return a DrugClass object or null
	 */
	public DrugClass<VirusT> getDRMDrugClass();

	/**
	 * Checks if the mutation is a TSM
	 *
	 * @return true if the mutation is a TSM
	 */
	public boolean isTSM();

	/**
	 * Retrieves drug class of the mutation if it's a TSM
	 * 
	 * @return a DrugClass object or null
	 */
	public DrugClass<VirusT> getTSMDrugClass();

	/**
	 * Gets gene of this mutation.
	 *
	 * @return A Gene
	 */
	public Gene<VirusT> getGene();

	/**
	 * Gets string gene of this mutation.
	 * 
	 * @return String
	 */
	public String getAbstractGene();

	/**
	 * Gets the reference AA at this position.
	 *
	 * Note: this function returns a String. If you want Character, use `getRefChar()`.
	 *
	 * @return A single character string
	 */
	public String getReference();
	
	/**
	 * Gets the reference AA at this position.
	 *
	 * Note: this function returns a Character. If you want String, use `getReference()`.
	 *
	 * @return A single character string
	 */
	public Character getRefChar();

	/**
	 * Gets the position number of this mutation.
	 *
	 * @return An integer &gt;= 1
	 */
	public int getPosition();

	/**
	 * Gets the GenePosition object of this mutation.
	 *
	 * @return A GenePosition object
	 */
	public GenePosition<VirusT> getGenePosition();

	/**
	 * Gets the amino acids string of this mutation.
	 *
	 * @return A String of amino acids (include insertion/deletion)
	 */
	public String getDisplayAAs();

	/**
	 * Gets a set of amino acid characters.
	 *
	 * - "_" represents an insertion;
	 * - "-" represents a deletion;
	 * - "*" represents a stop codon.
	 *
	 * @return A Set&lt;Character&gt; of amino acids/insertion/deletion/stop codons
	 */
	public Set<Character> getDisplayAAChars();

	/**
	 * Gets original amino acids string of this mutation.
	 *
	 * @return A String of amino acids (include insertion/deletion)
	 */
	public String getAAs();
	
	/**
	 * Gets unusual amino acids string of this mutation.
	 *
	 * @return A String of unusual amino acids
	 */
	public String getUnusualAAs();

	/**
	 * Gets the original set of amino acid characters.
	 *
	 * - "_" represents an insertion;
	 * - "-" represents a deletion;
	 * - "*" represents a stop codon.
	 *
	 * @return A Set&lt;Character&gt; of amino acids/insertion/deletion/stop codons
	 */
	public Set<Character> getAAChars();

	/**
	 * Gets a set of single amino acid mutations.
	 *
	 * @return A Set&lt;Mutation&gt;
	 */
	public Set<Mutation<VirusT>> split();

	/**
	 * Gets the triplet (codon) of this mutation.
	 *
	 * @return A string contains three characters (ACGT...)
	 */
	public String getTriplet();

	/**
	 * Gets the inserted NAs of this mutation.
	 *
	 * @return A string contains multiples of three characters (ACGT...)
	 */
	public String getInsertedNAs();

	/**
	 * Checks if the mutation contains an insertion.
	 *
	 * @return true if the mutation contains an insertion.
	 */
	public boolean isInsertion();

	/**
	 * Checks if the mutation contains a deletion.
	 *
	 * @return true if the mutation contains a deletion.
	 */
	public boolean isDeletion();

	/**
	 * Checks if the mutation contains an insertion or a deletion.
	 *
	 * @return true if the mutation contains an insertion or a deletion.
	 */
	public boolean isIndel();

	/**
	 * Checks if the Mutation object is a mixture of
	 * multiple amino acids / insertion / deletion / stop codon
	 *
	 * @return true if the mutation is a mixture
	 */
	public boolean isMixture();

	/**
	 * Checks if the Mutation object is a mixture with
	 * subtype B reference amino acid presents
	 *
	 * @return true if AAs contains subtype B reference
	 */
	public boolean hasReference();

	/**
	 * Checks if the mutation contains a stop codon
	 *
	 * @return true if the mutation contains a stop codon
	 */
	public boolean hasStop();

	/**
	 * Checks if the mutation is considered unusual
	 *
	 * @return true if the mutation is considered unusual
	 */
	public boolean isUnusual();
	
	/**
	 * Checks if the mutation is considered unusual
	 *  
	 * @param aaPcnts
	 * @return true if the mutation is considered unusual
	 */
	public boolean isUnusual(AminoAcidPercents<VirusT> aaPcnts);

	/**
	 * Checks if the mutation is considered unusual
	 *  
	 * @param treatment
	 * @param subtype
	 * @return true if the mutation is considered unusual
	 */
	public boolean isUnusual(String treatment, String subtype);
	
	/**
	 * Checks if the mutation is an SDRM mutation
	 *
	 * @return true if the mutation is an SDRM mutation
	 */
	public boolean isSDRM();

	/**
	 * Retrieves drug class of the mutation if it's a SDRM
	 * 
	 * @return a DrugClass object or null
	 */
	public DrugClass<VirusT> getSDRMDrugClass();

	/**
	 * Checks if the mutation codon contains BDHVN
	 *
	 * @return true if the codon contains BDHVN
	 */
	public boolean hasBDHVN();

	/**
	 * Checks if the mutation is highly ambiguous
	 *
	 * @return true if the mutation is highly ambiguous
	 */
	public boolean isAmbiguous();

	/**
	 * Checks if the mutation is highly ambiguous without BDHVN
	 *
	 * @return true if the mutation is highly ambiguous without BDHVN
	 */
	public boolean isAmbiguousWithoutBDHVN();
	
	/**
	 * Checks if the mutation is considered APOBEC-mediated (non-DRM)
	 *
	 * @return true if the mutation is considered APOBEC-mediated
	 */
	public boolean isApobecMutation();

	/**
	 * Checks if the mutation is a DRM considered APOBEC-mediated
	 *
	 * @return true if the mutation is a DRM considered APBOEC-mediated
	 */
	public boolean isApobecDRM();

	/**
	 * Gets the highest AA prevalence of this mutation.
	 *
	 * @return A double number of the prevalence (max: 100)
	 */
	public double getHighestMutPrevalence();
	
	/**
	 * Gets all prevalence result of this mutation position.
	 * 
	 * @return A list of MutationPrevalence&lt;VirusT&gt; object
	 */
	public List<MutationPrevalence<VirusT>> getPrevalences();

	/**
	 * Re-ordering AAs to place subtype B
	 * reference (if presented) at the first
	 *
	 * @return a string of AAs
	 */
	public String getAAsWithRefFirst();

	/**
	 * Gets AAs without subtype B reference (if presented)
	 *
	 * @return a string of AAs
	 */
	public String getAAsWithoutReference ();

	/**
	 * Retrieve the primary mutation type of current mutation.
	 *
	 * @return a MutType object
	 */
	public MutationType<VirusT> getPrimaryType();

	/**
	 * Retrieve all mutation types of current mutation.
	 *
	 * The reason we return a List here is because some mixture,
	 * for example, RT215SWY has multiple types.
	 *
	 * @return a List of MutType objects
	 */
	public List<MutationType<VirusT>> getTypes();

	public boolean equals(Object o);

	public int hashCode();

	/**
	 * Standard method to convert this mutation to a string.
	 *
	 * @return a String to represent this mutation.
	 */
	public String toString();

	/**
	 * Method to convert this mutation to a short version string.
	 *
	 * @return String
	 */
	public String getShortText();

	/**
	 * Sort by gene, position and aas.
	 *
	 * The order implemented by this method is similar
	 * to the order implemented by this SQL query:
	 *   SELECT * FROM mutations ORDER BY gene, pos, aas;
	 *
	 */
	public int compareTo (Mutation<VirusT> mut);

	/**
	 * Compares two mutations to determine if they share a non-reference amino acid
	 *
	 * Reference and stop codon are not be responsible for a match.
	 *
	 * @param queryMut Mutation
	 * 
	 * @return true if the Mutation and queryMut share at least one non-reference amino acid
	 */
	public boolean containsSharedAA(Mutation<VirusT> queryMut);

	/**
	 * Compares with given queryAAChars to determine if they shared at least an amino acid
	 *
	 * Reference and stop codon are not be responsible for a match.
	 *
	 * @param queryAAChars		Amino acid Characters set		
	 * @param ignoreRefOrStops  Ignore reference amino acid and stop codon
	 * 
	 * @return true if the Mutation and queryMut share at least one amino acid
	 */
	public boolean containsSharedAA(Set<Character> queryAAChars, boolean ignoreRefOrStops);

	/**
	 * The ASI format consists of an optional reference aa and a position followed by one or more
	 * upper case amino acids. Insertions are represented by 'i', deletions by 'd', and stops by 'Z;
	 * 
	 * @return ASI format
	 */
	public String getASIFormat();

	/**
	 * In HIVDB_Rules, insertions are denoted by '#' and deletions by '~'
	 * This differs from Insertion, _, and i and from Deletion, '-', and d
	 * 
	 * @return HIVDBformat
	 */
	public String getHIVDBFormat();

	/**
	 * If the insertion is known, report it out in the following format T69S_SS
	 * If the insertion is not known report it out as T69Insertion
	 * Report deletions as T69Deletion
	 * If there is a mixture that contains the reference aa, move the ref to
	 * the beginning of the mixture.
	 * Report the reference before the position (i.e. M184V)
	 * 
	 * @return Human formatted mutation
	 */
	public String getHumanFormat();

	public String getShortHumanFormat();

	/**
	 * Similar to getHumanFormat() except that the preceding reference is removed.
	 * 
	 * @return Human format mutation presentation without reference 
	 */
	public String getHumanFormatWithoutLeadingRef();

	public String getHumanFormatWithGene();
	public String getHumanFormatWithAbstractGene();
	
	public Collection<BoundComment<VirusT>> getComments();
	
}
