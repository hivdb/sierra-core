/*

    Copyright (C) 2019 Stanford HIVDB team

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

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.builder.HashCodeBuilder;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.CodonModifier;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.utilities.Json;

public class Gene<VirusT extends Virus<VirusT>> implements Comparable<Gene<VirusT>> {

	public static final Character WILDCARD = '.';
	public static final Pattern TRIM_WILDCARD_PATTERN = Pattern.compile("(^\\.+|\\.+$)");

	private final String name;
	private final Integer ordinal;
	private final String strain;
	private final String abstractGene;
	private final String refSequence;
	private final List<CodonModifier<VirusT>> codonModifiers;
	private final List<String> mutationTypes;
	private final Integer nucaminoMinNumOfAA;

	private transient VirusT virusInstance;
	private transient boolean codonModifiersAttached;
	private transient Set<MutationType<VirusT>> mutationTypeObjs;
	private transient Set<DrugClass<VirusT>> drugClasses;

	public static <VirusT extends Virus<VirusT>> Map<String, Gene<VirusT>> loadJson(String raw, VirusT virusIns) {
		Map<String, Gene<VirusT>> genes = new LinkedHashMap<>();
		List<Gene<VirusT>> geneList = Json.loads(
			raw, new TypeToken<List<Gene<VirusT>>>(){});
		for (Gene<VirusT> gene : geneList) {
			gene.virusInstance = virusIns;
			genes.put(gene.getName(), gene);
		}
		return Collections.unmodifiableMap(genes);
	}

	public VirusT getVirusInstance() {
		return virusInstance;
	}

	private Gene(
		String name, Integer ordinal, String strain, String abstractGene,
		String refSequence, List<String> mutationTypes,
		List<CodonModifier<VirusT>> codonModifiers,
		Integer nucaminoMinNumOfAA) {
		this.name = name;
		this.ordinal = ordinal;
		this.strain = strain;
		this.abstractGene = abstractGene;
		this.refSequence = refSequence;
		this.mutationTypes = mutationTypes;
		this.codonModifiers = codonModifiers;
		this.nucaminoMinNumOfAA = nucaminoMinNumOfAA;
	}
	
	public String name() {
		return name;
	}
	
	public String getName() {
		return name;
	}
	
	@Deprecated
	public int getLength() {
		// deprecated, use getAASize instead
		return refSequence.length();
	}
	
	public int getAASize() {
		return refSequence.length();
	}
	
	public int getNASize() {
		return refSequence.length() * 3;
	}
	
	public int getNucaminoMinNumOfAA() {
		return nucaminoMinNumOfAA;
	}

	public Character getRefChar(int pos) {
		return refSequence.charAt(pos - 1);
	}

	public String getRefSequence(int pos, int length) {
		return refSequence.substring(pos - 1, pos - 1 + length);
	}

	public String getRefSequence() {
		return refSequence;
	}

	public Strain<VirusT> getStrain() {
		return virusInstance.getStrain(strain);
	}

	public String getAbstractGene() {
		return abstractGene;
	}
	
	public List<CodonModifier<VirusT>> getTargetCodonModifiers(Strain<VirusT> targetStrain) {
		if (!codonModifiersAttached) {
			codonModifiers
			.stream()
			.forEach(cm -> cm.setVirusInstance(virusInstance));
			codonModifiersAttached = true;
		}
		return (
			codonModifiers
			.stream()
			.filter(cm -> cm.getTargetStrain() == targetStrain)
			.collect(Collectors.toList())
		);
	}

	public Set<DrugClass<VirusT>> getDrugClasses() {
		if (drugClasses == null) {
			Strain<VirusT> strain = getStrain();
			Set<DrugClass<VirusT>> drugClasses = (
				virusInstance
				.getDrugClasses()
				.stream()
				.filter(dc -> (
					dc.getStrains().contains(strain) &&
					dc.getAbstractGene().equals(abstractGene)
				))
				.collect(Collectors.toCollection(LinkedHashSet::new))
			);
			this.drugClasses = Collections.unmodifiableSet(drugClasses);
		}
		return drugClasses;
	}
	
	public Set<MutationType<VirusT>> getMutationTypes() {
		if (mutationTypeObjs == null) {
			Set<MutationType<VirusT>> mutationTypeObjs = (
				mutationTypes.stream()
				.map(mtype -> virusInstance.getMutationType(mtype))
				.collect(Collectors.toCollection(LinkedHashSet::new))
			);
			this.mutationTypeObjs = Collections.unmodifiableSet(mutationTypeObjs);
		}
		return mutationTypeObjs;
	}

	/**
	 * Apply codon modifiers for amino acid sequence
	 *
	 * The refSequence Positions between strains of same "Virus" can be
	 * varied. For example, for RT and IN of HIV-2 virus, there are
	 * deletions and insertions at non-DRM positions comparing to HXB
	 * reference.
	 * 
	 * This function allows to "re-fit" the sequences from one strain
	 * into another strain's numbering system. For example, ROD/EHO to
	 * HXB2, in order to make the sequence be compatible with HIV-1
	 * dedicated data structure.
	 *
	 * @param aaseq
	 * @param firstAA
	 * @param lastAA
	 * @param targetStrain
	 * @return String of adjusted AA alignment in full gene length
	 */
	public String applyCodonModifiersForAASeq(
		String aaseq, int firstAA, int lastAA, Strain<VirusT> targetStrain
	) {
		int numPrefixAAs = firstAA - 1;
		int numSuffixAAs = this.refSequence.length() - aaseq.length() - numPrefixAAs;
		aaseq =
			StringUtils.repeat(WILDCARD, numPrefixAAs) +
			aaseq +
			StringUtils.repeat(WILDCARD, numSuffixAAs);
		if (codonModifiers.size() == 0) {
			return aaseq;
		}
		List<CodonModifier<VirusT>> targetCodonModifiers = getTargetCodonModifiers(targetStrain);
		if (targetCodonModifiers.size() == 0) {
			return aaseq;
		}
		for (CodonModifier<VirusT> cm : targetCodonModifiers) {
			int pos = cm.getPosition();
			Integer insertAfter = cm.getInsertAfter();
			Integer deleteAfter = cm.getDeleteAfter();
			if (insertAfter != null) {
				if (insertAfter <= 0) {
					throw new RuntimeException(String.format(
						"unable to add %s AA(s) to aaseq", insertAfter
					));
				}
				// refSequence deletion, add AA(s) to aaseq
				aaseq =
					aaseq.substring(0, pos) +
					StringUtils.repeat(WILDCARD, insertAfter) +
					aaseq.substring(pos);
			}
			else {
				// refSequence insertion, delete AA(s) from aaseq
				if (deleteAfter > 0) {
					aaseq =
						aaseq.substring(0, pos) +
						aaseq.substring(pos + deleteAfter);

				}
				else if (deleteAfter == 0) {
					// delete anything after
					aaseq = aaseq.substring(0, pos);
				}
				else {
					throw new RuntimeException(String.format(
						"unable to remove %s AA(s) from aaseq", deleteAfter
					));
				}
			}
		}
		return aaseq;
	}

	/**
	 * Apply codon modifiers for nucleotide sequence
	 *
	 * The refSequence Positions between strains of same "Virus" can be
	 * varied. For example, for RT and IN of HIV-2 virus, there are
	 * deletions and insertions at non-DRM positions comparing to HXB
	 * reference.
	 * 
	 * This function allows to "re-fit" the sequences from one strain
	 * into another strain's numbering system. For example, ROD/EHO to
	 * HXB2, in order to make the sequence be compatible with HIV-1
	 * dedicated data structure.
	 *
	 * @param naseq
	 * @param firstAA
	 * @param lastAA
	 * @param targetStrain
	 * @return String of adjusted NA alignment in full gene length
	 */
	public String applyCodonModifiersForNASeq(
		String naseq, int firstAA, int lastAA, Strain<VirusT> targetStrain
	) {
		int numPrefixNAs = (firstAA - 1) * 3;
		int numSuffixNAs = getNASize() - naseq.length() - numPrefixNAs;
		naseq =
			StringUtils.repeat(WILDCARD, numPrefixNAs) +
			naseq +
			StringUtils.repeat(WILDCARD, numSuffixNAs);
		if (codonModifiers.size() == 0) {
			return naseq;
		}
		List<CodonModifier<VirusT>> targetCodonModifiers = getTargetCodonModifiers(targetStrain);
		if (targetCodonModifiers.size() == 0) {
			return naseq;
		}
		for (CodonModifier<VirusT> cm : targetCodonModifiers) {
			int pos = cm.getPosition();
			Integer insertAfter = cm.getInsertAfter();
			Integer deleteAfter = cm.getDeleteAfter();
			if (insertAfter != null) {
				if (insertAfter <= 0) {
					throw new RuntimeException(String.format(
						"unable to add %s codon(s) to naseq", insertAfter
					));
				}
				// refSequence deletion, add codon(s) to naseq
				naseq =
					naseq.substring(0, pos * 3) +
					StringUtils.repeat(WILDCARD, 3 * insertAfter) +
					naseq.substring(pos * 3);
			}
			else {
				// refSequence insertion, delete codon(s) from naseq
				pos = -pos;
				if (deleteAfter > 0) {
					naseq =
						naseq.substring(0, pos * 3) +
						naseq.substring((pos + deleteAfter) * 3);
				}
				else if (deleteAfter == 0) {
					// delete anything after
					naseq = naseq.substring(0, pos * 3);
				}
				else {
					throw new RuntimeException(String.format(
						"unable to remove %s codon(s) from naseq", deleteAfter
					));
				}
			}
		}
		return naseq;
	}

	public Collection<GenePosition<VirusT>> getGenePositionsBetween(int startPos, int endPos) {
		Set<GenePosition<VirusT>> genePositions = new LinkedHashSet<>();
		startPos = Math.max(startPos, 1);
		endPos = Math.min(endPos, getLength());
		for (int pos = startPos; pos <= endPos; pos ++) {
			genePositions.add(new GenePosition<>(this, pos));
		}
		return genePositions;
	}

	@Override
	public String toString() {
		return name();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// genes are singletons
		return false;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(46548647, 769543361)
			.append(name)
			.toHashCode();
	}

	@Override
	public int compareTo(Gene<VirusT> o) {
		return ordinal.compareTo(o.ordinal);
	}

}
