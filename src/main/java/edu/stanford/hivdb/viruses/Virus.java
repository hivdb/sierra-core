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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import edu.stanford.hivdb.comments.ConditionalComments;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.drugs.DrugResistanceAlgorithm;
import edu.stanford.hivdb.genotypes.Genotype;
import edu.stanford.hivdb.genotypes.GenotypeReference;
import edu.stanford.hivdb.genotypes.Genotyper;
import edu.stanford.hivdb.mutations.AminoAcidPercent;
import edu.stanford.hivdb.mutations.AminoAcidPercents;
import edu.stanford.hivdb.mutations.CodonPercents;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationPrevalence;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.MutationTypePair;
import edu.stanford.hivdb.mutations.MutationsValidator;
import edu.stanford.hivdb.seqreads.SequenceReads;
import edu.stanford.hivdb.seqreads.SequenceReadsValidator;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.sequences.SequenceValidator;
import edu.stanford.hivdb.utilities.ValidationResult;

public interface Virus<VirusT extends Virus<VirusT>> {

	public final static Map<String, Virus<?>> singletons = new HashMap<>();
	public final static Map<String, Virus<?>> singletonByClassName = new HashMap<>();
	public final static Map<String, SequenceValidator<?>> sequenceValidators = new HashMap<>();
	public final static Map<String, MutationsValidator<?>> mutationsValidators = new HashMap<>();;
	public final static Map<String, SequenceReadsValidator<?>> sequenceReadsValidators = new HashMap<>();

	public static void registerInstance(Virus<?> singleton) {
		singletons.put(singleton.getName(), singleton);
		singletonByClassName.put(singleton.getClass().getName(), singleton);
	}

	public static Virus<?> getInstance(String name) {
		return singletons.get(name);
	}

	public static <VirusT extends Virus<VirusT>> VirusT getInstance(Class<VirusT> klass) {
		return klass.cast(singletonByClassName.get(klass.getName()));
	}

	public String getName();

	public Collection<Strain<VirusT>> getStrains();

	public default Strain<VirusT> getStrain(String name) {
		for (Strain<VirusT> s : getStrains()) {
			if (s.getName().equals(name)) {
				return s;
			}
		}
		throw new NullPointerException(String.format("strain \"%s\" not found", name));
	}

	public Collection<Gene<VirusT>> getGenes(Strain<VirusT> strain);

	public Gene<VirusT> getGene(String name);

	public default Collection<String> getAbstractGenes() {
		Set<String> absGenes = new LinkedHashSet<>();
		for (Strain<VirusT> strain : getStrains()) {
			for (Gene<VirusT> gene : strain.getGenes()) {
				absGenes.add(gene.getAbstractGene());
			}
		}
		return absGenes;
	}

	public Collection<DrugClass<VirusT>> getDrugClasses();

	public Map<String, DrugClass<VirusT>> getDrugClassSynonymMap();

	public default DrugClass<VirusT> getDrugClass(String name) {
		return getDrugClassSynonymMap().get(name);
	}

	public Collection<Drug<VirusT>> getDrugs();

	public Map<String, Drug<VirusT>> getDrugSynonymMap();

	public default Collection<Drug<VirusT>> getDrugs(DrugClass<VirusT> drugClass) {
		return (
			getDrugs()
			.stream()
			.filter(d -> d.getDrugClass() == drugClass)
			.collect(Collectors.toList())
		);
	}

	public default Drug<VirusT> getDrug(String name) {
		return getDrugSynonymMap().get(name);
	}

	public Collection<Genotype<VirusT>> getGenotypes();

	public default Genotype<VirusT> getGenotype(String name) {
		for (Genotype<VirusT> gt : getGenotypes()) {
			if (gt.getIndexName().equals(name)) {
				return gt;
			}
		}
		throw new NullPointerException(String.format("Genotype \"%s\" not found", name));
	}

	public Genotype<VirusT> getGenotypeUnknown();

	public Genotyper<VirusT> getGenotyper();

	public List<GenotypeReference<VirusT>> getGenotypeReferences();

	public Strain<VirusT> getGenotypingCompatibleStrain();

	public Gene<VirusT> extractMutationGene(String mutText);

	public Mutation<VirusT> parseMutationString(Gene<VirusT> defaultGene, String mutText);

	public Mutation<VirusT> parseMutationString(String mutText);

	public MutationSet<VirusT> newMutationSet(String formattedMuts);

	public MutationSet<VirusT> newMutationSet(Collection<String> formattedMuts);

	public MutationSet<VirusT> newMutationSet(Gene<VirusT> defaultGene, String formattedMuts);

	public MutationSet<VirusT> newMutationSet(Gene<VirusT> defaultGene, Collection<String> formattedMuts);

	public Map<DrugClass<VirusT>, MutationSet<VirusT>> getDrugResistMutations();

	public default MutationSet<VirusT> getDrugResistMutations(DrugClass<VirusT> drugClass) {
		return getDrugResistMutations().get(drugClass);
	}

	public Map<DrugClass<VirusT>, MutationSet<VirusT>> getSurveilDrugResistMutations();

	public default MutationSet<VirusT> getSurveilDrugResistMutations(DrugClass<VirusT> drugClass) {
		return getSurveilDrugResistMutations().get(drugClass);

	}

	public Map<DrugClass<VirusT>, MutationSet<VirusT>> getRxSelectedMutations();

	public default MutationSet<VirusT> getRxSelectedMutations(DrugClass<VirusT> drugClass) {
		return getRxSelectedMutations().get(drugClass);
	}

	public MutationSet<VirusT> getApobecMutations();

	public MutationSet<VirusT> getApobecDRMs();

	public AminoAcidPercents<VirusT> getAminoAcidPercents(Strain<VirusT> strain, String treatment, String subtype);

	public CodonPercents<VirusT> getCodonPercents(Strain<VirusT> strain, String treatment, String subtype);

	public Collection<MutationType<VirusT>> getMutationTypes();

	public default MutationType<VirusT> getMutationType(String mutTypeText) {
		for (MutationType<VirusT> mtype : getMutationTypes()) {
			if (mtype.getName().equals(mutTypeText)) {
				return mtype;
			}
		}
		throw new NullPointerException(String.format("Mutation type \"%s\" not found", mutTypeText));
	}

	public default MutationType<VirusT> getOtherMutationType() {
		return getMutationType("Other");
	}

	public Collection<MutationTypePair<VirusT>> getMutationTypePairs();

	public default Collection<MutationTypePair<VirusT>> getMutationTypePairs(Strain<VirusT> strain) {
		return (
			getMutationTypePairs()
			.stream()
			.filter(mtp -> mtp.getStrain() == strain)
			.collect(Collectors.toList())
		);
	}

	public default Map<Gene<VirusT>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<VirusT> strain) {
		List<String> mainSubtypes = getMainSubtypes(strain);
		Map<Gene<VirusT>, Map<String, Integer[]>> aaPcntsNumPatients = new LinkedHashMap<>();
		String[] treatments = new String[] {"naive", "art"};
		for (Gene<VirusT> gene : strain.getGenes()) {
			Map<String, Integer[]> geneNumPatients = new LinkedHashMap<>();
			for (String subtype : mainSubtypes) {
				geneNumPatients.put(subtype, new Integer[] {0, 0});
				for (int i = 0; i < 2; i ++) {
					String rx = treatments[i];
					AminoAcidPercents<VirusT> aaPcnts = getAminoAcidPercents(strain, rx, subtype);
					int max = aaPcnts.get().stream().mapToInt(ap -> ap.getTotal()).max().getAsInt();
					geneNumPatients.get(subtype)[i] = max;
				}
			}
			aaPcntsNumPatients.put(gene, geneNumPatients);
		}
		return aaPcntsNumPatients;
	}

	public List<String> getMainSubtypes(Strain<VirusT> strain);

	public Collection<DrugResistanceAlgorithm<VirusT>> getDrugResistAlgorithms();

	public default Collection<DrugResistanceAlgorithm<VirusT>> getDrugResistAlgorithms(Strain<VirusT> strain) {
		return (
			getDrugResistAlgorithms()
			.stream()
			.filter(dra -> dra.getStrain() == strain)
			.collect(Collectors.toList())
		);
	}

	public default Collection<DrugResistanceAlgorithm<VirusT>> getDrugResistAlgorithmsByFamily(String family) {
		return (
			getDrugResistAlgorithms()
			.stream()
			.filter(dra -> dra.getFamily().equals(family))
			.collect(Collectors.toList())
		);
	}

	public default DrugResistanceAlgorithm<VirusT> getDrugResistAlgorithm(String name) {
		return (
			getDrugResistAlgorithms()
			.stream()
			.filter(dra -> dra.getName().equals(name))
			.findFirst()
			.get()
		);
	}

	public default DrugResistanceAlgorithm<VirusT> getDrugResistAlgorithm(String family, String version) {
		return (
			getDrugResistAlgorithms()
			.stream()
			.filter(dra -> dra.getFamily().equals(family) && dra.getVersion().equals(version))
			.findFirst()
			.get()
		);
	}

	public default Collection<DrugResistanceAlgorithm<VirusT>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
		return (
			algorithmNames.stream()
			.map(name -> getDrugResistAlgorithm(name))
			.collect(Collectors.toList())
		);
	}

	public default DrugResistanceAlgorithm<VirusT> getLatestDrugResistAlgorithm(String family) {
		return (
			getDrugResistAlgorithmsByFamily(family)
			.stream()
			.sorted((a, b) -> - a.getVersion().compareTo(b.getVersion()))
			.findFirst()
			.get()
		);
	}

	public default List<MutationPrevalence<VirusT>> getMutationPrevalence(GenePosition<VirusT> genePos) {
		Strain<VirusT> strain = genePos.getStrain();
		List<String> mainSubtypes = getMainSubtypes(strain);
		List<MutationPrevalence<VirusT>> mutPrevs = new ArrayList<>();
		for (String subtype : mainSubtypes) {
			Map<Character, List<AminoAcidPercent<VirusT>>> naiveSubtypeAAPcnts = new TreeMap<>();
			Map<Character, List<AminoAcidPercent<VirusT>>> artSubtypeAAPcnts = new TreeMap<>();
			for (String rx : new String[] {"naive", "art"}) {
				Map<Character, List<AminoAcidPercent<VirusT>>> subtypeAAPcnts;
				if (rx.equals("naive")) {
					subtypeAAPcnts = naiveSubtypeAAPcnts;
				}
				else {
					subtypeAAPcnts = artSubtypeAAPcnts;
				}
				for (
					AminoAcidPercent<VirusT> aaPcnt :
					getAminoAcidPercents(strain, rx, subtype).get(genePos)
				) {
					char aa = aaPcnt.getAA();
					char ref = aaPcnt.getRefChar();
					if (aa == ref) {
						continue;
					}
					if (!subtypeAAPcnts.containsKey(aa)) {
						subtypeAAPcnts.put(aa, new ArrayList<>());
					}
					subtypeAAPcnts.get(aa).add(aaPcnt);
				}
			}
			for (char aa : naiveSubtypeAAPcnts.keySet()) {
				List<AminoAcidPercent<VirusT>> naiveRxAAPcnts = naiveSubtypeAAPcnts.get(aa);
				AminoAcidPercent<VirusT> naiveAAPcnt = naiveRxAAPcnts.get(0);
				List<AminoAcidPercent<VirusT>> artRxAAPcnts = artSubtypeAAPcnts.get(aa);
				AminoAcidPercent<VirusT> artAAPcnt = artRxAAPcnts.get(0);
				if (naiveAAPcnt.getPercent() < 0.0001 && artAAPcnt.getPercent() < 0.0001) {
					// AA prevalence must be greater than 0.01%
					continue;
				}
				mutPrevs.add(new MutationPrevalence<>(subtype, naiveAAPcnt, artAAPcnt));
			}
		}
		if (mutPrevs.stream().anyMatch(mp -> (
			mp.getFrequencyNaive() > 0 || mp.getFrequencyTreated() > 0 ||
			mp.getPercentageNaive() >= 0.0001 || mp.getPercentageTreated() >= 0.0001))
		) {
			return Collections.unmodifiableList(mutPrevs);
		}
		return Collections.emptyList();
	}

	public ConditionalComments<VirusT> getConditionalComments();

	public default void registerSequenceValidator(SequenceValidator<VirusT> validator) {
		sequenceValidators.put(this.getClass().getName(), validator);
	}

	public default List<ValidationResult> validateSequence(AlignedSequence<VirusT> alignedSeq) {
		String className = this.getClass().getName();
		if (sequenceValidators.containsKey(className)) {
			@SuppressWarnings("unchecked")
			SequenceValidator<VirusT> validator = (SequenceValidator<VirusT>) sequenceValidators.get(this.getClass().getName());
			return validator.validate(alignedSeq);
		}
		else {
			return Collections.emptyList();
		}
	}

	public default void registerMutationsValidator(MutationsValidator<VirusT> validator) {
		mutationsValidators.put(this.getClass().getName(), validator);
	}

	public default List<ValidationResult> validateMutations(MutationSet<VirusT> mutations) {
		String className = this.getClass().getName();
		if (mutationsValidators.containsKey(className)) {
			@SuppressWarnings("unchecked")
			MutationsValidator<VirusT> validator = (MutationsValidator<VirusT>) mutationsValidators.get(this.getClass().getName());
			return validator.validate(mutations);
		}
		else {
			return Collections.emptyList();
		}
	}

	public default void registerSequenceReadsValidator(SequenceReadsValidator<VirusT> validator) {
		sequenceReadsValidators.put(this.getClass().getName(), validator);
	}

	public default List<ValidationResult> validateSequenceReads(SequenceReads<VirusT> seqReads) {
		String className = this.getClass().getName();
		if (sequenceReadsValidators.containsKey(className)) {
			@SuppressWarnings("unchecked")
			SequenceReadsValidator<VirusT> validator = (SequenceReadsValidator<VirusT>) sequenceReadsValidators.get(this.getClass().getName());
			return validator.validate(seqReads);
		}
		else {
			return Collections.emptyList();
		}
	}

}
