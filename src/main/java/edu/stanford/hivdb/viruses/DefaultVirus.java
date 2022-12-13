package edu.stanford.hivdb.viruses;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import edu.stanford.hivdb.comments.ConditionalComments;
import edu.stanford.hivdb.drugresistance.algorithm.DrugResistanceAlgorithm;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.genotypes.Genotype;
import edu.stanford.hivdb.genotypes.GenotypeReference;
import edu.stanford.hivdb.genotypes.Genotyper;
import edu.stanford.hivdb.mutations.AminoAcidPercents;
import edu.stanford.hivdb.mutations.CodonPercents;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationPrevalence;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.MutationTypePair;
import edu.stanford.hivdb.seqreads.SequenceReadsAssembler;
import edu.stanford.hivdb.sequences.AlignmentConfig;
import edu.stanford.hivdb.sequences.SequenceAssembler;

public abstract class DefaultVirus<VirusT extends DefaultVirus<VirusT>> implements Virus<VirusT> {

	private final VirusDataLoader<VirusT> dataLoader;

	protected abstract VirusDataLoader<VirusT> getVirusDataLoader();

	public DefaultVirus() {
		this.dataLoader = getVirusDataLoader();
	}

	@Override
	public String getName() {
		return dataLoader.getName();
	}

	@Override
	public Strain<VirusT> getMainStrain() {
		return dataLoader.getMainStrain();
	}

	@Override
	public Collection<Strain<VirusT>> getStrains() {
		return dataLoader.getStrains();
	}

	@Override
	public Strain<VirusT> getStrain(String name) {
		return dataLoader.getStrain(name);
	}

	@Override
	public Collection<Gene<VirusT>> getGenes(Strain<VirusT> strain) {
		return dataLoader.getGenes(strain);
	}

	@Override
	public Gene<VirusT> getGene(String name) {
		return dataLoader.getGene(name);
	}

	@Override
	public Collection<DrugClass<VirusT>> getDrugClasses() {
		return dataLoader.getDrugClasses();
	}

	@Override
	public Map<String, DrugClass<VirusT>> getDrugClassSynonymMap() {
		return dataLoader.getDrugClassSynonymMap();
	}

	@Override
	public DrugClass<VirusT> getDrugClass(String name) {
		return dataLoader.getDrugClass(name);
	}

	@Override
	public Collection<Drug<VirusT>> getDrugs() {
		return dataLoader.getDrugs();
	}

	@Override
	public Map<String, Drug<VirusT>> getDrugSynonymMap() {
		return dataLoader.getDrugSynonymMap();
	}

	@Override
	public Collection<DrugResistanceAlgorithm<VirusT>> getDrugResistAlgorithms() {
		return dataLoader.getDrugResistAlgorithms();
	}

	@Override
	public Collection<DrugResistanceAlgorithm<VirusT>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
		return dataLoader.getDrugResistAlgorithms(algorithmNames);
	}


	@Override
	public DrugResistanceAlgorithm<VirusT> getDrugResistAlgorithm(String name) {
		return dataLoader.getDrugResistAlgorithm(name);
	}

	@Override
	public DrugResistanceAlgorithm<VirusT> getDrugResistAlgorithm(String family, String version) {
		return dataLoader.getDrugResistAlgorithm(family, version);
	}

	@Override
	public Gene<VirusT> extractMutationGene(String mutText) {
		return dataLoader.extractMutationGene(mutText);
	}

	@Override
	public Mutation<VirusT> parseMutationString(Gene<VirusT> defaultGene, String mutText) {
		return dataLoader.parseMutationString(defaultGene, mutText);
	}

	@Override
	public Mutation<VirusT> parseMutationString(String mutText) {
		return dataLoader.parseMutationString(mutText);
	}

	@Override
	public MutationSet<VirusT> newMutationSet(String formattedMuts) {
		return dataLoader.newMutationSet(formattedMuts);
	}

	@Override
	public MutationSet<VirusT> newMutationSet(Collection<String> formattedMuts) {
		return dataLoader.newMutationSet(formattedMuts);
	}

	@Override
	public MutationSet<VirusT> newMutationSet(Gene<VirusT> defaultGene, String formattedMuts) {
		return dataLoader.newMutationSet(defaultGene, formattedMuts);
	}

	@Override
	public MutationSet<VirusT>	newMutationSet(Gene<VirusT> defaultGene, Collection<String> formattedMuts) {
		return dataLoader.newMutationSet(defaultGene, formattedMuts);
	}

	@Override
	public Map<DrugClass<VirusT>, MutationSet<VirusT>> getDrugResistMutations() {
		return dataLoader.getDrugResistMutations();
	}

	@Override
	public Map<DrugClass<VirusT>, MutationSet<VirusT>> getSurveilDrugResistMutations() {
		return dataLoader.getSurveilDrugResistMutations();
	}

	@Override
	public Map<DrugClass<VirusT>, MutationSet<VirusT>> getRxSelectedMutations() {
		return dataLoader.getRxSelectedMutations();
	}

	@Override
	public MutationSet<VirusT> getApobecMutations() {
		return dataLoader.getApobecMutations();
	}

	@Override
	public MutationSet<VirusT> getApobecDRMs() {
		return dataLoader.getApobecDRMs();
	}

	@Override
	public Collection<MutationType<VirusT>> getMutationTypes() {
		return dataLoader.getMutationTypes();
	}

	@Override
	public MutationType<VirusT> getMutationType(String mutTypeText) {
		return dataLoader.getMutationType(mutTypeText);
	}

	@Override
	public Collection<MutationTypePair<VirusT>> getMutationTypePairs() {
		return dataLoader.getMutationTypePairs();
	}

	@Override
	public AminoAcidPercents<VirusT> getAminoAcidPercents(Strain<VirusT> strain, String treatment, String subtype) {
		return dataLoader.getAminoAcidPercents(strain, treatment, subtype);
	}

	@Override
	public CodonPercents<VirusT> getCodonPercents(Strain<VirusT> strain, String treatment, String subtype) {
		return dataLoader.getCodonPercents(strain, treatment, subtype);
	}

	@Override
	public List<MutationPrevalence<VirusT>> getMutationPrevalence(GenePosition<VirusT> genePos) {
		return dataLoader.getMutationPrevalence(genePos);
	}

	@Override
	public ConditionalComments<VirusT> getConditionalComments() {
		return dataLoader.getConditionalComments();
	}

	@Override
	public List<String> getMainSubtypes(Strain<VirusT> strain) {
		return dataLoader.getMainSubtypes(strain);
	}

	@Override
	public Map<Gene<VirusT>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<VirusT> strain) {
		return dataLoader.getNumPatientsForAAPercents(strain);
	}

	@Override
	public Collection<Genotype<VirusT>> getGenotypes() {
		return dataLoader.getGenotypes();
	}

	@Override
	public Genotype<VirusT> getGenotype(String name) {
		return dataLoader.getGenotype(name);
	}

	@Override
	public Genotype<VirusT> getGenotypeUnknown() {
		return dataLoader.getGenotypeUnknown();
	}

	@Override
	public List<GenotypeReference<VirusT>> getGenotypeReferences() {
		return dataLoader.getGenotypeReferences();
	}

	@Override
	public Genotyper<VirusT> getGenotyper() {
		return dataLoader.getGenotyper();
	}

	@Override
	public Double getGenotypeUnknownThreshold() {
		return 0.01;
	}

	@Override
	public Double getGenotypeMaxFallbackToSecondaryDistanceDiff() {
		return 0.001;
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// VirusT instance is a singleton
		return false;
	}

	@Override
	public AlignmentConfig<VirusT> getAlignmentConfig() {
		return dataLoader.getAlignmentConfig();
	}

	@Override
	public SequenceReadsAssembler<VirusT> getSequenceReadsAssembler(Strain<VirusT> strain) {
		return dataLoader.getSeqReadsAssemblers().get(strain);
	}

	@Override
	public SequenceAssembler<VirusT> getSequenceAssembler(Strain<VirusT> strain) {
		return dataLoader.getSequenceAssemblers().get(strain);
	}

	@Override
	public String getGeneDisplay(Gene<VirusT> gene) {
		return getGeneDisplay(gene.getAbstractGene());
	}

	@Override
	public String getGeneDisplay(String geneName) {
		return geneName;
	}

}
