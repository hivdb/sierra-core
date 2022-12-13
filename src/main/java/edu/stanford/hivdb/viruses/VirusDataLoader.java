/*

    Copyright (C) 2020 Stanford HIVDB team

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

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;

import com.google.common.collect.Sets;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.comments.ConditionalComments;
import edu.stanford.hivdb.drugresistance.algorithm.DrugResistanceAlgorithm;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.genotypes.Genotype;
import edu.stanford.hivdb.genotypes.GenotypeReference;
import edu.stanford.hivdb.genotypes.Genotyper;
import edu.stanford.hivdb.mutations.AAMutation;
import edu.stanford.hivdb.mutations.AminoAcidPercents;
import edu.stanford.hivdb.mutations.CodonPercents;
import edu.stanford.hivdb.mutations.CodonMutation;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationPrevalence;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.MutationTypePair;
import edu.stanford.hivdb.seqreads.SequenceReadsAssembler;
import edu.stanford.hivdb.sequences.AlignmentConfig;
import edu.stanford.hivdb.sequences.SequenceAssembler;
import edu.stanford.hivdb.utilities.AAUtils;
import edu.stanford.hivdb.utilities.AssertUtils;
import edu.stanford.hivdb.utilities.Json;

public abstract class VirusDataLoader<T extends Virus<T>> {

	private static final Pattern NON_ASI_AA_PATTERN = Pattern.compile(
		"^([AC-IK-NP-TV-Z.*]+(?:[#_]?[AC-IK-NP-TV-Z.*]+)?|[id_#~-]|[iI]ns(?:ertion)?|[dD]el(?:etion)?)$"
	);

	protected static String loadResource(String resPath, ClassLoader classLoader, String absentDefault) {
		if (resPath != null && resPath.toLowerCase().startsWith("https://")) {
			try {
				URL url = new URL(resPath);
				HttpURLConnection conn = (HttpURLConnection) url.openConnection();
				conn.setRequestProperty("Accept-Encoding", "gzip, deflate");
				InputStream stream;
    	  if ("gzip".equals(conn.getContentEncoding())) {
    	     stream = new GZIPInputStream(conn.getInputStream());
    	  }
    	  else {
    	     stream = conn.getInputStream();
    	  }
				return IOUtils.toString(stream, StandardCharsets.UTF_8);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		else if (resPath != null) {
			try (
				InputStream stream = classLoader
					.getResourceAsStream(resPath);
			) {
				return IOUtils.toString(stream, StandardCharsets.UTF_8);
			} catch (IOException|NullPointerException e) {
				throw new ExceptionInInitializerError(
					String.format("Invalid resource name (%s)", resPath)
				);
			}
		}
		else if (absentDefault != null) {
			return absentDefault;
		}
		else {
			throw new ExceptionInInitializerError(
				"resPath and absentDefault are both NULL"
			);
		}
	}

	private static final String URL_ALLOW_PURGING_CACHE = "https://s3-us-west-2.amazonaws.com/cms.hivdb.org/sierra-sars2-allow-purging-cache";

	private final T virus;
	private final Pattern mutationPattern;
	private final String virusName;
	private final String mainStrainText;
	private final String strainsResPath;
	private final String genesResPath;
	private final String drugClassesResPath;
	private final String drugsResPath;
	private final String drmsResPath;
	private final String sdrmsResPath;
	private final String tsmsResPath;
	private final String apobecsResPath;
	private final String apobecDRMsResPath;
	private final String aaPcntsResPath;
	private final String codonPcntsResPath;
	private final String mutTypesResPath;
	private final String mutTypePairsResPath;
	private final String mainSubtypesResPath;
	private final String genotypeReferencesResPath;
	private final String genotypesRespath;
	private final String algorithmsIndexPath;
	private final String algorithmsResPath;
	private final String condCommentsResPath;
	private final String alignConfigResPath;
	private final String assemblyConfigResPath;

	private transient ClassLoader classLoader;
	private transient Map<String, Strain<T>> strains;
	private transient Map<String, Gene<T>> genes;
	private transient Map<String, DrugClass<T>> drugClasses;
	private transient Map<String, Drug<T>> drugs;
	private transient Map<DrugClass<T>, MutationSet<T>> drugResistMutations;
	private transient Map<DrugClass<T>, MutationSet<T>> surveilDrugResistMuts;
	private transient Map<DrugClass<T>, MutationSet<T>> rxSelectedMutations;
	private transient MutationSet<T> apobecMutations;
	private transient MutationSet<T> apobecDRMs;
	private transient Map<String, AminoAcidPercents<T>> aminoAcidPcnts = new HashMap<>();
	private transient Map<String, CodonPercents<T>> codonPcnts = new HashMap<>();
	private transient Map<String, MutationType<T>> mutationTypes;
	private transient List<MutationTypePair<T>> mutationTypePairs;
	private transient Map<Strain<T>, List<String>> mainSubtypes;
	private transient Map<GenePosition<T>, List<MutationPrevalence<T>>> mutPrevalenceMap = new HashMap<>();
	private transient Map<Strain<T>, Map<Gene<T>, Map<String, Integer[]>>> allAAPcntsNumPatients = new HashMap<>();
	private transient Map<String, Genotype<T>> allGenotypes;
	private transient List<GenotypeReference<T>> allGenotypeReferences;
	private transient Genotyper<T> genotyper;
	private transient List<DrugResistanceAlgorithm<T>> drugResistAlgs;
	private transient Map<String, DrugResistanceAlgorithm<T>> drugResistAlgLookup;
	private transient ConditionalComments<T> condComments;
	private transient AlignmentConfig<T> alignmentConfig;
	private transient Map<Strain<T>, SequenceReadsAssembler<T>> seqReadsAssemblers;
	private transient Map<Strain<T>, SequenceAssembler<T>> sequenceAssemblers;

	public VirusDataLoader(T virus) {
		this.virus = virus;
		this.mutationPattern = Pattern.compile(
			String.format(
				"^\\s*" +
				"(__ASI__)?((?i:%s))?[:_-]?" +
				"([AC-IK-NP-TV-Y])?" +
				"(\\d{1,4})" +
				"([AC-IK-NP-TV-Zid.*]+(?:[#_]?[AC-IK-NP-TV-Z.*]+)?|[id_#~-]|[iI]ns(?:ertion)?|[dD]el(?:etion)?)" +
				"(?::([ACGTRYMWSKBDHVN-]{3})?)?" +
				"\\s*$",
				getGenePattern()
			)
		);

		this.classLoader = virus.getClass().getClassLoader();
		this.virusName = getVirusName();
		this.mainStrainText = getMainStrainText();
		this.strainsResPath = getStrainsResPath();
		this.genesResPath = getGenesResPath();
		this.drugClassesResPath = getDrugClassesResPath();
		this.drugsResPath = getDrugsResPath();
		this.drmsResPath = getDRMsResPath();
		this.sdrmsResPath = getSDRMsResPath();
		this.tsmsResPath = getTSMsResPath();
		this.apobecsResPath = getApobecsResPath();
		this.apobecDRMsResPath = getApobecDRMsResPath();
		this.aaPcntsResPath = getAAPcntsResPath();
		this.codonPcntsResPath = getCodonPcntsResPath();
		this.mutTypesResPath = getMutTypesResPath();
		this.mutTypePairsResPath = getMutTypePairsResPath();
		this.mainSubtypesResPath = getMainSubtypesResPath();
		this.genotypeReferencesResPath = getGenotypeReferencesResPath();
		this.genotypesRespath = getGenotypesRespath();
		this.algorithmsIndexPath = getAlgorithmsIndexPath();
		this.algorithmsResPath = getAlgorithmsResPath();
		this.condCommentsResPath = getCondCommentsResPath();
		this.alignConfigResPath = getAlignConfigResPath();
		this.assemblyConfigResPath = getAssemblyConfigResPath();
	}

	/**
	 * Retrieve a string of RegExp for matching all valid gene names.
	 * The defined pattern should not use any group capturing "(...)".
	 *
	 * @return string
	 */
	protected abstract String getGenePattern();

	/**
	 * Retrieve the virus name, e.g. HIV2
	 *
	 * @return string
	 */
	protected abstract String getVirusName();

	/**
	 * Retrieve a string of main strain name, e.g. HIV2A
	 *
	 * @return string
	 */
	protected abstract String getMainStrainText();

	/**
	 * Retrieve resource path of strains
	 *
	 * @return string
	 */
	protected abstract String getStrainsResPath();
	protected abstract String getGenesResPath();
	protected abstract String getDrugClassesResPath();
	protected abstract String getDrugsResPath();
	protected abstract String getDRMsResPath();
	protected abstract String getSDRMsResPath();
	protected abstract String getTSMsResPath();
	protected abstract String getApobecsResPath();
	protected abstract String getApobecDRMsResPath();
	protected abstract String getAAPcntsResPath();
	protected abstract String getCodonPcntsResPath();
	protected abstract String getMutTypesResPath();
	protected abstract String getMutTypePairsResPath();
	protected abstract String getMainSubtypesResPath();
	protected abstract String getGenotypeReferencesResPath();
	protected abstract String getGenotypesRespath();
	protected abstract String getAlgorithmsIndexPath();
	protected abstract String getAlgorithmsResPath();
	protected abstract String getCondCommentsResPath();
	protected abstract String getAlignConfigResPath();
	protected abstract String getAssemblyConfigResPath();

	private MutationSet<T> loadMutationSetFromRes(String resPath, Collection<Strain<T>> strains) {
		String raw = loadResource(resPath, classLoader, "[]");
		return (
			strains.stream()
			.map(strain -> MutationSet.loadJson(raw, geneText -> strain.getGene(geneText)))
			.reduce(new MutationSet<>(), (acc, other) -> acc.mergesWith(other))
		);
	}

	private Map<DrugClass<T>, MutationSet<T>> loadMutationSetByDrugClassFromRes(String resPath, Collection<Strain<T>> strains) {
		Map<DrugClass<T>, MutationSet<T>> mutationsMap = new LinkedHashMap<>();
		String raw = loadResource(resPath, classLoader, "{}");

		Map<String, List<Map<String, ?>>> muts = Json.loads(
			raw, new TypeToken<Map<String, List<Map<String, ?>>>>(){});
		for (String drugClassText : muts.keySet()) {
			DrugClass<T> drugClass = getDrugClass(drugClassText);
			mutationsMap.put(
				drugClass,
				strains.stream()
				.map(strain -> MutationSet.loadJsonMap(
					muts.get(drugClassText),
					geneText -> strain.getGene(geneText)
				))
				.reduce(new MutationSet<>(), (acc, other) -> acc.mergesWith(other))
			);
		}
		return Collections.unmodifiableMap(mutationsMap);
	}

	private void initCondComments() {
		String raw = loadResource(condCommentsResPath, classLoader, "[]");
		this.condComments = new ConditionalComments<>(raw, virus);
	}

	private void initMainSubtypes() {
		String raw = loadResource(mainSubtypesResPath, classLoader, "{}");
		Map<String, List<String>> subtypes = Json.loads(raw, new TypeToken<Map<String, List<String>>>() {});
		Map<Strain<T>, List<String>> mainSubtypes = new LinkedHashMap<>();
		for (Map.Entry<String, List<String>> entry : subtypes.entrySet()) {
			mainSubtypes.put(
				getStrain(entry.getKey()),
				Collections.unmodifiableList(entry.getValue()));
		}
		this.mainSubtypes = Collections.unmodifiableMap(mainSubtypes);
	}

	private void initMutationTypes() {
		String raw = loadResource(mutTypesResPath, classLoader, null);
		mutationTypes = MutationType.loadJson(raw, virus);
	}

	private void initMutationTypePairs() {
		String raw = loadResource(mutTypePairsResPath, classLoader, "[]");
		mutationTypePairs = MutationTypePair.loadJson(raw, virus);
	}

	private void initStrains() {
		String raw = loadResource(strainsResPath, classLoader, null);
		this.strains = Strain.loadJson(raw, virus);
	}

	private void initGenes() {
		String raw = loadResource(genesResPath, classLoader, null);
		this.genes = Gene.loadJson(raw, virus);
	}

	private void initDrugClasses() {
		String raw = loadResource(drugClassesResPath, classLoader, "[]");
		this.drugClasses = DrugClass.loadJson(raw, virus);
	}

	private void initDrugs() {
		String raw = loadResource(drugsResPath, classLoader, "[]");
		this.drugs = Drug.loadJson(raw, virus);
	}

	private void initDrugResistAlgs() {
		String raw = loadResource(algorithmsIndexPath, classLoader, "[]");
		Map<String, List<List<String>>> algs = Json.loads(
			raw, new TypeToken<Map<String, List<List<String>>>>(){}.getType());
		List<DrugResistanceAlgorithm<T>> algList = new ArrayList<>();
		Map<String, DrugResistanceAlgorithm<T>> algMap = new LinkedHashMap<>();
		for (String family : algs.keySet()) {
			for (List<String> algData : algs.get(family)) {
				String version = algData.get(0);
				String publishDate = algData.get(1);
				String name = String.format("%s_%s", family, version);
				String xmlText = loadResource(
					String.format(algorithmsResPath, family, version),
					classLoader,
					null
				);
				DrugResistanceAlgorithm<T> alg = new DrugResistanceAlgorithm<>(
					name, family, version, publishDate, virus, xmlText);
				algList.add(alg);
				algMap.put(name, alg);
				algMap.put(alg.getEnumCompatName(), alg);
			}
		}
		this.drugResistAlgs = Collections.unmodifiableList(algList);
		this.drugResistAlgLookup = Collections.unmodifiableMap(algMap);
	}

	private void initGenotypes() {
		String raw = loadResource(genotypesRespath, classLoader, "{}");
		this.allGenotypes = Genotype.loadJson(raw, virus);
	}

	private void initGenotypeReferences() {
		String raw = loadResource(genotypeReferencesResPath, classLoader, "[]");
		this.allGenotypeReferences = GenotypeReference.loadJson(raw, virus);
	}

	private void initDrugResistMutations() {
		drugResistMutations = loadMutationSetByDrugClassFromRes(drmsResPath, getStrains());
	}

	private void initSurveilDrugResistMuts() {
		surveilDrugResistMuts = loadMutationSetByDrugClassFromRes(sdrmsResPath, getStrains());
	}

	private void initApobecMutations() {
		apobecMutations = loadMutationSetFromRes(apobecsResPath, getStrains());
	}

	private void initApobecDRMs() {
		apobecDRMs = loadMutationSetFromRes(apobecDRMsResPath, getStrains());
	}

	private void initRxSelectedMutations() {
		this.rxSelectedMutations = loadMutationSetByDrugClassFromRes(tsmsResPath, getStrains());
	}

	public String getName() {
		return virusName;
	}

	public Strain<T> getMainStrain() {
		return getStrain(mainStrainText);
	}

	public Collection<Strain<T>> getStrains() {
		if (strains == null) {
			initStrains();
		}
		return strains.values();
	}


	public Strain<T> getStrain(String name) {
		if (strains == null) {
			initStrains();
		}
		return AssertUtils.notNull(
			strains.get(name),
			"Strain \"%s\" not found", name
		);
	}


	public Collection<Gene<T>> getGenes(Strain<T> strain) {
		if (genes == null) {
			initGenes();
		}
		return (
			genes.values()
			.stream()
			.filter(gene -> gene.getStrain() == strain)
			.collect(Collectors.toCollection(LinkedHashSet::new))
		);
	}


	public Gene<T> getGene(String name) {
		if (genes == null) {
			initGenes();
		}
		return AssertUtils.notNull(
			genes.get(name),
			"Gene \"%s\" not found", name
		);
	}


	public Collection<DrugClass<T>> getDrugClasses() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return Sets.newLinkedHashSet(drugClasses.values());
	}


	public Map<String, DrugClass<T>> getDrugClassSynonymMap() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses;
	}


	public DrugClass<T> getDrugClass(String name) {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses.get(name);
	}


	public Collection<Drug<T>> getDrugs() {
		if (drugs == null) {
			initDrugs();
		}
		return Sets.newTreeSet(drugs.values());
	}


	public Map<String, Drug<T>> getDrugSynonymMap() {
		if (drugs == null) {
			initDrugs();
		}
		return drugs;
	}


	public Collection<DrugResistanceAlgorithm<T>> getDrugResistAlgorithms() {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return drugResistAlgs;
	}


	public Collection<DrugResistanceAlgorithm<T>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return (
			algorithmNames.stream()
			.map(name -> drugResistAlgLookup.get(name))
			.collect(Collectors.toList())
		);
	}



	public DrugResistanceAlgorithm<T> getDrugResistAlgorithm(String name) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return AssertUtils.notNull(
			drugResistAlgLookup.get(name),
			"Unable to locate algorithm %s", name
		);
	}


	public DrugResistanceAlgorithm<T> getDrugResistAlgorithm(String family, String version) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return AssertUtils.notNull(
			drugResistAlgLookup.get(String.format("%s_%s", family, version)),
			"Unable to locate algorithm %s_%s", family, version
		);
	}

	/**
	 * Extracts gene from mutText string
	 * @param mutText
	 * @return a Gene enum object
	 */

	public Gene<T> extractMutationGene(String mutText) {
		Gene<T> gene = null;
		Matcher m = mutationPattern.matcher(mutText);
		if (m.matches()) {
			try {
				gene = getGene(mainStrainText + m.group(2).toUpperCase());
			} catch (NullPointerException e) {
				throw new Mutation.InvalidMutationException(
					"Gene is not specified and also not found in the " +
					"given text: " + mutText + ". The correct format " +
					"for an input mutation string is, for example, " +
					"RT:215Y.", e);
			}
		}
		return gene;
	}

	/**
	 * Converts gene and mutText string into a Mutation object
	 * mutText may or may not have a preceding consensus
	 * @param gene, mutText
	 * @return a Mutation object
	 */

	public Mutation<T> parseMutationString(Gene<T> defaultGene, String mutText) {
		Matcher m = mutationPattern.matcher(mutText);
		Mutation<T> mut = null;
		if (m.matches()) {
			Gene<T> gene;
			boolean isASI = m.group(1) == null ? false : m.group(1).equals("__ASI__");
			try {
				gene = getGene(mainStrainText + m.group(2).toUpperCase());
			} catch (NullPointerException e) {
				if (defaultGene == null) {
					throw new Mutation.InvalidMutationException(
						"Gene is not specified and also not found in the " +
						"given text: " + mutText + ". The correct format " +
						"for an input mutation string is, for example, " +
						"RT:215Y.", e);
				}
				else {
					gene = defaultGene;
				}
			}
			int pos = Integer.parseInt(m.group(4));
			String aas = m.group(5);
			if (!isASI && !NON_ASI_AA_PATTERN.matcher(aas).matches()) {
				throw new Mutation.InvalidMutationException(
					"Tried to parse mutation string using invalid parameters: " + mutText);
			}
			aas = AAUtils.normalizeAAs(aas);
			String triplet = m.group(6);
			if (triplet == null) triplet = "";
			if (isASI) {
				mut = new AAMutation<>(gene, pos, aas.toCharArray());
			}
			else {
				mut = new CodonMutation<>(gene, pos, aas, triplet, "", 0xff);
			}
		} else {
			throw new Mutation.InvalidMutationException(
				"Tried to parse mutation string using invalid parameters: " + mutText);
		}
		return mut;
	}


	public Mutation<T> parseMutationString(String mutText) {
		return parseMutationString(null, mutText);
	}


	public MutationSet<T> newMutationSet(String formattedMuts) {
		return newMutationSet(null, formattedMuts);
	}


	public MutationSet<T> newMutationSet(Collection<String> formattedMuts) {
		return newMutationSet(null, formattedMuts);
	}

	/**
	 * Parse the given string then create mutation set.
	 *
	 * Supported delimiters:
	 * 	- space ( )
	 * 	- tabulation (\t)
	 * 	- new line (\n)
	 * 	- carriage return (\r)
	 * 	- dot (.)
	 * 	- comma (,)
	 * 	- semicolon (;)
	 * 	- plus (+)
	 *
	 * Supported single mutation formats:
	 * 	- with consensus (P1X)
	 * 	- without consensus (52R)
	 * 	- lowercase indels (69i or 44d)
	 * 	- dash/underscore indels (69_XX or 44-)
	 * 	- "hivdb" indels (69#XX or 44~)
	 * 	- word indels (69Insertion, 44Deletion)
	 * 	- stop codon (122*)
	 *
	 * All duplicated mutations are removed before returning.
	 *
	 * @param defaultGene
	 * @param fomattedMuts
	 * @return A list of Mutation objects
	 */

	public MutationSet<T> newMutationSet(Gene<T> defaultGene, String formattedMuts) {
		if (formattedMuts == null) {
			return new MutationSet<>();
		}
		return newMutationSet(
			defaultGene,
			Arrays.asList(formattedMuts.split("[\\s,;+\\.]+"))
		);
	}


	public MutationSet<T> newMutationSet(Gene<T> defaultGene, Collection<String> formattedMuts) {
		return MutationSet.parseString(
			defaultGene, formattedMuts, (gene, mStr) -> parseMutationString(gene, mStr));
	}


	public Map<DrugClass<T>, MutationSet<T>> getDrugResistMutations() {
		if (drugResistMutations == null) {
			initDrugResistMutations();
		}
		return drugResistMutations;
	}


	public Map<DrugClass<T>, MutationSet<T>> getSurveilDrugResistMutations() {
		if (surveilDrugResistMuts == null) {
			initSurveilDrugResistMuts();
		}
		return surveilDrugResistMuts;
	}


	public Map<DrugClass<T>, MutationSet<T>> getRxSelectedMutations() {
		if (rxSelectedMutations == null) {
			initRxSelectedMutations();
		}
		return rxSelectedMutations;
	}


	public MutationSet<T> getApobecMutations() {
		if (apobecMutations == null) {
			initApobecMutations();
		}
		return apobecMutations;
	}


	public MutationSet<T> getApobecDRMs() {
		if (apobecDRMs == null) {
			initApobecDRMs();
		}
		return apobecDRMs;
	}


	public Collection<MutationType<T>> getMutationTypes() {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.values();
	}


	public MutationType<T> getMutationType(String mutTypeText) {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.get(mutTypeText);
	}


	public Collection<MutationTypePair<T>> getMutationTypePairs() {
		if (mutationTypePairs == null) {
			initMutationTypePairs();
		}
		return mutationTypePairs;
	}

	/**
	 * Get an AminoAcidPercents instance
	 *
	 * @param treatment "all", "naive" or "art"
	 * @param subtype "all", "A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG", "other"
	 */

	public AminoAcidPercents<T> getAminoAcidPercents(Strain<T> strain, String treatment, String subtype) {
		String resourceName = String.format(aaPcntsResPath, treatment, subtype);
		String resourceKey = String.format("%s::%s", resourceName, strain.getName());
		if (!aminoAcidPcnts.containsKey(resourceKey)) {
			aminoAcidPcnts.put(resourceKey, new AminoAcidPercents<>(resourceName, virus, strain));
			// Example of empty Instance:
			// aminoAcidPcnts.put(resourceName, AminoAcidPercents.newEmptyInstance());
		}
		return aminoAcidPcnts.get(resourceKey);
	}

	/**
	 * Get a CodonPercents instance
	 *
	 * @param treatment "all", "naive" or "art"
	 * @param subtype "all", "A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG"
	 */

	public CodonPercents<T> getCodonPercents(Strain<T> strain, String treatment, String subtype) {
		String resourceName = String.format(codonPcntsResPath, treatment, subtype);
		if (!codonPcnts.containsKey(resourceName)) {
			codonPcnts.put(resourceName, new CodonPercents<>(resourceName, virus, strain));
			// Example of emptyInstance:
			// codonPcnts.put(resourceName, CodonPercents.newEmptyInstance());
		}
		return codonPcnts.get(resourceName);
	}


	public List<MutationPrevalence<T>> getMutationPrevalence(GenePosition<T> genePos) {
		if (!mutPrevalenceMap.containsKey(genePos)) {
			mutPrevalenceMap.put(genePos, virus.defaultGetMutationPrevalence(genePos));
		}
		return new ArrayList<>(mutPrevalenceMap.get(genePos));
	}


	public ConditionalComments<T> getConditionalComments() {
		if (condComments == null) {
			initCondComments();
		}
		return condComments;
	}


	public List<String> getMainSubtypes(Strain<T> strain) {
		if (mainSubtypes == null) {
			initMainSubtypes();
		}
		return mainSubtypes.getOrDefault(strain, Collections.emptyList());
	}


	public Map<Gene<T>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<T> strain) {
		if (!allAAPcntsNumPatients.containsKey(strain)) {
			allAAPcntsNumPatients.put(strain, virus.defaultGetNumPatientsForAAPercents(strain));
		}
		return allAAPcntsNumPatients.get(strain);
	}


	public Collection<Genotype<T>> getGenotypes() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.values();
	}


	public Genotype<T> getGenotype(String name) {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get(name);
	}


	public Genotype<T> getGenotypeUnknown() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get("U");
	}


	public List<GenotypeReference<T>> getGenotypeReferences() {
		if (allGenotypeReferences == null) {
			initGenotypeReferences();
		}
		return allGenotypeReferences;
	}


	public Genotyper<T> getGenotyper() {
		if (genotyper == null) {
			genotyper = new Genotyper<>(virus);
		}
		return genotyper;
	}

	public AlignmentConfig<T> getAlignmentConfig() {
		if (alignmentConfig == null) {
			String raw = loadResource(alignConfigResPath, classLoader, null);
			alignmentConfig = AlignmentConfig.loadJson(raw, virus);
		}
		return alignmentConfig;
	}

	public Map<Strain<T>, SequenceReadsAssembler<T>> getSeqReadsAssemblers() {
		if (this.seqReadsAssemblers == null) {
			String raw = loadResource(assemblyConfigResPath, classLoader, null);
			Map<Strain<T>, SequenceReadsAssembler<T>> seqReadsAssemblers = SequenceReadsAssembler.loadJson(
				raw,
				virus
			);
			this.seqReadsAssemblers = seqReadsAssemblers;
		}
		return this.seqReadsAssemblers;
	}

	public Map<Strain<T>, SequenceAssembler<T>> getSequenceAssemblers() {
		if (sequenceAssemblers == null) {
			String raw = loadResource(assemblyConfigResPath, classLoader, null);
			Map<Strain<T>, SequenceAssembler<T>> sequenceAssemblers = SequenceAssembler.loadJson(
				raw,
				virus);
			this.sequenceAssemblers = sequenceAssemblers;
		}
		return sequenceAssemblers;
	}

	public Boolean purgeCache() {
		// check if purging cache is allowed
		try {
			URL url = new URL(URL_ALLOW_PURGING_CACHE);
			HttpURLConnection conn = (HttpURLConnection) url.openConnection();
			conn.connect();
			int code = conn.getResponseCode();
			if (code != 200) {
				return false;
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		drugResistMutations = null;

		return true;
	}

}
