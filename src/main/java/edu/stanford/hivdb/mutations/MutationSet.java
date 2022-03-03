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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.base.Predicates;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

public class MutationSet<VirusT extends Virus<VirusT>> extends TreeSet<Mutation<VirusT>> implements Comparable<MutationSet<VirusT>> {

	private static final String MESSAGE_ON_WRITE =
		"Modification to MutationSet is not allowed.";

	private static final long serialVersionUID = -1835692164622753L;

	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> loadJson(
		String raw, Function<String, Gene<VirusT>> getGene
	) {
		List<Map<String, ?>> muts = Json.loads(
			raw, new TypeToken<List<Map<String, ?>>>(){});
		return loadJsonMap(muts, getGene);
	}

	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> loadJsonMap(
		List<Map<String, ?>> muts, Function<String, Gene<VirusT>> getGene
	) {
		return new MutationSet<VirusT>(
			muts.stream()
			.map(m -> new AAMutation<VirusT>(
				getGene.apply((String) m.get("gene")),
				((Double) m.get("position")).intValue(),
				((String) m.get("aa")).toCharArray()
			))
			.collect(Collectors.toList()));
	}

	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> parseString(
		VirusT virusInstance,
		String formattedMuts
	) {
		return parseString(
			null, formattedMuts,
			(gene, mutStr) -> virusInstance.parseMutationString(gene, mutStr));
	}

	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> parseString(
		VirusT virusInstance,
		Collection<String> formattedMuts
	) {
		return parseString(
			null, formattedMuts,
			(gene, mutStr) -> virusInstance.parseMutationString(gene, mutStr));
	}

	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> parseString(
		Gene<VirusT> defaultGene,
		Collection<String> formattedMuts
	) {
		return parseString(
			defaultGene, formattedMuts, null
		);
	}
	
	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> parseString(
			Gene<VirusT> defaultGene,
			String formattedMuts
		) {
			if (formattedMuts == null) {
				return new MutationSet<>();
			}
			return parseString(
				defaultGene,
				Arrays.asList(formattedMuts.split("[\\s,;+\\.]+")),
				null
			);
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
	 * @param <VirusT>			A virus subtype
	 * @param defaultGene		Gene
	 * @param formattedMuts		Formatted mutations separated by comma or colon or space or dot
	 * @param mutationParser	Mutation parsing function
	 * @return 					A list of Mutation objects
	 */
	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> parseString(
		Gene<VirusT> defaultGene,
		String formattedMuts,
		BiFunction<Gene<VirusT>, String, Mutation<VirusT>> mutationParser
	) {
		if (formattedMuts == null) {
			return new MutationSet<>();
		}
		return parseString(
			defaultGene,
			Arrays.asList(formattedMuts.split("[\\s,;+\\.]+")),
			mutationParser
		);
	}

	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> parseString(
		Gene<VirusT> defaultGene,
		Collection<String> formattedMuts,
		BiFunction<Gene<VirusT>, String, Mutation<VirusT>> mutationParser
	) {
		if (defaultGene == null && mutationParser == null) {
			throw new NullPointerException("defaultGene and mutationParser can not be both null");
		}
		BiFunction<Gene<VirusT>, String, Mutation<VirusT>> localMutParser;
		if (mutationParser == null) {
			localMutParser = (
				(gene, mutStr) -> defaultGene.getVirusInstance().parseMutationString(gene, mutStr)
			);
		}
		else {
			localMutParser = mutationParser;
		}
		return new MutationSet<>(
			formattedMuts
			.stream()
			.filter(mStr -> mStr.length() > 0)
			.map(mStr -> localMutParser.apply(defaultGene, mStr))
			.collect(Collectors.toList())
		);
	}

	private transient Map<GenePosition<VirusT>, Mutation<VirusT>> genePositionMap;

	private transient Map<Mutation<VirusT>, List<MutationPrevalence<VirusT>>> prevalences;

	public MutationSet(Collection<Mutation<VirusT>> mutations) {
		genePositionMap = new TreeMap<>();
		_addAll(mutations);
		// prevent any future changes to genePositionMap
		genePositionMap = Collections.unmodifiableMap(genePositionMap);
	}

	@SafeVarargs
	public MutationSet(Mutation<VirusT>... mutations) {
		this(Arrays.asList(mutations));
	}

	public MutationSet<VirusT> displayAmbiguities() {
		List<Mutation<VirusT>> tmpMuts = new ArrayList<>();
		for (Mutation<VirusT> mut : this) {
			tmpMuts.add(new AAMutation<>(
				mut.getGene(), mut.getPosition(),
				mut.getAAChars(), 0xff));
		}
		return new MutationSet<>(tmpMuts);
	}

	// private addAll method for internal usage
	private boolean _addAll(Collection<Mutation<VirusT>> muts) {
		List<Boolean> rr = muts
			.stream().map(mut -> _add(mut)).collect(Collectors.toList());
		return rr.stream().anyMatch(r -> r);
	}

	// private add method for internal usage
	private boolean _add(Mutation<VirusT> mut) {
		if (mut == null) {
			return false;
		}
		GenePosition<VirusT> gp = mut.getGenePosition();
		Mutation<VirusT> origMut = genePositionMap.getOrDefault(gp, null);
		if (mut.equals(origMut) || mut.getReference().equals(mut.getAAs())) {
			return false;
		}
		if (origMut != null) {
			mut = origMut.mergesWith(mut.getAAChars());
			super.remove(origMut);
		}
		super.add(mut);
		genePositionMap.put(gp, mut);
		return true;
	}

	// Begin of all write methods
	@Override
	@Deprecated
	public boolean addAll(Collection<? extends Mutation<VirusT>> muts) {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}

	@Override
	@Deprecated
	public boolean add(Mutation<VirusT> mut) {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}

	@Override
	@Deprecated
	public void clear() {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}

	@Override
	@Deprecated
	public boolean removeAll(Collection<?> muts) {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}

	@Override
	@Deprecated
	public boolean retainAll(Collection<?> muts) {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}

	@Override
	@Deprecated
	public boolean remove(Object m) {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}

	@Override
	@Deprecated
	public Mutation<VirusT> pollFirst() {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}

	@Override
	@Deprecated
	public Mutation<VirusT> pollLast() {
		throw new UnsupportedOperationException(MESSAGE_ON_WRITE);
	}
	
	// End of all write methods

	/**
	 * Merges with another MutationSet and returns the result.
	 *
	 * Example:
	 *   self = new MutationSet(RT, "48VER");
	 *   another = new MutationSet(RT, "48A,48L,36E");
	 *   assertEquals(
	 *     new MutationSet(RT, "48AELRV,36E"),
	 *     self.mergesWith(another));
	 *
	 * @param another	Mutations
	 * @return 			A new MutationSet object contains matched mutations
	 */
	public MutationSet<VirusT> mergesWith(Collection<Mutation<VirusT>> another) {
		List<Mutation<VirusT>> newList = new ArrayList<>(this);
		newList.addAll(another);
		return new MutationSet<>(newList);
	}

	/**
	 * Just a shortcut to mergesWith(Collection).
	 * 
	 * @param mutation 	Mutation
	 * @return			Mutation set
	 */
	public MutationSet<VirusT> mergesWith(Mutation<VirusT> mutation) {
		return mergesWith(new MutationSet<>(mutation));
	}

	/**
	 * Intersects with another MutationSet and returns the result.
	 *
	 * Example:
	 *   self = new MutationSet(RT, "48VER");
	 *   another = new MutationSet(RT, "48E,48AR,36E");
	 *   assertEquals(
	 *     new MutationSet(RT, "48ER"),
	 *     self.intersectsWith(another));
	 *
	 * @param another	Mutations
	 * @return 			A new MutationSet object contains matched mutations
	 */
	public MutationSet<VirusT> intersectsWith(Collection<Mutation<VirusT>> another) {
		Set<GenePosition<VirusT>> gpKeys = new HashSet<>(genePositionMap.keySet());
		Set<GenePosition<VirusT>> gpKeysAnother;
		if (another instanceof MutationSet) {
			gpKeysAnother = ((MutationSet<VirusT>) another).genePositionMap.keySet();
		} else {
			gpKeysAnother = another
				.stream()
				.map(m -> m.getGenePosition())
				.collect(Collectors.toSet());
		}
		gpKeys.retainAll(gpKeysAnother);
		List<Mutation<VirusT>> mutations = new ArrayList<>();
		for (GenePosition<VirusT> gp : gpKeys) {
			Mutation<VirusT> thisMut = genePositionMap.get(gp);
			Set<Character> otherAAChars = new TreeSet<>();
			if (another instanceof MutationSet) {
				otherAAChars.addAll(
					((MutationSet<VirusT>) another).genePositionMap.get(gp).getAAChars()
				);
			}
			else {
				for (Mutation<VirusT> mut : another) {
					if (mut.getGenePosition().equals(gp)) {
						otherAAChars.addAll(mut.getAAChars());
					}
				}
			}
			mutations.add(thisMut.intersectsWith(otherAAChars));
		}
		return new MutationSet<>(mutations);
	}

	/**
	 * Just a shortcut to intersectsWith(Collection).
	 * 
	 * @param mutation 	Mutation
	 * @return			Mutationset
	 */
	public MutationSet<VirusT> intersectsWith(Mutation<VirusT> mutation) {
		return intersectsWith(new MutationSet<>(mutation));
	}

	/**
	 * New MutationSet with elements in this but not in another.
	 * 
	 * @param another	Mutations
	 * @return			Mutation set
	 */
	public MutationSet<VirusT> subtractsBy(Collection<Mutation<VirusT>> another) {
		MutationSet<VirusT> anotherSet;
		if (another instanceof MutationSet) {
			anotherSet = (MutationSet<VirusT>) another;
		}
		else {
			anotherSet = new MutationSet<>(another);
		}
		List<Mutation<VirusT>> mutations = new ArrayList<>();
		for (GenePosition<VirusT> gp : genePositionMap.keySet()) {
			Mutation<VirusT> thisMut = genePositionMap.get(gp);
			Mutation<VirusT> anotherMut = anotherSet.genePositionMap.get(gp);
			if (anotherMut == null) {
				mutations.add(thisMut);
			}
			else {
				mutations.add(thisMut.subtractsBy(anotherMut.getAAChars()));
			}
		}
		return new MutationSet<>(mutations);
	}

	public MutationSet<VirusT> subtractsBy(Mutation<VirusT> mutations) {
		return subtractsBy(new MutationSet<>(mutations));
	}

	public MutationSet<VirusT> filterBy(Predicate<Mutation<VirusT>> predicate) {
		return new MutationSet<>(
			this
			.getSplitted()
			.stream()
			.filter(predicate)
			.collect(Collectors.toList())
		);
	}

	public MutationSet<VirusT> filterByNoSplit(Predicate<Mutation<VirusT>> predicate) {
		return new MutationSet<>(
			this
			.stream()
			.filter(predicate)
			.collect(Collectors.toList())
		);
	}
	
	public Long count() {
		return countIf(Predicates.alwaysTrue());
	}
	
	public Long countIf(Predicate<Mutation<VirusT>> predicate) {
		List<Mutation<VirusT>> muts = this
			.getSplitted()
			.stream()
			.filter(predicate)
			.sorted(Mutation::compare)
			.collect(Collectors.toList());
		long counter = 0;
		int mutSize = muts.size();
		for (int i = 0; i < mutSize; i ++) {
			if (i == 0) {
				counter += 1;
				continue;
			}
			Mutation<VirusT> mut = muts.get(i);
			if (!mut.isDeletion()) {
				counter += 1;
				continue;
			}
			Mutation<VirusT> prevMut = muts.get(i - 1);
			if (!prevMut.isDeletion()) {
				counter += 1;
				continue;
			}
			if (prevMut.getPosition() + 1 < mut.getPosition()) {
				counter += 1;
				continue;
			}
		}
		return counter;
	}

	public <T> Map<T, MutationSet<VirusT>> filterAndGroupBy(
			Predicate<Mutation<VirusT>> predicate, Function<Mutation<VirusT>, T> function) {
		Map<T, List<Mutation<VirusT>>> tmp = this
			.stream()
			.filter(predicate)
			.collect(Collectors.groupingBy(function));
		Map<T, MutationSet<VirusT>> r = new TreeMap<>();
		for (Map.Entry<T, List<Mutation<VirusT>>> e : tmp.entrySet()) {
			r.put(e.getKey(), new MutationSet<>(e.getValue()));
		}
		return r;
	}

	public <T> Map<T, MutationSet<VirusT>> groupBy(Function<Mutation<VirusT>, T> function) {
		Map<T, List<Mutation<VirusT>>> tmp = this
			.stream()
			.collect(Collectors.groupingBy(function));
		Map<T, MutationSet<VirusT>> r = new TreeMap<>();
		for (Map.Entry<T, List<Mutation<VirusT>>> e : tmp.entrySet()) {
			r.put(e.getKey(), new MutationSet<>(e.getValue()));
		}
		return r;
	}

	/**
	 * Filter mutations by gene then groups the result by their primary MutType.
	 *
	 * @param gene	Gene
	 * @return 		Map&lt;MutType, MutationSet&gt;
	 */
	public SortedMap<MutationType<VirusT>, MutationSet<VirusT>> groupByMutType(final Gene<VirusT> gene) {
		SortedMap<MutationType<VirusT>, MutationSet<VirusT>> r = new TreeMap<>(
			filterAndGroupBy(
				mut -> mut.getGene() == gene,
				mut -> mut.getPrimaryType()
			)
		);
		Set<MutationType<VirusT>> mutTypes = gene.getMutationTypes();
		if (r.size() < mutTypes.size()) {
			for (MutationType<VirusT> mt : mutTypes) {
				r.putIfAbsent(mt, new MutationSet<>());
			}
		}
		return r;
	}

	/**
	 * Filter mutations by given mutation type.
	 *
	 * @param mutType	Mutation type
	 * @return 			MutationSet
	 */
	public MutationSet<VirusT> getByMutType(final MutationType<VirusT> mutType) {
		return filterBy(mut -> mut.getPrimaryType() == mutType);
	}

	/**
	 * Groups mutations by their gene.
	 *
	 * @return Map&gt;Gene, MutationSet&gt;
	 */
	public Map<Gene<VirusT>, MutationSet<VirusT>> groupByGene() {
		return groupBy(mut -> mut.getGene());
	}

	/**
	 * Returns only mutations of a specific gene.
	 * @param gene	Gene
	 * @return 		A MutationSet contains all mutations of the given gene
	 */
	public MutationSet<VirusT> getGeneMutations(Gene<VirusT> gene) {
		return filterBy(mut -> mut.getGene() == gene);
	}

	/**
	 * Returns only mutations of a specific gene.
	 * @param gene	Gene
	 * @return 		A MutationSet contains all mutations of the given gene
	 */
	public MutationSet<VirusT> getGeneMutationsNoSplit(Gene<VirusT> gene) {
		return filterByNoSplit(mut -> mut.getGene() == gene);
	}

	/**
	 * Returns only mutations which are insertions.
	 * @return A new MutationSet instance
	 */
	public MutationSet<VirusT> getInsertions() {
		return filterBy(mut -> mut.isInsertion());
	}

	/**
	 * Returns only mutations which are deletions.
	 * @return A new MutationSet instance
	 */
	public MutationSet<VirusT> getDeletions() {
		return filterBy(mut -> mut.isDeletion());
	}

	/**
	 * Returns only mutations which contain stop codons.
	 * @return A new MutationSet instance
	 */
	public MutationSet<VirusT> getStopCodons() {
		return filterBy(Mutation::hasStop);
	}

	/**
	 * Returns only mutations which contain ambiguous codons.
	 * @return A new MutationSet instance
	 */
	public MutationSet<VirusT> getAmbiguousCodons() {
		return filterBy(mut -> mut.isAmbiguous());
	}

	/**
	 * Returns mutations that contain at least an unusual AA.
	 * @return A new MutationSet instance
	 */
	public MutationSet<VirusT> getUnusualMutations() {
		return filterBy(mut -> mut.isUnusual());
	}

	/**
	 * Creates a map of mutations and their prevalence in HIVDB.
	 * When a mutation has more than one amino acid, the prevalence of the
	 * most prevalent mutation is returned.
	 *
	 * @return Map&lt;Mutation, Double&gt; mutPrevalences
	 */
	public Map<Mutation<VirusT>, Double> getHighestMutPrevalences() {
		Map<Mutation<VirusT>, Double> r = new LinkedHashMap<>();
		for (Mutation<VirusT> mut : this) {
			r.put(mut, mut.getHighestMutPrevalence());
		}
		return r;
	}

	/**
	 * Returns only mutations which are at DR positions.
	 * @return A new MutationSet instance
	 */
	public MutationSet<VirusT> getAtDRPMutations() {
		return filterBy(mut -> mut.isAtDrugResistancePosition());
	}

	/**
	 * Looks up mutation types for each given mutations.
	 *
	 * @param muts	Mutation set
	 * @return 			Map&lt;Mutation, List&lt;MutationType&gt;&gt; Map of types for each mutation.
	 */
	public Map<Mutation<VirusT>, List<MutationType<VirusT>>> getMutationTypes(MutationSet<VirusT> muts) {
		Map<Mutation<VirusT>, List<MutationType<VirusT>>> r = new LinkedHashMap<>();
		for (Mutation<VirusT> mut: muts) {
			List<MutationType<VirusT>> types = mut.getTypes();
			if (types.size() > 0) {
				r.put(mut, types);
			}
		}
		return r;
	}
	
	public MutationSet<VirusT> getTSMs() {
		return filterBy(mut -> mut.isTSM());
	}
	
	public MutationSet<VirusT> getTSMs(DrugClass<VirusT> drugClass) {
		return filterBy(mut -> mut.getTSMDrugClass() == drugClass);
	}


	/**
	 * Returns only mutations which are DRMs.
	 * @return A new MutationSet instance
	 */
	public MutationSet<VirusT> getDRMs() {
		return filterBy(mut -> mut.isDRM());
	}

	/**
	 * Returns only mutations which are DRMs.
	 * 
	 * @param drugClass Drugclass
	 * @return 			A new MutationSet instance
	 */
	public MutationSet<VirusT> getDRMs(DrugClass<VirusT> drugClass) {
		return filterBy(mut -> mut.getDRMDrugClass() == drugClass);
	}

	public MutationSet<VirusT> getSDRMs() {
		return filterBy(mut -> mut.isSDRM());
	}

	/**
	 * Returns only mutations which are SDRMs.
	 * 
	 * @param drugClass	Drug class
	 * @return			Drug class related SDRMs
	 */
	public MutationSet<VirusT> getSDRMs(DrugClass<VirusT> drugClass) {
		return filterBy(mut -> mut.getSDRMDrugClass() == drugClass);
	}
	
	public MutationSet<VirusT> getApobecMutations() {
		return filterBy(mut -> mut.isApobecMutation());
	}
	
	public MutationSet<VirusT> getApobecDRMs() {
		return filterBy(mut -> mut.isApobecDRM());
	}

	/** Returns a mutation at specified gene position.
	 *
	 * @param gp a GenePosition object
	 * @return The matched mutation
	 */
	public Mutation<VirusT> get(GenePosition<VirusT> gp) {
		return genePositionMap.get(gp);
	}

	/** Returns a mutation at specified gene position.
	 *
	 * @param gene	Gene
	 * @param pos	Position
	 * @return 		The matched mutation
	 */
	public Mutation<VirusT> get(Gene<VirusT> gene, int pos) {
		GenePosition<VirusT> gp = new GenePosition<VirusT>(gene, pos);
		return genePositionMap.get(gp);
	}

	/** Returns a set of non-mixture mutations for all mutations.
	 *
	 * @return The mutation set
	 */
	public Set<Mutation<VirusT>> getSplitted() {
		Set<Mutation<VirusT>> splittedMuts = new TreeSet<>();
		for(Mutation<VirusT> mut : genePositionMap.values()) {
			splittedMuts.addAll(mut.split());
		}
		return splittedMuts;
	}

	/** Returns list of mutation positions
	 *
	 * @return a list of mutation positions
	 */
	public Set<GenePosition<VirusT>> getPositions() {
		return new TreeSet<>(genePositionMap.keySet());
	}

	/** Check if the given position is an insertion
	 *
	 * @param gene		Gene
	 * @param pos		Position
	 * @return Boolean	Has Insertion at this position
	 */
	public boolean hasInsertionAt(Gene<VirusT> gene, int pos) {
		Mutation<VirusT> mut = get(gene, pos);
		return mut == null ? false : mut.isInsertion();
	}

	/** Check if the given position is a deletion
	 *
	 * @param gene		Gene
	 * @param pos		Position
	 * @return Boolean	Has deletion at this position
	 */
	public boolean hasDeletionAt(Gene<VirusT> gene, int pos) {
		Mutation<VirusT> mut = get(gene, pos);
		return mut == null ? false : mut.isDeletion();
	}

	/**
	 * Returns true if contains any mutation that its AAs shared with mutRef.
	 *
	 * @param anotherMut	Another mutation
	 * @return				Has shared amino acid mutation
	 */
	public boolean hasSharedAAMutation(Mutation<VirusT> anotherMut) {
		return hasSharedAAMutation(anotherMut, true);
	}
	
	/**
	 * Returns true if contains any mutation that its AAs shared with mutRef.
	 *
	 * @param anotherMut		Another Mutation
	 * @param ignoreRefOrStops	Ignore reference mutation or stop mutation
	 * @return					Has shared Amino acid mutation
	 */
	public boolean hasSharedAAMutation(Mutation<VirusT> anotherMut, boolean ignoreRefOrStops) {
		GenePosition<VirusT> gp = anotherMut.getGenePosition();
		Mutation<VirusT> selfMut = genePositionMap.get(gp);
		if (selfMut == null) {
			return false;
		}
		return selfMut.containsSharedAA(anotherMut.getAAChars(), ignoreRefOrStops);
	}

	public Map<Mutation<VirusT>, List<MutationPrevalence<VirusT>>> getPrevalences() {
		if (this.prevalences == null) {
			Map<Mutation<VirusT>, List<MutationPrevalence<VirusT>>> prevalences = new TreeMap<>();
			for (Mutation<VirusT> mut : this) {
				prevalences.put(mut, mut.getPrevalences());
			}
			this.prevalences = Collections.unmodifiableMap(prevalences);
		}
		return this.prevalences;
	}

	public String join(
			CharSequence delimiter,
			Function<Mutation<VirusT>, CharSequence> mutationToString) {
		if (this.isEmpty()) {
			// to be compatible with the old output
			return "None";
		}
		return this
			.stream()
			.map(mutationToString)
			.collect(Collectors.joining(delimiter));
	}

	public String join(
			Character delimiter,
			Function<Mutation<VirusT>, CharSequence> mutationToString) {
		return join("" + delimiter, mutationToString);
	}

	public String join(CharSequence delimiter) {
		return join(delimiter, Mutation<VirusT>::getHumanFormat);
	}

	public String join(Character delimiter) {
		return join("" + delimiter, Mutation<VirusT>::getHumanFormat);
	}

	public String join(Function<Mutation<VirusT>, CharSequence> mutationToString) {
		return join(",", mutationToString);
	}

	public String join() {
		return join(",", Mutation<VirusT>::getHumanFormat);
	}

	public List<String> toStringList(
			Function<Mutation<VirusT>, CharSequence> mutationToString) {
		return this
			.stream()
			.map(mut -> mutationToString.apply(mut).toString())
			.collect(Collectors.toList());
	}

	public List<String> toStringList() {
		return toStringList(Mutation<VirusT>::getHumanFormat);
	}

	public List<String> toASIFormat() {
		return toStringList(Mutation<VirusT>::getASIFormat);
	}

	@Override
	public int compareTo(MutationSet<VirusT> o) {
		List<Mutation<VirusT>> mutsLeft = new ArrayList<>(getSplitted());
		List<Mutation<VirusT>> mutsRight = new ArrayList<>(o.getSplitted());
		int minSize = Math.min(mutsLeft.size(), mutsRight.size());
		for (int i = 0; i < minSize; i ++) {
			Mutation<VirusT> m1 = mutsLeft.get(i);
			Mutation<VirusT> m2 = mutsRight.get(i);
			int cmp = m1.compareTo(m2);
			if (cmp != 0) {
				return cmp;
			}
		}
		return mutsLeft.size() - mutsRight.size();
	}
}
