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

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.io.IOUtils;

import com.google.common.collect.Lists;
import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;


/**
 * There are two public methods: getHighestMutPrevalence, unusualMutations
 *
 */
public class AminoAcidPercents<VirusT extends Virus<VirusT>> {

	final static protected Gson gson = new Gson();

	final private List<AminoAcidPercent<VirusT>> aminoAcidPcnts;
	final private Map<GenePosition<VirusT>, Map<Character, AminoAcidPercent<VirusT>>> aminoAcidPcntMap = new HashMap<>();


	/**
	 * AminoAcidPercents initializer
	 *
	 * @param resourceName	Amino acid percent file name
	 * @param virusInstance Virus instance
	 * @param strain		Virus strain
	 *
	 */
	public AminoAcidPercents(String resourceName, VirusT virusInstance, Strain<VirusT> strain) {

		try (
			InputStream stream = this
				.getClass().getClassLoader()
				.getResourceAsStream(resourceName);
		) {
			String raw = IOUtils.toString(stream, StandardCharsets.UTF_8);
			List<Map<String, ?>> aminoAcidPcntData = gson.fromJson(
				raw, new TypeToken<List<Map<String, ?>>>(){}.getType());
			List<AminoAcidPercent<VirusT>> aminoAcidPcnts = (
				aminoAcidPcntData.stream()
				.map(aaPcnt -> new AminoAcidPercent<>(
					strain.getGene((String) aaPcnt.get("gene")),
					((Double) aaPcnt.get("position")).intValue(),
					((String) aaPcnt.get("aa")).charAt(0),
					(Double) aaPcnt.get("percent"),
					((Double) aaPcnt.get("count")).intValue(),
					((Double) aaPcnt.get("total")).intValue(),
					(String) aaPcnt.get("reason"),
					(Boolean) aaPcnt.get("isUnusual")
				))
				.collect(Collectors.toList())
			);
			this.aminoAcidPcnts = Collections.unmodifiableList(aminoAcidPcnts);
		} catch (IOException|NullPointerException e) {
			throw new ExceptionInInitializerError(
				String.format("Invalid resource name (%s)", resourceName)
			);
		}

		for (AminoAcidPercent<VirusT> aaPcnt : aminoAcidPcnts) {
			GenePosition<VirusT> gp = aaPcnt.getGenePosition();
			aminoAcidPcntMap.putIfAbsent(gp, new LinkedHashMap<>());
			aminoAcidPcntMap.get(gp).put(aaPcnt.getAA(), aaPcnt);
		}
	}
	
	public static <VirusT extends Virus<VirusT>> AminoAcidPercents<VirusT> newEmptyInstance() {
		return new AminoAcidPercents<>();
	}
	
	/**
	 * Initialize an empty AminoAcidPercents
	 * 
	 */
	private AminoAcidPercents() {
		aminoAcidPcnts = Collections.emptyList();
	}

	public List<AminoAcidPercent<VirusT>> get() {
		// make a copy in case of any modification
		return new ArrayList<>(aminoAcidPcnts);
	}

	public List<AminoAcidPercent<VirusT>> get(Gene<VirusT> gene) {
		return (aminoAcidPcnts
				.stream().filter(aap -> aap.getGene().equals(gene))
				.collect(Collectors.toList()));
	}

	public List<AminoAcidPercent<VirusT>> get(Gene<VirusT> gene, int pos) {
		return new ArrayList<>(
			aminoAcidPcntMap.get(new GenePosition<>(gene, pos))
			.values());
	}

	public List<AminoAcidPercent<VirusT>> get(GenePosition<VirusT> genePos) {
		return new ArrayList<>(
			aminoAcidPcntMap.getOrDefault(genePos, Collections.emptyMap())
			.values());
	}

	public AminoAcidPercent<VirusT> get(Gene<VirusT> gene, int pos, char aa) {
		return get(new GenePosition<>(gene, pos), aa);
	}

	public AminoAcidPercent<VirusT> get(GenePosition<VirusT> genePos, char aa) {
		Map<Character, AminoAcidPercent<VirusT>> posAAPcnts = aminoAcidPcntMap
			.getOrDefault(genePos, Collections.emptyMap());
		if (posAAPcnts.isEmpty()) {
			return null;
		}
		AminoAcidPercent<VirusT> aaPcnt = posAAPcnts.get(aa);
		if (aaPcnt != null) {
			return aaPcnt;
		}
		else {
			int total = Lists.newArrayList(posAAPcnts.values()).get(0).getTotal();
			return new AminoAcidPercent<>(
				genePos.getGene(), genePos.getPosition(),
				aa, .0, 0, total, "PCNT", true);
		}
	}

	/**
	 * Returns the highest amino acid prevalence associated with each of
	 * the AA in a mixture.
	 *
	 * @param gene	Gene
	 * @param pos	Position
	 * cons 	consensus at the position
	 * @param mixture Mixture
	 *
	 * @return Double highest amino acid prevalence
	 */
	public Double getHighestAAPercentValue(
		Gene<VirusT> gene, int pos, /* char cons,*/ String mixture
	) {
		Double pcntVal = 0.0;
		GenePosition<VirusT> gpos = new GenePosition<>(gene, pos);

		for (char aa : mixture.toCharArray()) {
			/* if (aa == cons || aa == '*') {
				// ignore consensus and stop codon
				continue;
			} */
			double aaPcntVal = aminoAcidPcntMap.get(gpos).get(aa).getPercent();
			pcntVal = Math.max(pcntVal, aaPcntVal);
		}
		return pcntVal;
	}
	
	/**
	 * Finds unusual AAs from given AAs
	 * 
	 * @param gene
	 * @param pos
	 * @param aas
	 * @return Set
	 */
	public Set<Character> getUnusualAAs(Gene<VirusT> gene, int pos, Collection<Character> aas) {
		GenePosition<VirusT> gpos = new GenePosition<VirusT>(gene, pos);
		return aas.stream()
			.filter(aa -> {
				AminoAcidPercent<VirusT> aaPcnt = get(gpos, aa);
				return aaPcnt == null || aaPcnt.isUnusual();
			})
			.collect(Collectors.toCollection(TreeSet::new));		
	}

	/**
	 * Returns true if the given mutation contains any unusual AA
	 *
	 * @param gene	Gene
	 * @param pos 	Position
	 * @param aas	Amino acids
	 * 
	 * @return true if contains unusual AA
	 */
	public Boolean containsUnusualAA(Gene<VirusT> gene, int pos, String aas) {
		GenePosition<VirusT> gpos = new GenePosition<VirusT>(gene, pos);
		for (char aa : aas.toCharArray()) {
			AminoAcidPercent<VirusT> aaPcnt = get(gpos, aa);
			if (aaPcnt == null || aaPcnt.isUnusual()) {
				return true;
			}
		}
		return false;
	}

}
