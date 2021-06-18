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
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.io.IOUtils;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;


/**
 * There are two public methods: getHighestMutPrevalence, unusualMutations
 *
 */
public class CodonPercents<VirusT extends Virus<VirusT>> {

	final static protected Gson gson = new Gson();

	final protected List<CodonPercent<VirusT>> codonPcnts;
	final private Map<GenePosition<VirusT>, Map<String, CodonPercent<VirusT>>> codonPcntMap = new HashMap<>();
	
	/**
	 * CodonPercents initializer
	 *
	 * @param resourceName	Codon percent resource file name
	 * @param virusInstance Virus instance
	 * @param strain		Virus strain
	 */
	public CodonPercents(String resourceName, VirusT virusInstance, Strain<VirusT> strain) {

		try (
			InputStream stream = this
				.getClass().getClassLoader()
				.getResourceAsStream(resourceName);
		) {
			String raw = IOUtils.toString(stream, StandardCharsets.UTF_8);
			List<Map<String, ?>> codonPcntData = gson.fromJson(
				raw, new TypeToken<List<Map<String, ?>>>(){}.getType());
			List<CodonPercent<VirusT>> codonPcnts = (
				codonPcntData.stream()
				.map(aaPcnt -> new CodonPercent<>(
					strain.getGene((String) aaPcnt.get("gene")),
					((Double) aaPcnt.get("position")).intValue(),
					(String) aaPcnt.get("codon"),
					(Double) aaPcnt.get("percent"),
					((Double) aaPcnt.get("count")).intValue(),
					((Double) aaPcnt.get("total")).intValue()
				))
				.collect(Collectors.toList())
			);
			this.codonPcnts = Collections.unmodifiableList(codonPcnts);
		} catch (IOException|NullPointerException e) {
			throw new ExceptionInInitializerError(
				String.format("Invalid resource name (%s)", resourceName)
			);
		}

		for (CodonPercent<VirusT> cdPcnt : codonPcnts) {
			GenePosition<VirusT> gp = cdPcnt.getGenePosition();
			codonPcntMap.putIfAbsent(gp, new LinkedHashMap<>());
			codonPcntMap.get(gp).put(cdPcnt.getCodon(), cdPcnt);
		}
	}
	
	public static <VirusT extends Virus<VirusT>> CodonPercents<VirusT> newEmptyInstance() {
		return new CodonPercents<>();
	}
	
	/**
	 * Initialize an empty CodonPercents
	 * 
	 */
	private CodonPercents() {
		codonPcnts = Collections.emptyList();
	}

	public List<CodonPercent<VirusT>> get() {
		// make a copy in case of any modification
		return new ArrayList<>(codonPcnts);
	}

	public List<CodonPercent<VirusT>> get(Gene<VirusT> gene) {
		return (codonPcnts
				.stream().filter(cdp -> cdp.getGene().equals(gene))
				.collect(Collectors.toList()));
	}

	public List<CodonPercent<VirusT>> get(Gene<VirusT> gene, int pos) {
		return new ArrayList<>(
			codonPcntMap.getOrDefault(new GenePosition<VirusT>(gene, pos), Collections.emptyMap())
			.values());
	}

	public CodonPercent<VirusT> get(Gene<VirusT> gene, int pos, String codon) {
		Map<String, CodonPercent<VirusT>> posCodons =
			codonPcntMap.getOrDefault(new GenePosition<VirusT>(gene, pos), Collections.emptyMap());
		if (posCodons.containsKey(codon)) {
			return posCodons.get(codon);
		}
		else if (posCodons.isEmpty()) {
			throw new IllegalArgumentException(
				String.format("Argument 'pos' is out of range: %d", pos));
		}
		else if (codon.matches("^(ins|del)$")) {
			int total = posCodons.values().iterator().next().getTotal();
			CodonPercent<VirusT> posCodon = new CodonPercent<>(gene, pos, codon, .0, 0, total);
			posCodons.put(codon, posCodon);
			return posCodon;
		}
		else if (codon.matches("^[ACGT]{3}$")) {
			int total = posCodons.values().iterator().next().getTotal();
			CodonPercent<VirusT> posCodon = new CodonPercent<>(gene, pos, codon, .0, 0, total);
			posCodons.put(codon, posCodon);
			return posCodon;
		}
		else {
			throw new IllegalArgumentException(
				String.format("Invalid argument codon \"%s\" at %s%d", codon, gene, pos));
		}
	}

	/**
	 * Returns the highest codon prevalence associated with each of
	 * the codon in a mixture.
	 *
	 * @param gene			Gene
	 * @param pos			Position
	 * @param codonMixture	Mixture codon list
	 *
	 * @return Double highest amino acid prevalence
	 */
	public Double getHighestCodonPercentValue(
		Gene<VirusT> gene, int pos, String... codonMixture
	) {
		Double pcntVal = 0.0;

		for (String cd : codonMixture) {
			double cdPcntVal = get(gene, pos, cd).getPercent();
			pcntVal = Math.max(pcntVal, cdPcntVal);
		}
		return pcntVal;
	}

}
