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
public class AminoAcidPercents<VirusT extends Virus<VirusT>> {

	final static protected Gson gson = new Gson();

	final private List<AminoAcidPercent<VirusT>> aminoAcidPcnts;
	final private Map<GenePosition<VirusT>, Map<Character, AminoAcidPercent<VirusT>>> aminoAcidPcntMap = new HashMap<>();


	/**
	 * AminoAcidPercents initializer
	 *
	 * @param resourceName
	 * @param virusInstance
	 * @param strain
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
		return (
			aminoAcidPcntMap
			.getOrDefault(genePos, Collections.emptyMap())
			.get(aa));
	}

	/**
	 * Returns the highest amino acid prevalence associated with each of
	 * the AA in a mixture.
	 *
	 * @param gene
	 * @param pos
	 * @param cons consensus at the position
	 * @param mixture
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
	 * Returns true if the given mutation contains any unusual AA
	 *
	 * @param gene
	 * @param pos
	 * @param aas
	 * @return true if contains unusual AA
	 */
	public Boolean containsUnusualAA(Gene<VirusT> gene, int pos, String aas) {
		GenePosition<VirusT> gpos = new GenePosition<VirusT>(gene, pos);
		for (char aa : aas.toCharArray()) {
			AminoAcidPercent<VirusT> aaPcnt = get(gpos, aa);
			if (aaPcnt != null && aaPcnt.isUnusual()) {
				return true;
			}
		}
		return false;
	}

}
