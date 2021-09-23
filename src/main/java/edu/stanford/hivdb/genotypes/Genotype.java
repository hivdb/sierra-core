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
package edu.stanford.hivdb.genotypes;

import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Virus;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.math.BigDecimal;
import java.math.RoundingMode;


public class Genotype<VirusT extends Virus<VirusT>> {

	public static class RegionalGenotype<VirusT extends Virus<VirusT>> {
		private Double proportion;
		private Genotype<VirusT> genotype;

		public RegionalGenotype(
			Genotype<VirusT> genotype, double proportion
		) {
			this.proportion = proportion;
			this.genotype = genotype;
		}

		public Genotype<VirusT> getGenotype() {
			return genotype;
		}

		public Double getProportion() {
			return proportion;
		}

		@Override
		public String toString() {
			int scale = 2;
			if (proportion + 1e-8 > 1.0) {
				scale = 0;
			}
			else if (proportion > 0.1) {
				scale = 1;
			}
			BigDecimal pcnt = new BigDecimal(proportion * 100);
			pcnt = pcnt.setScale(scale, RoundingMode.HALF_UP);
			return genotype + " (" + pcnt.toPlainString() + "%)";
		}

	}


	private static final Double MIN_PRIMARY_REGIONAL_GENOTYPE_PROPORTION = 0.9; // 90%
	// private static Map<String, HIVGenotype> genotypes = new LinkedHashMap<>();
	private String name;
	private Boolean isSimpleCRF;
	private String displayName;
	private String parentGenotypes;
	private GenotypeClassificationLevel classificationLevel;
	private Double distanceUpperLimit;
	private List<CRFRegion> regions;
	
	private transient VirusT virusInstance;

	private static class CRFRegion {
		String genotypeName;
		int start;
		int end;
	}

	public static <VirusT extends Virus<VirusT>> Map<String, Genotype<VirusT>> loadJson(String raw, VirusT virusIns) {
		// InputStream json = HIVGenotype.class.getClassLoader().getResourceAsStream("HIVGenotypes.json");
		Map<String, Genotype<VirusT>> genotypes = new LinkedHashMap<>();
		genotypes = Json.loads(
			raw, new TypeToken<Map<String, Genotype<VirusT>>>(){});
		for (Genotype<VirusT> genotype : genotypes.values()) {
			genotype.virusInstance = virusIns;
		}
		return Collections.unmodifiableMap(genotypes);
	}

	private Genotype() {}

	/*public static HIVGenotype getInstance(String name) {
		HIVGenotype genotype = genotypes.get(name);
		if (genotype == null) {
			throw new NullPointerException(String.format("Cannot find subtype for %s.", name));
		}
		return genotype;
	}*/

	/** get genotype based given boundaries of partial sequence
	 *
	 * A recombination (CRFs) of other genotypes can be divided
	 * into regions by certain breakpoints. Therefore, a partial
	 * sequence can be conditionally considered as one of the
	 * recombinant genotypes, if all of its nucleotides were
	 * presented in such a single region.
	 *
	 * @param firstNA the start position of the partial sequence
	 * @param lastNA the end position
	 * @return Genotype object; returns itself when no better
	 * matches were found
	 */
	public RegionalGenotype<VirusT> getPrimaryRegionalGenotype(int firstNA, int lastNA) {
		List<RegionalGenotype<VirusT>> results = getRegionalGenotypes(firstNA, lastNA);
		RegionalGenotype<VirusT> primary = results
			.stream()
			.sorted((r1, r2) -> r2.proportion.compareTo(r1.proportion))
			.findFirst().get();
		if (primary.proportion >= getMinPrimaryRegionalGenotypeProportion()) {
			return primary;
		}
		return new RegionalGenotype<>(this, 1.0);
	}

	public Double getMinPrimaryRegionalGenotypeProportion() {
		return MIN_PRIMARY_REGIONAL_GENOTYPE_PROPORTION;
	}

	/** get genotype index name
	 *
	 * @return string
	 */
	public String getIndexName() {
		return name;
	}

	/** get genotype display name
	 *
	 * @return string
	 */
	public String getDisplayName() {
		return displayName;
	}

	/** get parent genotypes
	 *
	 * Sub-subtype and Some CRF belong to one or more parent
	 * genotypes. Sometime when the distance between the
	 * submitted and reference sequence was too high to report
	 * the original genotype itself, the parent genotypes were
	 * reported instead. This method returns these parent
	 * genotypes if them exist. Otherwise, a null pointer should
	 * be returned.
	 *
	 * @return List of Genotype object
	 */
	public List<Genotype<VirusT>> getParentGenotypes() {
		if (parentGenotypes != null) {
			return Arrays
				.stream(StringUtils.split(parentGenotypes, '|'))
				.map(n -> virusInstance.getGenotype(n))
				.collect(Collectors.toList());
		}
		else {
			return null;
		}
	}

	/** check distance
	 *
	 * Accepts a double value of the distance between submitted
	 * sequence and the reference. Returns false when distance
	 * is too high for this genotype.
	 *
	 * @param distance	Distance between submitted sequence and reference sequence.
	 * @return			Is distance less than upper limit
	 */
	public Boolean checkDistance(double distance) {
		// System.out.printf("name: %s; distance: %f; distanceUL: %f\n", name, distance, distanceUpperLimit);
		return distance < distanceUpperLimit;
	}

	public List<RegionalGenotype<VirusT>> getRegionalGenotypes(int firstNA, int lastNA) {
		Map<Genotype<VirusT>, Double> mapResults = new LinkedHashMap<>();
		double length = lastNA - firstNA;
		if (isSimpleCRF) {
			for (CRFRegion region : regions) {
				if (lastNA >= region.start && firstNA <= region.end) {
					// intersected
					int start = firstNA > region.start ? firstNA : region.start;
					int end = lastNA < region.end ? lastNA : region.end;
					Genotype<VirusT> genotype = virusInstance.getGenotype(region.genotypeName);
					double proportion = mapResults.getOrDefault(genotype, 0.0);
					proportion += (end - start) / length;
					mapResults.put(genotype, proportion);
				}
			}
		}
		else {
			mapResults.put(this, 1.0);
		}
		return mapResults
			.entrySet().stream()
			.map(e -> new RegionalGenotype<>(e.getKey(), e.getValue()))
			.collect(Collectors.toList());
	}

	/** get if the current genotype is treated as a recombination form
	 *
	 * @return boolean
	 */
	public Boolean hasParentGenotypes() {
		return parentGenotypes != null;
	}

	@Override
	public String toString() {
		return getDisplayName();
	}

	public GenotypeClassificationLevel getClassificationLevel() {
		return classificationLevel;
	}

}
