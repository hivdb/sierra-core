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

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.CharMatcher;

import edu.stanford.hivdb.genotypes.Genotype.RegionalGenotype;
import edu.stanford.hivdb.mutations.StrainModifier;
import edu.stanford.hivdb.viruses.Virus;

public class BoundGenotype<VirusT extends Virus<VirusT>> {

	private static final Double MAX_FALLBACK_TO_SECONDARY_DISTANCE_DIFF = 0.01; // 1%
	private String sequence;
	private GenotypeReference<VirusT> reference;
	private int firstNA;
	private int lastNA;
	private List<Integer> discordanceList;
	private double distance;
	private final VirusT virusInstance;

	protected BoundGenotype(
		GenotypeReference<VirusT> reference,
		String sequence, int firstNA, int lastNA,
		List<Integer> discordanceList, VirusT virusIns
	) {
		this.reference = reference;
		this.sequence = sequence;
		this.firstNA = firstNA;
		this.lastNA = lastNA;
		this.discordanceList = discordanceList;
		int numWildcards = CharMatcher.is(StrainModifier.WILDCARD).countIn(sequence);
		this.distance = (
				(double) discordanceList.size() /
				(lastNA - firstNA + 1 - numWildcards));
		this.virusInstance = virusIns;
	}

	/** get the sequence tested
	 *
	 * @return String
	 */
	public String getSequence() {
		return sequence;
	};

	/** get first compared NA position
	 *
	 * @return the first compared NA position
	 */
	public int getFirstNA() {
		return firstNA;
	}

	/** get last compared NA position
	 *
	 * @return the last compared NA position
	 */
	public int getLastNA() {
		return lastNA;
	}

	/** get original distance
	 *
	 * @return distance (double, &lt; 1.0)
	 */
	public Double getDistance() {
		return distance;
	}

	/** get reference
	 *
	 * @return GenotypeReference object
	 */
	public GenotypeReference<VirusT> getReference() {
		return reference;
	}

	public String getReferenceAccession() {
		return reference.getAccession();
	}
	
	public String getReferenceCountry() {
		return reference.getCountry();
	}
	
	public Integer getReferenceYear() {
		return reference.getYear();
	}

	/** get genotype
	 *
	 * @return Genotype object
	 */
	public Genotype<VirusT> getGenotype() {
		return reference.getGenotype();
	}
	
	public Genotype<VirusT> getSubtype() {
		return reference.getGenotype();
	}

	/** get discordance positions
	 *
	 * @return an array of position numbers
	 */
	public List<Integer> getDiscordanceList() {
		return discordanceList;
	}

	/** get distance in percent form
	 *
	 * @return String
	 */
	public String getDistancePcnt() {
		int scale = 2; // 0.00% ~ 9.99%
		if (distance + 1e-8 > 1.0) {
			// 100%
			scale = 0;
		}
		else if (distance > 0.1) {
			// > 10.0%
			scale = 1;
		}
		BigDecimal pcnt = new BigDecimal(distance * 100);
		pcnt = pcnt.setScale(scale, RoundingMode.HALF_UP);
		return pcnt.toPlainString() + "%";
	}

	/** get parent / regional genotype(s)
	 *
	 * @return List of Genotype object
	 */
	public List<Genotype<VirusT>> getDisplayGenotypes() {
		List<Genotype<VirusT>> displayGenotypes = new ArrayList<>();
		if (shouldDisplayUnknown()) {
			displayGenotypes.add(virusInstance.getGenotypeUnknown());
		}
		else {
			Genotype<VirusT> origGenotype = getGenotype();
			Genotype<VirusT> regionalGenotype = getPrimaryRegionalGenotype().getGenotype();
			if (!checkDistance() && regionalGenotype == origGenotype &&
					origGenotype.hasParentGenotypes()) {
				// distance is too far and no regional Genotype found
				return origGenotype.getParentGenotypes();
			}
			displayGenotypes.add(regionalGenotype);
		}
		return displayGenotypes;
	}
	
	public List<Genotype<VirusT>> getDisplaySubtypes() {
		return getDisplayGenotypes();
	}

	/** get display string
	 *
	 * The display result normally contains the name of
	 * current genotype and the distance percent between
	 * the following parentheses. For example: "B (1.4%)".
	 *
	 * There are two exceptions:
	 *
	 * 1. The genotype is a recombination of other types.
	 *    In this case, the given sequence can be
	 *    conditionally considered as one of the
	 *    recombinant genotypes. For example, the genotype
	 *    is CRF51_01B, the display string can be "B (2.5%)"
	 *    if the given sequence only contained PR and RT.
	 *
	 * 2. The genotype is a sub-subtype of another subtype.
	 *    In this case the subtype will be showed first.
	 *    The sub-subtype will be between the following
	 *    parentheses. For example, (previous) genotype
	 *    A-FSU will be always displayed like
	 *    "A (A_FSU) (3.1%)".
	 *
	 * @return a human-friendly string
	 */
	public String getDisplay() {
		StringBuffer buf = new StringBuffer();
		buf.append(getDisplayWithoutDistance());
		if (!buf.toString().equals("Unknown")) {
			buf.append(" (");
			buf.append(getDistancePcnt());
			buf.append(")");
		}
		return buf.toString();
	}

	public String getDisplayWithoutDistance() {
		StringBuffer buf = new StringBuffer();
		for (Genotype<VirusT> genotype : getDisplayGenotypes()) {
			buf.append(genotype.getDisplayName());
			buf.append(" + ");
		}
		buf.setLength(buf.length() - 3);
		return buf.toString();
	}

	/** get genotype based given boundaries of partial sequence
	 *
	 * A recombination (CRFs) of other genotypes can be divided
	 * into regions by certain breakpoints. Therefore, a partial
	 * sequence can be conditionally considered as one of the
	 * recombinant genotypes, if all of its nucleotides were
	 * presented in such a single region.
	 *
	 * @return Genotype object; returns itself when no better
	 * matches were found
	 */
	public RegionalGenotype<VirusT> getPrimaryRegionalGenotype() {
		return getGenotype()
			.getPrimaryRegionalGenotype(getFirstNA(), getLastNA());
	}

	public List<RegionalGenotype<VirusT>> getRegionalGenotypes() {
		return getGenotype()
			.getRegionalGenotypes(getFirstNA(), getLastNA());
	}

	/** check distance
	 *
	 * Check the distance between current testing sequence
	 * and the reference. Returns false when distance is too
	 * large for this genotype.
	 *
	 * @return Boolean
	 */
	public Boolean checkDistance() {
		return getGenotype().checkDistance(getDistance());
	}

	/** check if "Unknown" should be reported
	 *
	 * A universal distance cut-off (11%) is applied to all
	 * genotypes. The 11% cut-off percent was estimated by
	 * analyzing the overall distances of ~11,000 sequences.
	 *
	 * @return Boolean
	 */
	public Boolean shouldDisplayUnknown() {
		return getDistance() > 0.11;
	}

	public List<Genotype<VirusT>> getParentGenotypes() {
		return getGenotype().getParentGenotypes();
	}

	public Boolean shouldFallbackTo(BoundGenotype<VirusT> fallback) {
		if (checkDistance() ||
				fallback.getDistance() - getDistance() >
				MAX_FALLBACK_TO_SECONDARY_DISTANCE_DIFF) {

			return false;
		}
		return true;
	}

	@Override
	public String toString() {
		return getDisplay();
	}

}
