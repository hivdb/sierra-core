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

import java.util.Collections;
import java.util.List;

import com.google.gson.reflect.TypeToken;
import com.neovisionaries.i18n.CountryCode;

import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Virus;

public class GenotypeReference<VirusT extends Virus<VirusT>> {

	private String genotypeName;
	private String country;
	private String authorYear;
	private Integer year;
	private String accession;
	private Integer firstNA;
	private Integer lastNA;
	private String sequence;
	private transient VirusT virusInstance;

	public static <VirusT extends Virus<VirusT>> List<GenotypeReference<VirusT>> loadJson(String raw, VirusT virusIns) {
		//InputStream json = (
		//	HIVGenotypeReference.class.getClassLoader()
		//	.getResourceAsStream("HIVGenotypeReferences.json"));
		List<GenotypeReference<VirusT>> references = Json.loads(
				raw,
			    new TypeToken<List<GenotypeReference<VirusT>>>(){});
		for (GenotypeReference<VirusT> ref : references) {
			ref.virusInstance = virusIns;
		}
		return Collections.unmodifiableList(references);
	}

	/** get BoundGenotype object of given sequence
	 *
	 * @param sequence a string of DNA sequence
	 * @param firstNA starting position of the given sequence
	 * @param lastNA ending position of the given sequence
	 * @param discordanceList discordance from the comparison
	 * @return BoundGenotype object contained the comparison result
	 */
	public BoundGenotype<VirusT> getBoundGenotype(
			String sequence, int firstNA, int lastNA,
			List<Integer> discordanceList) {
		firstNA = Math.max(firstNA, getFirstNA());
		lastNA = Math.min(lastNA, getLastNA());
		return new BoundGenotype<VirusT>(
			this, sequence, firstNA, lastNA, discordanceList, virusInstance);
	}

	/** a getter to get current reference's start position in HXB2
	 *
	 * @return integer
	 */
	public Integer getFirstNA() {
		return firstNA;
	}

	/** a getter to get current reference's end position in HXB2
	 *
	 * @return integer
	 */
	public Integer getLastNA() {
		return lastNA;
	}

	/** a getter to get current reference's genotype
	 *
	 * @return a Genotype object that indicates the corresponding genotype
	 */
	public Genotype<VirusT> getGenotype() {
		return virusInstance.getGenotype(genotypeName);
	}

	public String getCountry() {
		CountryCode countryCode = CountryCode.getByCode(country);
		if (countryCode == null) {
			return country;
		}
		else {
			return countryCode.getName();
		}
	}

	public String getAuthorYear() {
		return authorYear;
	}

	public Integer getYear() {
		return year != null ? year : null;
	}

	/** a getter to get current reference's genbank accession id
	 *
	 * @return string
	 */
	public String getAccession() {
		return accession;
	}

	/** a getter to get current reference's sequence
	 *
	 * @return string
	 */
	public String getSequence() {
		return sequence;
	}

	@Override
	public String toString() {
		return getAccession() + " (" + genotypeName + ")";
	}
}
