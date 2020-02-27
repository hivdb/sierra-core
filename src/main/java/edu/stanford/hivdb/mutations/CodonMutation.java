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

import java.util.Map;

import org.apache.commons.lang3.StringUtils;

import edu.stanford.hivdb.utilities.AAUtils;
import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

public class CodonMutation<VirusT extends Virus<VirusT>> extends AAMutation<VirusT> {

	private final String aas;
	private final String triplet;
	private final String insertedNAs;
	
	private static char[] calcAACharArray(String aas) {
		aas = AAUtils.normalizeAAs(aas);
		if (aas.contains("_")) {
			aas = "_";
		}
		return aas.toCharArray();
	}

	public static <VirusT extends Virus<VirusT>> CodonMutation<VirusT> fromNucAminoMutation(Gene<VirusT> gene, int aaStart, Map<?, ?> mut) {
		int pos = ((Double) mut.get("Position")).intValue() - aaStart + 1;

		String codon = "";
		String insertedCodon = "";
		boolean isInsertion = (Boolean) mut.get("IsInsertion");
		boolean isDeletion = (Boolean) mut.get("IsDeletion");

		StringBuilder aas = new StringBuilder();
		if (isDeletion) {
			aas.append('-');
		}
		else {
			codon = (String) mut.get("CodonText");
			codon = codon.replace(' ', '-');
			// The length of `CodonText` from NucAmino always equals to 3
			aas.append(CodonUtils.translateNATriplet(codon));
			if (isInsertion) {
				aas.append('_');
				insertedCodon = (String) mut.get("InsertedCodonsText");
				aas.append(CodonUtils.simpleTranslate(insertedCodon));
			}
		}
		return new CodonMutation<>(gene, pos, aas.toString(), codon, insertedCodon);
	}
	
	/**
	 *
	 * @param gene
	 * @param pos
	 * @param aas
	 * @param triplet
	 * @param insertedNAs
	 * @param maxDisplayAAs
	 */
	public CodonMutation(
		Gene<VirusT> gene, int position, String aas,
		String triplet, String insertedNAs,
		int maxDisplayAAs
	) {
		super(gene, position, calcAACharArray(aas), maxDisplayAAs);
		this.aas = AAUtils.normalizeAAs(aas);
		this.triplet = triplet.toUpperCase();
		this.insertedNAs = insertedNAs;
	}
	
	/**
	 *
	 * @param gene
	 * @param pos
	 * @param aas
	 * @param triplet
	 * @param insertedNAs
	 */
	public CodonMutation(
		Gene<VirusT> gene, int position, String aas,
		String triplet, String insertedNAs
	) {
		this(
			gene, position, aas, triplet, insertedNAs,
			AAMutation.DEFAULT_MAX_DISPLAY_AAS
		);
	}

	public CodonMutation(Gene<VirusT> gene, int position, String aas, String triplet) {
		this(gene, position, aas, triplet, "");
	}

	public CodonMutation(Gene<VirusT> gene, int position, String aas) {
		this(gene, position, aas, "", "");
	}

	public CodonMutation(Gene<VirusT> gene, int position, Character aa) {
		this(gene, position, "" + aa, "", "");
	}

	@Override
	@Deprecated
	public Mutation<VirusT> mergesWith(Mutation<VirusT> another) {
		if (isIndel() || another.isIndel()) {
			throw new UnsupportedOperationException(String.format(
				"Can not merge indel mutations (%s with %s)",
				toString(), another.toString()));
		}
		return super.mergesWith(another);
	}

	@Override
	public boolean isUnsequenced() {
		// "NNN", "NN-", "NNG" should be considered as unsequenced region
		return !isInsertion() &&
			StringUtils.countMatches(triplet.replace('-', 'N'), "N") > 1;
	}

	@Override
	public String getDisplayAAs() {
		String[] splited = aas.split("_", 2);
		if (splited[0].length() > getMaxDisplayAAs()) {
			splited[0] = "X";
		}
		return StringUtils.join(splited, '_');
	}

	@Override
	public String getAAs() { return aas; }

	@Override
	public String getTriplet() { return triplet; }

	@Override
	public String getInsertedNAs() { return insertedNAs; }

	@Override
	public boolean hasBDHVN() {
		// TODO: what if BDHVN doesn't affect the amimo acid?
		return triplet.matches(".*[BDHVN].*");
	}

	@Override
	public String getAAsWithRefFirst() {
		String[] aas = getDisplayAAs().split("_", 2);
		String ref = getReference();

		if (aas[0].contains(ref)) {
			aas[0] = ref + aas[0].replaceAll(ref, "");
		}
		return String.join("_", aas);
	}

	@Override
	public String getAAsWithoutReference() {
		String[] aas = getDisplayAAs().split("_", 2);
		String ref = getReference();
		aas[0] = aas[0].replace(ref, "");
		return String.join("_", aas);
	}

}
