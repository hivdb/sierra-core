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

import java.util.LinkedHashMap;
import java.util.Map;

import edu.stanford.hivdb.utilities.CodonUtils;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;

public class CodonReads<VirusT extends Virus<VirusT>> implements WithGene<VirusT> {
	private final String codon;
	private final long reads;
	private final long totalReads;
	private final Gene<VirusT> gene;
	private final int position;
	private final double proportion;
	private transient Mutation<VirusT> mutation;
	private transient Character aminoAcid;
	private transient Double aaPcnt;
	private transient Double codonPcnt;

	
	public static String normalizeCodon(String codon) {
		// Tolerant spaces, commas, colons, semicolons and dashes
		return codon.replaceAll("[ ,:;-]", "");
	}

	public CodonReads(
		final Gene<VirusT> gene, final int position,
		final String codon, final long reads, final long totalReads
	) {
		this.gene = gene;
		this.position = position;
		this.codon = normalizeCodon(codon);
		this.reads = reads;
		this.totalReads = totalReads;
		this.proportion = (double) reads / totalReads;
	}
	
	public String getCodon() { return codon; }
	public Long getReads() { return reads; }
	public Long getTotalReads() { return totalReads; }

	@Override
	public Strain<VirusT> getStrain() { return gene.getStrain(); }
	
	@Override
	public Gene<VirusT> getGene() { return gene; }
	
	@Override
	public String getAbstractGene() { return gene.getAbstractGene(); }
	
	public Character getRefAminoAcid() {
		return gene.getRefChar(position);
	}
	
	public Character getAminoAcid() {
		if (aminoAcid == null) {
			if (!codon.matches("^[ACGT]*$")) {
				// do not allow ambiguous codes
				aminoAcid = 'X';
			}
			if (codon.length() > 5) {
				aminoAcid = '_';  // insertion
			}
			else if (codon.length() < 3) {
				aminoAcid = '-';  // deletion
			}
			else {
				String aminoAcids = CodonUtils.translateNATriplet(codon.substring(0, 3));
				if (aminoAcids.length() > 1) {
					// Ambiguous codon should not happen in NGS codons
					aminoAcid = 'X';
				}
				else {
					aminoAcid = aminoAcids.charAt(0);
				}
			}
		}
		return aminoAcid;
	}

	private Mutation<VirusT> getMutation() {
		if (mutation == null) {
			mutation = new AAMutation<>(gene, position, getAminoAcid());
		}
		return mutation;
	}
	
	public Double getProportion() {
		return this.proportion;
	}
	
	public Double getCodonPercent() {
		if (this.codonPcnt == null) {
			Double codonPcntNumber;
			CodonPercent<VirusT> codonPcnt;
			String cleanedCodon = codon;
			if (codon.length() > 5) {
				cleanedCodon = "ins";
			}
			else if (codon.length() < 3) {
				cleanedCodon = "del";
			}
			try {
				codonPcnt = (
					gene.getVirusInstance()
					.getCodonPercents(gene.getStrain(), "all", "all")
					.get(gene, position, cleanedCodon)
				);
				if (codonPcnt == null) {
					codonPcntNumber = .0;
				}
				else {
					codonPcntNumber = codonPcnt.getPercent();
				}
			}
			catch (IllegalArgumentException e) {
				codonPcntNumber = .0;
			}
			this.codonPcnt = codonPcntNumber;
		}
		return this.codonPcnt;
	}
	
	public Double getAAPercent() {
		if (aaPcnt == null) {
			aaPcnt = (
				gene.getVirusInstance()
				.getAminoAcidPercents(gene.getStrain(), "all", "all")
				.get(gene, position, getAminoAcid())
				.getPercent()
			);
				
		}
		return aaPcnt;
	}
	
	public boolean isReference() {
		return getAminoAcid() == gene.getRefChar(position);
	}
	
	public boolean isApobecMutation() {
		return isReference() ? false : getMutation().isApobecMutation();
	}
	
	public boolean isApobecDRM() {
		return isReference() ? false : getMutation().isApobecDRM();
	}
	
	public boolean hasStop() {
		return isReference() ? false : getMutation().hasStop();
	}
	
	public boolean isUnusual() {
		return isReference() ? false : getMutation().isUnusual();
	}
	
	public boolean isUnusualByCodon() {
		return getCodonPercent() < 0.0001;
	}
	
	public boolean isDRM() {
		return isReference() ? false : getMutation().isDRM();
	}
	
	/**
	 * Get Extended map
	 * 
	 * @return Map&lt;String, Object&gt;
	 */
	public Map<String, Object> extMap() {
		Map<String, Object> result = new LinkedHashMap<>();
		result.put("codon", getCodon());
		result.put("reads", getReads());
		result.put("refAminoAcid", getRefAminoAcid());
		result.put("aminoAcid", getAminoAcid());
		result.put("proportion", getProportion());
		result.put("aaPercent", getAAPercent());
		result.put("codonPercent", getCodonPercent());
		result.put("isApobecMutation", isApobecMutation());
		result.put("isApobecDRM", isApobecDRM());
		result.put("isUnusual", isUnusual());
		result.put("isDRM", isDRM());
		return result;
	}

}