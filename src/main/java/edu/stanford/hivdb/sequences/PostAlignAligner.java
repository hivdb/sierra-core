/*

	Copyright (C) 2017-2020 Stanford HIVDB team

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

package edu.stanford.hivdb.sequences;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.io.IOUtils;

import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executor;
import java.util.concurrent.Executors;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.CodonMutation;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.utilities.FastaUtils;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

/**
 * Receives unaligned sequence objects; uses the Go program
 * NucAmino to align the input sequences and returns a list of
 * aligned sequences.
 *
 * The aligned sequence is empty when the amount of POL
 * sites below MIN_NUM_OF_SITES (50). Certain gene (PR, RT, or
 * IN) is ignored when the the length of gene sequence
 * below MIN_NUM_OF_SITES_PER_GENE.
 *
 */
public class PostAlignAligner<VirusT extends Virus<VirusT>> implements Aligner<VirusT> {
	
	private final VirusT virusInstance;
	// private final AlignmentConfig<VirusT> alignmentConfig;
	
	private final String EXECUTABLE;
	private final Map<Gene<VirusT>, Double> MIN_MATCH_PCNT;
	private final Map<Gene<VirusT>, Double> SEQUENCE_SHRINKAGE_WINDOW;
	private final Map<Gene<VirusT>, Double> SEQUENCE_SHRINKAGE_CUTOFF_PCNT;
	private final Map<Gene<VirusT>, File> REF_SEQUENCE;
	private final Map<Gene<VirusT>, List<String>> POST_PROCESSORS;
	private final Executor executor = Executors.newFixedThreadPool(4);
	
	protected PostAlignAligner(VirusT virusIns) {
		this.virusInstance = virusIns;
		String executable = System.getenv("POSTALIGN_PROGRAM");
		if (executable == null || executable == "") {
			// use "postalign" as default program path
			executable = "postalign";
		}
		EXECUTABLE = executable;
		
		AlignmentConfig<VirusT> alignConfig = virusIns.getAlignmentConfig();
		MIN_MATCH_PCNT = Collections.unmodifiableMap(alignConfig.getConfigField("minMatchPcnt"));
		SEQUENCE_SHRINKAGE_WINDOW = Collections.unmodifiableMap(alignConfig.getConfigField("seqShrinkageWindow"));
		SEQUENCE_SHRINKAGE_CUTOFF_PCNT = Collections.unmodifiableMap(alignConfig.getConfigField("seqShrinkageCutoffPcnt"));

		Map<Gene<VirusT>, String> refSeqText = alignConfig.getConfigField("refSequence");
		REF_SEQUENCE = Collections.unmodifiableMap(
			refSeqText.entrySet().stream()
			.collect(Collectors.toMap(
				e -> e.getKey(),
				e -> makeReferenceFASTA(e.getKey(), e.getValue())
			))
		);
		
		POST_PROCESSORS = Collections.unmodifiableMap(alignConfig.getConfigField("postProcessors"));
	}

	private File makeReferenceFASTA(Gene<VirusT> gene, String refSeq) {
		File referenceFasta;
		try {
			referenceFasta = File.createTempFile("postalign-ref-" + gene.getName(), ".fas");
			BufferedWriter bw = new BufferedWriter(new FileWriter(referenceFasta));
			bw.write(">Ref_" + gene.getName() + "\n" + refSeq + "\n");
			bw.close();
			referenceFasta.setReadOnly();
			referenceFasta.deleteOnExit();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		return referenceFasta;
	}
	
	/**
	 * Receives a sequence and aligns it to each gene by NucAmino.
	 *
	 * @param sequence	Sequence waiting to be aligned.
	 * @return 			an AlignedSequence object
	 */
	@Override
	public AlignedSequence<VirusT> align(Sequence sequence) {
		List<Sequence> seqs = new ArrayList<>();
		seqs.add(sequence);
		List<AlignedSequence<VirusT>> result = parallelAlign(seqs);
		if (result.isEmpty()) {
			return null;
		}
		return result.get(0);
	}

	/**
	 * Receives set of sequences and aligns them to each gene in parallel by NucAmino.
	 *
	 * @param sequences		Sequence list waiting to be aligned
	 * @return 				list of AlignedSequence objects
	 */
	@Override
	public List<AlignedSequence<VirusT>> parallelAlign(Collection<Sequence> sequences) {
		return parallelAlign(sequences, false);
	}

	/**
	 * Uses PostAlign to align HIV sequences.
	 *  
	 * @param sequences
	 * @return
	 */
	protected Map<Strain<VirusT>, Map<Gene<VirusT>, List<Map<String, ?>>>> execute(Collection<Sequence> sequences) {
		Map<Gene<VirusT>, CompletableFuture<List<Map<String, ?>>>> futures = new TreeMap<>();
		
		for (Gene<VirusT> gene : REF_SEQUENCE.keySet()) {
			File refSeqFile = REF_SEQUENCE.get(gene);
			List<String> cmd = Lists.newArrayList(
				EXECUTABLE,
				"-i", "-",
				"-o", "-",
				"-f", "MINIMAP2",
				"-r", refSeqFile.getAbsolutePath()
			);
			cmd.addAll(POST_PROCESSORS.get(gene));
			cmd.add("save-json");
			CompletableFuture<List<Map<String, ?>>> future = CompletableFuture.supplyAsync(() -> {
				List<Map<String, ?>> jsonObjs = new ArrayList<>();
				Iterable<List<Sequence>> partialSets = Iterables.partition(sequences, 100);
				for (List<Sequence> partialSet : partialSets) {
					try {
						Process proc = Runtime.getRuntime().exec(cmd.toArray(String[]::new));
						OutputStream stdin = proc.getOutputStream();
						FastaUtils.writeStream(partialSet, stdin);
						String jsonString = IOUtils.toString(proc.getInputStream(), "UTF-8");
						String errString = IOUtils.toString(proc.getErrorStream(), "UTF-8");
						if (errString.length() > 0) {
							throw new RuntimeException(errString);
						}
						jsonObjs.addAll(Json.loads(jsonString, new TypeToken<List<Map<String, ?>>>(){}.getType()));
						proc.waitFor();
					} catch (InterruptedException e) {
						throw new RuntimeException(e);
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
				return jsonObjs;
			}, executor);
			futures.put(gene, future);
		}

		CompletableFuture.allOf(futures.values().toArray(new CompletableFuture<?>[0])).join();
		
		Map<Gene<VirusT>, List<Map<String, ?>>> results = new TreeMap<>();
		for (Gene<VirusT> gene : futures.keySet()) {
			CompletableFuture<List<Map<String, ?>>> future = futures.get(gene);
			try {
				results.put(gene, future.get());
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			} catch (ExecutionException e) {
				throw new RuntimeException(e);
			}
		}
		return (
			results
			.entrySet()
			.stream()
			.collect(Collectors.toMap(
				e -> e.getKey().getStrain(),
				e -> {
					Map<Gene<VirusT>, List<Map<String, ?>>> result = new TreeMap<>();
					result.put(e.getKey(), e.getValue());
					return result;
				},
				(a, b) -> {a.putAll(b); return a;},
				TreeMap::new
			))
		);
	}

	private Map<Sequence, AlignedSequence<VirusT>> selectBestAlignments(
		List<AlignedSequence<VirusT>> newAlignments,
		Map<Sequence, AlignedSequence<VirusT>> knownAlignments
	) {
		for (AlignedSequence<VirusT> alignedSeq : newAlignments) {
			Sequence inputSeq = alignedSeq.getInputSequence();
			if (knownAlignments.containsKey(inputSeq)) {
				if (alignedSeq.isEmpty()) {
					// no overwrite
					continue;
				}
				AlignedSequence<VirusT> knownAlignedSeq = knownAlignments.get(inputSeq);
				// if (knownAlignedSeq.getAvailableGenes().size() < alignedSeq.getAvailableGenes().size()) {
				// 	knownAlignments.put(inputSeq, alignedSeq);
				// }
				if (knownAlignedSeq.getNumMatchedNAs() < alignedSeq.getNumMatchedNAs()) {
					knownAlignments.put(inputSeq, alignedSeq);
				}
			}
			else {
				knownAlignments.put(inputSeq, alignedSeq);
			}
		}
		return knownAlignments;
	}

	private List<AlignedSequence<VirusT>> parallelAlign(Collection<Sequence> sequences, boolean reversingSequence) {
		Map<Sequence, Map<Strain<VirusT>, StringBuilder>> errors = new LinkedHashMap<>();
		Collection<Sequence> preparedSeqs = sequences;
		if (reversingSequence) {
			preparedSeqs = preparedSeqs.stream()
				.map(s -> s.reverseCompliment())
				.collect(Collectors.toList());
		}
		Map<Strain<VirusT>, Map<Gene<VirusT>, List<Map<String, ?>>>> allJsonObjs;
		
		allJsonObjs = execute(preparedSeqs);
		
		Map<Sequence, AlignedSequence<VirusT>> results = new LinkedHashMap<>();
		for (Strain<VirusT> strain : allJsonObjs.keySet()) {
			List<AlignedSequence<VirusT>> alignedSeqs = processCommandOutput(
				strain, sequences, allJsonObjs.get(strain),
				reversingSequence, errors
			);
			results = selectBestAlignments(alignedSeqs, results);
		}
		if (!reversingSequence && !errors.isEmpty()) {
			// a second run for reverse complement

			int numStrains = virusInstance.getStrains().size();
			List<Sequence> errorSeqs = errors
				.entrySet().stream()
				.filter(e -> e.getValue().size() == numStrains)
				.map(e -> e.getKey())
				.collect(Collectors.toList());
			if (!errorSeqs.isEmpty()) {
				List<AlignedSequence<VirusT>> reversedResults = parallelAlign(errorSeqs, true);
				results = selectBestAlignments(reversedResults, results);
			}
		}
		return Lists.newArrayList(results.values());
	}

	private AlignedGeneSeq<VirusT> geneSeqFromReport(
			Sequence sequence, Gene<VirusT> gene, Map<?, ?> report,
			boolean sequenceReversed) {
		int geneLength = gene.getAASize();
		int firstAA = Math.max(1, ((Double) report.get("FirstAA")).intValue());
		int lastAA = Math.max(geneLength, ((Double) report.get("LastAA")).intValue());
		int aaSize = Math.max(0, lastAA - firstAA + 1);
		final int minNumOfSites = gene.getNucaminoMinNumOfAA();
		if (aaSize < minNumOfSites) {
			throw new MisAlignedException(String.format(
				"Alignment of gene %s was discarded " +
				"since the length of alignment was too short (< %d).",
				gene, minNumOfSites
			), aaSize == 0);
		}

		List<?> polAlignedSites = (List<?>) report.get("AlignedSites");
		List<AlignedSite> alignedSites = polAlignedSites.stream()
			.map(m -> (Map<?, ?>) m)
			.map(m -> new AlignedSite(
				((Double) m.get("PosAA")).intValue(),
				((Double) m.get("PosNA")).intValue(),
				((Double) m.get("LengthNA")).intValue()
			))
			.collect(Collectors.toList());

		int firstNA = alignedSites.get(0).getPosNA();
		AlignedSite lastSite = alignedSites.get(alignedSites.size() - 1);
		int lastNA = lastSite.getPosNA() - 1 + lastSite.getLengthNA();

		List<?> polMutations = (List<?>) report.get("Mutations");
		List<Mutation<VirusT>> mutations = polMutations.stream()
			.map(m -> (Map<?, ?>) m)
			.map(m -> CodonMutation.fromNucAminoMutation(gene, 1, m))
			.collect(Collectors.toList());

		List<?> polFrameShifts = (List<?>) report.get("FrameShifts");
		List<FrameShift<VirusT>> frameShifts = polFrameShifts.stream()
			.map(fs -> (Map<?, ?>) fs)
			.map(fs -> FrameShift.fromNucAminoFrameShift(gene, 1, fs))
			.collect(Collectors.toList());

		int[] trimDels = trimGaps(sequence, firstAA, lastAA, mutations, frameShifts);
		int trimDelsLeft = trimDels[0];
		int trimDelsRight = trimDels[1];

		AlignedGeneSeq<VirusT> geneSeq = new AlignedGeneSeq<>(
			sequence, gene,
			firstAA + trimDelsLeft,
			lastAA - trimDelsRight,
			firstNA + trimDelsLeft * 3,
			lastNA - trimDelsRight * 3,
			alignedSites, mutations, frameShifts, 0, 0, sequenceReversed);
		if (geneSeq.getMatchPcnt() < MIN_MATCH_PCNT.get(gene)) {
			throw new MisAlignedException(String.format(
				"Alignment of gene %s was discarded " +
				"since the discordance rate was too high (%.1f%% > %d%%).",
				gene, 100 - geneSeq.getMatchPcnt(), 100 - MIN_MATCH_PCNT.get(gene)
			), false);
		}

		int[] trimUUs = trimLowQualities(
			gene, sequence,
			geneSeq.getFirstAA(),
			geneSeq.getLastAA(),
			geneSeq.getMutations(),
			geneSeq.getFrameShifts());
		int trimUUsLeft = trimUUs[0];
		int trimUUsRight = trimUUs[1];
		if (trimUUsLeft > 0 || trimUUsRight > 0) {
			geneSeq = new AlignedGeneSeq<>(
				sequence, gene,
				geneSeq.getFirstAA() + trimUUsLeft,
				geneSeq.getLastAA() - trimUUsRight,
				geneSeq.getFirstNA() + trimUUsLeft * 3,
				geneSeq.getLastNA() - trimUUsRight * 3,
				geneSeq.getAlignedSites(),
				geneSeq.getMutations(),
				geneSeq.getFrameShifts(), trimUUsLeft, trimUUsRight, sequenceReversed);
		}

		if (geneSeq.getSize() < minNumOfSites) {
			throw new MisAlignedException(String.format(
				"Alignment of gene %s was discarded " +
				"since the length of alignment (%d) was too short (< %d).",
				gene, aaSize, minNumOfSites
			), false);
		}

		return geneSeq;
	}

	/**
	 * Remove only deletion/NNNs from the beginning and the end of alignment
	 *
	 * @param sequence
	 * @param firstAA
	 * @param lastAA
	 * @param mutations
	 * @param frameShifts
	 * @return
	 */
	private int[] trimGaps(
		Sequence sequence, int firstAA, int lastAA,
		Collection<Mutation<VirusT>> mutations,
		Collection<FrameShift<VirusT>> frameShifts
	) {
		int trimLeft = 0;
		int trimRight = 0;
		int proteinSize = lastAA - firstAA + 1;
		List<Boolean> gapSites = new ArrayList<>(Collections.nCopies(proteinSize, false));
		for (Mutation<VirusT> mut : mutations) {
			int idx = mut.getPosition() - firstAA;
			if (mut.isDeletion() || mut.isUnsequenced()) {
				gapSites.set(idx, true);
			}
		}
		// remove initial deletions
		for (int idx=0; idx < proteinSize; idx ++) {
			if (!gapSites.get(idx)) {
				if (idx > trimLeft) {
					trimLeft = idx;
				}
				break;
			}
		}
		// remove trailing deletions
		for (int idx=proteinSize-1; idx > -1; idx --) {
			if (!gapSites.get(idx)) {
				if (proteinSize - idx - 1 > trimRight) {
					trimRight = proteinSize - idx - 1;
				}
				break;
			}
		}
		return new int[]{trimLeft, trimRight};
	}

	/**
	 *  Input sequence may contain non-POL NAs in the beginning and the end, e.g. AF442565, KF134931.
	 *
	 * Following code did two things:
	 *
	 *   1. Remove large (length > SEQUENCE_TRIM_SITES_CUTOFF) low quality pieces from gene sequence
	 *   2. Keep small (length <= SEQUENCE_TRIM_SITES_CUTOFF) low quality pieces from gene sequence
	 *
	 * A site is considered "low quality" if it:
	 *   - is unusual mutation;
	 *   - has "X" in aas; or
	 *   - has stop codon
	 */
	private int[] trimLowQualities(
		Gene<VirusT> gene,
		Sequence sequence,
		int firstAA,
		int lastAA,
		Collection<Mutation<VirusT>> mutations,
		Collection<FrameShift<VirusT>> frameShifts
	) {
		int badPcnt;
		int trimLeft = 0;
		int trimRight = 0;
		int problemSites = 0;
		int sinceLastBadQuality = 0;
		int proteinSize = lastAA - firstAA + 1;
		List<Integer> candidates = new ArrayList<>();
		List<Boolean> invalidSites = new ArrayList<>(Collections.nCopies(proteinSize, false));
		for (Mutation<VirusT> mut : mutations) {
			int idx = mut.getPosition() - firstAA;
			if (!mut.isUnsequenced() && (
					mut.isUnusual()
					|| mut.getDisplayAAs().equals("X") || mut.isApobecMutation() || mut.hasStop())) {
				invalidSites.set(idx, true);
			}
		}
		for (FrameShift<VirusT> fs : frameShifts) {
			int idx = fs.getPosition() - firstAA;
			invalidSites.set(idx,  true);
		}
		// forward scan for trimming left
		for (int idx=0; idx < proteinSize; idx ++) {
			if (sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW.get(gene)) {
				break;
			} else if (invalidSites.get(idx)) {
				problemSites ++;
				trimLeft = idx + 1;
				badPcnt = trimLeft > 0 ? problemSites * 100 / trimLeft : 0;
				if (badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT.get(gene)) {
					candidates.add(trimLeft);
				}
				sinceLastBadQuality = 0;
			} else {
				sinceLastBadQuality ++;
			}
		}
		trimLeft = candidates.size() > 0 ? candidates.get(candidates.size() - 1) : 0;
		candidates.clear();
		// backward scan for trimming right
		problemSites = 0;
		sinceLastBadQuality = 0;
		for (int idx=proteinSize-1; idx > -1; idx --) {
			if (sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW.get(gene)) {
				break;
			} else if (invalidSites.get(idx)) {
				problemSites ++;
				trimRight = proteinSize - idx;
				badPcnt = trimRight > 0 ? problemSites * 100 / trimRight : 0;
				if (badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT.get(gene)) {
					candidates.add(trimRight);
				}
				sinceLastBadQuality = 0;
			} else {
				sinceLastBadQuality ++;
			}
		}
		trimRight = candidates.size() > 0 ? candidates.get(candidates.size() - 1) : 0;
		return new int[]{trimLeft, trimRight};
	}

	/**
	 * Process the JSON output of NucAmino
	 * @param sequences - input unaligned sequences
	 * @param jsonString - string output of NucAmino
	 * @param errors - map of failed sequences and the errors
	 * @return list of AlignedSequence for all input sequences
	 */
	private List<AlignedSequence<VirusT>> processCommandOutput(
			Strain<VirusT> strain,
			Collection<Sequence> sequences,
			Map<Gene<VirusT>, List<Map<String, ?>>> jsonObjs,
			boolean sequenceReversed,
			Map<Sequence, Map<Strain<VirusT>,
			StringBuilder>> errors) {

		int numAlignments = 0;
		for (List<Map<String, ?>> alnResults : jsonObjs.values()) {
			numAlignments = Math.max(numAlignments, ((List<?>) alnResults).size());
		}
		List<AlignedSequence<VirusT>> alignedSequences = new ArrayList<>();
		Map<String, Sequence> sequenceMap = sequences.stream()
			.collect(Collectors.toMap(seq -> seq.getHeader(), seq -> seq));
		
		for (int idx = 0; idx < numAlignments; idx ++) {
			Sequence sequence = null;
			Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> alignedGeneSeqs = new TreeMap<>();
			Map<Gene<VirusT>, String> discardedGenes = new LinkedHashMap<>();

			for (Gene<VirusT> gene : jsonObjs.keySet()) {
				Map<String, ?> result = jsonObjs.get(gene).get(idx);
				// TODO: should we use hash key to prevent name conflict?
				String name = (String) result.get("Name");
				sequence = sequenceMap.get(name);
				Map<?, ?> report = (Map<?, ?>) result.get("Report");
				String error = (String) result.get("Error");
				if (!error.isEmpty()) {
					errors.putIfAbsent(sequence, new TreeMap<>());
					errors.get(sequence).putIfAbsent(strain, new StringBuilder());
					errors.get(sequence).get(strain).append(error);
				} else {
					try {
						alignedGeneSeqs.put(gene, geneSeqFromReport(sequence, gene, report, sequenceReversed));
					} catch (MisAlignedException e) {
						if (!e.isSuppressible()) {
							discardedGenes.put(gene, e.getMessage());
						}
					}
				}
			}
			
			if (sequence == null) {
				throw new RuntimeException("Nucamino returns malformed results.");
			}
			
			if (alignedGeneSeqs.isEmpty()) {
				errors.putIfAbsent(sequence, new TreeMap<>());
				errors.get(sequence).putIfAbsent(strain, new StringBuilder());
				errors.get(sequence).get(strain).append("No aligned results were found.");
			}
			alignedSequences.add(
				new AlignedSequence<>(
					strain, sequence, alignedGeneSeqs,
					discardedGenes, sequenceReversed)
			);
		}
		
		return alignedSequences;

		// AlignmentExtension extResult = new AlignmentExtension(
		// 	sequence, gene, firstAA, lastAA, firstNA, lastNA,
		// 	alignResult.getAlignedNAs(),
		// 	alignResult.getControlLine(),
		// 	alignResult.getAATripletLine());
	}

}
