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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Pair;

import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executor;
import java.util.concurrent.Executors;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

import com.amazonaws.services.lambda.AWSLambda;
import com.amazonaws.services.lambda.AWSLambdaClientBuilder;
import com.amazonaws.services.lambda.model.InvokeRequest;
import com.amazonaws.services.lambda.model.InvokeResult;
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
public class NucAminoAligner<VirusT extends Virus<VirusT>> {
	private final int MIN_MATCH_PCNT = 60;
	private final int SEQUENCE_SHRINKAGE_WINDOW = 15;
	private final int SEQUENCE_SHRINKAGE_CUTOFF_PCNT = 30;
	private final Executor executor = Executors.newFixedThreadPool(20);
	private final Map<Strain<VirusT>, String[]> NUCAMINO_LOCAL_COMMANDS;
	
	private final VirusT virusInstance;

	private final static Map<String, NucAminoAligner<?>> singletons = new HashMap<>(); 
	
	private static class MisAlignedException extends IllegalArgumentException {
		/**
		 *
		 */
		private static final long serialVersionUID = 46495128315347L;
		private final boolean suppressible;

		public MisAlignedException(String message, boolean suppressible) {
			super(message);
			this.suppressible = suppressible;
		}

		public boolean isSuppressible() { return suppressible; }
	}

	private static String getJoinedNucaminoGenes(Strain<?> strain) {
		return strain.getNucaminoGeneMap()
			.keySet().stream()
			.collect(Collectors.joining(","));
	}
	
	@SuppressWarnings("unchecked")
	public static <VirusT extends Virus<VirusT>> NucAminoAligner<VirusT> getInstance(VirusT virusIns) {
		String className = virusIns.getClass().getName();
		if (!singletons.containsKey(className)) {
			singletons.put(className, new NucAminoAligner<>(virusIns));
		}
		return (NucAminoAligner<VirusT>) singletons.get(className);
	}

	private NucAminoAligner(VirusT virusIns) {
		this.virusInstance = virusIns;
		String executable = System.getenv("NUCAMINO_PROGRAM");
		if (executable == null) {
			// use "nucamino" as default program path
			executable = "nucamino";
		}
		Map<Strain<VirusT>, String[]> nucaminoCommands = new HashMap<>();
		for (Strain<VirusT> strain : virusIns.getStrains()) {
			File profileFile = strain.makeNucaminoProfileFile();
			nucaminoCommands.put(
				strain,
				new String[] {
					/* Command */
					executable,	 	// path to nucamino binary
					"align-with",	// sub-command: use custom alignment profile
					profileFile.getAbsolutePath(),	// specify profile path
					getJoinedNucaminoGenes(strain),	// specify gene to align against
					/* Flags */
					"-q", 			// quiet mode
					"-f", "json", 	// return output format as json
				}
			);
		}
		
		NUCAMINO_LOCAL_COMMANDS = Collections.unmodifiableMap(nucaminoCommands);
	
	}
	
	/**
	 * Receives a sequence and aligns it to each HIV gene by NucAmino.
	 *
	 * @param sequence	Sequence waiting to be aligned.
	 * @return 			an AlignedSequence object
	 */
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
	 * Receives set of sequences and aligns them to each HIV gene in parallel by NucAmino.
	 *
	 * @param sequences		Sequence list waiting to be aligned
	 * @return 				list of AlignedSequence objects
	 */
	public List<AlignedSequence<VirusT>> parallelAlign(Collection<Sequence> sequences) {
		return parallelAlign(sequences, false);
	}

	/**
	 * Uses locally installed NucAmino to align HIV sequences.
	 *  
	 * @param sequences
	 * @return
	 */
	private Map<Strain<VirusT>, List<String>> localNucamino(Collection<Sequence> sequences) {
		Map<Strain<VirusT>, CompletableFuture<List<String>>> futures = new TreeMap<>();
		
		for (Strain<VirusT> strain : virusInstance.getStrains()) {
			String[] cmd = NUCAMINO_LOCAL_COMMANDS.get(strain);
			CompletableFuture<List<String>> future = CompletableFuture.supplyAsync(() -> {
				List<String> jsonStrings = new ArrayList<>();
				Iterable<List<Sequence>> partialSets = Iterables.partition(sequences, 10);
				for (List<Sequence> partialSet : partialSets) {
					try {
						Process proc = Runtime.getRuntime().exec(cmd);
						OutputStream stdin = proc.getOutputStream();
						BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
						// System.out.println(partialSet);
						FastaUtils.writeStream(partialSet, stdin);
						jsonStrings.add(stdout.lines().collect(Collectors.joining()));
						stdout.close();
						proc.waitFor();
					} catch (InterruptedException e) {
						throw new RuntimeException(e);
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
				return jsonStrings;
			}, executor);
			futures.put(strain, future);
		}

		CompletableFuture.allOf(futures.values().toArray(new CompletableFuture<?>[0])).join();
		
		Map<Strain<VirusT>, List<String>> results = new TreeMap<>();
		for (Strain<VirusT> strain : virusInstance.getStrains()) {
			CompletableFuture<List<String>> future = futures.get(strain);
			try {
				results.put(strain, future.get());
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			} catch (ExecutionException e) {
				throw new RuntimeException(e);
			}
		}
		return results;
	}

	private Map<Strain<VirusT>, List<String>> awsNucamino(Collection<Sequence> sequences, String awsFuncAndQual) {
		List<CompletableFuture<Pair<Strain<VirusT>, String>>> futures = new ArrayList<>();
		String[] funcAndQual = awsFuncAndQual.split(":");

		Iterable<List<Sequence>> partialSets = Iterables.partition(sequences, 5);
		for (List<Sequence> partialSet : partialSets) {
			for (Strain<VirusT> strain : virusInstance.getStrains()) {
				Map<String, String> payload = new HashMap<>();
				payload.put("profile_yaml", strain.makeNucaminoProfileString());
				payload.put("genes", getJoinedNucaminoGenes(strain));
				payload.put("fasta", FastaUtils.writeString(partialSet));
				String payloadText = Json.dumps(payload);
				AWSLambda client = AWSLambdaClientBuilder.standard().build();
				InvokeRequest request = new InvokeRequest()
					.withFunctionName(funcAndQual[0])
					.withPayload(payloadText)
					.withQualifier(funcAndQual[1]);
				CompletableFuture<Pair<Strain<VirusT>, String>> future = CompletableFuture.supplyAsync(() -> {
					InvokeResult response = client.invoke(request);
					ByteBuffer respPayload = response.getPayload();
					return Pair.of(
						strain,
						new String(respPayload.array(), Charset.forName("UTF-8"))
					);
				}, executor);
				futures.add(future);
			}
		}

		CompletableFuture.allOf(futures.toArray(new CompletableFuture<?>[0])).join();
		Map<Strain<VirusT>, List<String>> results = new TreeMap<>();
		for (CompletableFuture<Pair<Strain<VirusT>, String>> future : futures) {
			Pair<Strain<VirusT>, String> result;
			try {
				result = future.get();
				Strain<VirusT> strain = result.getLeft();
				String jsonString = result.getRight();
				if (!results.containsKey(strain)) {
					results.put(strain, new ArrayList<>());
				}
				results.get(strain).add(jsonString);
			} catch (InterruptedException | ExecutionException e) {
				throw new RuntimeException(e);
			}
		}
		return results;
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
		Map<Strain<VirusT>, List<String>> jsonStrings;
		
		String awsFunc = System.getenv("NUCAMINO_AWS_LAMBDA");
		if (awsFunc == null || awsFunc.equals("")) {
			jsonStrings = localNucamino(preparedSeqs);
		} else {
			jsonStrings = awsNucamino(preparedSeqs, awsFunc);
		}
		
		Map<Sequence, AlignedSequence<VirusT>> results = new LinkedHashMap<>();
		for (Strain<VirusT> strain : jsonStrings.keySet()) {
			for (String jsonString : jsonStrings.get(strain)) {
				List<AlignedSequence<VirusT>> alignedSeqs = processCommandOutput(
					strain, sequences, jsonString,
					reversingSequence, errors
				);
				results = selectBestAlignments(alignedSeqs, results);
			}
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
		if (geneSeq.getMatchPcnt() < MIN_MATCH_PCNT) {
			throw new MisAlignedException(String.format(
				"Alignment of gene %s was discarded " +
				"since the discordance rate was too high (%.1f%% > %d%%).",
				gene, 100 - geneSeq.getMatchPcnt(), 100 - MIN_MATCH_PCNT
			), false);
		}

		int[] trimUUs = trimLowQualities(
			sequence,
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
		Sequence sequence, int firstAA, int lastAA,
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
			if (sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW) {
				break;
			} else if (invalidSites.get(idx)) {
				problemSites ++;
				trimLeft = idx + 1;
				badPcnt = trimLeft > 0 ? problemSites * 100 / trimLeft : 0;
				if (badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT) {
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
			if (sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW) {
				break;
			} else if (invalidSites.get(idx)) {
				problemSites ++;
				trimRight = proteinSize - idx;
				badPcnt = trimRight > 0 ? problemSites * 100 / trimRight : 0;
				if (badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT) {
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
			Strain<VirusT> strain, Collection<Sequence> sequences, String jsonString,
			boolean sequenceReversed, Map<Sequence, Map<Strain<VirusT>, StringBuilder>> errors) {

		Map<?, ?> jsonObj = Json.loads(
			jsonString, new TypeToken<Map<?, ?>>(){}.getType());
		int numAlignments = 0;
		for (Object alnResults : jsonObj.values()) {
			numAlignments = Math.max(numAlignments, ((List<?>) alnResults).size());
		}
		List<AlignedSequence<VirusT>> alignedSequences = new ArrayList<>();
		Map<String, Sequence> sequenceMap = sequences.stream()
			.collect(Collectors.toMap(seq -> seq.getHeader(), seq -> seq));
		Map<String, Gene<VirusT>> nucaminoGeneMap = strain.getNucaminoGeneMap();
		for (int idx = 0; idx < numAlignments; idx ++) {
			Sequence sequence = null;
			Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> alignedGeneSeqs = new TreeMap<>();
			Map<Gene<VirusT>, String> discardedGenes = new LinkedHashMap<>();

			for (Map.Entry<String, Gene<VirusT>> entry : nucaminoGeneMap.entrySet()) {
				String geneText = entry.getKey();
				Gene<VirusT> gene = entry.getValue();
				Map<?, ?> result = ((Map<?, ?>) ((List<?>) jsonObj.get(geneText)).get(idx));
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
