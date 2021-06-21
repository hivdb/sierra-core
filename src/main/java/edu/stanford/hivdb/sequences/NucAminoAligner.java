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

import com.amazonaws.services.lambda.AWSLambda;
import com.amazonaws.services.lambda.AWSLambdaClientBuilder;
import com.amazonaws.services.lambda.model.InvokeRequest;
import com.amazonaws.services.lambda.model.InvokeResult;
import com.google.gson.reflect.TypeToken;

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
public class NucAminoAligner<VirusT extends Virus<VirusT>> implements Aligner<VirusT> {
	private final double MIN_MATCH_PCNT = 60;
	private final double SEQUENCE_SHRINKAGE_WINDOW = 15;
	private final double SEQUENCE_SHRINKAGE_CUTOFF_PCNT = 30;
	private final Executor executor = Executors.newFixedThreadPool(20);
	private final Map<Strain<VirusT>, String[]> NUCAMINO_LOCAL_COMMANDS;
	
	private final VirusT virusInstance;

	private static String getJoinedNucaminoGenes(Strain<?> strain) {
		return strain.getNucaminoGeneMap()
			.keySet().stream()
			.collect(Collectors.joining(","));
	}
	
	protected NucAminoAligner(VirusT virusIns) {
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
	
	@Override
	public VirusT getVirusInstance() { return virusInstance; }
	
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
	
	@Override
	public Map<Sequence, AlignedSequence<VirusT>> commandParallelAlign(
		Collection<Sequence> sequences,
		boolean reversingSequence,
		Map<Sequence, Map<Strain<VirusT>, StringBuilder>> errors
	) {
		Map<Strain<VirusT>, List<String>> jsonStrings;
		
		String awsFunc = System.getenv("NUCAMINO_AWS_LAMBDA");
		if (awsFunc == null || awsFunc.equals("")) {
			jsonStrings = localNucamino(sequences);
		} else {
			jsonStrings = awsNucamino(sequences, awsFunc);
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
		return results;
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
					Map<Gene<VirusT>, Double> minMatchPcnt = new HashMap<>();
					minMatchPcnt.put(gene, MIN_MATCH_PCNT);
					Map<Gene<VirusT>, Double> seqShrinkWindow = new HashMap<>();
					seqShrinkWindow.put(gene, SEQUENCE_SHRINKAGE_WINDOW);
					Map<Gene<VirusT>, Double> seqShrinkCutoff = new HashMap<>();
					seqShrinkWindow.put(gene, SEQUENCE_SHRINKAGE_CUTOFF_PCNT);
					try {
						alignedGeneSeqs.put(
							gene,
							geneSeqFromReport(
								sequence,
								gene,
								report,
								sequenceReversed,
								minMatchPcnt.get(gene),
								seqShrinkWindow.get(gene),
								seqShrinkCutoff.get(gene),
								0 // minNumOfSites
							));
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
