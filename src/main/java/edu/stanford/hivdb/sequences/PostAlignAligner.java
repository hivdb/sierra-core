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

	@Override
	public VirusT getVirusInstance() { return virusInstance; }
	
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

	@Override
	public Map<Sequence, AlignedSequence<VirusT>> commandParallelAlign(
		Collection<Sequence> sequences,
		boolean reversingSequence,
		Map<Sequence, Map<Strain<VirusT>, StringBuilder>> errors
	) {
		Map<Strain<VirusT>, Map<Gene<VirusT>, List<Map<String, ?>>>> allJsonObjs;
		
		allJsonObjs = execute(sequences);
		
		Map<Sequence, AlignedSequence<VirusT>> results = new LinkedHashMap<>();
		for (Strain<VirusT> strain : allJsonObjs.keySet()) {
			List<AlignedSequence<VirusT>> alignedSeqs = processCommandOutput(
				strain, sequences, allJsonObjs.get(strain),
				reversingSequence, errors
			);
			results = selectBestAlignments(alignedSeqs, results);
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
						alignedGeneSeqs.put(gene, geneSeqFromReport(
							sequence,
							gene,
							report,
							sequenceReversed,
							MIN_MATCH_PCNT,
							SEQUENCE_SHRINKAGE_WINDOW,
							SEQUENCE_SHRINKAGE_CUTOFF_PCNT
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
