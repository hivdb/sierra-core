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
import org.apache.commons.lang3.tuple.Pair;

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
	private final Map<String, Double> MIN_MATCH_PCNT;
	private final Map<String, Double> SEQUENCE_SHRINKAGE_WINDOW;
	private final Map<String, Double> SEQUENCE_SHRINKAGE_CUTOFF_PCNT;
	private final Map<String, Double> MIN_NUM_OF_AA;
	private final Map<String, File> REF_SEQUENCE;
	private final Map<String, List<String>> POST_PROCESSORS;
	private final Map<String, String> FROM_FRAGMENT;
	private final Map<String, Gene<VirusT>> GENE;
	private final Map<String, List<Pair<Long, Long>>> REF_RANGES;
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
		MIN_NUM_OF_AA = Collections.unmodifiableMap(alignConfig.getConfigField("minNumOfAA"));
		FROM_FRAGMENT = Collections.unmodifiableMap(alignConfig.getConfigField("fromFragment"));
		GENE = Collections.unmodifiableMap(alignConfig.getConfigField("gene"));
		REF_RANGES = Collections.unmodifiableMap(alignConfig.getConfigField("refRanges"));

		Map<String, String> refSeqText = alignConfig.getConfigField("refSequence");
		REF_SEQUENCE = Collections.unmodifiableMap(
			refSeqText.entrySet().stream()
			.filter(e -> e.getValue() != null)
			.collect(Collectors.toMap(
				e -> e.getKey(),
				e -> makeReferenceFASTA(e.getKey(), e.getValue())
			))
		);

		Map<String, List<String>> postProcessors = alignConfig.getConfigField("postProcessors");
		POST_PROCESSORS = Collections.unmodifiableMap(
			postProcessors.entrySet()
			.stream()
			.collect(
				Collectors.toMap(
					e -> e.getKey(),
					e -> {
						List<String> cmds = e.getValue();
						if (cmds == null) {
							return Collections.emptyList();
						}
						return cmds
							.stream()
							.filter(n -> !n.startsWith("#"))
							.collect(Collectors.toList());
					}
				)
			)
		);
	}


	private File makeReferenceFASTA(String fragmentName, String refSeq) {
		File referenceFasta;
		try {
			referenceFasta = File.createTempFile("postalign-ref-" + fragmentName, ".fas");
			BufferedWriter bw = new BufferedWriter(new FileWriter(referenceFasta));
			bw.write(">Ref_" + fragmentName + "\n" + refSeq + "\n");
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
	protected Map<String, List<Map<String, ?>>> execute(Collection<Sequence> sequences) {
		Map<String, CompletableFuture<List<Map<String, ?>>>> futures = new TreeMap<>();

		for (String refFragmentName : REF_SEQUENCE.keySet()) {
			File refSeqFile = REF_SEQUENCE.get(refFragmentName);
			List<String> cmd = Lists.newArrayList(
				EXECUTABLE,
				"-i", "-",
				"-o", "-",
				"-f", "MINIMAP2",
				"-r", refSeqFile.getAbsolutePath()
			);
			cmd.addAll(POST_PROCESSORS.get(refFragmentName));
			cmd.add("save-json");
			for (String fragmentName : FROM_FRAGMENT.keySet()) {
				String fromFragment = FROM_FRAGMENT.get(fragmentName);
				if (fromFragment != null && fromFragment.equals(refFragmentName)) {
					cmd.add(fragmentName);
					for (Pair<Long, Long> range : REF_RANGES.get(fragmentName)) {
						cmd.add(range.getLeft().toString());
						cmd.add(range.getRight().toString());
					}
				}
			}
			CompletableFuture<List<Map<String, ?>>> future = CompletableFuture.supplyAsync(() -> {
				List<Map<String, ?>> jsonObjs = new ArrayList<>();
				Iterable<List<Sequence>> partialSets = Iterables.partition(sequences, 100);
				for (List<Sequence> partialSet : partialSets) {
					try {
						Process proc = Runtime.getRuntime().exec(cmd.stream().toArray(String[]::new));
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
			futures.put(refFragmentName, future);
		}

		CompletableFuture.allOf(futures.values().toArray(new CompletableFuture<?>[0])).join();

		Map<String, List<Map<String, ?>>> results = new TreeMap<>();
		for (String refFragmentName : futures.keySet()) {
			CompletableFuture<List<Map<String, ?>>> future = futures.get(refFragmentName);
			try {
				results.put(refFragmentName, future.get());
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			} catch (ExecutionException e) {
				throw new RuntimeException(e);
			}
		}
		return results;
	}

	@Override
	public Map<Sequence, AlignedSequence<VirusT>> commandParallelAlign(
		List<Sequence> sequences,
		boolean reversingSequence,
		Map<Sequence, Map<Strain<VirusT>, StringBuilder>> errors
	) {
		Map<String, List<Map<String, ?>>> allJsonObjs = execute(sequences);
		Map<Sequence, AlignedSequence<VirusT>> results = new LinkedHashMap<>();
		Map<Strain<VirusT>, Map<Gene<VirusT>, String>> geneFragmentLookup = (
			GENE.entrySet().stream()
			.filter(e -> e.getValue() != null)
			.collect(Collectors.toMap(
				e -> e.getValue().getStrain(),
				e -> new TreeMap<>(Map.of(e.getValue(), e.getKey())),
				(a, b) -> {
					b.forEach((k, v) -> a.merge(k, v, (v1, v2) -> {
						throw new RuntimeException(String.format(
							"Same `geneName` %s is referred in multiple fragments",
							k.name()
						));
					}));
					return a;
				},
				TreeMap::new
			))
		);

		for (Strain<VirusT> strain : geneFragmentLookup.keySet()) {
			Map<Gene<VirusT>, String> geneLookup = geneFragmentLookup.get(strain);
			List<AlignedSequence<VirusT>> alignedSeqs = processCommandOutput(
				strain, sequences, allJsonObjs, geneLookup,
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
			List<Sequence> sequences,
			Map<String, List<Map<String, ?>>> jsonObjs,
			Map<Gene<VirusT>, String> geneLookup,
			boolean sequenceReversed,
			Map<Sequence, Map<Strain<VirusT>,
			StringBuilder>> errors) {

		int numSequences = sequences.size();
		List<AlignedSequence<VirusT>> alignedSequences = new ArrayList<>();

		for (int idx = 0; idx < numSequences; idx ++) {
			Sequence sequence = sequences.get(idx);
			Map<Gene<VirusT>, AlignedGeneSeq<VirusT>> alignedGeneSeqs = new TreeMap<>();
			Map<Gene<VirusT>, String> discardedGenes = new LinkedHashMap<>();

			for (Gene<VirusT> gene : geneLookup.keySet()) {
				String fragmentName = geneLookup.get(gene);
				String refFragmentName = FROM_FRAGMENT.get(fragmentName);
				Map<String, ?> result = jsonObjs.get(refFragmentName).get(idx);

				Map<?, ?> geneResult = ((List<?>) result.get("GeneReports"))
					.stream()
					.filter(map -> (
						((String) ((Map<?, ?>) map).get("Gene")).equals(fragmentName)
					))
					.map(map -> (Map<?, ?>) map)
					.findFirst()
					.get();
				Map<?, ?> report = (Map<?, ?>) geneResult.get("Report");
				String error = (String) geneResult.get("Error");
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
							MIN_MATCH_PCNT.get(fragmentName),
							SEQUENCE_SHRINKAGE_WINDOW.get(fragmentName),
							SEQUENCE_SHRINKAGE_CUTOFF_PCNT.get(fragmentName),
							MIN_NUM_OF_AA.get(fragmentName).intValue()
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
