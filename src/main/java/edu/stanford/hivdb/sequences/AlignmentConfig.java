package edu.stanford.hivdb.sequences;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Pair;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

public class AlignmentConfig<VirusT extends Virus<VirusT>> {

	private final AlignmentMethod method;
	private final Map<String, Map<String, ?>> fragmentConfig;

	public static <VirusT extends Virus<VirusT>> AlignmentConfig<VirusT> loadJson(String raw, VirusT virusIns) {
		Map<String, ?> configMap = Json.loads(raw, new TypeToken<LinkedHashMap<String, ?>>(){});
		AlignmentMethod method = AlignmentMethod.valueOf((String) configMap.get("method"));
		@SuppressWarnings("unchecked")
		List<Map<String, Object>> fragmentConfigList = (List<Map<String, Object>>) configMap.get("fragmentConfig");
		Map<String, Map<String, ?>> fragmentConfig = new TreeMap<>();
		for (Map<String, Object> config : fragmentConfigList) {
			String fragmentName = (String) config.get("fragmentName");

			String geneName = (String) config.getOrDefault("geneName", null);
			if (geneName != null) {
				// add Gene object
				Gene<VirusT> gene = virusIns.getGene(geneName);
				config.put("gene", gene);
			}

			List<?> refRangesRaw = (List<?>) config.getOrDefault("refRanges", null);
			try {
				if (refRangesRaw != null) {
					List<Pair<Long, Long>> refRanges = (
						refRangesRaw
						.stream()
						.map(
							pair -> {
								Long left = ((Double) ((List<?>) pair).get(0)).longValue();
								Long right = ((Double) ((List<?>) pair).get(1)).longValue();
								if (left < 1) {
									throw new IllegalArgumentException(
										String.format(
											"Parameter `refRanges` contains start position < 1 (fragment %s)",
											fragmentName
										)
									);
								}
								if (right - 2 < left) {
									throw new IllegalArgumentException(
										String.format(
											"Parameter `refRanges` contains a range %d-%d is too short to include at least one codon (fragment %s)",
											left, right, fragmentName
										)
									);
								}
								return Pair.of(left, right);
							}
						)
						.collect(Collectors.toList())
					);
					config.put("refRanges", refRanges);
				}
			}
			catch (NullPointerException | IndexOutOfBoundsException exc) {
					throw new IllegalArgumentException(
						String.format(
							"Parameter `refRanges` is illegal (fragment %s)", fragmentName
						)
					);

			}
			fragmentConfig.put(fragmentName, config);
		}
		return new AlignmentConfig<>(method, fragmentConfig);
	}

	private AlignmentConfig(AlignmentMethod method, Map<String, Map<String, ?>> fragmentConfig) {
		this.method = method;
		this.fragmentConfig = fragmentConfig;
	}

	public AlignmentMethod getMethod() {
		return method;
	}

	@SuppressWarnings("unchecked")
	public <T> Map<String, T> getConfigField(String fieldName) {
		Map<String, T> resultMap = new TreeMap<>();

		for (String fragmentName : fragmentConfig.keySet()) {
			Object val = fragmentConfig.get(fragmentName).getOrDefault(fieldName, null);
			resultMap.put(fragmentName, (T) val);
		}
		return resultMap;
	}

}
