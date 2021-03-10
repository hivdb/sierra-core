package edu.stanford.hivdb.sequences;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

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
			
			Double refStart = (Double) config.getOrDefault("refStart", null);
			Double refEnd = (Double) config.getOrDefault("refEnd", null);
			if (refStart != null || refEnd != null) {
				// convert refStart and refEnd from Double to Long
				try {
					config.put("refStart", Long.valueOf(refStart.longValue()));
					config.put("refEnd", Long.valueOf(refEnd.longValue()));
				}
				catch (NullPointerException exc) {
					throw new IllegalArgumentException(
						String.format(
							"Parameters `refStart` and `refEnd` must be both presented for fragment %s", fragmentName
						)
					);
				}
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
			resultMap.put(fragmentName, (T) fragmentConfig.get(fragmentName).getOrDefault(fieldName, null));
		}
		return resultMap;
	}

}
