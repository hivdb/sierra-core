package edu.stanford.hivdb.sequences;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;

public class AlignmentConfig<VirusT extends Virus<VirusT>> {

	private final AlignmentMethod method;
	private final Map<Gene<VirusT>, Map<String, ?>> geneConfig;
	
	public static <VirusT extends Virus<VirusT>> AlignmentConfig<VirusT> loadJson(String raw, VirusT virusIns) {
		Map<String, ?> configMap = Json.loads(raw, new TypeToken<Map<String, ?>>(){});
		AlignmentMethod method = AlignmentMethod.valueOf((String) configMap.get("method"));
		@SuppressWarnings("unchecked")
		List<Map<String, ?>> geneConfigList = (List<Map<String, ?>>) configMap.get("geneConfig");
		Map<Gene<VirusT>, Map<String, ?>> geneConfig = new TreeMap<>();
		for (Map<String, ?> config : geneConfigList) {
			Gene<VirusT> gene = virusIns.getGene((String) config.get("geneName"));
			geneConfig.put(gene, config);
		}
		return new AlignmentConfig<>(method, geneConfig);
	}
	
	private AlignmentConfig(AlignmentMethod method, Map<Gene<VirusT>, Map<String, ?>> geneConfig) {
		this.method = method;
		this.geneConfig = geneConfig;
	}
	
	public AlignmentMethod getMethod() {
		return method;
	}
	
	@SuppressWarnings("unchecked")
	public <T> Map<Gene<VirusT>, T> getConfigField(String fieldName) {
		Map<Gene<VirusT>, T> resultMap = new TreeMap<>();
		for (Gene<VirusT> gene : geneConfig.keySet()) {
			resultMap.put(gene, (T) geneConfig.get(gene).get(fieldName));
		}
		return resultMap;
	}

}
