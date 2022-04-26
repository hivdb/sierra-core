package edu.stanford.hivdb.viruses;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.AssemblyRegion.AssemblyRegionType;

public abstract class Assembler<
	VirusT extends Virus<VirusT>,
	GeneSeqT extends WithGene<VirusT>,
	RegionT extends AssemblyRegion<VirusT, GeneSeqT, RegionT>,
	T
> {
	
	public static <
		VirusT extends Virus<VirusT>,
		GeneSeqT extends WithGene<VirusT>,
		RegionT extends AssemblyRegion<VirusT, GeneSeqT, RegionT>,
		T extends Assembler<VirusT, GeneSeqT, RegionT, T>
	>
	Map<Strain<VirusT>, T> loadJson(
		Class<T> myClass,
		Class<RegionT> regionClass,
		String raw,
		VirusT virusIns) {
		Map<String, List<Map<String, Object>>> rawObj = Json.loads(
			raw,
			new com.google.gson.reflect.TypeToken<Map<String, List<Map<String, Object>>>>(){}
		);
		Map<Strain<VirusT>, T> results = new TreeMap<>();
		for (Map.Entry<String, List<Map<String, Object>>> entry : rawObj.entrySet()) {
			Strain<VirusT> strain = virusIns.getStrain(entry.getKey());
			List<RegionT> regions = new ArrayList<>();
			for (Map<String, Object> rawRegion : entry.getValue()) {
				Double refStart = (Double) rawRegion.get("refStart");
				Double refEnd = (Double) rawRegion.get("refEnd");
				String geneName = (String) rawRegion.get("geneName");
				List<?> trim = (List<?>) rawRegion.get("trim");
				List<Long> trimLong = new ArrayList<>();
				if (trim != null) {
					for (Object trimOne : trim) {
						if (trimOne instanceof List) {
							List<?> trimRange = ((List<?>) trimOne);
							if (trimRange.size() == 1) {
								trimLong.add(((Double) trimRange.get(0)).longValue());
							}
							else if (trimRange.size() > 1) {
								long start = ((Double) trimRange.get(0)).longValue();
								long end = ((Double) trimRange.get(1)).longValue();
								for (long pos = start; pos <= end; pos ++) {
									trimLong.add(pos);
								}
							}
						}
						else if (trimOne instanceof Double) {
							trimLong.add(((Double) trimOne).longValue());
						}
					}
				}
				regions.add(RegionT.of(
					regionClass,
					(String) rawRegion.get("name"),
					AssemblyRegionType.valueOf((String) rawRegion.get("type")),
					refStart == null ? null : refStart.longValue(),
					refEnd == null ? null : refEnd.longValue(),
					geneName == null ? null : virusIns.getGene(geneName),
					trim == null ? null : trimLong
				));
			}
			T assembler = T.of(myClass, regions);
			results.put(strain, assembler);
		}
		return results;
	}
	
	public static <
		VirusT extends Virus<VirusT>,
		GeneSeqT extends WithGene<VirusT>,
		RegionT extends AssemblyRegion<VirusT, GeneSeqT, RegionT>,
		T extends Assembler<VirusT, GeneSeqT, RegionT, T>
	>
	T of(Class<T> myClass, List<RegionT> regions) {
		
		
		try {
			T instance = myClass
				.getDeclaredConstructor(List.class)
				.newInstance(regions);
			return instance;
		} catch (
			InstantiationException
			| IllegalAccessException
			| IllegalArgumentException
			| InvocationTargetException
			| NoSuchMethodException
			| SecurityException e
		) {
			throw new RuntimeException(e);
		}
	}
	
	private final List<RegionT> regions;
	
	protected Assembler(List<RegionT> regions) {
		this.regions = Collections.unmodifiableList(regions);
	}
	
	public List<RegionT> getRegions() { return regions; }
	
	public abstract String assemble(
		final Map<Gene<VirusT>, GeneSeqT> alignedGeneSeqs,
		final Map<String, UntranslatedRegion> utrLookup,
		final boolean skipIns
	);
	
}