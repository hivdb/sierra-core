package edu.stanford.hivdb.drugs;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.fstrf.stanfordAsiInterpreter.resistance.ASIParsingException;
import org.fstrf.stanfordAsiInterpreter.resistance.xml.XmlAsiTransformer;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class DrugResistanceAlgorithm<VirusT extends Virus<VirusT>> {
	
	final private String name;
	final private String family;
	final private String version;
	final private String publishDate;
	final private Strain<VirusT> strain;
	final private String xmlText;
	final private transient String originalLevelText;
	final private transient String originalLevelSIR;
	final private transient Map<Gene<VirusT>, org.fstrf.stanfordAsiInterpreter.resistance.definition.Gene> geneMap;

	private static <VirusT extends Virus<VirusT>> Map<Gene<VirusT>, org.fstrf.stanfordAsiInterpreter.resistance.definition.Gene> initGeneMap(
		String xmlText, Strain<VirusT> strain
	) {
		InputStream resource = new ByteArrayInputStream(xmlText.getBytes());
		Map<?, ?> geneMap;
		XmlAsiTransformer transformer = new XmlAsiTransformer(true);
		try {
			geneMap = transformer.transform(resource);
		} catch (Exception e) {
			throw new ExceptionInInitializerError(e);
		}
		Set<String> absGenes = (
			strain.getGenes().stream()
			.map(gene -> gene.getAbstractGene())
			.collect(Collectors.toSet())
		);
		return Collections.unmodifiableMap(geneMap
			.entrySet()
			.stream()
			// removes genes supported by algorithm but not by Virus implementation
			.filter(e -> absGenes.contains(e.getKey()))
			.collect(Collectors.toMap(
				e -> strain.getGene((String) e.getKey()),
				e -> (org.fstrf.stanfordAsiInterpreter.resistance.definition.Gene) e.getValue()
			))
		);
	}
	
	private static Map<?, ?> getAlgorithmInfo(String xmlText) {
		InputStream resource = new ByteArrayInputStream(xmlText.getBytes());
		XmlAsiTransformer transformer = new XmlAsiTransformer(true);
		try {
			Map<?, ?> algInfo = (Map<?, ?>) transformer.getAlgorithmInfo(resource);
			return algInfo;
		} catch (ASIParsingException e) {
			throw new ExceptionInInitializerError(e);
		}
	}
	
	public DrugResistanceAlgorithm(Strain<VirusT> strain, String xmlText) {
		this(null, null, null, null, strain, xmlText);
	}
	
	public DrugResistanceAlgorithm(String name, Strain<VirusT> strain, String xmlText) {
		this(name, null, null, null, strain, xmlText);
	}
	
	public DrugResistanceAlgorithm(
		String name, String family, String version,
		String publishDate, Strain<VirusT> strain, String xmlText
	) {
		Map<?, ?> algInfo = getAlgorithmInfo(xmlText);
		Map<?, ?> algNVD = (Map<?, ?>) algInfo.get("ALGNAME_ALGVERSION_ALGDATE");
		Map<?, ?> originalLevel = (Map<?, ?>) algInfo.get("ORDER1_ORIGINAL_SIR");
		this.name = name == null ? String.format("%s_%s", algNVD.get("ALGNAME"), algNVD.get("ALGVERSION")) : name;
		this.family = family == null ? (String) algNVD.get("ALGNAME") : family;
		this.version = version == null ? (String) algNVD.get("ALGVERSION") : version;
		this.publishDate = publishDate == null ? (String) algNVD.get("ALGDATE") : publishDate;
		this.originalLevelText = (String) originalLevel.get("ORIGINAL");
		this.originalLevelSIR = (String) originalLevel.get("SIR");
		this.strain = strain;
		this.xmlText = xmlText;
		this.geneMap = initGeneMap(xmlText, strain);
	}
	
	public String getName() {
		return name;
	}
	
	public String getDisplay() {
		return String.format("%s %s", family, version);
	}
	
	public String name() {
		return name;
	}
	
	public String getFamily() {
		return family;
	}
	
	public String getVersion() {
		return version;
	}
	
	public String getPublishDate() {
		return publishDate;
	}
	
	public String getOriginalLevelText() {
		return originalLevelText;
	}
	
	public String getOriginalLevelSIR() {
		return originalLevelSIR;
	}
	
	public Strain<VirusT> getStrain() {
		return strain;
	}
	
	public org.fstrf.stanfordAsiInterpreter.resistance.definition.Gene getASIGene(Gene<VirusT> gene) {
		return geneMap.get(gene);
	}
	
	public String getXMLText() {
		return xmlText;
	}
	
	public String getEnumCompatName() {
		String name = getName().replaceAll("[^_0-9A-Za-z-]", "_");
		name = name.replace("-stanford", "stanford").replace('-', 'p');
		if (name.matches("^\\d")) {
			name = "_" + name;
		}
		return name;
	}
	
	@Override
	public String toString() {
		return getName();
	}

}
