package edu.stanford.hivdb.viruses;

public class DefaultVirusGraphQLExtension implements VirusGraphQLExtension {
	private static DefaultVirusGraphQLExtension singleton = new DefaultVirusGraphQLExtension();

	private DefaultVirusGraphQLExtension() {}
	
	public static DefaultVirusGraphQLExtension getInstance() { return singleton; }
	
}