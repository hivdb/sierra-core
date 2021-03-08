package edu.stanford.hivdb.viruses;

import graphql.schema.GraphQLCodeRegistry;
import graphql.schema.GraphQLObjectType;

public interface VirusGraphQLExtension {

	public default GraphQLCodeRegistry getExtendedCodonRegistry() {
		return GraphQLCodeRegistry.newCodeRegistry().build();
	}
	
	public default GraphQLObjectType.Builder extendObjectBuilder(String objectName, GraphQLObjectType.Builder builder) {
		return builder;
	}
	
}