package edu.stanford.hivdb.utilities;

import org.apache.commons.lang3.concurrent.Computable;
import org.apache.commons.lang3.concurrent.Memoizer;

public class SimpleMemoizer<O> extends Memoizer<String, O> {
	
	public SimpleMemoizer(Computable<String, O> computable) {
		super(computable);
	}
	
	public O get(String name) {
		try {
			return super.compute(name);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
	}

}
