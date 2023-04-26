package edu.stanford.hivdb.sequences;

import org.apache.commons.lang3.builder.HashCodeBuilder;

public class AlignmentMessage implements Comparable<AlignmentMessage> { 
	
	public enum AlignmentMessageLevel {
		INFO,
		WARNING,
		ERROR
	}
	
	private final AlignmentMessageLevel level;
	private final String message;
	
	public AlignmentMessage(String level, String message) {
		this.level = AlignmentMessageLevel.valueOf(level);
		this.message = message;
	}
	
	public AlignmentMessageLevel getLevel() {
		return level;
	}
	
	public String getMessage() {
		return message;
	}

	@Override
	public String toString() {
		return String.format("%s: %s", level, message);
	}
	
	@Override
	public int hashCode() {
		return new HashCodeBuilder(467547, 323223465)
			.append(level)
			.append(message)
			.toHashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (this == o) { return true; }
		if (o == null) { return false; }
		if (!(o instanceof AlignmentMessage)) { return false; }
		AlignmentMessage other = (AlignmentMessage) o;
		return toString().equals(other.toString());
	}

	@Override
	public int compareTo(AlignmentMessage o) {
		if (level != o.level) {
			return level.compareTo(o.level);
		}
		if (!message.equals(o.message)) {
			return message.compareTo(o.message);
		}
		return 0;
	}

}
