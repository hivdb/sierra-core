package edu.stanford.hivdb.utilities;

public interface ValidationMessage {
	
	public ValidationLevel getLevel();
	public String getTemplate();
	
	
	public default ValidationResult format(Object... args) {
		String message = String.format(getTemplate(), args);
		return new ValidationResult(getLevel(), message);
	}

	public default ValidationResult formatWithLevel(ValidationLevel overrideLevel, Object... args) {
		String message = String.format(getTemplate(), args);
		return new ValidationResult(overrideLevel, message);
	}

}
