package edu.stanford.hivdb.utilities;

public class AssertUtils {
	
	/**
	 * Assert if input value is null
	 * 
	 * @param <T>
	 * @param value
	 * @param msg
	 * @return
	 */
	public static <T> T notNull(T value, String msg, Object... args) {
		if (value == null) {
			throw new IllegalArgumentException(String.format(msg, args));
		}
		return value;
	}
	
	public static void isTrue(Boolean value, String msg, Object... args) {
		if (!value) {
			throw new IllegalArgumentException(String.format(msg, args));
		}
	}
	
}
