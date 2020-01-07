package edu.stanford.hivdb.drugresistance.algorithm;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

public class AsiDrugComparableResult {
	public final SIREnum SIR;
	public final String Interpretation;
	public final String Explanation;

	public AsiDrugComparableResult(
			String SIR, String Interpretation, String Explanation) {
		this.SIR = SIREnum.valueOf(SIR);
		this.Interpretation = Interpretation;
		this.Explanation = Explanation;
	}

	public String getInterpretation() { return Interpretation; }
	public String getExplanation() { return Explanation; }

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		if (o == null) { return false; }
		if (!(o instanceof AsiDrugComparableResult)) { return false;}
		AsiDrugComparableResult a = (AsiDrugComparableResult) o;

		return new EqualsBuilder()
			.append(SIR, a.SIR)
			.append(Interpretation, a.Interpretation)
			.append(Explanation, a.Explanation)
			.isEquals();
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(84535, 467345)
			.append(SIR)
			.append(Interpretation)
			.append(Explanation)
			.toHashCode();
	}

	@Override
	public String toString() {
		return String.format(
			"AsiDrugComparableResult(%s, %s, %s)",
			"\"" + SIR + "\"",
			"\"" + Interpretation.replace("\"", "\\\"") + "\"",
			"\"" + Explanation.replace("\"", "\\\"") + "\"");
	}
}