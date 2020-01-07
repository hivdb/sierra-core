/*

    Copyright (C) 2019 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.mutations;

import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class CodonModifier<VirusT extends Virus<VirusT>> {
	
	private final int position;
	private Integer insertAfter;
	private Integer deleteAfter;
	private String targetStrain;

	private transient VirusT virusInstance;
	
	private CodonModifier(
		int position, int insertAfter,
		int deleteAfter, String targetStrain
	) {
		this.position = position;
		this.insertAfter = insertAfter;
		this.deleteAfter = deleteAfter;
		this.targetStrain = targetStrain;
	}
	
	public void setVirusInstance(VirusT virusInstance) {
		this.virusInstance = virusInstance;
	}
	
	private void checkVirusInstance() {
		if (virusInstance == null) {
			throw new ExceptionInInitializerError(
				"Object not properly initialzed: virusInstance is empty"
			);
		}
	}
	
	public Strain<VirusT> getTargetStrain() {
		checkVirusInstance();
		return virusInstance.getStrain(targetStrain);
	}
	
	public int getPosition() {
		return position;
	}
	
	public Integer getInsertAfter() {
		return insertAfter;
	}
	
	public Integer getDeleteAfter() {
		return deleteAfter;
	}

}
