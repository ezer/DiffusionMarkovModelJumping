package edu.cam.fgrip;

import java.util.List;


public abstract class ProbabilityModel {
	protected Parameters params;
	
	public ProbabilityModel(Parameters p){
		params=p;
	}
	
	public abstract double forwardProbability(int oldID, int tfID, List<TranscriptionFactor> tfs, double totalTF);
	
	public abstract double backwardProbability(int oldID, int tfID, List<TranscriptionFactor> tfs);
	
	public abstract double getSlideProbability(int tfIDstart, int tfIDstop, int oldID, List<TranscriptionFactor> tfs);
	
	//TODO move this to some other class
	protected int numTFsBasedOnState(int oldID, int tfID, List<TranscriptionFactor> tfs) {
		int usedUp=0;
		int id=oldID;
		int i=0;
		while(id>0){
			if((id&1)>0 && tfs.get(i).hasSameNameAs(tfs.get(tfID))){
				
				usedUp++;
			}
			id=id >>> 1;
			i++;
		}
		return usedUp;
	}

	

	
}
