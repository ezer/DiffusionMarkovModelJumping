package edu.cam.fgrip;

import java.util.ArrayList;
import java.util.List;


public class SwitchSlideModel extends ProbabilityModel{
	
	public SwitchSlideModel(Parameters p) {
		super(p);
	}
	
	@Override
	public double forwardProbability(int oldID, int tfID,
			List<TranscriptionFactor> tfs, double totalTF) {
	
	
		int usedUp = numTFsBasedOnState(oldID, tfID, tfs);
		double rightmostBound=0;
		double leftmostBound=0;
		double s_l_real=params.get("s_l_real");
		
	
		
		
		rightmostBound = Math.min((int)(s_l_real/2), getRightBound(oldID, tfID, tfs));
		leftmostBound = Math.min((int)(s_l_real/2), getLeftBound(oldID, tfID, tfs));
		/*
		if(oldID==3){
			System.out.println("&&&&&&&&&&&&&&&&&&&&&&&&&&&");
			System.out.println("forward of 3: "+oldID+" adding "+tfID+" totalTF="+totalTF);
			System.out.println("rightbound: "+getRightBound(oldID, tfID, tfs)+" sliding:"+(int)(s_l_real/2));
			System.out.println("lefboundbound: "+getLeftBound(oldID, tfID, tfs)+" sliding:"+(int)(s_l_real/2));
			System.out.println("model Length overall:"+(tfs.get(0).getCoveringTo(tfs.get(tfs.size()-1))+2*params.get("s_l_obs")));
		    System.out.println("composed of: "+tfs.get(0).getCoveringTo(tfs.get(tfs.size()-1))+" "+2*params.get("s_l_obs"));
		    System.out.println("binding input: "+tfs.get(tfID).getStart()+" "+(leftmostBound+rightmostBound)+" "+usedUp+" "+totalTF);
		    System.out.println("&&&&&&&&&&&&&&&&&&&&&&&&&&&"); 
		}
		*/
		//System.out.println("^^^^^^^^^ leftmostbound: "+leftmostBound+"   rightmostBound: "+rightmostBound);
		double modelLength=tfs.get(0).getCoveringTo(tfs.get(tfs.size()-1))+2*params.get("s_l_obs");
		//System.out.println("oldID: "+oldID+ " TFIF "+tfID+" modelLength: "+modelLength+" usedUp: "+usedUp+" totalTF: "+totalTF);
		
		return propensityBindingUniform(tfs.get(tfID), leftmostBound+rightmostBound, usedUp, totalTF, (int)modelLength, tfs);
	}

	private double getLeftBound(int oldID, int tfID, List<TranscriptionFactor> tfs) {
		double leftmostBound;
		int id;
		int i;
		i=tfID+1;
		id=oldID>>>i;
		//System.out.println("Starts with: oldID: "+oldID+" tfID "+tfID+"iter i "+i+" id "+id);
		while(i<tfs.size() && ((id&1)==0) && !tfs.get(i).hasSameNameAs(tfs.get(tfID))){
			i++;
			id=id>>>1;
			//System.out.println("iter i "+i+" id "+id);
		}
		
		double maxSlide=params.get("s_l_obs")/2.0;
		if(i>=tfs.size())
			leftmostBound=maxSlide;
		else if((id&1)==0 && tfs.get(i).hasSameNameAs(tfs.get(tfID))){ //cluster
			leftmostBound= Math.min(maxSlide, tfs.get(tfID).getDistanceTo(tfs.get(i))/2.0);
		}
		else{ //barrier
			leftmostBound= Math.min(maxSlide, tfs.get(i).getDistanceBetween(tfs.get(tfID)));
		}
		return leftmostBound;
	}

	private double getRightBound(int oldID, int tfID, List<TranscriptionFactor> tfs) {
	
		double maxSlide=params.get("s_l_obs")/2.0;
		
		
		double rightmostBound;
		int id=oldID;
		int i=tfID-1;
		while(i>=0 && (((id>>>i)&1)==0 && !tfs.get(i).hasSameNameAs(tfs.get(tfID))) ){
			i--;
			id=oldID;
		}
		
	
		
		if(i<0)
			rightmostBound=maxSlide;
		else if(((id>>>i)&1)==0 && tfs.get(i).hasSameNameAs(tfs.get(tfID))){ //cluster
			rightmostBound= Math.min(maxSlide, tfs.get(tfID).getDistanceTo(tfs.get(i))/2.0);
		}
		else{ //barrier
			rightmostBound= Math.min(maxSlide,tfs.get(tfID).getDistanceBetween(tfs.get(i)));
		}
		return rightmostBound;
	}

	//TODO


	
	
	//TODO clean this up: hack to test new functiom
	private int getLeftBoundWithoutCluster(int oldID, int tfID, List<TranscriptionFactor> tfs) {
		int leftmostBound;
		int id;
		int i;
		i=tfID+1;
		id=oldID>>>i;
		//System.out.println("Starts with: oldID: "+oldID+" tfID "+tfID+"iter i "+i+" id "+id);
		while(i<tfs.size() && ((id&1)==0) ){
			i++;
			id=id>>>1;
			//System.out.println("iter i "+i+" id "+id);
		}
		
		double maxSlide=params.get("s_l_obs")/2;
		if(i>=tfs.size())
			leftmostBound=(int)maxSlide;
		else{
			leftmostBound=(int) Math.min(maxSlide, tfs.get(i).getDistanceBetween(tfs.get(tfID)));
		}
		return leftmostBound;
	}

	private int getRightBoundWithoutCluster(int oldID, int tfID, List<TranscriptionFactor> tfs) {
		int rightmostBound;
		int id=oldID;
		int i=tfID-1;
		while(i>=0 && (((id>>>i)&1)==0 ) ){
			i--;
			id=oldID;
		}
		
		double maxSlide=params.get("s_l_obs")/2;
		
		if(i<0)
			rightmostBound=(int)maxSlide;
		else{
			rightmostBound=(int) Math.min(maxSlide,tfs.get(tfID).getDistanceBetween(tfs.get(i)));
		}
		return rightmostBound;
	}

	
	
	
	
	
	
	@Override
	public double backwardProbability(int oldID, int tfID,
			List<TranscriptionFactor> tfs) {
		//TODO
		double s_l_real=params.get("s_l_real");
		double rightmostBound = Math.min((int)(s_l_real/2), getRightBound(oldID, tfID, tfs));
		double leftmostBound = Math.min((int)(s_l_real/2), getLeftBound(oldID, tfID, tfs));
		
		//int rightmostBound = //getRightBoundWithoutCluster(oldID, tfID, tfs);
		//int leftmostBound = //getLeftBoundWithoutCluster(oldID, tfID, tfs);
		//List<Integer> clusters=countClusters(rightmostBound, leftmostBound, tfID, tfs);
		//List<Integer> clusters=new ArrayList<Integer>();
			//System.out.println("for state: "+oldID+" ridding tf: "+tfID+" distance: "+(rightmostBound+leftmostBound));
		double expectedVisitsOrigin=tfs.get(tfID).visitsToOrigin(rightmostBound+leftmostBound);//tfs.get(tfID).visitsToOrigin(rightmostBound+leftmostBound-tfs.get(tfID).getSize());
			//System.out.println("expected visits origin: "+expectedVisitsOrigin+" clusters size: "+clusters.size());
		double expectedVisitsNotOrigin = tfs.get(tfID).slides()-expectedVisitsOrigin;//expectedVisitsOrigin*(1+clusters.size());
		//TODO (only works if bindng sites have same strength)	
		double expectedTimeInOrigin=tfs.get(tfID).timeInBindingSite()*expectedVisitsOrigin;
		double expectedTime=expectedTimeInOrigin + tfs.get(tfID).timeInBackground() * expectedVisitsNotOrigin;
		
		//check:
		/*
		double estimatedClockResets=1;
		for(int i: clusters){
			
			System.out.println("timeToSlide: "+slideProb(tfID, i, oldID, tfs)+"  expectedTime: "+expectedTime);
			//TODO
			estimatedClockResets=estimatedClockResets+1.0/slideProb(tfID, i, oldID, tfs);
		}
		System.out.println("clooooock resets: "+estimatedClockResets);
		*/
		//System.out.println("expectedTimeInOrigin: "+expectedTimeInOrigin);
		//System.out.println("visits: "+expectedVisitsOrigin+", "+expectedVisitsNotOrigin+", time in binding: "+tfs.get(tfID).timeInBindingSite() + ", background: " + tfs.get(tfID).timeInBackground());
		//System.out.println("expected time: "+expectedTime);
		double k=(1.0/(expectedTime));
		return k;
	}
/*
	private double getBackwardTimeBeforeClusterAdjustment(int oldID, int tfID,
			List<TranscriptionFactor> tfs){
		int rightmostBound = getRightBoundWithoutCluster(oldID, tfID, tfs);
		int leftmostBound = getLeftBoundWithoutCluster(oldID, tfID, tfs);
		List<Integer> clusters=countClusters(rightmostBound, leftmostBound, tfID, tfs);
		//List<Integer> clusters=new ArrayList<Integer>();
			System.out.println("for state: "+oldID+" ridding tf: "+tfID+" distance: "+(rightmostBound+leftmostBound));
		double expectedVisitsOrigin=tfs.get(tfID).visitsToOrigin(rightmostBound+leftmostBound);//tfs.get(tfID).visitsToOrigin(rightmostBound+leftmostBound-tfs.get(tfID).getSize());
			System.out.println("expected visits origin: "+expectedVisitsOrigin+" clusters size: "+clusters.size());
		double expectedVisitsNotOrigin = tfs.get(tfID).slides()-expectedVisitsOrigin*(1+clusters.size());
			
		double expectedTimeInOrigin=tfs.get(tfID).timeInBindingSite()*expectedVisitsOrigin*(1+clusters.size());
		double expectedTime=expectedTimeInOrigin + tfs.get(tfID).timeInBackground() * expectedVisitsNotOrigin;
		
		return expectedTime;
	}
	*/
	private List<Integer> countClusters(double rightmostBound, double leftmostBound, int tfID,
			List<TranscriptionFactor> tfs) {
		
		List<Integer> clusters=new ArrayList<Integer>();
		int i=tfID-1;
		while(i>=0 && tfs.get(i).getDistanceBetween(tfs.get(tfID))<=rightmostBound){
			if(tfs.get(i).hasSameNameAs(tfs.get(tfID))){
				clusters.add(i);
			}
			i--;
		}
		i=tfID+1;
		while(i<tfs.size() && tfs.get(i).getDistanceBetween(tfs.get(tfID))<=leftmostBound){
			if(tfs.get(i).hasSameNameAs(tfs.get(tfID))){
				clusters.add(i);
			}
			i++;
		}
		return clusters;
	}

	
	private double slideProb(int tfIDstart, int tfIDstop, int oldID,
			List<TranscriptionFactor> tfs){
		
		int distance=tfs.get(tfIDstart).getDistanceTo(tfs.get(tfIDstop));
		//System.out.println("***************************************  "+distance);
		double N=2*(distance)*(distance);
		return (N+0.0)/tfs.get(tfIDstart).slides();
	}
	
	@Override
	public double getSlideProbability(int tfIDstart, int tfIDstop, int oldID,
			List<TranscriptionFactor> tfs) {	
		/*
		int rightmostBound = getRightBoundWithoutCluster(oldID, tfIDstart, tfs);
		int leftmostBound = getLeftBoundWithoutCluster(oldID, tfIDstart, tfs);
		int rightWithCluster = getRightBound(oldID, tfIDstart, tfs);
		int leftWithCluster = getLeftBound(oldID, tfIDstart, tfs);
		int a=rightWithCluster+leftWithCluster;
		int b=rightmostBound+leftmostBound-a;
		double s=tfs.get(tfIDstart).slides();
		double timeInSite=(s/(a+b))*tfs.get(tfIDstart).timeInBindingSite();
		double expectedTime=timeInSite+(s*a/(a+b)-s/(a+b))*tfs.get(tfIDstart).timeInBackground();
		double k=1.0/expectedTime-backwardProbability(oldID, tfIDstart, tfs);
		*/
		
		double k=0;
		double distance=tfs.get(tfIDstart).getDistanceTo(tfs.get(tfIDstop));
		//System.out.println("***************************************  "+distance);
		double N=2*(distance)*(distance);
		double passesThroughOrigin=N/(distance);
		double t_bp=tfs.get(tfIDstart).timeInBackground();
		k=1.0/((N-passesThroughOrigin)*t_bp+passesThroughOrigin*tfs.get(tfIDstart).timeInBindingSite());
		return k;
		//double slideP=slideProb(tfIDstart, tfIDstop, oldID, tfs);
		//double slideT=slideP*getBackwardTimeBeforeClusterAdjustment(oldID, tfIDstart, tfs);
		//return 1/slideT;  //about 1000x off
	}
	
	
	private double timeOfBinding3D(TranscriptionFactor tf, int usedUpTFs, double totalTF, int modelLength, List<TranscriptionFactor> tfs){	
		if(tf.getConcentration()==0)
			return 0;
		
		/*
		 * propensity of binding 3D summary:
		 * either calculate k_a or just use it
		 * 
		 */
		
		if(params.contains("k_a")){
			return 1.0/params.get("k_a");
		}
		
		
		double t_r=params.get("t_r");
		double nonCognateTF=params.get("nonCognateTF");
		double TFlength=params.get("TFlength");
		double DNAlength=params.get("DNAlength");
		double f=tf.f();
		
		double TF_free=tf.freeTFs(usedUpTFs);
		double TF_bound=tf.getConcentration()-TF_free;

		double A_ratio=1-((totalTF+nonCognateTF)*f*TFlength)/DNAlength;
		double k_a=1/t_r*(TF_bound/TF_free)*1.0/A_ratio;
		//TODO
		k_a=1800;
		double lambda=tf.getLambda(tfs, modelLength);
		//System.out.println("f: "+f);
		double k_a_adj=k_a*(lambda-f*lambda)/(1-f*lambda);
		//TODO
		//k_a_adj=1.31;
		//System.out.println("TF_free: "+TF_free+" TF_bound: "+TF_bound+" A_ratio: "+A_ratio+" k_a: "+k_a+" lambda: "+lambda+" k_a_adj: "+k_a_adj);
		return 1/k_a_adj;
	}
	
	public double propensityBindingUniform(TranscriptionFactor tf, double d, int usedUpTfs, double totalTF, int modelLength, List<TranscriptionFactor> tfs){
		
		double DNAlength=params.get("DNAlength");
		double M=DNAlength;
		double t_3D=timeOfBinding3D(tf, usedUpTfs, totalTF, modelLength, tfs);
		double t_1D=tf.slides()*tf.timeInBackground();
		double t_binding=((M/(d))*(t_1D+t_3D));	
		//System.out.println("t_binding: "+t_binding+" t_3D "+t_3D+" t_1D "+t_1D+" tf.getConcentration: "+tf.getConcentration()+" usedUp: "+usedUpTfs);
		return (1.0/t_binding)*(tf.getConcentration()-usedUpTfs);	
	}
	
}
