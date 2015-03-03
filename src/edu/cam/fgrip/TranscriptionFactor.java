package edu.cam.fgrip;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


public class TranscriptionFactor implements Comparable<TranscriptionFactor>{
	
		private String myName;
		private int myStart; 
		private int myStop;
		private double myWeight;
		private double myConcentration;
		private double myLambda=-1;
		
	
		
		private Parameters params;
		private static Map<String, Double> t_0s;
		private static Map<String, Double> fs;
	
	public TranscriptionFactor(String name, int start, int stop, double weight, double concentration, Parameters p){
		myName=name;
		myStart=start;
		myStop=stop;
		myWeight=weight;
		myConcentration=concentration;
		params=p;
		
		if(t_0s==null){
			try{
				
				Scanner in=new Scanner(new File("src/resources/t_0"));
				t_0s=new HashMap<String, Double>();
				fs=new HashMap<String, Double>();
				while(in.hasNext()){
					String tf=in.next().toUpperCase();
					Double val=in.nextDouble();
					Double fVal=in.nextDouble();
					t_0s.put(tf, val);
					fs.put(tf, fVal);
				}
			}catch(Throwable t){
				System.out.println("t_0 file not found");
			}
		}
	}
	
	public boolean hasSameNameAs(TranscriptionFactor tf2){
		return myName.equals(tf2.myName);
	}
	
	public double getLambda(List<TranscriptionFactor> tfs, int distance){
		if(myLambda>0)
			return myLambda;
		else{
			int covered= spaceCoveredByTF(tfs, distance);
			System.out.println("covered: "+covered);
			System.out.println("uncpveered: "+(distance-covered));
			System.out.println("time in binding: "+timeInBindingSite());
			System.out.println("time in background: "+timeInBackground());
			double timeInRegion=(covered*timeInBindingSite()+(distance-covered)*timeInBackground());
			System.out.println("timeInRegion: "+timeInRegion+" denominator "+(timeInRegion+(params.get("DNAlength")-distance)*timeInBackground()));
			myLambda= timeInRegion/(timeInRegion+(4600000-distance)*timeInBackground());
			myLambda= 0.00434;
			return myLambda;
		}
	}
	
	public double f(){
		
		if(fs.containsKey(myName.toUpperCase())){
		
		return fs.get(myName.toUpperCase());
		}
		return params.get("f");
	}
	
	private int spaceCoveredByTF(List<TranscriptionFactor> tfs, int distance){
		return tfs.size();
	}
	
	public String toString(){
		return myName+": "+myStart+", "+myStop;
	}
	
	@Override
	public int compareTo(TranscriptionFactor tf) {
		return myStart-tf.myStart;
	}
	
	public boolean overlapsWith(TranscriptionFactor tf){
		return !((myStart <= tf.myStart && myStop < tf.myStart) || (tf.myStart <= myStart && tf.myStop < myStart));
	}
	
	public double timeInBindingSite(){
	
		
		/*
		 * Two strategies to calculate this
		 * 
		 * Case 1: each TF has different t_0 value stored in myWeight
		 * Time in binding site=t_0
		 * 
		 * Case 2: each TF has a different weight and the same t_0 value
		 * time in binding site = t_0*exp(-2*weight)
		 */
		
		if(params.contains("e_star")){
			return myWeight;
		}
		
		if(myWeight<=0)
			return timeInBackground();
		
		//double avgExponent=Math.exp(2);
		
		double t_0=0;
		if(t_0s.containsKey(myName.toUpperCase())){
			t_0=t_0s.get(myName.toUpperCase());
		}else{
			t_0=Math.pow(10, -6);//((t_r/slides())/avgExponent);
		}
		double beta=params.get("beta");
		return t_0*Math.exp(beta*myWeight);
	}
	
	public double timeInBackground(){
		
		/*
		 * Two strategies to calculate this:
		 * 
		 * Case 1: each TF has its own t_0 in myWeight
		 * time in background=t_0*Math.exp(e_star)
		 * 
		 * Case 2: each TF has a different weight and the same t_0 value
		 * time in background=t_r
		 */
		if(params.contains("e_star")){
			return myWeight*Math.exp(-params.get("e_star"));
		}
		double t_r=params.get("t_r");
		return t_r/slides();
	}
	
	public double freeTFs(int usedUpTfs){
		double f=params.get("f");
		double TF_nonspecific=f*(myConcentration-usedUpTfs);
		return myConcentration-usedUpTfs-TF_nonspecific;
	}
	
	public double visitsToOrigin(double slideWindow){
		//System.out.println("*******Slide Window"+slideWindow+", "+slides());
		return slides()/slideWindow;
	}
	
	public double propensityOfDissociation(){
		return propensityOfDissociation(1, 1);
	}
	
	public double propensityOfDissociation(double a, double b){
		return a*Math.exp(-b*myWeight);
	}
	
	public double backgroundPropensityOfDissociation(double a, double b){
		return a*Math.exp(-b);
	}

	public double getConcentration() {
		return myConcentration;
	}

	public double slides() {
		double s_l_obs=params.get("s_l_obs");
		return (s_l_obs)*(s_l_obs)/2.0;
	}

	public int getDistanceTo(TranscriptionFactor other) {
		return Math.abs(myStart-other.myStart);	
	}
	
	public int getDistanceBetween(TranscriptionFactor other){
		return Math.max(other.myStart-myStop, myStart-other.myStop);
	}

	public double getCoveringTo(TranscriptionFactor other) {
		return other.myStop-myStart+1;
	}

	public int getSize() {
		
		return myStop-myStart+1;
	}

	public String getName() {
		
		return myName;
	}

	public int getStart() {
		
		return myStart;
	}

	public int getEnd() {
		
		return myStop;
	}

	
		
	
	
	
}
