package edu.cam.fgrip;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


public class MarkovChainGenerator {

	private List<TranscriptionFactor> myTFs;
	private ProbabilityModel myModel;
	private double totalMolecules=0;
	private Parameters myParams;
	public static int numTFs;
	public static JumpingEvent[][] jumpingTable;
    public static int numRows=10000;
	
	public void inputJumpingEvents(String file){
		try {
			Scanner in=new Scanner(new File(file));
			String row=in.nextLine();
			String[] vals = row.split("\\s+");
			//System.out.println("vals length: "+vals.length);
			jumpingTable=new JumpingEvent[vals.length/2][numRows];
			for(int i=1; i<vals.length; i+=2){
				if(Double.parseDouble(vals[i])<0 || Double.parseDouble(vals[i])>99){
					jumpingTable[i/2][1]=null;
				}else{
				jumpingTable[i/2][1]=new JumpingEvent(Double.parseDouble(vals[i]), (int) (Double.parseDouble(vals[i+1])-1));
			    }
			}
			int index=1;
			while(index<numRows){
				row=in.nextLine();
				vals = row.split("\\s+");
				for(int i=1; i<vals.length; i+=2){
					if(Double.parseDouble(vals[i])<0 || Double.parseDouble(vals[i])>99){
						jumpingTable[i/2][index]=null;
					}else{
					jumpingTable[i/2][index]=new JumpingEvent(Double.parseDouble(vals[i]), (int) (Double.parseDouble(vals[i+1])-1));
				    }
				}
			index++;	
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//System.out.println("out"+jumpingTable);
	}
	
	public MarkovChain buildMarkovChainFromStream(String stream, int numTokens, ProbabilityModel m, String paramFile){
		myModel=m;
		myParams=new Parameters(paramFile);
		parseTFs(stream, true);
		List<State> stateGraph=linkState(myTFs.size(), numTokens); 
		/*
		for(int i=stateGraph.size()-1; i>=0; i--){
			State s=stateGraph.get(i);
			if(s.numTransitions()==0){
				stateGraph.remove(s);
			}
		}*/
		//TODO I'll probably end up either pushing linkState up or factoring out all the TF-centered stuff
		MarkovChain mc=new MarkovChain(stateGraph);
		mc.addTFs(myTFs);
		return mc;
	}
	
	
	public MarkovChain buildMarkovChain(String file, String jumpFile, int numStates, int numTokens, ProbabilityModel m, String paramFile) {
		inputJumpingEvents(jumpFile);
		
		myModel=m;
		myParams=new Parameters(paramFile);
		numTFs=numStates;
		//System.out.println("params read in");
		parseTFs(file, false);
		//System.out.println("tfs parsed: "+numStates);
		List<State> stateGraph=linkState(numStates, numTokens); 
		//System.out.println("I am printing the markov chain now");
		//for(State s: stateGraph){
		//	System.out.println(s);
		//}
		//System.out.println("The end");
		/*
		for(int i=stateGraph.size()-1; i>=0; i--){
			State s=stateGraph.get(i);
			if(s.numTransitions()==0){
				stateGraph.remove(s);
			}
		}*/
		//TODO I'll probably end up either pushing linkState up or factoring out all the TF-centered stuff
		MarkovChain mc=new MarkovChain(stateGraph);
		for(State s: stateGraph){
			s.setMarkovChain(mc);
		}
		mc.addTFs(myTFs);
		return mc;
	}
	
	//TODO right now all the Tokens start off in the empty state, but I could change that to allow the tokens start off in stateID x
	/*
	 * Builds the connected State network in the Markov chain
	 * 
	 */
	public List<State> linkState(int numStates, int numTokens){
		List<State> transitions=new ArrayList<State>();
		transitions.add(new State(numTokens, 0, numStates));
		for(int i=1; i<Math.pow(2, numStates); i++){
			transitions.add(new State(0, i, numStates));
		}
		
		//add connections:
		for(int i=0; i<Math.pow(2, numStates); i++){
			for(int j=0; j<numStates; j++){	
				int temp=(int) Math.pow(2, j);
				//each state is addressed by its index and the pattern of 1s and 0s corresponds to a TF being present or not
				//if xor of the state with 2^j is not 0, then that means the jth TF is absent. 
				//two transitions must be made: one linking state i to state i+2^j and one going backwards
				
				
				
				if( (i&temp) == 0  ){	
					if( !containsSwitch(i, j)) {
						transitions.get(i).addState(transitions.get((int)(i+Math.pow(2, j))), myModel.forwardProbability(i, j, myTFs, totalMolecules));
						transitions.get((int)(i+Math.pow(2, j))).addState(transitions.get(i), myModel.backwardProbability(i, j, myTFs));
						}
				}	
			}
		}
		addClusterTransitions(numStates, transitions);
		printStateGraph(transitions);
		return transitions;
	}
	
	

	private void printStateGraph(List<State> transitions) {
		//for(State s: transitions){
		//	System.out.println(s);
		//}
	}

	private void addClusterTransitions(int n, List<State> states) {
		Collections.sort(myTFs);
		for(int i=1; i<myTFs.size(); i++){
			//while the TFs are closer than s_l_real
			int j=1;
			while(i>=j && myTFs.get(i-j).getDistanceTo(myTFs.get(i))<myParams.get("s_l_real")){
				//System.out.println("compare "+(i-j)+", "+i);
			
				if(myTFs.get(i-j).hasSameNameAs(myTFs.get(i))){
					//System.out.println("found cluster");
					for(State s: states){
						int id=s.getID();
						int v= id;
						int v2=id;
						v=v >>> (i-j);
						v2 = v2 >>> i;
						if((v&1)==1 && (v2&1)==0){
							
							boolean valid=true;
							for(int k=i-j; k<i; k++){
								v=v >>> 1;
								if((v&1)==1){
									valid=false;
									break;
								}
							}
							if(valid){
								int nextState=id-(int)Math.pow(2, i-j)+(int)Math.pow(2, i);
								//System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^ If you add: "+i+" to state "+(id-(int)Math.pow(2, i-j))+" then..."+((int) Math.log(states.size()+1)+1));
								//System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^ If you add: "+(i-j)+" to state "+(id-(int)Math.pow(2, i-j))+" then...");
								
								if(!containsSwitch(id-(int)Math.pow(2, i-j), i) && !containsSwitch(id-(int)Math.pow(2, i-j), i-j)){
									//System.out.println("you get no switch so link: "+nextState+" and "+id);
									double p= myModel.getSlideProbability(i-j, i, id, myTFs);
									//System.out.println("cluster transitions: "+id+" "+nextState+": "+p);
									states.get(id).addState(states.get(nextState), p);
									states.get(nextState).addState(states.get(id), p);
								}
							}
						}
					}
				}
				j=j+1;
			}
			
			/**
			if(myTFs.get(i-1).hasSameNameAs(myTFs.get(i))){
				System.out.println("found cluster");
				for(State s: states){
					int id=s.getID();
					int v= id;
					v=v >>> (i-1);
					if((v&1)==1 && (v&2)==0){
						int nextState=id+(int)Math.pow(2, i-1);
						double p= myModel.getSlideProbability(i-1, i, myTFs);
						System.out.println("cluster transitions: "+i+" "+nextState+": "+p);
						states.get(i).addState(states.get(nextState), p);
						states.get(nextState).addState(states.get(i), p);
					}
				}
			}
			*/
		}
		
	}

	

	
	public void parseTFs(String input, boolean stream){
		myTFs=new ArrayList<TranscriptionFactor>();
		Map<String, Double> conc=new HashMap<String, Double>();
		try {
			Scanner in;
			if(!stream){
			in=new Scanner(new File(input));
			}else{
			in=new Scanner(input);	
			}
			while(in.hasNext()){
				String name=in.next();
				int start=in.nextInt();
				int stop=in.nextInt();
				double weight=in.nextDouble();
				double concentration=in.nextDouble();
				//double x=in.nextDouble();
				//double y=in.nextDouble();
				//double z=in.nextDouble();
				if(!conc.containsKey(name)){
					conc.put(name, concentration);
				}
				//System.out.println(name+" "+start+" "+stop+" "+weight+" "+concentration);
				//System.out.println("x,y,z:"+x+" "+y+" "+z);
				myTFs.add(new TranscriptionFactor(name, start, stop, weight, concentration, myParams));
			}
			
			Collections.sort(myTFs);
			//System.out.println(myTFs.size());
			for(Double val: conc.values()){
				totalMolecules+=val;
			}
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public TranscriptionFactor get(int i){
		return myTFs.get(i);
	}
	
	
	public boolean containsSwitch(int state, int tf){
		
		int j=tf;
		int i=state;
		
		if((state&(int)Math.pow(2, j))!=0){
			return true;
		}
		
		//if(!(!(j-1>=0 && (i&((int)Math.pow(2, j-1)))!=0 && isSwitch(j-1, j)) && (!(j+1<numStates && (i&((int)Math.pow(2, j+1)))!=0 && isSwitch(j+1, j))))){ 
		//	return true;
		//}
		//moving left
		j=tf;
		while(j-1>=0 && isSwitch(j-1, tf)){
			if( (state & ((int)Math.pow(2, j-1)))!=0 ){
				return true;
			}
			j=j-1;
		}
		//moving right
		j=tf;
		while(j+1<myTFs.size() && isSwitch(tf, j+1)){
			if( (state & ((int)Math.pow(2, j+1)))!=0 ){
				return true;
			}
			j=j+1;
		}
		
		return false;
	}
	
	public boolean isSwitch(int a, int b) {
		if(a>=0 && a<myTFs.size() && b>=0 && b<myTFs.size())
		return get(a).overlapsWith(myTFs.get(b));
		
		return false;
	}

	public List<TranscriptionFactor> getTFs() {
		if(myTFs==null){
			System.out.println("aaaah!");
		}
		return myTFs;
	}
}
