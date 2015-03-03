
package edu.cam.fgrip;


import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EmptyStackException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.TreeMap;


public class State{
		private static Random r=new Random();
		private static int tokensInSystem=0; //keeps track of the total number of tokens that are in this state machine
		
		private List<Token> myTokens; //keeps all tokens that are currently in this state
		private Map<State, Double> myTransitions; //keeps all the states that I can transition to in a single reaction, mapped to the propensity of that reaction
		private MarkovChain mc;
		
		private double myPropensitySum=0; //to speed up the program, I'm maintaining the sum of the propensities of going from this state to any other state
		public int myID; //id number, the bits that are 1 are those that have a TF there at that position
		private State mySink=null; //a sink state is a state is a state that has no transitions out of it.  Tokens go there when the simulation is done.
		private double myMaxTime=0; //maximum amount of time that the simulation can run
		
		private int myTFsTotal;
		public 	class Token{
			public Token(int id){
				tokenID=id;
			}
			
			public String toString(){
				return ""+tokenID+" which is "+recentlyMoved;
			}
			
			
			private PriorityQueue<JumpingEvent> jumpingQueue=new PriorityQueue<JumpingEvent>();
			
			private int tokenID;
			//TODO if we base the simulation just on time and NOT number of steps then we can remove this and make the simulation run faster
			private boolean recentlyMoved=false;		
			private double prevTime=0; //time previous reaction started
			private double time=0; //time previous reaction completed
			private int visitedTFBits=0; //the TFs that were visited at least once.  For instance, the first bit is 1 if the first TF was visited and 0 otherwise.
			private int visitedStateBits=0; //the states that were visited at least once.  For instance, the first bit is 1 if the first state was visited and 0 otherwise.
			
			//the reason these are maps rather than arrays is that 1. we don't know the final size of the arrays and 2. they might be very sparse in complex promoters
			//TODO I might make the total number of states in the system a static variable and then make these into arrays, but I haven't decided yet if that is better or worse
			private Map<Integer, Double> timesOfFirstTFOccupancy=new TreeMap<Integer,Double>(); //TFs that were visited mapped to the time that they were first visited
			private Map<Integer, Double> timesAtEachConfig=new TreeMap<Integer, Double>(); //The states that were visited mapped to the total time that state was occupied
			private Map<Integer, Double> timeOfFirstConfigOccupancy=new TreeMap<Integer, Double>(); //The states that were visited mapped to the time that each was first occupied
			private TreeMap<Double, Integer> completePath=new TreeMap<Double, Integer>();
		}
		
		public String toString(){
			String s="";
			s=s+myID+": \n";
			for(State s2: myTransitions.keySet()){
				s=s+"->"+s2.myID+": "+myTransitions.get(s2)+"\n";
			}
			return s;
		}
		
		
		public void setMarkovChain(MarkovChain m){
			mc=m;
		}
		public State getState(int id){
			return mc.getState(id);
		}
		
		public int numTransitions(){
			return myTransitions.size();
		}
		public State(int numberTokens, int stateID, int numTFs){
			myTFsTotal=numTFs;
			myID=stateID;
			myTokens=new ArrayList<Token>();
			for(int i=0; i<numberTokens; i++){
				tokensInSystem++;
				myTokens.add(new Token(tokensInSystem));
			}
			myTransitions=new HashMap<State, Double>();
			
		}
		public State(int numberTokens, int stateID, int numTFs, MarkovChain m){
			myTFsTotal=numTFs;
			myID=stateID;
			myTokens=new ArrayList<Token>();
			for(int i=0; i<numberTokens; i++){
				tokensInSystem++;
				myTokens.add(new Token(tokensInSystem));
			}
			myTransitions=new HashMap<State, Double>();
			mc=m;
		}
		
		public static int extraStates(State first, State second){
			return first.myID & (first.myID ^ second.myID);
		}
		
		
		public static int[][][] getTimeCourseTransitions(double time, double delta, State sink){
			int maxState=0;
			for(Token t: sink.getTokens()){
				for(int i: t.completePath.values()){
					if(i>maxState){
						maxState=i;
					}
				}
			}
			
			int maxBin=(int) (time/delta); 
			
			int[][][] times=new int[maxBin][maxState+1][maxState+1]; 
			
			
			for(Token token: sink.getTokens()){
				int prevBin=-1;
				int prevState=0;
				for(Double t: token.completePath.keySet()){
					int bin=(int)(t/delta);
					int state=token.completePath.get(t);
					int tempPrev=prevState;
					while(bin>=prevBin){
						prevBin++;
						prevState=state;
						if(maxBin>prevBin){
							times[prevBin][tempPrev][state]+=1;
						}
					}
				
				}
			}

			return times;
		}
		
		
		public static int[][] getTimeCourse(double time, double delta, State sink){
			
			int maxState=0;
			for(Token t: sink.getTokens()){
				for(int i: t.completePath.values()){
					if(i>maxState){
						maxState=i;
					}
				}
			}
			
			int maxBin=(int) (time/delta); 
			
			int[][] times=new int[maxState+1][maxBin]; 
			
			for(Token token: sink.getTokens()){
				int prevBin=-1;
				for(Double t: token.completePath.keySet()){
					int bin=(int)(t/delta);
					while(bin>=prevBin){
						prevBin++;
						int state=token.completePath.get(t);
						if(maxBin>prevBin){
							times[state][prevBin]+=1;
						}
					}
				
				}
			}
			
			return times;
			
		}
		
		
		
		
		/*
		 * This method iterates through all the tokens in this state and updates their state.  
		 */
		public void updateState(){
			
			for(int i=myTokens.size()-1; i>=0; i--){
				
				
				double selection=r.nextDouble();
			
				if(!myTokens.get(i).recentlyMoved){
					//System.out.println("Token "+i+" is not recently moved in State:"+myID);
					//call Gillespie
					State nextState=this;
					State nextStateGillespie = Gillespie.getNextReaction(selection*myPropensitySum, myTransitions); 
					//update time
					Token t=myTokens.remove(i);
					t.prevTime=t.time;
					
					//get next time according to Gillespie
					double nextTimeGillespie=t.time+Gillespie.computeNextReactionTime(myPropensitySum, r);
					//System.out.println("     time according to Gillespie:"+nextTimeGillespie);
					
					//get next diffusion reaction according to the jumping Queue
					 JumpingEvent jump=null;
					   //  System.out.println("       A peek at the top jumpingEvent: "+t.jumpingQueue.peek());
					     //iterate through priority queue until you get a diffusion reaction that is possible
					     while(t.jumpingQueue.peek()!=null && t.jumpingQueue.peek().time<nextTimeGillespie){
					    	// System.out.println("                        "+t.jumpingQueue.peek());
					    	 JumpingEvent jumpAttempt=t.jumpingQueue.poll();
					    	 if(jumpAttempt.canHappen(this) && jumpAttempt.time>t.prevTime){
					    	//	 System.out.println("     can happen");
					    		 jump=jumpAttempt;
					    		 break;
					    	 }else{
					    		 if(jumpAttempt.time>t.prevTime){
					    	//	 System.out.println("     adding new event&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
				             //    System.out.println(jumpAttempt.time);
					    		//for each impossible diffusion reaction, re-calculate diffusion and put into the queue
				                 jumpAttempt=jumpAttempt.nextEvent(jumpAttempt.time);
				                 if(jumpAttempt!=null){
				                 t.jumpingQueue.add(jumpAttempt);	
				                 }
					    		 }
					    	 }
					     }
					     			
					//if Gillespie happens first,
					     if(jump==null){
					    	// System.out.println("                               Gillespie event happened first");
					//calculate nextState the standard way
					    	 nextState=nextStateGillespie;
					    	// System.out.println("      next State gillespie: "+nextState.myID);
					    	 t.time=nextTimeGillespie;
					    	 
					    	 //if next reaction is dissociation, add jumping event
					    	 if(isDissociation(this, nextState)){
					    		// System.out.println("it is a dissociation event "+this.myID+" to "+nextState.myID);
					    		 JumpingEvent jp=JumpingEvent.generateNewJumpingEvent(this, nextState, t.time);
					    		// System.out.println(" so it can diffuse "+jp);
					    		 if(jp!=null){
					    		//	 System.out.println("*****************************************");
					    		 t.jumpingQueue.add(jp);
					    		 }
					    	 }
					    	 
					     }else{
					    	// System.out.println("jumping event happened first");
					    	 
					//if jumping Queue happens first,
					//grab nextState from the Jumping event
					     nextState=jump.nextState(this);
					    // System.out.println("Jump description: "+jump+" next State "+nextState.myID+"  at time "+jump.time);
					     t.time=jump.time;
					     }
					
					
					//save data
					
					
					//System.out.println("current state is "+ myID);
					//System.out.println("next reaction will be to "+nextState.myID);
					
					//update time in state
					if(!t.timesAtEachConfig.containsKey(myID)){
						t.timesAtEachConfig.put(myID, 0.0);
					}
					
					//update time at this config
					//System.out.println("state: "+myID+" prevTime here:"+t.timesAtEachConfig.get(myID)+" current time:"+t.time+" prevTime:"+t.prevTime);
					if((t.time-t.prevTime)<0){
						System.out.println("negative time? What? "+t.time+" "+t.prevTime);
						throw new EmptyStackException();
					}
					t.timesAtEachConfig.put(myID, t.timesAtEachConfig.get(myID)+(t.time-t.prevTime)); 
					
					
					
					//add to Path:
					t.completePath.put(t.time, nextState.myID);
					
					//add to next state
					if(mySink!=null && t.time>myMaxTime){
						//***added recently
						
						//update time at this config
						t.timesAtEachConfig.put(myID, t.timesAtEachConfig.get(myID)-(t.time-myMaxTime)); 
                   //     System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"+myMaxTime);
						//***
						mySink.addToken(t, myID);
					}else{
						nextState.addToken(t, myID);
					}
				
				}else{
					myTokens.get(i).recentlyMoved=false;
				}
			}
			
		}
		
		private boolean isDissociation(State state, State nextState) {
		//	System.out.println("             is dissociation: "+state.myID+" "+nextState.myID);
		//	System.out.println("                 "+extraStates(state, nextState)+" "+extraStates(nextState, state));
		//	System.out.println(extraStates(state, nextState)>0 && extraStates(nextState, state)==0);
			return extraStates(state, nextState)>0 && extraStates(nextState, state)==0;
		}

		/*
		 * adds state to transition map
		 */
		public void addState(State state, double probability) {
			myTransitions.put(state, probability);	
			myPropensitySum+=probability;
		}
		
		/*
		 * adds token to this state and, if new TFs or new states were visited, records the time
		 */
		public void addToken(Token t, int prevState){
			
			if(prevState<myID){
				t.recentlyMoved=true;
			}
			myTokens.add(t);
			
			//this is where I update the "first time seen array"
				//1. find which bit has not been visited yet, if any
				//2. add the current token time to the map 
				//3. update which bits have been visited
			int state=myID;
			int visited=t.visitedTFBits;		
			int newbits=(state&visited)^state;
			Map<Integer, Double> m=t.timesOfFirstTFOccupancy;
			int i=0;
			while(newbits>0){
				if((newbits&1)>0){
			
					m.put(i, t.time);
				}
				
				newbits = newbits >>> 1;
				i++;
				
			}
			t.visitedTFBits=visited|state;
			
			//do the same for the first config occupancy map
			int visitedStates=t.visitedStateBits;
			int newStateBit=1<<myID;
			if( ((newStateBit|visitedStates)-visitedStates)!=0 ){
				t.timeOfFirstConfigOccupancy.put(myID, t.time);
			}
			t.visitedStateBits=newStateBit|visitedStates;
				
		}
		
		public void printTokenData(){
			for(Token t: myTokens){
				for(Integer i: t.timesOfFirstTFOccupancy.keySet()){
					System.out.print(t.timesOfFirstTFOccupancy.get(i)+"\t");
				}
				System.out.println();
			}
		}
		
		public static String getFirstOccupancyString(List<State> states){
			
			if(states==null)
				return "states is null";
			
			List<Token> allTokens=new ArrayList<Token>();
			for (State s: states){
				allTokens.addAll(s.myTokens);
			}
			//TODO writing to a file should be extracted to a separate method instead of repeated thrice, but its a little more difficult because you cycle over a different method of t each time, so I'm ignoring this for now
			//save bitTimes (of first reach)
			String time="Occupancy<br />";
			
				for(Token t: allTokens){
					
					int j=0;
					for(Integer i: t.timesOfFirstTFOccupancy.keySet()){
						while(i>j){
							time=time+"\t";
							j++;
						}
						time=time+t.timesOfFirstTFOccupancy.get(i)+"\t";
						j=i+1;
					}
					time=time+"<br />";
					
					j=0;
						
					
				}
				
			return time;
			
		}
		
public static String getFirstStateOccupancyString(List<State> states){
			
			if(states==null)
				return "states is null";
			
			List<Token> allTokens=new ArrayList<Token>();
			for (State s: states){
				allTokens.addAll(s.myTokens);
			}
			//TODO writing to a file should be extracted to a separate method instead of repeated thrice, but its a little more difficult because you cycle over a different method of t each time, so I'm ignoring this for now
			//save bitTimes (of first reach)
			String time="State Occupancy<br />";
			
				for(Token t: allTokens){
					
					int j=0;
					for(Integer i: t.timeOfFirstConfigOccupancy.keySet()){
						while(i>j){
							time=time+"\t";
							j++;
						}
						time=time+t.timeOfFirstConfigOccupancy.get(i)+"\t";
						j=i+1;
					}
					time=time+"<br />";
					
					j=0;
						
					
				}
				
			return time;
			
		}
		
		public static void saveTokenData(List<State> states, String basename){
			
			List<Token> allTokens=new ArrayList<Token>();
			int tfsTotal=0;
			for (State s: states){
				allTokens.addAll(s.myTokens);
				if(s.myTFsTotal>tfsTotal){
					tfsTotal=s.myTFsTotal;
					
				}
			}
			//TODO writing to a file should be extracted to a separate method instead of repeated thrice, but its a little more difficult because you cycle over a different method of t each time, so I'm ignoring this for now
			//save bitTimes (of first reach)
			try {
				FileWriter outFirstOccurance=new FileWriter(basename+"_firstoccurance");
				FileWriter firstStateOccurance=new FileWriter(basename+"_firstStateOccurance");
				FileWriter outOccupancy=new FileWriter(basename+"_occupancy");
				
			
				for(Token t: allTokens){
					
					int j=0;
					for(Integer i: t.timesOfFirstTFOccupancy.keySet()){
						while(i>j){
							outFirstOccurance.write("\t");
							j++;
						}
						outFirstOccurance.write(t.timesOfFirstTFOccupancy.get(i)+"\t");
						j=i+1;
					}
					outFirstOccurance.write("\n");
					
					j=0;
					
					for(Integer i: t.timeOfFirstConfigOccupancy.keySet()){
							
						while(i>j){
							
							firstStateOccurance.write("\t");
							j++;
						}
						firstStateOccurance.write(t.timeOfFirstConfigOccupancy.get(i)+"\t");
						j=i+1;
						
					}
					firstStateOccurance.write("\n");
					
					j=0;
					for(Integer i: t.timesAtEachConfig.keySet()){
						while(i>j){
							outOccupancy.write("0\t");
							j++;
						}
						outOccupancy.write(t.timesAtEachConfig.get(i)+"\t");
						j=i+1;
					}
					while(j<=Math.pow(2, MarkovChainGenerator.numTFs)){
						outOccupancy.write("0\t");
						j++;
					}
					outOccupancy.write("\n");
					
					
				}
				
				outFirstOccurance.close();
				firstStateOccurance.close();
				outOccupancy.close();
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}		
		
		
		public List<Token> getTokens(){
			return myTokens;
		}
	
		public int countTokens(){
			return myTokens.size();
		}

		public int getID() {
			return myID;
		}

		public void addSink(State sink, double time) {
			mySink=sink;
			myMaxTime=time;
		}


		public Map<State, Double> removeTransitions() {
			Map<State, Double> temp=myTransitions;
			myTransitions=new TreeMap<State, Double>();
			return temp;
		}
		
		public void addTransitions(Map<State, Double> t) {
			myTransitions=t;
		}
		
	}
