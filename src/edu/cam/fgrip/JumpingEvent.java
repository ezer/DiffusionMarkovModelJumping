package edu.cam.fgrip;

import java.util.Random;


public class JumpingEvent implements Comparable<JumpingEvent>{

	//remember time is the length of time for a diffusion to happen, not an absolute time in the simulation
	public double time;
    private int myTF;
	private static Random r=new Random();
	
	

	public static JumpingEvent generateNewJumpingEvent(State first, State second, double cellularTime) {
		//calculate which dissociation happened
		int dissociation=State.extraStates(first, second);
		//look up next event, if any, from the lookup table
		return nextEvent(dissociation, cellularTime);
	}
	
	public JumpingEvent(double t, int tf) {
		myTF=tf;
		time=t;
	}

	public boolean canHappen(State current) {
		//look at whether state has tf position open
		
		return ((int)(Math.pow(2, myTF)) & current.myID)==0;
	}
    public String toString(){
    	return "Jump "+time+" "+myTF;
    }
	
	private static JumpingEvent nextEvent(int dissociation, double cellularTime){
		//System.out.println("jumping table 1:"+dissociation+" "+(int) Math.floor((Math.log(dissociation)+0.1)/Math.log(2)));
		//System.out.println("jumping table 2:"+r.nextInt(MarkovChainGenerator.numRows));
		JumpingEvent jp=MarkovChainGenerator.jumpingTable[(int) Math.floor((Math.log(dissociation)+0.1)/Math.log(2))][r.nextInt(MarkovChainGenerator.numRows)]; 
		if(jp==null){
			return null;
		}
		//System.out.println("Selected Table Event: "+jp);
		jp.time=jp.time+cellularTime; 
		//System.out.println("Next Event: "+jp);
		return jp;
	}
	public JumpingEvent nextEvent(double cellularTime) {
		// this gives you the next event assuming that the current diffusion event fails
		return nextEvent((int)(Math.pow(2, myTF)), cellularTime);
	}

	public State nextState(State currentState) {
		//get ID of next state
		int id=currentState.myID | (int)(Math.pow(2, myTF));
		//System.out.println("next state function outputs: "+currentState.myID+" to "+id+" via "+myTF);
		//get next state reference, via State
		return currentState.getState(id);
	}


	@Override
	public int compareTo(JumpingEvent o) {
		
		return (int)(100*(time-o.time));
	}

}
