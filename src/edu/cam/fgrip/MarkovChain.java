package edu.cam.fgrip;

import java.util.*;

public class MarkovChain {

	List<State> stateList= new ArrayList<State>();
	List<TranscriptionFactor> tfs;
	
	public MarkovChain(List<State> stateGraph) {
		stateList=stateGraph;
	}
	
	public void addTFs(List<TranscriptionFactor> t){
		tfs=t;
	}
	public void addState(State state){
		state.setMarkovChain(this);
		stateList.add(state);
	}
	
	public State getState(int id){
		return stateList.get(id);
	}
	
	public int[][] runTimeCourse(double time, double delta, String fileName){
		List<State> sink=runSimulationUntil(time, fileName);
		return State.getTimeCourse(time, delta, sink.get(0));
		
	}
	
	public int[][][] getTimeCourseTransitions(double time, double delta, String fileName){
		List<State> sink=runSimulationUntil(time, fileName);
		return State.getTimeCourseTransitions(time, delta, sink.get(0));
	}
	
	public int[][] runTimeCourse(double time, double delta, String fileName, List<State> sink){
		
		return State.getTimeCourse(time, delta, sink.get(0));
		
	}
	
	public int[][][] getTimeCourseTransitions(double time, double delta, String fileName, List<State> sink){
		
		return State.getTimeCourseTransitions(time, delta, sink.get(0));
	}
	
	public List<State> runSimulationUntil(double time, String str){
		State sink=new State(0, -1, 0, this);
		int totalTokens=0;
		for(State s: stateList){
			s.addSink(sink, time);
			totalTokens+=s.countTokens();
		}
		
		     List<State> sinks=runSimulationUntilSink(sink, totalTokens, str);
		
			
		return sinks;
	
	}
	
	
	public ArrayList<String> getJSONs(int time, int lapse, String filename){
		List<State> sink=runSimulationUntil(time, filename);
		int[][][] test2=getTimeCourseTransitions(time, lapse, filename, sink);
		JSONtest json=new JSONtest();
		ArrayList<String> jsons= json.generateGraphJSON(test2, tfs);
		test2=null;
		int[][] test1=runTimeCourse(time, lapse, filename, sink);
		jsons.add(json.generateAreaJSON(test1, lapse));
		String out=State.getFirstOccupancyString(sink);
		//System.out.println(out.substring(2, 100));
		jsons.add(out);
		return jsons;
	}
	
	

	
	public void runSimulationUntilStateReached(int stateID, String filebasename){
		State sink=null;
		int totalTokens=0;
		for(State s: stateList){
			if(s.getID()==stateID){
				sink=s;
			}
			totalTokens+=s.countTokens();
		}
		
		Map<State, Double> temp=sink.removeTransitions();
		
		runSimulationUntilSink(sink, totalTokens, filebasename);
		
		sink.addTransitions(temp);
	}

	
	public void runSimulationUntilTFReached(int TF, String filebasename){
		State sink=new State(0, -1, 0, this);
		
		//find all states that have that TF
		int tfbit=1<<TF;
		int totalTokens=0;
		Map<State, Map<State, Double>> temp=new HashMap<State, Map<State, Double>>();
		for(State s: stateList){
			totalTokens+=s.countTokens();
			if((tfbit&s.getID())!=0){
				temp.put(s, s.removeTransitions());
				s.addState(sink, 1);
			}
		}
		
		runSimulationUntilSink(sink, totalTokens, filebasename);
		
		for(State s: temp.keySet()){
			s.addTransitions(temp.get(s));
		}
	}
	
	private List<State> runSimulationUntilSink(State sink, int totalTokens, String filebasename) {
		
		while(sink.countTokens()<totalTokens){
			
			for(State s: stateList){		
				s.updateState();
			}
		
		}
		
			sink.printTokenData();
			ArrayList<State> sinkList=new ArrayList<State>();
			sinkList.add(sink);
		
			
		    State.saveTokenData(sinkList, filebasename);
		    return sinkList;
	}

	public void writeProportionFiles(String occupancyFile, List<TranscriptionFactor> tfs){
		/*
		 * Time in each configuration = product of 1/s_l for each TF that is bound
		 */
		//rightmostBound = Math.min((int)(s_l_real/2), getRightBound(oldID, tfID, tfs));
		//leftmostBound = Math.min((int)(s_l_real/2), getLeftBound(oldID, tfID, tfs));
		
		
		
	}
		
	
}
