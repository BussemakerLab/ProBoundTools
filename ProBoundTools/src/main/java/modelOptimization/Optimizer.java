package modelOptimization;

import java.util.Formatter;

import org.json.JSONObject;

import modelComponents.ModelComponent;
import proBoundTools.Array;
import proBoundTools.Misc;

abstract public class Optimizer {
	
	
	public Optimizer() {
	}
	
	abstract public boolean optimize(LossFunction l) throws Exception;  //Optimizes a loss function, returns true if converged
	
	abstract public JSONObject saveToJSON_settings();
	
	protected void printStep(int iterations, int calls, double likelihood, 
			double distance, double ... params) {
		Formatter fmt = new Formatter();
		
		System.out.printf("   %7d      %7d   ", iterations, calls);
		fmt.format("%18.18s   %18.18s", String.format("%10.15f", likelihood), 
				String.format("%10.15f", distance));
		System.out.print(fmt);
		fmt = new Formatter();
		for (int i=0; i<params.length; i++) {
			fmt.format("   %18.18s", String.format("%10.15f", params[i]));
		}
		System.out.print(fmt+"\n");
		fmt.close();
	}
	

	
	public class Optimizer_MaxIteractionException extends Exception
	{
		private static final long serialVersionUID  = 5215251214L;
		// Parameterless Constructor
		public Optimizer_MaxIteractionException() {}

		// Constructor that accepts a message
		public Optimizer_MaxIteractionException(String message)
		{
			super("Maximum Optimizer Iteration Exceeded: "+message);
		}
		
	}
	
	public class Optimizer_LineSearchException extends Exception
	{
		private static final long serialVersionUID  = 351373532623L;
		// Parameterless Constructor
		public Optimizer_LineSearchException() {}

		// Constructor that accepts a message
		public Optimizer_LineSearchException(String message)
		{
			super("Line search failure: "+message);
		}
	}
	
}
