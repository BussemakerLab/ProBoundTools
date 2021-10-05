package modelOptimization;

abstract public class LossFunction {
	
	public int nParameters;      //Number of parameters
	public double[] parameters;  // Current value of parameter. Set by setParameters(double[]);
	public double value;         // Current value of the loss function. Computed by updateValue()
	public double[] gradient;    // Current value of the gradient of the loss function. Computed by updateGradient();
	public boolean isNan;        // Indicates that the evaluation has failed. Indicates that an evaluation has failed.
	public double[] lastReported;// Contains the the last reported parameters.
	
	protected boolean computeVariance = false; //Indicates if the variance of the value and gradient should be computed
	public double valueVariance;               //Variance of the function value
	
	abstract public void lossFunction_setParameters(double[] in);  // Sets the parameters
	abstract public void lossFunction_updateValue();               // Computes the value of the loss function
	abstract public void lossFunction_reportPosition();            // Is called after the optimizer takes as step. Used for logging
	abstract public void lossFunction_nextBatch(int nData);        // Gets a new batch of nData data points to use for fuction/gradient evaluation.
	abstract public long lossFunction_getDataSize();	           // Returns the number of datapoints.
	abstract public int  lossFunction_getBatchSize();	           // Returns the number of datapoints in the batch.
	
}
