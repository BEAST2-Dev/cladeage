package beast.math.distributions;

import javax.swing.JProgressBar;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;

import beast.app.ca.CAPanel;
import beast.app.ca.CladeAgeProbabilities;
import beast.core.Description;
import beast.core.Input;

@Description("Distribution based on fossil information")
public class FossilCalibration extends ParametricDistribution {
	public Input<Double> minOccuranceAgeInput = new Input<Double>("minOccuranceAge", CAPanel.OCCURRENCE_AGE_HELP, 0.0);
	public Input<Double> maxOccuranceAgeInput = new Input<Double>("maxOccuranceAge", CAPanel.OCCURRENCE_AGE_HELP, 0.0);

	public Input<Double> minDivRateInput = new Input<Double>("minDivRate", CAPanel.DIV_RATE_HELP, 0.01);
	public Input<Double> maxDivRateInput = new Input<Double>("maxDivRate", CAPanel.DIV_RATE_HELP, 0.01);
	
	public Input<Double> minTurnoverRateInput = new Input<Double>("minTurnoverRate", CAPanel.TURNOVER_RATE_HELP, 0.1);
	public Input<Double> maxTurnoverRateInput = new Input<Double>("maxTurnoverRate", CAPanel.TURNOVER_RATE_HELP, 0.1);
	
	public Input<Double> minSamplingRateInput = new Input<Double>("minSamplingRate", CAPanel.SAMPLING_RATE_HELP, 0.01);
	public Input<Double> maxSamplingRateInput = new Input<Double>("maxSamplingRate", CAPanel.SAMPLING_RATE_HELP, 0.01);
	
	public Input<Double> minSamplingGapInput = new Input<Double>("minSamplingGap", CAPanel.SAMPLING_GAP_HELP, 0.0);
	public Input<Double> maxSamplingGapInput = new Input<Double>("maxSamplingGap", CAPanel.SAMPLING_GAP_HELP, 0.0);
	
	public Input<Integer> NumberOfTreeSimulationsInput = new Input<Integer>("numberOfTreeSimulations", CAPanel.NR_SIMULATIONS_HELP, 1000);
	public Input<Integer> MaxNrOfBranchesInput = new Input<Integer>("maxNrOfBranches", CAPanel.MAX_NR_TREES_HELP, 100000);
	public Input<Integer> SamplingReplicatesPerTreeInput = new Input<Integer>("samplingReplicatesPerTree", CAPanel.REPS_PER_TREE_HELP, 10);
	
    ContinuousDistribution m_dist = null;

	private double minOccuranceAge;
	private double minDivRate;
	private double minTurnoverRate;
	private double minSamplingRate;
	private double minSamplingGap;

	private double maxOccuranceAge;
	private double maxDivRate;
	private double maxTurnoverRate;
	private double maxSamplingRate;
	private double maxSamplingGap;
	private int NumberOfTreeSimulations;
	private int MaxNrOfBranches;
	private int SamplingReplicatesPerTree;


	@Override
	public void initAndValidate() throws Exception {
		minOccuranceAge = minOccuranceAgeInput.get();
		minDivRate = minDivRateInput.get();
		minTurnoverRate = minTurnoverRateInput.get();
		minSamplingRate = minSamplingRateInput.get();
		minSamplingGap = minSamplingGapInput.get();

		maxOccuranceAge = maxOccuranceAgeInput.get();
		maxDivRate = maxDivRateInput.get();
		maxTurnoverRate = maxTurnoverRateInput.get();
		maxSamplingRate = maxSamplingRateInput.get();
		maxSamplingGap = maxSamplingGapInput.get();

		NumberOfTreeSimulations = NumberOfTreeSimulationsInput.get();
		MaxNrOfBranches = MaxNrOfBranchesInput.get();
		SamplingReplicatesPerTree = SamplingReplicatesPerTreeInput.get();
		
		//TODO: deal with offset
	}
    
	@Override
	public Distribution getDistribution() {
		if (m_dist == null) {
			updateDistribution();
		}
		return m_dist;
	}



	private void updateDistribution() {
		CladeAgeProbabilities probs = new CladeAgeProbabilities();
		probs.bd_simulate(
				minOccuranceAge, maxOccuranceAge,
				minDivRate, maxDivRate,
				minTurnoverRate, maxTurnoverRate,
				minSamplingRate, maxSamplingRate,
				minSamplingGap, maxSamplingGap,
				NumberOfTreeSimulations, MaxNrOfBranches, SamplingReplicatesPerTree, new JProgressBar());
		m_dist = probs.fitExponential();
	}

}
