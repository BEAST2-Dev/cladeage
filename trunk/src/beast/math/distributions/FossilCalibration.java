package beast.math.distributions;



import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.app.ca.CAPanel;
import beast.app.ca.CladeAgeProbabilities;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;


@Description("Distribution based on fossil information")
public class FossilCalibration extends ParametricDistribution {
	public Input<RealParameter> minOccuranceAgeInput = new Input<RealParameter>("minOccuranceAge", CAPanel.OCCURRENCE_AGE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxOccuranceAgeInput = new Input<RealParameter>("maxOccuranceAge", CAPanel.OCCURRENCE_AGE_HELP, Validate.REQUIRED);

	public Input<RealParameter> minDivRateInput = new Input<RealParameter>("minDivRate", CAPanel.DIV_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxDivRateInput = new Input<RealParameter>("maxDivRate", CAPanel.DIV_RATE_HELP, Validate.REQUIRED);
	
	public Input<RealParameter> minTurnoverRateInput = new Input<RealParameter>("minTurnoverRate", CAPanel.TURNOVER_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxTurnoverRateInput = new Input<RealParameter>("maxTurnoverRate", CAPanel.TURNOVER_RATE_HELP, Validate.REQUIRED);
	
	public Input<RealParameter> minSamplingRateInput = new Input<RealParameter>("minSamplingRate", CAPanel.SAMPLING_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxSamplingRateInput = new Input<RealParameter>("maxSamplingRate", CAPanel.SAMPLING_RATE_HELP, Validate.REQUIRED);
	
	public Input<RealParameter> minSamplingGapInput = new Input<RealParameter>("minSamplingGap", CAPanel.SAMPLING_GAP_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxSamplingGapInput = new Input<RealParameter>("maxSamplingGap", CAPanel.SAMPLING_GAP_HELP, Validate.REQUIRED);
	
	public Input<Integer> NumberOfTreeSimulationsInput = new Input<Integer>("numberOfTreeSimulations", CAPanel.NR_SIMULATIONS_HELP, 10000);
	public Input<Integer> MaxNrOfBranchesInput = new Input<Integer>("maxNrOfBranches", CAPanel.MAX_NR_TREES_HELP, 100000);
	public Input<Integer> SamplingReplicatesPerTreeInput = new Input<Integer>("samplingReplicatesPerTree", CAPanel.REPS_PER_TREE_HELP, 10);
	
//	public boolean canSetMinOccuranceAge(Object o) throws Exception {
//        final Double value = ((RealParameter) o).getValue(); 
//        if (maxOccuranceAgeInput.get() != null) {
//        	return (value <= maxOccuranceAgeInput.get().getValue());
//        } 
//        return true;
//		
//	}
//	public boolean canSetMaxOccuranceAge(Object o) throws Exception {
//        final Double value = ((RealParameter) o).getValue(); 
//        if (minOccuranceAgeInput.get() != null) {
//        	return (value >= minOccuranceAgeInput.get().getValue());
//        }
//        return true;
//		
//	}
	
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
		minOccuranceAge = minOccuranceAgeInput.get().getValue();
		minDivRate = minDivRateInput.get().getValue();
		minTurnoverRate = minTurnoverRateInput.get().getValue();
		minSamplingRate = minSamplingRateInput.get().getValue();
		minSamplingGap = minSamplingGapInput.get().getValue();

		maxOccuranceAge = maxOccuranceAgeInput.get().getValue();
		maxDivRate = maxDivRateInput.get().getValue();
		maxTurnoverRate = maxTurnoverRateInput.get().getValue();
		maxSamplingRate = maxSamplingRateInput.get().getValue();
		maxSamplingGap = maxSamplingGapInput.get().getValue();

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
				NumberOfTreeSimulations, MaxNrOfBranches, SamplingReplicatesPerTree, null);
		
//		if (maxOccuranceAge == minOccuranceAge) {
//			m_dist = probs.fitExponential(null);
//			double rmsd = probs.getApprox_distribution_rmsd();
//			ContinuousDistribution tmp = probs.fitExpGamma(null);
//			if (probs.getApprox_distribution_rmsd() < rmsd) {
//				m_dist = tmp;
//			}
//		} else {
//			m_dist = probs.fitGamma(null);
//			double rmsd = probs.getApprox_distribution_rmsd();
//			ContinuousDistribution tmp = probs.fitLognormal(null);
//			if (probs.getApprox_distribution_rmsd() < rmsd) {
//				 m_dist = tmp;
//			}
//		}
		try {
			EmpiricalDistribution dist = new EmpiricalDistribution();
			dist.setup(probs.getAges(), probs.getProbabilities(), true);
			m_dist = dist;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

}
