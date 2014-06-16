package beast.math.distributions;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.app.ca.CAPanel;
import beast.app.ca.CladeAgeProbabilities;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.math.distributions.ParametricDistribution;


@Description("Distribution based on fossil information")
public class FossilCalibration extends ParametricDistribution {
	
	
    public static enum CladeAgeMethod {
    	empirical("empirical CladeAge"), 
    	fitted("fitted CladeAge"), 
    	fitted_star("fitted CladeAge*"), 
    	//fitted_rapid("fitted CladeAge rapid"), 
    	standardLogNormal("Lognormal"), 
    	standardGamma("Gamma"), 
    	standardExponential("Exponential"), 
    	standardNormal("Normal");

    	CladeAgeMethod(final String name) {
            this.ename = name;
        }

        public String toString() {
            return ename;
        }

        private final String ename;
    }

	
	public Input<CladeAgeMethod> cladeAgeMethodInput = new Input<FossilCalibration.CladeAgeMethod>("method", 
	"	1.) "+CladeAgeMethod.empirical+" calculates probabilities for 100 ages and returns EmpiricalCladeAgeDistribution " +
	"(this is the only method that allows a sampling gap) " +
	"" +
	"2.) "+CladeAgeMethod.fitted+" calculates probabilities for 100 ages and returns FittedCladeAgeDistribution " +
	"" +
	"3.) "+CladeAgeMethod.fitted_star+" calculates probabilities for 100 ages and returns a FittedCladeAgeDistribution based on only 4 ages " +
	"(this is only there to visualize the fit of the next method, run_fitted_cladeage_rapid) " +
	"" +
	//"4.) "+CladeAgeMethod.fitted_rapid+" calculates probabilities for 4 ages only, and also returns a FittedCladeAgeDistribution based on these 4 ages " +
	//"(this is supposed to be used repeatedly during the MCMC run with updated values for lambda and mu) " +
	//"" +
	"4.) "+CladeAgeMethod.standardLogNormal+", "+CladeAgeMethod.standardGamma+", "+CladeAgeMethod.standardExponential+", "+CladeAgeMethod.standardNormal+" standard calculates probabilities for 100 ages and returns LogNormalImpl, GammaDistributionImpl, ExponentialDistributionImpl, or NormalDistributionImpl, depending on its arguments."
	, CladeAgeMethod.empirical, CladeAgeMethod.values());

	
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
		
		//TODO: deal with offset
	}
    
	@Override
	public Distribution getDistribution() {
		if (m_dist == null) {
			updateEmpiricalCladeAgeDistribution();
		}
		return m_dist;
	}

	private void updateEmpiricalCladeAgeDistribution() {
		Log.info.println("Updating FossilCalibration " + getID());
		CladeAgeProbabilities probs = new CladeAgeProbabilities();
								
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
			switch (cladeAgeMethodInput.get()) {
			case empirical:
			m_dist = probs.run_empirical_cladeage(
					minOccuranceAge, maxOccuranceAge,
					minDivRate, maxDivRate,
					minTurnoverRate, maxTurnoverRate,
					minSamplingRate, maxSamplingRate,
					minSamplingGap, maxSamplingGap,
					null);
			break;
			case fitted:
				m_dist = probs.run_fitted_cladeage(minOccuranceAge, maxOccuranceAge, minDivRate, maxDivRate, minTurnoverRate, maxTurnoverRate, minSamplingRate, maxSamplingRate, null);
			break;
			case fitted_star:
				m_dist = probs.run_fitted_cladeage_star(minOccuranceAge, maxOccuranceAge, minDivRate, maxDivRate, minTurnoverRate, maxTurnoverRate, minSamplingRate, maxSamplingRate, null);
			break;
			//case fitted_rapid:
			//	m_dist = probs.run_fitted_cladeage_rapid(minOccuranceAge, maxOccuranceAge, minDivRate, maxDivRate, minTurnoverRate, maxTurnoverRate, minSamplingRate, maxSamplingRate);
			//break;
			case standardLogNormal:
				m_dist = probs.run_standard(minOccuranceAge, maxOccuranceAge, minDivRate, maxDivRate, minTurnoverRate, maxTurnoverRate, minSamplingRate, maxSamplingRate, 
						"LogNormal", null);
			break;
			case standardGamma:
				m_dist = probs.run_standard(minOccuranceAge, maxOccuranceAge, minDivRate, maxDivRate, minTurnoverRate, maxTurnoverRate, minSamplingRate, maxSamplingRate, 
						"Gamma", null);
				break;
			case standardExponential:
				m_dist = probs.run_standard(minOccuranceAge, maxOccuranceAge, minDivRate, maxDivRate, minTurnoverRate, maxTurnoverRate, minSamplingRate, maxSamplingRate, 
						"Exponential", null);
				break;
			case standardNormal:
				m_dist = probs.run_standard(minOccuranceAge, maxOccuranceAge, minDivRate, maxDivRate, minTurnoverRate, maxTurnoverRate, minSamplingRate, maxSamplingRate, 
						"Normal", null);
			break;
			}				
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

//	public String toString() {
//		String str = "";
//		try {
//			for (Input input : listInputs()) {
//				str += input.getName() +" = " + input.get() + "\n";
//			}
//		} catch (IllegalArgumentException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IllegalAccessException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		return str;
//	}
}
