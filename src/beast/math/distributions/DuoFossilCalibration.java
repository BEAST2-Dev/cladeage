package beast.math.distributions;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.app.ca.CAPanel;
import beast.app.ca.CladeAgeProbabilities;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.inference.distribution.ParametricDistribution;


@Description("Distribution based on the ages of two fossils")
public class DuoFossilCalibration extends ParametricDistribution implements Serializable {

	private static final long serialVersionUID = 1L;

	public static enum CladeAgeMethod {
    	duo("duo CladeAge"); 

    	CladeAgeMethod(final String name) {
            this.ename = name;
        }

        public String toString() {
            return ename;
        }

        private final String ename;
    }

	
	public Input<CladeAgeMethod> cladeAgeMethodInput = new Input<DuoFossilCalibration.CladeAgeMethod>("method", 
		CladeAgeMethod.duo+" calculates probabilities for 100 ages and returns CladeAgeDistribution.",
		CladeAgeMethod.duo,
		CladeAgeMethod.values());

	
	public Input<RealParameter> minYoungerOccuranceAgeInput = new Input<RealParameter>("minYoungerOccuranceAge", CAPanel.YOUNGER_OCCURRENCE_AGE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxYoungerOccuranceAgeInput = new Input<RealParameter>("maxYoungerOccuranceAge", CAPanel.YOUNGER_OCCURRENCE_AGE_HELP, Validate.REQUIRED);

	public Input<RealParameter> minOlderOccuranceAgeInput = new Input<RealParameter>("minOlderOccuranceAge", CAPanel.OLDER_OCCURRENCE_AGE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxOlderOccuranceAgeInput = new Input<RealParameter>("maxOlderOccuranceAge", CAPanel.OLDER_OCCURRENCE_AGE_HELP, Validate.REQUIRED);

	public Input<RealParameter> weightYoungerOccuranceInput = new Input<RealParameter>("weightYoungerOccurance", CAPanel.WEIGHT_YOUNGER_OCCURRENCE_HELP, Validate.REQUIRED);

	public Input<RealParameter> weightOlderOccuranceInput = new Input<RealParameter>("weightOlderOccurance", CAPanel.WEIGHT_OLDER_OCCURRENCE_HELP, Validate.REQUIRED);

	public Input<RealParameter> minDivRateInput = new Input<RealParameter>("minDivRate", CAPanel.DIV_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxDivRateInput = new Input<RealParameter>("maxDivRate", CAPanel.DIV_RATE_HELP, Validate.REQUIRED);
	
	public Input<RealParameter> minTurnoverRateInput = new Input<RealParameter>("minTurnoverRate", CAPanel.TURNOVER_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxTurnoverRateInput = new Input<RealParameter>("maxTurnoverRate", CAPanel.TURNOVER_RATE_HELP, Validate.REQUIRED);
	
	public Input<RealParameter> minSamplingRateInput = new Input<RealParameter>("minSamplingRate", CAPanel.SAMPLING_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxSamplingRateInput = new Input<RealParameter>("maxSamplingRate", CAPanel.SAMPLING_RATE_HELP, Validate.REQUIRED);
	
	
    ContinuousDistribution m_dist = null;

	private double minYoungerOccuranceAge;
	private double minOlderOccuranceAge;
	private double minDivRate;
	private double minTurnoverRate;
	private double minSamplingRate;

	private double maxYoungerOccuranceAge;
	private double maxOlderOccuranceAge;
	private double maxDivRate;
	private double maxTurnoverRate;
	private double maxSamplingRate;

	private double weightYoungerOccurance;
	private double weightOlderOccurance;

	@Override
	public void initAndValidate() {
		minYoungerOccuranceAge = minYoungerOccuranceAgeInput.get().getValue();
		minOlderOccuranceAge = minOlderOccuranceAgeInput.get().getValue();
		minDivRate = minDivRateInput.get().getValue();
		minTurnoverRate = minTurnoverRateInput.get().getValue();
		minSamplingRate = minSamplingRateInput.get().getValue();

		maxYoungerOccuranceAge = maxYoungerOccuranceAgeInput.get().getValue();
		maxOlderOccuranceAge = maxOlderOccuranceAgeInput.get().getValue();
		maxDivRate = maxDivRateInput.get().getValue();
		maxTurnoverRate = maxTurnoverRateInput.get().getValue();
		maxSamplingRate = maxSamplingRateInput.get().getValue();

		weightYoungerOccurance = weightYoungerOccuranceInput.get().getValue();
		weightOlderOccurance = weightOlderOccuranceInput.get().getValue();
	}
    
	@Override
	public Distribution getDistribution() {
		if (m_dist == null) {
			if (true) {
				updateCladeAgeDistribution();
			} else {
				File serialisedDist = new File(getID() + ".ser");
				if (serialisedDist.exists()) {
					try {
						FileInputStream fis = new FileInputStream(serialisedDist);
						ObjectInputStream in = new ObjectInputStream(fis);
						m_dist = (CladeAgeDistribution) in.readObject();
						in.close();
					} catch (Exception e) {
						e.printStackTrace();
						Log.err.println(e.getMessage());
						updateCladeAgeDistribution();
					}
				} else {
					updateCladeAgeDistribution();
					try {
						FileOutputStream fos = new FileOutputStream(serialisedDist);
						ObjectOutputStream out = new ObjectOutputStream(fos);
						out.writeObject(m_dist);
						out.close();
					} catch (IOException e) {
						e.printStackTrace();
						Log.err.println(e.getMessage());
					}
				}
			}
		}
		return m_dist;
	}

	private void updateCladeAgeDistribution() {
		Log.info.println("Updating FossilCalibration " + getID());
		CladeAgeProbabilities probs = new CladeAgeProbabilities();
								
		try {
			switch (cladeAgeMethodInput.get()) {
			case duo:
			m_dist = probs.run_duo_cladeage(
					minYoungerOccuranceAge, maxYoungerOccuranceAge, weightYoungerOccurance,
					minOlderOccuranceAge, maxOlderOccuranceAge, weightOlderOccurance,
					minDivRate, maxDivRate,
					minTurnoverRate, maxTurnoverRate,
					minSamplingRate, maxSamplingRate,
					null);
			break;
			}				
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
}