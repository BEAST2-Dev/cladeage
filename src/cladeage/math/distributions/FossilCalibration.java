package cladeage.math.distributions;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import cladeage.app.ca.CAPanel;
import cladeage.app.ca.CladeAgeProbabilities;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.inference.distribution.ParametricDistribution;


@Description("Distribution based the age of a single fossil")
public class FossilCalibration extends ParametricDistribution implements Serializable {

	private static final long serialVersionUID = 1L;

	public static enum CladeAgeMethod {
    	standard("standard CladeAge"); 

    	CladeAgeMethod(final String name) {
            this.ename = name;
        }

        public String toString() {
            return ename;
        }

        private final String ename;
    }

	
	public Input<CladeAgeMethod> cladeAgeMethodInput = new Input<FossilCalibration.CladeAgeMethod>("method", 
		CladeAgeMethod.standard+" calculates probabilities for 100 ages and returns CladeAgeDistribution ",
		CladeAgeMethod.standard,
		CladeAgeMethod.values());

	
	public Input<RealParameter> minOccuranceAgeInput = new Input<RealParameter>("minOccuranceAge", CAPanel.OCCURRENCE_AGE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxOccuranceAgeInput = new Input<RealParameter>("maxOccuranceAge", CAPanel.OCCURRENCE_AGE_HELP, Validate.REQUIRED);

	public Input<RealParameter> minDivRateInput = new Input<RealParameter>("minDivRate", CAPanel.DIV_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxDivRateInput = new Input<RealParameter>("maxDivRate", CAPanel.DIV_RATE_HELP, Validate.REQUIRED);
	
	public Input<RealParameter> minTurnoverRateInput = new Input<RealParameter>("minTurnoverRate", CAPanel.TURNOVER_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxTurnoverRateInput = new Input<RealParameter>("maxTurnoverRate", CAPanel.TURNOVER_RATE_HELP, Validate.REQUIRED);
	
	public Input<RealParameter> minSamplingRateInput = new Input<RealParameter>("minSamplingRate", CAPanel.SAMPLING_RATE_HELP, Validate.REQUIRED);
	public Input<RealParameter> maxSamplingRateInput = new Input<RealParameter>("maxSamplingRate", CAPanel.SAMPLING_RATE_HELP, Validate.REQUIRED);
	
	
    ContinuousDistribution m_dist = null;

	private double minOccuranceAge;
	private double minDivRate;
	private double minTurnoverRate;
	private double minSamplingRate;

	private double maxOccuranceAge;
	private double maxDivRate;
	private double maxTurnoverRate;
	private double maxSamplingRate;

	@Override
	public void initAndValidate() {
		minOccuranceAge = minOccuranceAgeInput.get().getValue();
		minDivRate = minDivRateInput.get().getValue();
		minTurnoverRate = minTurnoverRateInput.get().getValue();
		minSamplingRate = minSamplingRateInput.get().getValue();

		maxOccuranceAge = maxOccuranceAgeInput.get().getValue();
		maxDivRate = maxDivRateInput.get().getValue();
		maxTurnoverRate = maxTurnoverRateInput.get().getValue();
		maxSamplingRate = maxSamplingRateInput.get().getValue();
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
			case standard:
			m_dist = probs.run_standard_cladeage(
					minOccuranceAge, maxOccuranceAge,
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
