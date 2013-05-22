package beast.math.distributions;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("Prior based on information from the fossil record")
public class FossilPrior extends MRCAPrior {
	public Input<FossilCalibration> callibrationDistr = new Input<FossilCalibration>("fossilDistr", "", Validate.REQUIRED);

	public FossilPrior() {
		m_distInput.setRule(Validate.FORBIDDEN);
	}
	
	public void initAndValidate() throws Exception {
		m_distInput.setValue(callibrationDistr.get(), this);
		super.initAndValidate();
		m_distInput.setValue(null, this);
	};
}
