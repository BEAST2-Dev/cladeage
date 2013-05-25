package beast.math.distributions;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("Prior based on information from the fossil record")
public class FossilPrior extends MRCAPrior {
	public Input<FossilCalibration> callibrationDistr = new Input<FossilCalibration>("fossilDistr", "", Validate.REQUIRED);

	public FossilPrior() {
		m_distInput.setRule(Validate.OPTIONAL);
	}
	
	public void initAndValidate() throws Exception {
		// this only makes sense as a prior on the parent of a clade,
		// so tipsonly and useOriginate will be set to reflect this.
		m_bOnlyUseTipsInput.setValue(false, this);
		m_bUseOriginateInput.setValue(true, this);
		
		m_distInput.setValue(callibrationDistr.get(), this);
		super.initAndValidate();
		m_distInput.setValue(null, this);
	};
}
