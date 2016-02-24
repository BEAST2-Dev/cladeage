package beast.math.distributions;

import java.io.PrintStream;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.math.distributions.MRCAPrior;


@Description("Prior based on information from the fossil record")
@Citation("Michael Matschiner and Remco Bouckaert. A rough guide to CladeAge, 2013")
public class FossilPrior extends MRCAPrior {
	public Input<FossilCalibration> calibrationDistr = new Input<>("fossilDistr", "", Validate.REQUIRED);

	public FossilPrior() {
		distInput.setRule(Validate.OPTIONAL);
	}
	
	public void initAndValidate() {
		// this only makes sense as a prior on the parent of a clade,
		// so tipsonly and useOriginate will be set to reflect this.
		onlyUseTipsInput.setValue(false, this);
		useOriginateInput.setValue(true, this);
		
		distInput.setValue(calibrationDistr.get(), this);
		super.initAndValidate();
		distInput.setValue(null, this);
	};
	
	
    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
    	String id = (taxonsetInput.get() != null ? taxonsetInput.get().getID() : null);
    	if (id == null) {
    		id = getID();
    	}
        if (dist != null || isMonophyleticInput.get()) {
            out.print("logP(mrca(" + id + "))\t");
        }
        out.print("mrcatime(" + id + ")\t");
    }

}
