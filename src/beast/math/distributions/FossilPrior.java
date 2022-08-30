package beast.math.distributions;

import java.io.PrintStream;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.MRCAPrior;


@Description("Prior based on information from the fossil record")
@Citation("Matschiner M, Musilova Z, Barth JMI, Starostova Z, Salzburger W, Steel M, and\n" +
	"  Bouckaert R (2017) Bayesian phylogenetic estimation of clade ages supports\n" + 
	"  trans-Atlantic dispersal of cichlid fishes. Systematic Biology 66: 3–22.")
public class FossilPrior extends MRCAPrior {
	public Input<FossilCalibration> calibrationDistr = new Input<>("fossilDistr", "", Validate.REQUIRED);

	public FossilPrior() {
		distInput.setRule(Validate.OPTIONAL);
		isMonophyleticInput.setValue(true, this);
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
