package beast.math.distributions;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;


@Deprecated // use CompoundDistribution instead
@Description("Takes a collection of fossil priors and combines them into the compound of these distributions.")
public class CompoundFossilPrior extends Distribution {
    public Input<List<FossilPrior>> distributionsInput =
            new Input<>("fossilPrior", "individual priors based on the fossil record", new ArrayList<FossilPrior>());
    
    List<FossilPrior> distributions;

    @Override
    public void initAndValidate() {
    	distributions = distributionsInput.get();
        if (distributions.size() == 0) {
            logP = 0;
        }
    }

    @Override
    public double calculateLogP() {
        logP = 0;
        for (Distribution dists : distributions) {
            if (dists.isDirtyCalculation()) {
                logP += dists.calculateLogP();
            } else {
                logP += dists.getCurrentLogP();
            }
            if (Double.isInfinite(logP) || Double.isNaN(logP)) {
                return logP;
            }
        }
        return logP;
    }
    
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}
