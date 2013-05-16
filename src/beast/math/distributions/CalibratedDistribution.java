package beast.math.distributions;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

@Description("Distribution based on parameters estimated from a fossil record.")
public class CalibratedDistribution extends ParametricDistribution {
    public Input<RealParameter> alphaInput = new Input<RealParameter>("alphaMale", "first shape parameter, defaults to 1");
    public Input<RealParameter> betaInput = new Input<RealParameter>("beta", "the other shape parameter, defaults to 1");
    public Input<RealParameter> gammaInput = new Input<RealParameter>("gamma", "yet another shape parameter, defaults to 1");

    final org.apache.commons.math.distribution.GammaDistribution dist1 = new GammaDistributionImpl(1.0, 1.0);
    final org.apache.commons.math.distribution.ExponentialDistribution dist2 = new ExponentialDistributionImpl(1.0);;
    ContinuousDistribution m_dist;
    
    @Override
    public void initAndValidate() throws Exception {
    	refresh();
    	m_dist = new ContinuousDistribution() {
			
			@Override
			public double cumulativeProbability(double x0, double x1) throws MathException {
				return (dist1.cumulativeProbability(x0, x1) + dist2.cumulativeProbability(x0, x1)) / 2.0;
			}
			
			@Override
			public double cumulativeProbability(double x) throws MathException {
				return (dist1.cumulativeProbability(x) + dist2.cumulativeProbability(x))/ 2.0;
			}
			
			@Override
			public double logDensity(double x) {
				return (dist1.logDensity(x) + dist2.logDensity(x))/ 2.0;
			}
			
			@Override
			public double inverseCumulativeProbability(double p) throws MathException {
				return (dist1.inverseCumulativeProbability(p) + dist2.inverseCumulativeProbability(p))/ 2.0;
			}
			
			@Override
			public double density(double x) {
				return (dist1.density(x) + dist2.density(x))/ 2.0;
			}
		};
    }

	private void refresh() {
    	double alpha = alphaInput.get().getValue();
    	double beta = betaInput.get().getValue();
    	double gamma = gammaInput.get().getValue();
		dist1.setAlpha(alpha);
		dist1.setBeta(beta);
		dist2.setMean(gamma);
	}

	@Override
	public ContinuousDistribution getDistribution() {
        refresh();
		return m_dist;
	}

}

