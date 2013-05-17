package beast.math.distributions;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.core.Description;

@Description("Mixture of exponential and gamma distribution")
public class ExpGamma implements ContinuousDistribution {
    final org.apache.commons.math.distribution.ExponentialDistribution dist1 = new ExponentialDistributionImpl(1.0);
    final org.apache.commons.math.distribution.GammaDistribution dist2 = new GammaDistributionImpl(1.0, 1.0);
    
    
    public double getMean() {
    	return dist1.getMean();
    }
    public double getAlpha() {
    	return dist2.getAlpha();
    }
    public double getBeta() {
    	return dist2.getBeta();
    }
    public double getWeight() {
    	return weight;
    }
    
    
    // weight of the exponential, 1-weight will be used for gamma
    double weight = 0.5;
    
    public void setParameters(double weight, double mean, double alpha, double beta) {
    	this.weight = weight;
		dist1.setMean(mean);
		dist2.setAlpha(alpha);
		dist2.setBeta(beta);
    }
    
	@Override
	public double cumulativeProbability(double x0, double x1) throws MathException {
		return (weight * dist1.cumulativeProbability(x0, x1) + (1.0-weight) * dist2.cumulativeProbability(x0, x1)) / 2.0;
	}
	
	@Override
	public double cumulativeProbability(double x) throws MathException {
		return (weight * dist1.cumulativeProbability(x) + (1.0-weight) * dist2.cumulativeProbability(x))/ 2.0;
	}
	
	@Override
	public double logDensity(double x) {
		return (weight * dist1.logDensity(x) + (1.0-weight) * dist2.logDensity(x))/ 2.0;
	}
	
	@Override
	public double inverseCumulativeProbability(double p) throws MathException {
		return (weight * dist1.inverseCumulativeProbability(p) + (1.0-weight) * dist2.inverseCumulativeProbability(p))/ 2.0;
	}
	
	@Override
	public double density(double x) {
		return (weight * dist1.density(x) + (1.0-weight) * dist2.density(x))/ 2.0;
	}
	


}
