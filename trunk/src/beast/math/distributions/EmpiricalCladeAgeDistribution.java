package beast.math.distributions;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import beast.core.Description;

@Description("Continuous distribution defined by a set of empirically determined densities")
public class EmpiricalCladeAgeDistribution implements ContinuousDistribution {

	// The ages of points, sorted from lowest to highest.
	double [] ages;
	
	// The densities associated with ages.
	double [] densities;
	
	// A flag to indicate that the upper end is approximated by an exponential 
	// that fits through the last point and third last point
	boolean useExponentialTail = true;

	// The coefficients for exponential function f(x) = beta*alpha^x
	double alpha, beta;
	
	// The normaliser for the distribution.
	double normaliser = 0;
	
	public EmpiricalCladeAgeDistribution(double[] ages, double[] densities, boolean useExponentialTail) {
		this.ages = ages;
		this.densities = densities;
		this.useExponentialTail = useExponentialTail;
	}
	
	
	@Override
	public double cumulativeProbability(double age) throws MathException {		
		// XXX todo
		return 0.0;
	}

	@Override
	public double cumulativeProbability(double age0, double age1) throws MathException {
		return cumulativeProbability(age1) - cumulativeProbability(age0);
	}

	@Override
	public double inverseCumulativeProbability(double p) throws MathException {
		// XXX todo
		return 0.0;
	}
	
	@Override
	public double density(double x) {
		// XXX todo
		return 0.0;
	}
	
	@Override
	public double logDensity(double x) {
		return Math.log(density(x));
	}

}