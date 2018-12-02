package beast.math.distributions;

import java.awt.geom.Point2D.Double;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;

import beast.core.Description;

@Description("Continuous distribution defined by a set of densities determined through simulation")
public class CladeAgeDistribution implements ContinuousDistribution, Serializable {

	private static final long serialVersionUID = 1L;

	// The ages of points, sorted from lowest to highest.
	double [] ages;
	
	// The densities associated with ages.
	double [] densities;
	
	// Cumulative density associated with point x.
	double [] cummulative;

	// A flag to indicate that the upper end is approximated by an exponential 
	// that fits through the last point and third last point.
	boolean useExponentialTail = true;

	// The coefficients for exponential function f(x) = beta*alpha^x.
	double alpha, beta;
	
	// The normaliser for the distribution.
	double normaliser = 0;
	
	public CladeAgeDistribution(double[] ages, double[] densities, boolean useExponentialTail) throws Exception {
		if (ages.length != densities.length) {
			throw new Exception("Arrays 'ages' and 'densitites' should be of same length!");
		}
		
		// Sort points by ages.
		java.awt.geom.Point2D.Double [] points = new java.awt.geom.Point2D.Double[ages.length];
		for (int i = 0; i < points.length; i++) {
			points[i] = new java.awt.geom.Point2D.Double();
			points[i].x = ages[i];
			points[i].y = (java.lang.Double.isNaN(densities[i]) ? 0: densities[i]);
		}
		Arrays.sort(points, new Comparator<java.awt.geom.Point2D.Double>() {
			@Override
			public int compare(Double o1, Double o2) {
				if (o1.x > o2.x) {
					return 1;
				}
				if (o1.x < o2.x) {
					return -1;
				}
				return 0;
			}
		});
		this.ages = new double[points.length];
		this.densities = new double[points.length];
		for (int i = 0; i < points.length; i++) {
			this.ages[i] = points[i].x;
			this.densities[i] = points[i].y;
		}

		// Calculate parameters of the exponential tail.
		if (useExponentialTail) {
			// Take the mean decrease of the last X points.
			double X = 1;
			double m = 0;
			for (int i = 0; i < X; i++) {
				m += this.densities[points.length - 1 - i]/this.densities[points.length - 2 - i];
			}
			m /= X;
			double d = this.ages[points.length - 1] - this.ages[points.length - 2];
			this.alpha =  Math.pow(m, 1/d);
			double x2 = this.ages[points.length - 2];
			double y2 = this.densities[points.length - 2];
			this.beta = y2 / Math.pow(alpha, x2);
		}

		// Normalise the distribution.
		this.normaliser = 0;
		this.normaliser += this.densities[0] * (this.ages[1] - this.ages[0]);
		for (int i = 1; i < points.length - 1; i++) {
			this.normaliser += this.ages[i] * (this.ages[i+1] - this.ages[i-1])/2.0;
		}
		this.normaliser += this.densities[this.ages.length-1] * (this.ages[ages.length-1] - this.ages[ages.length-2]);
		if (useExponentialTail) {
			 double expContribution = -this.beta * Math.pow(alpha, this.ages[points.length - 1])/Math.log(this.alpha);
			 this.normaliser += expContribution;
		}
		for (int i = 0; i < points.length; i++) {
			this.densities[i] = this.densities[i]/normaliser;
		}
		if (useExponentialTail) {
			this.beta = this.beta / this.normaliser;
		}

		// Calculate the cumulative density so that it can be returned with function 'cumulativeProbability()'
		this.cummulative = new double[points.length];
		this.cummulative[0] = densities[0] * (ages[1] - ages[0]);
		for (int i = 1; i < points.length - 1; i++) {
			this.cummulative[i] = this.cummulative[i-1] + this.densities[i] * (this.ages[i+1] - this.ages[i-1])/2.0;
		}
		this.cummulative[this.ages.length-1] = this.cummulative[this.ages.length-2] + this.densities[this.ages.length-1] * (this.ages[this.ages.length-1] - this.ages[this.ages.length-2]);

		// Memorize the boolean for the exponential tail. 
		this.useExponentialTail = useExponentialTail;
	}
	
	@Override
	public double cumulativeProbability(double age) throws MathException {		
		if (age < 0) {
			return 0;
		}
		if (age > this.ages[this.ages.length-1]) {
			if (useExponentialTail) {
				return 1 - beta * Math.pow(alpha, age);
			} else {
				return 1;
			}
		}
		int i = Arrays.binarySearch(this.ages, age);
		if (i >= 0) {
			return cummulative[i];
		}
		i = -i-2;
		double w = this.ages[i+1] - this.ages[i]; 
		double a = 1.0 -  (age - this.ages[i])/w;
		double cumulative = this.cummulative[i] * a + this.cummulative[i+1] * (1.0 - a); 
		return cumulative;
	}

	@Override
	public double cumulativeProbability(double age0, double age1) throws MathException {
		return cumulativeProbability(age1) - cumulativeProbability(age0);
	}

	@Override
	public double inverseCumulativeProbability(double p) throws MathException {
		if (p < 0 || p > 1) {
			throw new MathException();
		}
		if (p > cummulative[cummulative.length-1]) {
			if (useExponentialTail) {
				// int_q^\infty beta*alpha^x dx = 1-p
				// solve integral
				// beta*alpha^x/log alpha |_q^\infty = alpha^x = 1-p
				// alpha^x = (1-p) * log alpha / beta
				// x = log((1-p) * log alpha / beta)
				return (Math.log(1.0 - p) + Math.log(Math.abs(Math.log(alpha))) - Math.log(beta))/Math.log(alpha);
			} else {
				return this.ages[cummulative.length-1];
			}
		}
		int i = Arrays.binarySearch(cummulative, p);
		if (i >= 0) {
			return this.ages[i];
		}
		i = -i-2;
		if (i < 0) {
			return ages[0];
		}
		double w = cummulative[i+1] - cummulative[i]; 
		double a = 1.0 -  (p - cummulative[i])/w;
		double age = this.ages[i] * a + this.ages[i+1] * (1.0 - a); 
		return age;
	}
	
	@Override
	public double density(double age) {
		if (age < this.ages[0]) {
			return 0;
		}
		if (age > this.ages[this.ages.length-1]) {
			if (useExponentialTail) {
				return beta * Math.pow(alpha, age);
			} else {
				return 0;
			}
		}
		int i = Arrays.binarySearch(this.ages, age);
		if (i >= 0) {
			return densities[i];
		}
		i = -i-2;
		if (i < 0) {
			return 0;
		}
		double w = this.ages[i+1] - this.ages[i];
		double a = 1.0 -  (age - this.ages[i])/w;
		double density = this.densities[i] * a + this.densities[i+1] * (1-a); 
		return density;
	}
	
	@Override
	public double logDensity(double age) {
		return Math.log(density(age));
	}

}