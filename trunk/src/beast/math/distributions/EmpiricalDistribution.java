package beast.math.distributions;

import java.awt.geom.Point2D.Double;
import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;

import beast.core.Description;

@Description("Continuous distribution defined by a set of density points")
public class EmpiricalDistribution implements ContinuousDistribution {
	
	// x - coordinate of points, sorted from lowest to highest
	double [] x;
	// density associated with point x
	double [] p;
	// cumulative density associated with point x
	double [] cummulative;
	
	// flag to indicate that the upper end is approximated by an exponential 
	// that fits through the last point and third last point
	boolean useExponentialTail = true;
	
	// coefficients for exponential function f(x) = beta.alpha^x
	double alpha, beta;

	double normaliser = 0;

	public void setup(double [] x_, double [] p_, boolean useExponentialTail) throws Exception {
		if (x_.length != p_.length) {
			throw new Exception("x and p should be of same length");
		}
		
		// sort points by x
		java.awt.geom.Point2D.Double [] points= new java.awt.geom.Point2D.Double[x_.length];
		for (int i = 0; i < points.length; i++) {
			points[i] = new java.awt.geom.Point2D.Double();
			points[i].x = x_[i];
			points[i].y = (java.lang.Double.isNaN(p_[i]) ? 0: p_[i]);
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
		x = new double[points.length];
		p = new double[points.length];
		for (int i = 0; i < points.length; i++) {
			x[i] = points[i].x;
			p[i] = points[i].y;
		}
		if (useExponentialTail) {
			// take the mean decrease of the last X points
			double X = 1;
			double m = 0;
			for (int i = 0; i < X; i++) {
				m += p[points.length - 1 - i]/p[points.length - 2 - i];
			}
			m /= X;
			double d = x[points.length - 1] - x[points.length - 2];
			alpha =  Math.pow(m, 1/d);

			double x2 = x[points.length - 2];
			double y2 = p[points.length - 2];
			beta = y2 / Math.pow(alpha, x2);
		}
		
		// normalise
		normaliser = 0;
		normaliser += p[0] * (x[1] - x[0]);
		for (int i = 1; i < points.length - 1; i++) {
			normaliser += p[i] * (x[i+1] - x[i-1])/2.0;
		}
		normaliser += p[p.length-1] * (x[p.length-1] - x[p.length-2]);
		if (useExponentialTail) {
			 double expContribution = -beta * Math.pow(alpha, x[points.length - 1])/Math.log(alpha);
			 normaliser += expContribution;
		}
		for (int i = 0; i < points.length; i++) {
			p[i] = p[i]/normaliser;
			p_[i] = p_[i]/normaliser;
		}
		if (useExponentialTail) {
			beta = beta / normaliser;
		}
				
		cummulative = new double[points.length];
		cummulative[0] = p[0] * (x[1] - x[0]);
		for (int i = 1; i < points.length - 1; i++) {
			cummulative[i] = cummulative[i-1] + p[i] * (x[i+1] - x[i-1])/2.0;
		}
		cummulative[p.length-1] = cummulative[p.length-2] + p[p.length-1] * (x[p.length-1] - x[p.length-2]);
		
		this.useExponentialTail = useExponentialTail;
	}

	@Override
	public double cumulativeProbability(double x) throws MathException {
		if (x < 0) {
			return 0;
		}
		if (x > this.x[this.x.length-1]) {
			if (useExponentialTail) {
				return 1 - beta * Math.pow(alpha, x);
			} else {
				return 1;
			}
		}
		int i = Arrays.binarySearch(this.x, x);
		if (i >= 0) {
			return cummulative[i];
		}
		i = -i-2;
		double w = this.x[i+1] - this.x[i]; 
		double a = 1.0 -  (x - this.x[i])/w;
		double c = this.cummulative[i] * a + this.cummulative[i+1] * (1.0 - a); 
		return c;
	}
	
	@Override
	public double cumulativeProbability(double x0, double x1) throws MathException {
		return cumulativeProbability(x1) - cumulativeProbability(x0);
	}
	
	@Override
	public double inverseCumulativeProbability(double p) throws MathException {
		if (p < 0 || p > 1) {
			throw new MathException();
		}
		if (p > cummulative[cummulative.length-1]) {
			if (useExponentialTail) {
				return Math.log((1.0 - p) / beta)/Math.log(alpha);
			} else {
				return x[cummulative.length-1];
			}
		}
		int i = Arrays.binarySearch(cummulative, p);
		if (i >= 0) {
			return x[i];
		}
		i = -i-2;
		double w = cummulative[i+1] - cummulative[i]; 
		double a = 1.0 -  (p - cummulative[i])/w;
		double x = this.x[i] * a + this.x[i+1] * (1.0 - a); 
		return x;
	}
	
	@Override
	public double density(double x) {
		if (x < this.x[0]) {
			return 0;
		}
		if (x > this.x[this.x.length-1]) {
			if (useExponentialTail) {
				return beta * Math.pow(alpha, x);
			} else {
				return 0;
			}
		}
		int i = Arrays.binarySearch(this.x, x);
		if (i >= 0) {
			return p[i];
		}
		i = -i-2;
		double w = this.x[i+1] - this.x[i];
		double a = 1.0 -  (x - this.x[i])/w;
		double p = this.p[i] * a + this.p[i+1] * (1-a); 
		return p;
	}
	
	@Override
	public double logDensity(double x) {
		return Math.log(density(x));
	}
}