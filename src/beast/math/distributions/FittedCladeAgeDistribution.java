package beast.math.distributions;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import beast.core.Description;

@Description("Continuous distribution fitted to set of empirically determined densities")
public class FittedCladeAgeDistribution implements ContinuousDistribution {
	
	// The first occurrence age minimum.
	double first_occurrence_age_min;
	
	// The first occurrence age maximum.
	double first_occurrence_age_max;
	
	// Parameter C of the fitted CladeAge distribution.
	double c;
	
	// Parameter S of the fitted CladeAge distribution.
	double s;
	
	// Parameter M of the fitted CladeAge distribution.
	double m;
	
	// Parameter W of the fitted CladeAge distribution.
	double w;
	
	// The root mean square deviation of the fitted CladeAge distribution.
	double rmsd;
	
	public FittedCladeAgeDistribution(double first_occurrence_age_min, double first_occurrence_age_max, double fittedCladeAgeC, double fittedCladeAgeS, double fittedCladeAgeM, double fittedCladeAgeW, double fittedCladeAgeRmsd) {
		 this.first_occurrence_age_min = first_occurrence_age_min;
		 this.first_occurrence_age_max = first_occurrence_age_max;
		 this.c = fittedCladeAgeC;
		 this.s = fittedCladeAgeS;
		 this.m = fittedCladeAgeM;
		 this.w = fittedCladeAgeW;
		 this.rmsd = fittedCladeAgeRmsd;
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
	public double inverseCumulativeProbability(double age) throws MathException {
		// XXX todo
		return 0.0;
	}
	
	@Override
	public double density(double age) {
		
		if (age <= this.first_occurrence_age_min) {
			return 0;
		} else {			
			if (this.first_occurrence_age_min == this.first_occurrence_age_max) {
				return (this.c/(age-first_occurrence_age_min+this.s)) * Math.exp( -Math.pow(Math.log(age-first_occurrence_age_min+this.s)-this.m,2)/this.w );
			} else {
				// Calculate the integral of the fitted CladeAge distribution pdf for an upper bound, which is just the same as this age.
				// The fitted CladeAge distribution pdf is c/(x+s) * exp(-(log(x+s)-m)**2/w), and the integral of this pdf is
				// -1/2 * c * sqrt(pi) * sqrt(w) * erf(-(log(x+s)-m)/w)
				// whereby erf is the error function (see http://en.wikipedia.org/wiki/Error_function)
				double erf = 0;
				double x_upper = age-this.first_occurrence_age_min;
				double erf_argument = -(Math.log(x_upper+this.s)-this.m) / Math.sqrt(this.w);
				// Calculate the error function using the Wikipedia approximation (which is from Numerical Recipes in Fortran 77).
				double erf_t = 0;
				if (erf_argument > 0) {
					erf_t = 1/(1.0+0.5*erf_argument);
				} else {
					erf_t = 1/(1.0-0.5*erf_argument);
				}
				double erf_polynomial = -Math.pow(erf_argument, 2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t, 2) + 0.09678418*Math.pow(erf_t, 3) - 0.18628806*Math.pow(erf_t, 4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
				double erf_tau = erf_t * Math.exp(erf_polynomial);
				if (erf_argument >= 0) {
					erf = 1 - erf_tau;
				} else {
					erf = erf_tau - 1;
				}
				double int_upper = -0.5 * this.c * Math.sqrt(Math.PI) * Math.sqrt(this.w) * erf;
				
				// Calculate the integral of the probabilities for a lower bound which is either first_occurrence_age_min, or this age minus the difference between first_occurrence_age_max and first_occurrence_age_min, depending on which of the two is greater.
				double x_lower = 0;
				if (age > this.first_occurrence_age_max) {
					x_lower = age - this.first_occurrence_age_max;
				}
				erf_argument = -(Math.log(x_lower+this.s)-this.m) / Math.sqrt(this.w);
				// Calculate the error function using the Wikipedia approximation (which is from Numerical Recipes in Fortran 77).
				if (erf_argument > 0) {
					erf_t = 1/(1.0+0.5*erf_argument);
				} else {
					erf_t = 1/(1.0-0.5*erf_argument);
				}
				erf_polynomial = -Math.pow(erf_argument, 2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t, 2) + 0.09678418*Math.pow(erf_t, 3) - 0.18628806*Math.pow(erf_t, 4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
				erf_tau = erf_t * Math.exp(erf_polynomial);
				if (erf_argument >= 0) {
					erf = 1 - erf_tau;
				} else {
					erf = erf_tau - 1;
				}
				double int_lower = -0.5 * this.c * Math.sqrt(Math.PI) * Math.sqrt(this.w) * erf;
				
				// The density for this age is the integral of the upper bound minus the integral of the lower bound, divided by the difference between first_occurrence_age_max and first_occurrence_age_min.
				return (int_upper-int_lower)/(this.first_occurrence_age_max-this.first_occurrence_age_min);
			}
		}
	}
	
	@Override
	public double logDensity(double age) {
		return Math.log(density(age));
	}

}