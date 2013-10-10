package beast.app.ca;

import java.util.ArrayList;

import javax.swing.JProgressBar;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.GammaDistributionImpl;

// Import CladeAge distributions.
import beast.math.distributions.EmpiricalCladeAgeDistribution;
import beast.math.distributions.FittedCladeAgeDistribution;

public class CladeAgeProbabilities {

	// Initialize objects that are to become variables of class instances.
	private int number_of_ages = 100;
	private double[] ages = new double[number_of_ages+1];
	private double[] fix_probabilities = new double[number_of_ages+1];
	private double[] int_probabilities = new double[number_of_ages+1];
	private boolean cancel = false;
	private int first_occurrence_age_min;
	private int first_occurrence_age_max;
	private double[] distribution_parameters = new double[8];
	private double distribution_rmsd;
	
	// Initialize constants for Lanczos approximation of the gamma function.
	private double lanczos_g = 7;
	private double[] lanczos_p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,771.32342877765313, -176.61502916214059, 12.507343278686905,-0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

	// Initialize constants for the Nelder-Mead Downhill Simplex (according to Wikipedia example).
	private double alph = 1.0;
	private double gamm = 2.0;
	private double beta = -0.5;
	private double delt  = 0.5;
	
	// Initialize the normalizer for simulated probabilities
	private double normaliser = 1.0;

	// Set the amount of progress for the progress indicator.
	final static int DELTA_PROGRESS = 10;


	public int getFirst_occurrence_age_min() {
		return first_occurrence_age_min;
	}

	public int getFirst_occurrence_age_max() {
		return first_occurrence_age_max;
	}

	public void setAges(double[] ages) {
		this.ages = ages;
	}

	public double[] getAges() {
		return ages;
	}

	public void setFix_probabilities(double[] fix_probabilities) {
		this.fix_probabilities = fix_probabilities;
	}

	public double[] getFix_probabilities() {
		return fix_probabilities;
	}

	public void setInt_probabilities(double[] int_probabilities) {
		this.int_probabilities = int_probabilities;
	}

	public double[] getInt_probabilities() {
		return int_probabilities;
	}

	public double[] getDistribution_parameters() {
		return distribution_parameters;
	}

	public double getDistribution_rmsd() {
		return distribution_rmsd;
	}

	public void setCancel() {
		cancel = true;
	}
	public boolean getCancel() {
		return cancel;
	}

	public EmpiricalCladeAgeDistribution run_empirical_cladeage(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, double sampling_gap_min, double sampling_gap_max, JProgressBar dpb) throws Exception {

		// Make sure cancel is set to false.
		cancel = false;

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		int_probabilities = new double[number_of_ages + 1];
		double[] raw_probabilities = new double[number_of_ages + 1];
		
		// Part 1: Calculate the maximum age, and 100 time points between fossil age and maximum age.
		
		// Roughly calculate a maximum tree duration so that the probability for the maximum tree duration is 0.001 * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		double ndr_mean = (ndr_max+ndr_min)/2.0;
		double psi_mean = (psi_max+psi_min)/2.0;
		
		// Initiate the simplex, find 2 vertices.
		// vertex0
		// Math.log(ndr_mean/psi_mean)/ndr_mean is the peak of p(t), as determined through its first derivative.
		double vertex0t = 0;
		if (Math.log(ndr_mean/psi_mean)/ndr_mean > 0) {
			vertex0t = Math.log(ndr_mean/psi_mean)/ndr_mean;
		}
		
		// vertex1
		double vertex1t = vertex0t + 10;
		
		// Prepare for the Nelder-Mead loop.
		boolean keepGoing = true;
		
		// Until converged do the loop.
		while (keepGoing == true) {
			
			// Calculate the y value of each vertex.
			// vertex0
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(0.001)), 2);

			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Find the best (=lowest) and the worst (=highest) y value, and the respective t values.
			double bestY = vertex0Y;
			double bestt = vertex0t;
			double worstY = vertex1Y;
			double worstt = vertex1t;
			if (vertex1Y < bestY) {
				bestY = vertex1Y;
				bestt = vertex1t;
				worstY = vertex0Y;
				worstt = vertex0t;
			}

			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
			double reflectiont = bestt + alph * (bestt - worstt);
			
			// Calculate the y value of the reflection.
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Figure out which values to use as replacement for the values of the worst vertex.
				double replacet = 0;
				if (reflectionY < extensionY) {
					replacet = reflectiont;
				} else {
					replacet = extensiont;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0t = replacet;
				} else if (vertex1Y == worstY) {
					vertex1t = replacet;
				}

			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractiont = bestt + beta * (bestt - worstt);
				
				// Calculate the y value of the contraction.
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0t = contractiont;
					} else if (vertex1Y == worstY) {
						vertex1t = contractiont;
					}
				
				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that delt = 0.5, the default setting)
					// vertex0
					vertex0t = bestt + delt  * (vertex0t - bestt);
					
					// vertex1
					vertex1t = bestt + delt  * (vertex1t - bestt);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY) {
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0t - vertex1t) > 0.0001) {
				keepGoing = true;
			} 

		} // while (keepGoing == true) {

		// Memorize the result as the max_simulation_age.
		double max_simulation_age = vertex0t + first_occurrence_age_min;

		// Determine the ages.
		ages[0] = first_occurrence_age_min;
		for (int i = 0; i < number_of_ages; i++) {
			ages[i+1] = first_occurrence_age_min + ((i+1)/(double) number_of_ages)*(max_simulation_age-first_occurrence_age_min);
		}

		// Part 2: Estimate clade age probability densities
		
		// Prepare the arrays raw_probabilities and successful simulations (those in which at least one taxon survives), and fill them with 0.
		int[] successful_simulations = new int[number_of_ages];
		
		while (successful_simulations[-1] < 10000) {
			if (cancel) {
				return null;
			}
			// Increment the progress indicator.
			if (dpb != null) {
				dpb.setValue(successful_simulations[-1]);
			}

			// Draw values for parameters ndr (net diversification rate, lambda - mu) and epsilon (turnover rate, mu/lambda) from uniform distributions, and calculate lambda and mu from them.
			double ndr = ndr_min + Math.random() * (ndr_max-ndr_min);
			double epsilon = epsilon_min + Math.random() * (epsilon_max-epsilon_min);
			double mu = (ndr*epsilon)/(1.0 - epsilon);
			double lambda = ndr + mu;

			// Draw a duration t for this simulation.
			double first_occurrence_age = first_occurrence_age_min + Math.random()*(first_occurrence_age_max-first_occurrence_age_min);
			
			// Draw a sampling gap.
			double sampling_gap = sampling_gap_min + Math.random()*(sampling_gap_max-sampling_gap_min);
			
			// Determine the maximum tree duration.
			double max_tree_duration = ages[-1] - first_occurrence_age;
			
			// Initiate arrays for branch origin and branch termination.
			ArrayList<Double> branch_origin = new ArrayList<Double>();
			ArrayList<Double> branch_termination = new ArrayList<Double>();
			ArrayList<Double> new_branch_origin = new ArrayList<Double>();
			ArrayList<Double> new_branch_termination = new ArrayList<Double>();
			new_branch_origin.add(0.0);

			// Start the tree generation loop.
			while (new_branch_origin.size() > 0) {

				// For each new origin, add a new termination.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					new_branch_termination.add(new_branch_origin.get(o) + (Math.log(Math.random())/(-(lambda+mu))));
				}

				// Add new origin and termination to the old collection.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					branch_origin.add(new_branch_origin.get(o));
				}
				for (int t = 0; t < new_branch_termination.size(); t++) {
					branch_termination.add(new_branch_termination.get(t));
				}
				
				// Empty the new origin array.
				new_branch_origin = new ArrayList<Double>();
				
				// For each new termination, add it to the new origin array if it is < max_tree_duration and rand < lambda/(lambda+mu) - this represents a speciation event.
				for (int t = 0; t < new_branch_termination.size(); t++) {										
					if (new_branch_termination.get(t) < max_tree_duration) {
						if (Math.random() < lambda/(lambda+mu)) {
							new_branch_origin.add(new_branch_termination.get(t));
							new_branch_origin.add(new_branch_termination.get(t));
						}
					}
				}
				
				// Empty the new termination array.
				new_branch_termination = new ArrayList<Double>();
				
			} // while (new_branch_origin.size() > 0)
			
			// Analyse the obtained set of branches.
			for (int i = ages.length-1; i >= 0; i--) {
				
				// For each age, determine the tree duration.
				double tree_duration = ages[i] - first_occurrence_age;
				
				// Trim the tree so that no branches are older than ages[i] - first_occurrence_age. While trimming, make sure the tree contains at least one extant species.
				if (tree_duration >= 0) {
					double sum_of_species_durations = 0.0;
					int number_of_extant_taxa = 0;
					for (int o = branch_origin.size()-1; o >= 0; o--) {
						boolean remove_this_branch = false;
						if (branch_origin.get(o) > tree_duration) {
							remove_this_branch = true;
						} else {
							if (branch_termination.get(o) >= tree_duration) {
								branch_termination.set(o,tree_duration);
								number_of_extant_taxa += 1;
							}
							if (branch_termination.get(o) > sampling_gap) {
								if (branch_origin.get(o) > sampling_gap) {
									sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
								} else {
									sum_of_species_durations += branch_termination.get(o)-sampling_gap;
								}
							}
						}
						if (remove_this_branch == true) {
							branch_origin.remove(o);
							branch_termination.remove(o);
						}
					} // for (int o = branch_origin.size()-1; o >= 0; o--)
					
					// If at least one species is extant at this age, get the sum of lineage durations, excluding the sampling gap.
					if (number_of_extant_taxa > 0) {
						// Draw a value for the sampling rate psi.
						for (int pp = 0; pp < 10; pp++) {
							double psi = psi_min + Math.random()*(psi_max-psi_min);
							if (tree_duration >= sampling_gap) {
								raw_probabilities[i] += (psi*Math.exp(-psi*sum_of_species_durations)*number_of_extant_taxa)/10.0;
							}
						}
						successful_simulations[i] += 1;
					} // if (extant_at_this_age == true)

				} else {
					successful_simulations[i] += 1;
				} // if (tree_duration > 0)
				
			} // for (int i = 0; i < ages.length; i++)
			
		} // while (successful_simulations[0] < bd_sample_size)
		
		for (int i = 0; i < raw_probabilities.length; i ++) {
			int_probabilities[i] = raw_probabilities[i]/(double) successful_simulations[i];
		}
		
		// Return the empirical CladeAge distribution.
		return new EmpiricalCladeAgeDistribution(ages, int_probabilities, true);

	} // public void run_empirical_cladeage(...)

	public FittedCladeAgeDistribution run_fitted_cladeage(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, JProgressBar dpb) {
		
		// Make sure cancel is set to false.
		cancel = false;

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		fix_probabilities = new double[number_of_ages + 1];
		int_probabilities = new double[number_of_ages + 1];
		double[] raw_probabilities = new double[number_of_ages + 1];

		// Part 1: Calculate the maximum age, and 100 time points between fossil age and maximum age.
		
		// Roughly calculate a maximum age so that the probability for the maximum age is 0.001 * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		double ndr_mean = (ndr_max+ndr_min)/2.0;
		double psi_mean = (psi_max+psi_min)/2.0;
		
		// Initiate the simplex, find 2 vertices.
		// vertex0
		// Math.log(ndr_mean/psi_mean)/ndr_mean is the peak of p(t), as determined through its first derivative.
		double vertex0t = 0;
		if (Math.log(ndr_mean/psi_mean)/ndr_mean > 0) {
			vertex0t = Math.log(ndr_mean/psi_mean)/ndr_mean;
		}
		
		// vertex1
		double vertex1t = vertex0t + 10;
		
		// Prepare for the Nelder-Mead loop.
		boolean keepGoing = true;
		
		// Until converged do the loop.
		while (keepGoing == true) {
			
			// Calculate the y value of each vertex.
			// vertex0
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(0.001)),2);
			
			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(0.001)),2);
			
			// Find the best (=lowest) and the worst (=highest) y value, and the respective t values.
			double bestY = vertex0Y;
			double bestt = vertex0t;
			double worstY = vertex1Y;
			double worstt = vertex1t;
			if (vertex1Y < bestY) {
				bestY = vertex1Y;
				bestt = vertex1t;
				worstY = vertex0Y;
				worstt = vertex0t;
			}
			
			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
			double reflectiont = bestt + alph * (bestt - worstt);
			
			// Calculate the y value of the reflection.
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(0.001)),2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(0.001)),2);
				
				// Figure out which values to use as replacement for the values of the worst vertex.
				double replacet = 0;
				if (reflectionY < extensionY) {
					replacet = reflectiont;
				} else {
					replacet = extensiont;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0t = replacet;
				} else {
					vertex1t = replacet;
				}
				
			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractiont = bestt + beta * (bestt - worstt);
				
				// Calculate the y value of the contraction.
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(0.001)),2);
				
				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0t = contractiont;
					} else {
						vertex1t = contractiont;
					}
					
				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
					// vertex0
					vertex0t = bestt + delt  * (vertex0t - bestt);
					
					// vertex1
					vertex1t = bestt + delt  * (vertex1t - bestt);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY)
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0t-vertex1t) > 0.001) {
				keepGoing = true;
			}
			
		} // while (keepGoing == true)
			
		// Memorize the result as the max_simulation_age
		double max_simulation_age = vertex0t + first_occurrence_age_min;
				
		// Determine the ages.
		ages[0] = first_occurrence_age_min;
		for (int i = 0; i < number_of_ages; i++) {
			ages[i+1] = first_occurrence_age_min + ((i+1)/(double) number_of_ages)*(max_simulation_age-first_occurrence_age_min);
		}

		// Part 2: Estimate clade age probability densities.
		
		// Prepare the arrays raw_probabilities and successful simulations (those in which at least one taxon survives), and fill them with 0.
		int[] successful_simulations = new int[number_of_ages];
		
		while (successful_simulations[-1] < 10000) {
			if (cancel) {
				return null;
			}
			// Update the progress bar.
			if (dpb != null) {
				dpb.setValue(successful_simulations[-1]);
			}

			// Draw values for parameters ndr (net diversification rate, lambda - mu) and epsilon (turnover rate, mu/lambda) from uniform distributions, and calculate lambda and mu from them.
			double ndr = ndr_min + Math.random() * (ndr_max-ndr_min);
			double epsilon = epsilon_min + Math.random() * (epsilon_max-epsilon_min);
			double mu = (ndr*epsilon)/(1.0 - epsilon);
			double lambda = ndr + mu;
			
			// Determine the maximum tree duration.
			double max_tree_duration = ages[-1] - first_occurrence_age_min;
			
			// Initiate arrays for branch origin and branch termination.
			ArrayList<Double> branch_origin = new ArrayList<Double>();
			ArrayList<Double> branch_termination = new ArrayList<Double>();
			ArrayList<Double> new_branch_origin = new ArrayList<Double>();
			ArrayList<Double> new_branch_termination = new ArrayList<Double>();
			new_branch_origin.add(0.0);

			// Start the tree generation loop.
			while (new_branch_origin.size() > 0) {

				// For each new origin, add a new termination.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					new_branch_termination.add(new_branch_origin.get(o) + (Math.log(Math.random())/(-(lambda+mu))));
				}

				// Add new origin and termination to the old collection.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					branch_origin.add(new_branch_origin.get(o));
				}
				for (int t = 0; t < new_branch_termination.size(); t++) {
					branch_termination.add(new_branch_termination.get(t));
				}
				
				// Empty the new origin array.
				new_branch_origin = new ArrayList<Double>();
				
				// For each new termination, add it to the new origin array if it is < max_tree_duration and rand < lambda/(lambda+mu) - this represents a speciation event.
				for (int t = 0; t < new_branch_termination.size(); t++) {										
					if (new_branch_termination.get(t) < max_tree_duration) {
						if (Math.random() < lambda/(lambda+mu)) {
							new_branch_origin.add(new_branch_termination.get(t));
							new_branch_origin.add(new_branch_termination.get(t));
						}
					}
				}
				
				// Empty the new termination array.
				new_branch_termination = new ArrayList<Double>();
				
			} // while (new_branch_origin.size() > 0)
			
			// Analyse the obtained set of branches.
			for (int i = ages.length-1; i >= 0; i--) {
				
				// For each age, determine the tree length.
				double tree_duration = ages[i] - first_occurrence_age_min;
				
				// Trim the tree so that no branches are older than ages[i] - first_occurrence_age. While trimming, make sure the tree contains at least one extant species.
				double sum_of_species_durations = 0.0;
				int number_of_extant_taxa = 0;
				for (int o = branch_origin.size()-1; o >= 0; o--) {
					boolean remove_this_branch = false;
					if (branch_origin.get(o) > tree_duration) {
						remove_this_branch = true;
					} else {
						if (branch_termination.get(o) >= tree_duration) {
							branch_termination.set(o,tree_duration);
							number_of_extant_taxa += 1;
						}
						sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
					}
					if (remove_this_branch == true) {
						branch_origin.remove(o);
						branch_termination.remove(o);
					}
				} // for (int o = branch_origin.size()-1; o >= 0; o--)
					
				// If at least one species is extant at this age, get the sum of lineage durations.
				if (number_of_extant_taxa > 0) {
					if (psi_min == psi_max) {
						raw_probabilities[i] += psi_min * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa;
					} else {
						if (sum_of_species_durations == 0) {
							raw_probabilities[i] += ((psi_max+psi_min)/2.0) * number_of_extant_taxa;
						} else {
							raw_probabilities[i] += ((-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_max*sum_of_species_durations) * number_of_extant_taxa * (1+psi_max*sum_of_species_durations))-(-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa * (1+psi_min*sum_of_species_durations)))/(psi_max-psi_min);
						}
					}
					successful_simulations[i] += 1;
				} // if (number_of_extant_taxa > 0)
				
			} // for (int i = 0; i < ages.length; i++)
			
		} // while (successful_simulations[-1] < 10000)
		
		if (cancel) {
			return null;
		}

		// Calculate the fix_probabilities (that's when first_occurrence_age_min == first_occurrence_age_max, and int_probabilities which are
		// different from fix_probabilities when this is not the case.		
		for (int i = 0; i < raw_probabilities.length; i++) {
			fix_probabilities[i] = raw_probabilities[i]/(double) successful_simulations[i];
		}
		
		// If first_occurrence_age_min == first_occurrence_age_max, fix_probabilities and int_probabilities are equal.
		if (first_occurrence_age_min == first_occurrence_age_max) {
			
			for (int i = 0; i < fix_probabilities.length; i++) {
				int_probabilities[i] = fix_probabilities[i];
			}
			
		// If not, int_probabilities have to be calculated from fix_probabilities.
		} else {
			
			int_probabilities[0] = 0.0;
			for (int a = 1; a < ages.length; a++) {
				double probability_mass = 0;
				// If the difference between first_occurrence_age_max and first_occurrence_age_min is smaller than half the difference
				// between two adjacent ages, the probability mass for a given age is simply the difference between first_occurrence_age_max
				// and first_occurrence_age_min times the probability of that age.
				if (first_occurrence_age_max-first_occurrence_age_min <= (ages[1]-ages[0])/2.0) {
					probability_mass += fix_probabilities[a] * (first_occurrence_age_max-first_occurrence_age_min);
				// If the difference between first_occurrence_age_max and first_occurrence_age_min is more than half the difference between
				// two adjacent ages (this means the integral extends over more than one segment), three steps are necessary.
				} else {
					// i.) Half of the last segment is added in all cases.
					probability_mass += fix_probabilities[a] * (ages[1]-ages[0])/2.0;
					
					// ii.) The first, or part of the first segment is added if the integral extends beyond or into it.
					if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[0]) {
						// Add half the first segment.
						probability_mass += fix_probabilities[0] * (ages[1]-ages[0])/2.0;
					} else if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[0] + (ages[1]-ages[0])/2.0) {
						// Add a proportion of the first segment.
						probability_mass += fix_probabilities[0] * ((ages[0]+(ages[1]-ages[0])/2.0) - (ages[a]-(first_occurrence_age_max-first_occurrence_age_min)));
					}
					
					// iii.) For each segment between the first and the last, the whole segment, or part of that segment are added if the integral
					// extends beyond or into it.
					for (int aa = 1; aa < a; a++) {
						if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) <= ages[aa] - (ages[1]-ages[0])/2.0) {
							// Add the full segment.
							probability_mass += fix_probabilities[aa] * (ages[1]-ages[0]);
						} else if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[aa] + (ages[1]-ages[0])/2.0) {
							// Add a proportion of this segment.
							probability_mass += fix_probabilities[aa] * ((ages[aa]+(ages[1]-ages[0])/2.0) - (ages[a]-(first_occurrence_age_max-first_occurrence_age_min)));
						}
					}
				}
				// The probability mass is divided by the difference between first_occurrence_age_max first_occurrence_age_min.
				int_probabilities[a] = probability_mass/(first_occurrence_age_max-first_occurrence_age_min);
			} // for (int a = 1; a < ages.length; a++) {
			
		} // if (first_occurrence_age_min == first_occurrence_age_max)

		// Part 3: distribution fitting

		// Set the number of Nelder Mead Downhill Simplex repetitions.
		int nmRepetitions = 5;
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] ss = new double[nmRepetitions];
		double[] ms = new double[nmRepetitions];
		double[] ws = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];

		for (int x = 0; x < nmRepetitions; x++) {
			
			if (cancel != false) {
				return null;
			}

			// Initiate the simplex, find 4 vertices.
			// vertex0
			double vertex0s = 10 + Math.random()*100;
			double vertex0m = 1 + Math.random()*5;
			double vertex0w = 0.01 + Math.random();
			double vertex0c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex0s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex0s)-vertex0m),2)/vertex0w));
			double vertex0Y = 0;
			
			// vertex1
			double vertex1s = 10 + Math.random()*100;
			double vertex1m = 1 + Math.random()*5;
			double vertex1w = 0.01 + Math.random();
			double vertex1c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex1s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex1s)-vertex1m),2)/vertex1w));
			double vertex1Y = 0;
			
			// vertex2
			double vertex2s = 10 + Math.random()*100;
			double vertex2m = 1 + Math.random()*5;
			double vertex2w = 0.01 + Math.random();
			double vertex2c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex2s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex2s)-vertex2m),2)/vertex2w));
			double vertex2Y = 0;
			
			// vertex3
			double vertex3s = 10 + Math.random()*100;
			double vertex3m = 1 + Math.random()*5;
			double vertex3w = 0.01 + Math.random();
			double vertex3c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex3s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex3s)-vertex3m),2)/vertex3w));
			double vertex3Y = 0;
			
			// Prepare for the Nelder-Mead loop.
			keepGoing = true;
			int stepCounter = 0;

			// Until converged do the loop.
			while (keepGoing == true) {
				
				if (cancel != false) {
					return null;
				}
				
				// Increment the stepCounter.
				stepCounter += 1;

				// Calculate the y value of each vertex.
				// vertex0
				vertex0Y = 0;
				for (int i = 0; i < ages.length; i++) {
					double temp = (1.0/(ages[i]-first_occurrence_age_min+vertex0s))*Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex0s)-vertex0m),2)/vertex0w);
					vertex0Y += Math.pow((fix_probabilities[i]-vertex0c*temp),2);
				}
				
				// vertex1
				vertex1Y = 0;
				for (int i = 0; i < ages.length; i++) {
					double temp = (1.0/(ages[i]-first_occurrence_age_min+vertex1s))*Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex1s)-vertex1m),2)/vertex1w);
					vertex1Y += Math.pow((fix_probabilities[i]-vertex1c*temp),2);
				}
				
				// vertex2
				vertex2Y = 0;
				for (int i = 0; i < ages.length; i++) {
					double temp = (1.0/(ages[i]-first_occurrence_age_min+vertex2s))*Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex2s)-vertex2m),2)/vertex2w);
					vertex2Y += Math.pow((fix_probabilities[i]-vertex2c*temp),2);
				}

				// vertex3
				vertex3Y = 0;
				for (int i = 0; i < ages.length; i++) {
					double temp = (1.0/(ages[i]-first_occurrence_age_min+vertex3s))*Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex3s)-vertex3m),2)/vertex3w);
					vertex3Y += Math.pow((fix_probabilities[i]-vertex3c*temp),2);
				}

				// Find the best (=lowest) y value.
				double bestY = vertex0Y;
				if (vertex1Y < bestY) {
					bestY = vertex1Y;
				}
				if (vertex2Y < bestY) {
					bestY = vertex2Y;
				}
				if (vertex3Y < bestY) {
					bestY = vertex3Y;
				}

				// Find the worst (=highest) y value.
				double worstY = vertex0Y;
				if (vertex1Y > worstY) {
					worstY = vertex1Y;
				}
				if (vertex2Y > worstY) {
					worstY = vertex2Y;
				}
				if (vertex3Y > worstY) {
					worstY = vertex3Y;
				}

				// Find the second-worst (=second-highest) y value.
				double secondWorstY = bestY;
				if (vertex0Y > secondWorstY && vertex0Y != worstY) {
					secondWorstY = vertex0Y;
				}
				if (vertex1Y > secondWorstY && vertex1Y != worstY) {
					secondWorstY = vertex1Y;
				}
				if (vertex2Y > secondWorstY && vertex2Y != worstY) {
					secondWorstY = vertex2Y;
				}
				if (vertex3Y > secondWorstY && vertex3Y != worstY) {
					secondWorstY = vertex3Y;
				}

				// Find the parameter values of the best vertex.
				double bests = 0;
				double bestm = 0;
				double bestw = 0;
				if (vertex0Y == bestY) {
					bests = vertex0s;
					bestm = vertex0m;
					bestw = vertex0w;
				} else if (vertex1Y == bestY) {
					bests = vertex1s;
					bestm = vertex1m;
					bestw = vertex1w;
				} else if (vertex2Y == bestY) {
					bests = vertex2s;
					bestm = vertex2m;
					bestw = vertex2w;
				} else if (vertex3Y == bestY) {
					bests = vertex3s;
					bestm = vertex3m;
					bestw = vertex3w;
				}

				// Find the parameter values of the worst vertex.
				double worsts = 0;
				double worstm  = 0;
				double worstw = 0;
				if (vertex0Y == worstY) {
                    worsts = vertex0s;
                    worstm = vertex0m;
                    worstw = vertex0w;
				} else if (vertex1Y == worstY) {
                    worsts = vertex1s;
                    worstm = vertex1m;
                    worstw = vertex1w;
				} else if (vertex2Y == worstY) {
                    worsts = vertex2s;
                    worstm = vertex2m;
                    worstw = vertex2w;
				} else if (vertex3Y == worstY) {
                    worsts = vertex3s;
                    worstm = vertex3m;
                    worstw = vertex3w;
				}

				// Calculate the sum of the parameters over all vertices.
                double sums = vertex0s + vertex1s + vertex2s + vertex3s;
                double summ = vertex0m + vertex1m + vertex2m + vertex3m;
                double sumw = vertex0w + vertex1w + vertex2w + vertex3w;

				// Calculate the parameter values of the centroid.
                double centroids = (sums - worsts)/3.0;
                double centroidm = (summ - worstm)/3.0;
                double centroidw = (sumw - worstw)/3.0;

				// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
                double reflections = centroids + alph * (centroids - worsts);
                if (reflections <= 0) {
                	reflections = 10 + Math.random()*100;                	
                }
                double reflectionm = centroidm + alph * (centroidm - worstm);
                if (reflectionm <= 0) {
                	reflectionm = 1 + Math.random()*5;                	
                }
                double reflectionw = centroidw + alph * (centroidw - worstw);
                if (reflectionw <= 0) {
                	reflectionw = 0.01 + Math.random();                	
                }
                double reflectionc = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+reflections))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+reflections)-reflectionm),2)/reflectionw));

                // Calculate the y value of the reflection.
				double reflectionY = 0;
				for (int i = 0; i < ages.length; i++) {
					double temp = (1.0/(ages[i]-first_occurrence_age_min+reflections))*Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+reflections)-reflectionm),2)/reflectionw);
					reflectionY += Math.pow((fix_probabilities[i]-reflectionc*temp),2);
				}

				// Consider the three cases:
				// i.)   reflection is better than all vertices.
				// ii.)  reflection is better than the second-worst vertex.
				// iii.) reflection is worse than or equally good as the second-worst vertex.

				// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
				if (reflectionY < bestY) {
					
					// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
					double extensions = centroids + gamm * (centroids - worsts);
					if (extensions <= 0) {
						extensions = 10 + Math.random()*100;
					}
					double extensionm = centroidm + gamm * (centroidm - worstm);
					if (extensionm <= 0) {
						extensionm = 1 + Math.random()*5;
					}
					double extensionw = centroidw + gamm * (centroidw - worstw);
					if (extensionw <= 0) {
						extensionw = 0.01 + Math.random();
					}
					double extensionc = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+extensions))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+extensions)-extensionm),2)/extensionw));

					// Calculate the y value of the extension.
					double extensionY = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/(ages[i]-first_occurrence_age_min+extensions))*Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+extensions)-extensionm),2)/extensionw);
						extensionY += Math.pow((fix_probabilities[i]-extensionc*temp),2);
					}

					// Figure out which values to use as replacement for the values of the worst vertex.
                    double replaces = 0;
                    double replacem = 0;
                    double replacew = 0;
                    double replacec = 0;
                    if (reflectionY < extensionY) {
                    	replaces = reflections;
                    	replacem = reflectionm;
                    	replacew = reflectionw;
                    	replacec = reflectionc;
                    } else {
                    	replaces = extensions;
                    	replacem = extensionm;
                    	replacew = extensionw;
                    	replacec = extensionc;
                    }
                    
                    // Replace the parameter values of the worst vertex with the replacement values.
                    if (vertex0Y == worstY) {
                    	vertex0s = replaces;
                        vertex0m = replacem;
                        vertex0w = replacew;
                        vertex0c = replacec;
                    } else if (vertex1Y == worstY) {
                        vertex1s = replaces;
                        vertex1m = replacem;
                        vertex1w = replacew;
                        vertex1c = replacec;
                    } else if (vertex2Y == worstY) {
                        vertex2s = replaces;
                        vertex2m = replacem;
                        vertex2w = replacew;
                        vertex2c = replacec;
                    } else if (vertex3Y == worstY) {
                        vertex3s = replaces;
                        vertex3m = replacem;
                        vertex3w = replacew;
                        vertex3c = replacec;
                    }

                // Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
				} else if (reflectionY < secondWorstY) {
                	
					if (vertex0Y == worstY) {
                        vertex0s = reflections;
                        vertex0m = reflectionm;
                        vertex0w = reflectionw;
                        vertex0c = reflectionc;
					} else if  (vertex1Y == worstY) {
                        vertex1s = reflections;
                        vertex1m = reflectionm;
                        vertex1w = reflectionw;
                        vertex1c = reflectionc;
					} else if (vertex2Y == worstY) {
                        vertex2s = reflections;
                        vertex2m = reflectionm;
                        vertex2w = reflectionw;
                        vertex2c = reflectionc;
					} else if (vertex3Y == worstY) {
                        vertex3s = reflections;
                        vertex3m = reflectionm;
                        vertex3w = reflectionw;
                        vertex3c = reflectionc;
					}
                	
				// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
                } else {
                	
                	// Calculate the contraction.
                    double contractions = centroids + beta * (centroids - worsts);
                    double contractionm = centroidm + beta * (centroidm - worstm);
                    double contractionw = centroidw + beta * (centroidw - worstw);
                    double contractionc = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+contractions))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+contractions)-contractionm),2)/contractionw));

                    // Calculate the y value of the contraction.
    				double contractionY = 0;
    				for (int i = 0; i < ages.length; i++) {
    					double temp = (1.0/(ages[i]-first_occurrence_age_min+contractions))*Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+contractions)-contractionm),2)/contractionw);
    					contractionY += Math.pow((fix_probabilities[i]-contractionc*temp),2);
    				}
    				
    				// Consider two subcases of case iii:
    				// iiia) The contraction is better than the worst vertex.
    				// iiib) The contraction is not better than the worst vertex.
    				
    				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
					if (contractionY < worstY) {
						
						// Replace the parameter values of the worst vertex with the contraction values.
						if (vertex0Y == worstY) {
                            vertex0s = contractions;
                            vertex0m = contractionm;
                            vertex0w = contractionw;
                            vertex0c = contractionc;
						} else if (vertex1Y == worstY) {
                            vertex1s = contractions;
                            vertex1m = contractionm;
                            vertex1w = contractionw;
                            vertex1c = contractionc;
						} else if (vertex2Y == worstY) {
                            vertex2s = contractions;
                            vertex2m = contractionm;
                            vertex2w = contractionw;
                            vertex2c = contractionc;
						} else if (vertex3Y == worstY) {
                            vertex3s = contractions;
                            vertex3m = contractionm;
                            vertex3w = contractionw;
                            vertex3c = contractionc;
						}

					// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
					} else {
						
						// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that delt  = 0.5, the default setting)
						// vertex0
                        vertex0s = bests + delt * (vertex0s - bests);
                        vertex0m = bestm + delt * (vertex0m - bestm);
                        vertex0w = bestw + delt * (vertex0w - bestw);
                        vertex0c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex0s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex0s)-vertex0m),2)/vertex0w));

						// vertex1
                        vertex1s = bests + delt * (vertex1s - bests);
                        vertex1m = bestm + delt * (vertex1m - bestm);
                        vertex1w = bestw + delt * (vertex1w - bestw);
                        vertex1c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex1s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex1s)-vertex1m),2)/vertex1w));
						
						// vertex2
                        vertex2s = bests + delt * (vertex2s - bests);
                        vertex2m = bestm + delt * (vertex2m - bestm);
                        vertex2w = bestw + delt * (vertex2w - bestw);
                        vertex2c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex2s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex2s)-vertex2m),2)/vertex2w));
						
						// vertex3
                        vertex3s = bests + delt * (vertex3s - bests);
                        vertex3m = bestm + delt * (vertex3m - bestm);
                        vertex3w = bestw + delt * (vertex3w - bestw);
                        vertex3c = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+vertex3s))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+vertex3s)-vertex3m),2)/vertex3w));

					} // if (contractionY < worstY)
				
				} // if (reflectionY < bestY)
				
				// Stop the loop when all parameter values are identical in the first 10 decimals.
				keepGoing = false;
				if (Math.abs(vertex0s - vertex1s) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0m - vertex1m) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0w - vertex1w) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0s - vertex2s) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0m - vertex2m) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0w - vertex2w) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0s - vertex3s) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0m - vertex3m) > 0.0000001) {
					keepGoing = true;
				} else if (Math.abs(vertex0w - vertex3w) > 0.0000001) {
					keepGoing = true;
				}
				if (stepCounter > 1000) {
					keepGoing = false;
				}

			} // while (keepGoing == true)

			// Memorize all parameters and y.
			ss[x] = vertex0s;
            ms[x] = vertex0m;
            ws[x] = vertex0w;
            ys[x] = vertex0Y;
			
		} // for (int x = 0; x < nmRepetitions; x++)
		
		// Find the best result among the nmReplicates replicates.
		int index = 0;
		for (int x = 1; x < nmRepetitions; x++) {
			if (ys[x] < ys[index]) {
				index = x;
			}
		}
        double fittedCladeAgeS = ss[index];
        double fittedCladeAgeM = ms[index];
        double fittedCladeAgeW = ws[index];
        //double fittedCladeAgeC = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+fittedCladeAgeS))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+fittedCladeAgeS)-fittedCladeAgeM),2)/fittedCladeAgeW));
        double fittedCladeAgeRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));
		
		// Variable cladeAgeC, as determined above is adjusted to provide good fit with the empirical probability values, however, it may
		// not be the right scale factor for a distribution, which must sum to 1. Therefore, the integral under the distribution is calculated
		// between fossil_age_min and Inf, and a scaling factor cladeAgeCcorr is determined so that the distribution sums to 1.
		double erf_argument = -(Math.log(fittedCladeAgeS)-fittedCladeAgeM) / Math.sqrt(fittedCladeAgeW);
		double erf_t = 1/(1.0+0.5*Math.abs(erf_argument));
		double erf_polynomial = -Math.pow(erf_argument,2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t,2) + 0.09678418*Math.pow(erf_t,3) - 0.18628806*Math.pow(erf_t,4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
		double erf_tau = erf_t * Math.exp(erf_polynomial);
		double erf = 0;
		if (erf_argument >= 0) {
			erf = 1 - erf_tau;
		} else {
			erf = erf_tau - 1;
		}
		double fittedCladeAgeCcorr = 1/(0.5 * Math.sqrt(Math.PI) * Math.sqrt(fittedCladeAgeW) * (1+erf));

		// Return the fitted CladeAge distribution.
		return new FittedCladeAgeDistribution(first_occurrence_age_min, first_occurrence_age_max, fittedCladeAgeCcorr, fittedCladeAgeS, fittedCladeAgeM, fittedCladeAgeW, fittedCladeAgeRmsd);
		
	} // public ContinuousDistribution run_fitted_cladeage(...)
	
	public FittedCladeAgeDistribution run_fitted_cladeage_star(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, JProgressBar dpb) {

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		fix_probabilities = new double[number_of_ages + 1];
		int_probabilities = new double[number_of_ages + 1];
		double[] raw_probabilities = new double[number_of_ages + 1];

		// Part 1: Calculate the maximum age, and 100 time points between fossil age and maximum age.
		
		// Roughly calculate a maximum tree duration so that the probability for the maximum tree duration is 0.001 * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		double ndr_mean = (ndr_max+ndr_min)/2.0;
		double psi_mean = (psi_max+psi_min)/2.0;
		
		// Initiate the simplex, find 2 vertices.
		// vertex0
		// Math.log(ndr_mean/psi_mean)/ndr_mean is the peak of p(t), as determined through its first derivative.
		double vertex0t = 0;
		if (Math.log(ndr_mean/psi_mean)/ndr_mean > 0) {
			vertex0t = Math.log(ndr_mean/psi_mean)/ndr_mean;
		}
		
		// vertex1
		double vertex1t = vertex0t + 10;
		
		// Prepare for the Nelder-Mead loop.
		boolean keepGoing = true;
		
		// Until converged do the loop.
		while (keepGoing == true) {
			
			// Calculate the y value of each vertex.
			// vertex0
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(0.001)), 2);

			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Find the best (=lowest) and the worst (=highest) y value.
			double bestY = vertex0Y;
			double bestt = vertex0t;
			double worstY = vertex1Y;
			double worstt = vertex1t;
			if (vertex1Y < bestY) {
				bestY = vertex1Y;
				bestt = vertex1t;
				worstY = vertex0Y;
				worstt = vertex0t;
			}

			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
			double reflectiont = bestt + alph * (bestt - worstt);
			
			// Calculate the y value of the reflection.
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Figure out which values to use as replacement for the values of the worst vertex.
				double replacet = 0;
				if (reflectionY < extensionY) {
					replacet = reflectiont;
				} else {
					replacet = extensiont;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0t = replacet;
				} else if (vertex1Y == worstY) {
					vertex1t = replacet;
				}

			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractiont = bestt + beta * (bestt - worstt);
				
				// Calculate the y value of the contraction.
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0t = contractiont;
					} else if (vertex1Y == worstY) {
						vertex1t = contractiont;
					}
				
				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
					// vertex0
					vertex0t = bestt + delt  * (vertex0t - bestt);
					
					// vertex1
					vertex1t = bestt + delt  * (vertex1t - bestt);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY) {
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0t - vertex1t) > 0.0001) {
				keepGoing = true;
			} 

		} // while (keepGoing == true) {

		// Memorize the result as the max_simulation_age.
		// Here, this is added to first_occurrence_age_min, however in part 3, first_occurrence_age_min is ignored.
		double max_simulation_age = vertex0t + first_occurrence_age_min;
		
		// Determine the ages.
		ages[0] = first_occurrence_age_min;
		for (int i = 0; i < number_of_ages; i++) {
			ages[i+1] = first_occurrence_age_min + ((i+1)/(double) number_of_ages)*(max_simulation_age-first_occurrence_age_min);
		}

		// Part 2: Estimate clade age probability densities.

		// Prepare the arrays raw_probabilities and successful simulations (those in which at least one taxon survives), and fill them with 0.
		int[] successful_simulations = new int[number_of_ages];
		
		while (successful_simulations[-1] < 10000) {

			if (cancel) {
				return null;
			}

			// Update the progress bar.
			if (dpb != null) {
				dpb.setValue(successful_simulations[-1]);
			}

			// Draw values for parameters ndr (net diversification rate, lambda - mu) and epsilon (turnover rate, mu/lambda) from uniform distributions, and calculate lambda and mu from them.
			double ndr = ndr_min + Math.random() * (ndr_max-ndr_min);
			double epsilon = epsilon_min + Math.random() * (epsilon_max-epsilon_min);
			double mu = (ndr*epsilon)/(1.0 - epsilon);
			double lambda = ndr + mu;
			
			// Determine the maximum tree duration.
			double max_tree_duration = ages[-1] - first_occurrence_age_min;
			
			// Initiate arrays for branch origin and branch termination.
			ArrayList<Double> branch_origin = new ArrayList<Double>();
			ArrayList<Double> branch_termination = new ArrayList<Double>();
			ArrayList<Double> new_branch_origin = new ArrayList<Double>();
			ArrayList<Double> new_branch_termination = new ArrayList<Double>();
			new_branch_origin.add(0.0);

			// Start the tree generation loop.
			while (new_branch_origin.size() > 0) {

				// For each new origin, add a new termination.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					new_branch_termination.add(new_branch_origin.get(o) + (Math.log(Math.random())/(-(lambda+mu))));
				}

				// Add new origin and termination to the old collection.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					branch_origin.add(new_branch_origin.get(o));
				}
				for (int t = 0; t < new_branch_termination.size(); t++) {
					branch_termination.add(new_branch_termination.get(t));
				}
				
				// Empty the new origin array.
				new_branch_origin = new ArrayList<Double>();
				
				// For each new termination, add it to the new origin array if it is < max_tree_duration and rand < lambda/(lambda+mu) - this represents a speciation event.
				for (int t = 0; t < new_branch_termination.size(); t++) {										
					if (new_branch_termination.get(t) < max_tree_duration) {
						if (Math.random() < lambda/(lambda+mu)) {
							new_branch_origin.add(new_branch_termination.get(t));
							new_branch_origin.add(new_branch_termination.get(t));
						}
					}
				}
				
				// Empty the new termination array.
				new_branch_termination = new ArrayList<Double>();
				
			} // while (new_branch_origin.size() > 0)
			
			// Analyse the obtained set of branches.
			for (int i = ages.length-1; i >= 0; i--) {
				
				// For each age, determine the tree length.
				double tree_duration = ages[i] - first_occurrence_age_min;
				
				// Trim the tree so that no branches are older than ages[i] - first_occurrence_age. While trimming, make sure the tree contains at least one extant species.
				double sum_of_species_durations = 0.0;
				int number_of_extant_taxa = 0;
				for (int o = branch_origin.size()-1; o >= 0; o--) {
					boolean remove_this_branch = false;
					if (branch_origin.get(o) > tree_duration) {
						remove_this_branch = true;
					} else {
						if (branch_termination.get(o) >= tree_duration) {
							branch_termination.set(o,tree_duration);
							number_of_extant_taxa += 1;
						}
						sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
					}
					if (remove_this_branch == true) {
						branch_origin.remove(o);
						branch_termination.remove(o);
					}
				} // for (int o = branch_origin.size()-1; o >= 0; o--)
					
				// If at least one species is extant at this age, get the sum of lineage durations.
				if (number_of_extant_taxa > 0) {
					if (psi_min == psi_max) {
						raw_probabilities[i] += psi_min * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa;
					} else {
						if (sum_of_species_durations == 0) {
							raw_probabilities[i] += ((psi_max+psi_min)/2.0) * number_of_extant_taxa;
						} else {
							raw_probabilities[i] += ((-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_max*sum_of_species_durations) * number_of_extant_taxa * (1+psi_max*sum_of_species_durations))-(-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa * (1+psi_min*sum_of_species_durations)))/(psi_max-psi_min);
						}
					}
					successful_simulations[i] += 1;
				} // if (number_of_extant_taxa > 0)
				
			} // for (int i = 0; i < ages.length; i++)
			
		} // while (successful_simulations[-1] < 10000)
		
		if (cancel) {
			return null;
		}

		// Calculate the fix_probabilities (that's when first_occurrence_age_min == first_occurrence_age_max, and int_probabilities which are
		// different from fix_probabilities when this is not the case.		
		for (int i = 0; i < raw_probabilities.length; i++) {
			fix_probabilities[i] = raw_probabilities[i]/(double) successful_simulations[i];
		}
		
		// If first_occurrence_age_min == first_occurrence_age_max, fix_probabilities and int_probabilities are equal.
		if (first_occurrence_age_min == first_occurrence_age_max) {
			
			for (int i = 0; i < fix_probabilities.length; i++) {
				int_probabilities[i] = fix_probabilities[i];
			}
			
		// If not, int_probabilities have to be calculated from fix_probabilities.
		} else {
			
			int_probabilities[0] = 0.0;
			for (int a = 1; a < ages.length; a++) {
				double probability_mass = 0;
				// If the difference between first_occurrence_age_max and first_occurrence_age_min is smaller than half the difference
				// between two adjacent ages, the probability mass for a given age is simply the difference between first_occurrence_age_max
				// and first_occurrence_age_min times the probability of that age.
				if (first_occurrence_age_max-first_occurrence_age_min <= (ages[1]-ages[0])/2.0) {
					probability_mass += fix_probabilities[a] * (first_occurrence_age_max-first_occurrence_age_min);
				// If the difference between first_occurrence_age_max and first_occurrence_age_min is more than half the difference between
				// two adjacent ages (this means the integral extends over more than one segment), three steps are necessary.
				} else {
					// i.) Half of the last segment is added in all cases.
					probability_mass += fix_probabilities[a] * (ages[1]-ages[0])/2.0;
					
					// ii.) The first, or part of the first segment is added if the integral extends beyond or into it.
					if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[0]) {
						// Add half the first segment.
						probability_mass += fix_probabilities[0] * (ages[1]-ages[0])/2.0;
					} else if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[0] + (ages[1]-ages[0])/2.0) {
						// Add a proportion of the first segment.
						probability_mass += fix_probabilities[0] * ((ages[0]+(ages[1]-ages[0])/2.0) - (ages[a]-(first_occurrence_age_max-first_occurrence_age_min)));
					}
					
					// iii.) For each segment between the first and the last, the whole segment, or part of that segment are added if the integral
					// extends beyond or into it.
					for (int aa = 1; aa < a; a++) {
						if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) <= ages[aa] - (ages[1]-ages[0])/2.0) {
							// Add the full segment.
							probability_mass += fix_probabilities[aa] * (ages[1]-ages[0]);
						} else if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[aa] + (ages[1]-ages[0])/2.0) {
							// Add a proportion of this segment.
							probability_mass += fix_probabilities[aa] * ((ages[aa]+(ages[1]-ages[0])/2.0) - (ages[a]-(first_occurrence_age_max-first_occurrence_age_min)));
						}
					}
				}
				// The probability mass is divided by the difference between first_occurrence_age_max first_occurrence_age_min.
				int_probabilities[a] = probability_mass/(first_occurrence_age_max-first_occurrence_age_min);
			} // for (int a = 1; a < ages.length; a++) {
			
		} // if (first_occurrence_age_min == first_occurrence_age_max)

		// Part 1 & 2 are only included to visualize the fit of the distribution that is analytically calculated in parts 3-5.
		// Part 3: Calculate the maximum tree duration, which will be used in Part 4 to determine only 4 time points
		
		if (cancel) {
			return null;
		}
		
		// Roughly calculate a maximum tree duration so that the probability for the maximum tree duration is 0.001 * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		ndr_mean = (ndr_max+ndr_min)/2.0;
		psi_mean = (psi_max+psi_min)/2.0;
		
		// Initiate the simplex, find 2 vertices.
		// vertex0
		// Math.log(ndr_mean/psi_mean)/ndr_mean is the peak of p(t), as determined through its first derivative.
		vertex0t = 0;
		if (Math.log(ndr_mean/psi_mean)/ndr_mean > 0) {
			vertex0t = Math.log(ndr_mean/psi_mean)/ndr_mean;
		}
		
		// vertex1
		vertex1t = vertex0t + 10;
		
		// Prepare for the Nelder-Mead loop.
		keepGoing = true;
		
		// Until converged do the loop.
		while (keepGoing == true) {
			
			// Calculate the y value of each vertex.
			// vertex0
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(0.001)), 2);

			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Find the best (=lowest) and the worst (=highest) y value, and the respective t value.
			double bestY = vertex0Y;
			double bestt = vertex0t;
			double worstY = vertex1Y;
			double worstt = vertex1t;
			if (vertex1Y < bestY) {
				bestY = vertex1Y;
				bestt = vertex1t;
				worstY = vertex0Y;
				worstt = vertex0t;
			}

			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
			double reflectiont = bestt + alph * (bestt - worstt);
			
			// Calculate the y value of the reflection.
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Figure out which values to use as replacement for the values of the worst vertex.
				double replacet = 0;
				if (reflectionY < extensionY) {
					replacet = reflectiont;
				} else {
					replacet = extensiont;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0t = replacet;
				} else if (vertex1Y == worstY) {
					vertex1t = replacet;
				}

			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractiont = bestt + beta * (bestt - worstt);
				
				// Calculate the y value of the contraction.
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0t = contractiont;
					} else if (vertex1Y == worstY) {
						vertex1t = contractiont;
					}
				
				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
					// vertex0
					vertex0t = bestt + delt  * (vertex0t - bestt);
					
					// vertex1
					vertex1t = bestt + delt  * (vertex1t - bestt);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY) {
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0t - vertex1t) > 0.0001) {
				keepGoing = true;
			} 

		} // while (keepGoing == true) {

		// Memorize the result as the max_simulation_age
		double max_tree_duration = vertex0t;

		// Part 4: Select four time points and estimate clade age probabilities for these time points.
		
		// Set the number of replicates for both ndr and epsilon variation.
		int replicates = 30;
		
		// Specify the values of ndr to be used for simulations.
		double[] ndrs = new double[replicates];
		if (ndr_min == ndr_max) {
			for (int x = 0; x < replicates; x++) {
				ndrs[x] = ndr_min; 
			}
		} else {
			for (int x = 0; x < replicates; x++) {
				ndrs[x] = ndr_min + (ndr_max-ndr_min)/(double) (2*replicates) + (ndr_max-ndr_min)*(x/(double) replicates); 
			}			
		}
		
		// Specify the values of epsilon to be used for simulations.
		double[] epsilons = new double[replicates];
		if (epsilon_min == epsilon_max) {
			for (int x = 0; x < replicates; x++) {
				epsilons[x] = epsilon_min; 
			}
		} else {
			for (int x = 0; x < replicates; x++) {
				epsilons[x] = epsilon_min + (epsilon_max-epsilon_min)/(double) (2*replicates) + (epsilon_max-epsilon_min)*(x/(double) replicates); 
			}			
		}
		
		// Initiate arrays xs and ys.
		double[] xs = new double[4];
		xs[0] = 0.0;
		xs[1] = 0.2*max_tree_duration;
		xs[2] = 0.4*max_tree_duration;
		xs[3] = 0.8*max_tree_duration;
		double[] ys = new double[4];
		ys[0] = (psi_min+psi_max)/2.0;
		
		// Start tree simulation loop for different values of x, ndr, and epsilon.
		for (int z = 1; z < 4; z++) {
			double x = xs[z];
			double y = 0;
			for (int zz = 0; zz < replicates; zz++) {
				double ndr = ndrs[zz];
				for (int zzz = 0; zzz < replicates; zzz++) {
					double epsilon = epsilons[zzz];
					
					// Calculate values of lambda and mu from the net diversification rate and the turnover.
					double mu = (ndr*epsilon)/(1.0 - epsilon);
					double lambda = ndr + mu;
					
					// Set success to false
					boolean success = false;
					
					// Start the tree generation loop with these settings for x, lambda, and mu.
					while (success == false) {
						
						// Initiate arrays for branch origin and branch termination.
						ArrayList<Double> branch_origin = new ArrayList<Double>();
						ArrayList<Double> branch_termination = new ArrayList<Double>();
						ArrayList<Double> new_branch_origin = new ArrayList<Double>();
						ArrayList<Double> new_branch_termination = new ArrayList<Double>();
						new_branch_origin.add(0.0);

						// Start the tree generation loop.
						while (new_branch_origin.size() > 0) {
							
							// For each new origin, add a new termination.
							for (int o = 0; o < new_branch_origin.size(); o++) {
								new_branch_termination.add(new_branch_origin.get(o) + (Math.log(Math.random())/(-(lambda+mu))));
							}

							// Add new origin and termination to the old collection.
							for (int o = 0; o < new_branch_origin.size(); o++) {
								branch_origin.add(new_branch_origin.get(o));
							}
							for (int t = 0; t < new_branch_termination.size(); t++) {
								branch_termination.add(new_branch_termination.get(t));
							}

							// Empty the new origin array.
							new_branch_origin = new ArrayList<Double>();

							// For each new termination, add it to the new origin array if it is < x and rand < lambda/(lambda+mu) - this represents a speciation event.
							for (int t = 0; t < new_branch_termination.size(); t++) {										
								if (new_branch_termination.get(t) < x) {
									if (Math.random() < lambda/(lambda+mu)) {
										new_branch_origin.add(new_branch_termination.get(t));
										new_branch_origin.add(new_branch_termination.get(t));
									}
								}
							}
							
							// Empty the new termination array.
							new_branch_termination = new ArrayList<Double>();
							
						} // while (new_branch_origin.size() > 0)
						
						// Calculate the sum of species durations and the number of extant taxa.
						double sum_of_species_durations = 0.0;
						int number_of_extant_taxa = 0;
						for (int o = branch_origin.size()-1; o >= 0; o--) {
							if (branch_termination.get(o) >= x) {
								branch_termination.set(o,x);
								number_of_extant_taxa += 1;
							}
							sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
						} // for (int o = branch_origin.size()-1; o >= 0; o--)

						// Calculate the probability density of first observing a fossil at this x.
						if (number_of_extant_taxa > 0) {
							if (psi_min == psi_max) {
								y += (psi_min * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa)/(double) Math.pow(replicates,2);
							} else {
								if (sum_of_species_durations == 0) {
									y += (((psi_max+psi_min)/2.0) * number_of_extant_taxa)/(double) Math.pow(replicates,2);
								} else {
									y += (((-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_max*sum_of_species_durations) * number_of_extant_taxa * (1+psi_max*sum_of_species_durations))-(-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa * (1+psi_min*sum_of_species_durations)))/(psi_max-psi_min))/(double) Math.pow(replicates,2);
								}
							}
							success = true;
						} // if (number_of_extant_taxa > 0)
						
					} // while (success == false)
										
				} // for (int zzz = 0; zzz < replicates; zzz++)
			} // for (int zz = 0; zz < replicates; zz++)
			
			ys[z] = y;
			
		} // for (int z = 1; z < 4; z++)

		// Part 5: Analytical distribution fitting.
		
		if (cancel) {
			return null;
		}

		// Now that probability densities have been calculated, the CladeAge distribution is analytically calculated to fit these.
        // S is the only parameter of the CladeAge distribution which can not be solved analytically, it has to be solved numerically with the following
        // Nelder-Mead Downhill Simplex.
		
		// For easier readability, copy values for x and y.
        double x1 = xs[1];
        double y0 = ys[0];
        double y1 = ys[1];
        double y2 = ys[2];
        double y3 = ys[3];

		// Initiate the simplex, find 2 vertices.
        // vertex0
        double vertex0s = 1.0;
        
        // vertex1
        double vertex1s = 10.0;
        
        // Prepare for the Nelder-Mead loop.
        keepGoing = true;
		
        // Until converged do the loop.
		while (keepGoing == true) {

			// Calculate the y value of each vertex.
			// vertex0
			double lhs_denom1denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(2*x1+vertex0s)),2);
			double lhs_denom1nom = Math.log((y2/y0)*((2*x1+vertex0s)/vertex0s));
			double lhs_denom2denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(x1+vertex0s)),2);
			double lhs_denom2nom = Math.log((y1/y0)*((x1+vertex0s)/vertex0s));
			double lhs_nom1denom = 2*Math.log((x1+vertex0s)/vertex0s);
			double lhs_nom1nom = lhs_denom2nom;
			double lhs_nom2denom = 2*Math.log((2*x1+vertex0s)/vertex0s);
			double lhs_nom2_nom = lhs_denom1nom;
			double lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
			double rhs_denom1denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(4*x1+vertex0s)),2);
			double rhs_denom1nom = Math.log((y3/y0)*((4*x1+vertex0s)/vertex0s));
			double rhs_denom2denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(x1+vertex0s)),2);
			double rhs_denom2nom = Math.log((y1/y0)*((x1+vertex0s)/vertex0s));
			double rhs_nom1denom = 2*Math.log((x1+vertex0s)/vertex0s);
			double rhs_nom1nom = rhs_denom2nom;
			double rhs_nom2denom = 2*Math.log((4*x1+vertex0s)/vertex0s);
			double rhs_nom2_nom = rhs_denom1nom;
			double rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
			double vertex0Y = Math.pow((lhs-rhs),2);

			// vertex1
			lhs_denom1denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(2*x1+vertex1s)),2);
			lhs_denom1nom = Math.log((y2/y0)*((2*x1+vertex1s)/vertex1s));
			lhs_denom2denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(x1+vertex1s)),2);
			lhs_denom2nom = Math.log((y1/y0)*((x1+vertex1s)/vertex1s));
			lhs_nom1denom = 2*Math.log((x1+vertex1s)/vertex1s);
			lhs_nom1nom = lhs_denom2nom;
			lhs_nom2denom = 2*Math.log((2*x1+vertex1s)/vertex1s);
			lhs_nom2_nom = lhs_denom1nom;
			lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
			rhs_denom1denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(4*x1+vertex1s)),2);
			rhs_denom1nom = Math.log((y3/y0)*((4*x1+vertex1s)/vertex1s));
			rhs_denom2denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(x1+vertex1s)),2);
			rhs_denom2nom = Math.log((y1/y0)*((x1+vertex1s)/vertex1s));
			rhs_nom1denom = 2*Math.log((x1+vertex1s)/vertex1s);
			rhs_nom1nom = rhs_denom2nom;
			rhs_nom2denom = 2*Math.log((4*x1+vertex1s)/vertex1s);
			rhs_nom2_nom = rhs_denom1nom;
			rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
			double vertex1Y = Math.pow((lhs-rhs),2);

			// Find the best (=lowest) y value.
			double bestY = vertex0Y;
			double worstY = vertex1Y;
			double bests = vertex0s;
			double worsts = vertex1s;
			if (vertex1Y < vertex0Y) {
				bestY = vertex1Y;
				worstY = vertex0Y;
				bests = vertex1s;
				worsts = vertex0s;
			}
			
			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
			double reflections = bests + alph * (bests - worsts);
			
			// Calculate the y value of the reflection.
			lhs_denom1denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(2*x1+reflections)),2);
			lhs_denom1nom = Math.log((y2/y0)*((2*x1+reflections)/reflections));
			lhs_denom2denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(x1+reflections)),2);
			lhs_denom2nom = Math.log((y1/y0)*((x1+reflections)/reflections));
			lhs_nom1denom = 2*Math.log((x1+reflections)/reflections);
			lhs_nom1nom = lhs_denom2nom;
			lhs_nom2denom = 2*Math.log((2*x1+reflections)/reflections);
			lhs_nom2_nom = lhs_denom1nom;
			lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
			rhs_denom1denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(4*x1+reflections)),2);
			rhs_denom1nom = Math.log((y3/y0)*((4*x1+reflections)/reflections));
			rhs_denom2denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(x1+reflections)),2);
			rhs_denom2nom = Math.log((y1/y0)*((x1+reflections)/reflections));
			rhs_nom1denom = 2*Math.log((x1+reflections)/reflections);
			rhs_nom1nom = rhs_denom2nom;
			rhs_nom2denom = 2*Math.log((4*x1+reflections)/reflections);
			rhs_nom2_nom = rhs_denom1nom;
			rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
			double reflectionY = Math.pow((lhs-rhs),2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
				double extensions = bests + gamm * (bests - worsts);
				
				// Calculate the y value of the extension.
				lhs_denom1denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(2*x1+extensions)),2);
				lhs_denom1nom = Math.log((y2/y0)*((2*x1+extensions)/extensions));
				lhs_denom2denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(x1+extensions)),2);
				lhs_denom2nom = Math.log((y1/y0)*((x1+extensions)/extensions));
				lhs_nom1denom = 2*Math.log((x1+extensions)/extensions);
				lhs_nom1nom = lhs_denom2nom;
				lhs_nom2denom = 2*Math.log((2*x1+extensions)/extensions);
				lhs_nom2_nom = lhs_denom1nom;
				lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
				rhs_denom1denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(4*x1+extensions)),2);
				rhs_denom1nom = Math.log((y3/y0)*((4*x1+extensions)/extensions));
				rhs_denom2denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(x1+extensions)),2);
				rhs_denom2nom = Math.log((y1/y0)*((x1+extensions)/extensions));
				rhs_nom1denom = 2*Math.log((x1+extensions)/extensions);
				rhs_nom1nom = rhs_denom2nom;
				rhs_nom2denom = 2*Math.log((4*x1+extensions)/extensions);
				rhs_nom2_nom = rhs_denom1nom;
				rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
				double extensionY = Math.pow((lhs-rhs),2);

				// Figure out which values to use as replacement for the values of the worst vertex.
				double replaces = 0;
				if (reflectionY < extensionY) {
					replaces = reflections;
				} else {
					replaces = extensions;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0s = replaces;
				} else {
					vertex1s = replaces;
				}
				
			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractions = bests + beta * (bests - worsts);
				
				// Calculate the y value of the contraction.
				lhs_denom1denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(2*x1+contractions)),2);
				lhs_denom1nom = Math.log((y2/y0)*((2*x1+contractions)/contractions));
				lhs_denom2denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(x1+contractions)),2);
				lhs_denom2nom = Math.log((y1/y0)*((x1+contractions)/contractions));
				lhs_nom1denom = 2*Math.log((x1+contractions)/contractions);
				lhs_nom1nom = lhs_denom2nom;
				lhs_nom2denom = 2*Math.log((2*x1+contractions)/contractions);
				lhs_nom2_nom = lhs_denom1nom;
				lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
				rhs_denom1denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(4*x1+contractions)),2);
				rhs_denom1nom = Math.log((y3/y0)*((4*x1+contractions)/contractions));
				rhs_denom2denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(x1+contractions)),2);
				rhs_denom2nom = Math.log((y1/y0)*((x1+contractions)/contractions));
				rhs_nom1denom = 2*Math.log((x1+contractions)/contractions);
				rhs_nom1nom = rhs_denom2nom;
				rhs_nom2denom = 2*Math.log((4*x1+contractions)/contractions);
				rhs_nom2_nom = rhs_denom1nom;
				rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
				double contractionY = Math.pow((lhs-rhs),2);

				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0s = contractions;
					} else {
						vertex1s = contractions;
					}

				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
					// vertex0
					vertex0s = bests + delt * (vertex0s - bests);
					
					// vertex1
					vertex1s = bests + delt * (vertex1s - bests);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY)
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0s - vertex1s) > 0.0001) {
				keepGoing = true;
			}
			 
		} // while (keepGoing == true)

		double fittedCladeAgeS = vertex0s;
		double m_denom1denom = Math.pow((Math.log(fittedCladeAgeS)),2) - Math.pow((Math.log(2*x1+fittedCladeAgeS)),2);
		double m_denom1nom = Math.log((y2/y0)*((2*x1+fittedCladeAgeS)/fittedCladeAgeS));
		double m_denom2denom = Math.pow((Math.log(fittedCladeAgeS)),2) - Math.pow((Math.log(x1+fittedCladeAgeS)),2);
		double m_denom2nom = Math.log((y1/y0)*((x1+fittedCladeAgeS)/fittedCladeAgeS));
		double m_nom1denom = 2*Math.log((x1+fittedCladeAgeS)/fittedCladeAgeS);
		double m_nom1nom = m_denom2nom;
		double m_nom2denom = 2*Math.log((2*x1+fittedCladeAgeS)/fittedCladeAgeS);
		double m_nom2_nom = m_denom1nom;
		double fittedCladeAgeM = ((m_denom1denom/m_denom1nom)-(m_denom2denom/m_denom2nom))/((m_nom1denom/m_nom1nom)-(m_nom2denom/m_nom2_nom));
		double fittedCladeAgeW = (Math.pow((Math.log(fittedCladeAgeS)-fittedCladeAgeM),2) - Math.pow((Math.log(x1+fittedCladeAgeS)-fittedCladeAgeM),2)) / Math.log((y1/y0)*((x1+fittedCladeAgeS)/fittedCladeAgeS));
		double fittedCladeAgeC = y0 * fittedCladeAgeS * Math.exp(Math.pow((Math.log(fittedCladeAgeS)-fittedCladeAgeM),2)/fittedCladeAgeW);
		// Calculate the root mean square deviation (in order to be fair with other distributions, do that based on the 101 original ages. This can be removed if the fast distribution calculation is required.
		double sd = 0;
		for (int a = 0; a < ages.length; a++) {
			sd += Math.pow((fix_probabilities[a]-((fittedCladeAgeC/(ages[a]-first_occurrence_age_min+fittedCladeAgeS)) * Math.exp( -Math.pow((Math.log(ages[a]-first_occurrence_age_min+fittedCladeAgeS)-fittedCladeAgeM),2)/fittedCladeAgeW ))),2);
		}
		double msd = sd/(double) ages.length;
		double fittedCladeAgeRmsd = Math.sqrt(msd);
		
		// Variable fittedCladeAgeC, as determined above has been adjusted to provide good fit with the empirical probability values, 
		// however, it may not be the right scale factor for a distribution, which must sum to 1. Therefore, the integral under
		// the distribution is calculated between fossil_age_min and Inf, and a scaling factor @fittedCladeAgeCcorr is determined so that
		// the distribution sums to 1.
		double erf_argument = -(Math.log(fittedCladeAgeS)-fittedCladeAgeM) / Math.sqrt(fittedCladeAgeW);
		double erf_t = 1/(1.0+0.5*Math.abs(erf_argument));
		double erf_polynomial = -Math.pow(erf_argument,2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t,2) + 0.09678418*Math.pow(erf_t,3) - 0.18628806*Math.pow(erf_t,4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
		double erf_tau = erf_t * Math.exp(erf_polynomial);
		double erf = 0;
		if (erf_argument >= 0) {
			erf = 1 - erf_tau;
		} else {
			erf = erf_tau - 1;
		}
		double fittedCladeAgeCcorr = 1/(0.5 * Math.sqrt(Math.PI) * Math.sqrt(fittedCladeAgeW) * (1+erf));
		
		// Return the fitted CladeAge distribution.
		return new FittedCladeAgeDistribution(first_occurrence_age_min, first_occurrence_age_max, fittedCladeAgeCcorr, fittedCladeAgeS, fittedCladeAgeM, fittedCladeAgeW, fittedCladeAgeRmsd);
		
	} // public fittedCladeAge run_fitted_cladeage_star(...)

	public FittedCladeAgeDistribution run_fitted_cladeage_rapid(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max) {

		// Part 1: Calculate the maximum tree duration, which will be used in Part 4 to determine only 4 time points.
		// This is equivalent to part 3 of function run_fitted_cladeage_star(...).
		
		// Roughly calculate a maximum tree duration so that the probability for the maximum tree duration is 0.001 * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		double ndr_mean = (ndr_max+ndr_min)/2.0;
		double psi_mean = (psi_max+psi_min)/2.0;
		
		// Initiate the simplex, find 2 vertices.
		// vertex0
		// Math.log(ndr_mean/psi_mean)/ndr_mean is the peak of p(t), as determined through its first derivative.
		double vertex0t = 0;
		if (Math.log(ndr_mean/psi_mean)/ndr_mean > 0) {
			vertex0t = Math.log(ndr_mean/psi_mean)/ndr_mean;
		}
		
		// vertex1
		double vertex1t = vertex0t + 10;
		
		// Prepare for the Nelder-Mead loop.
		boolean keepGoing = true;
		
		// Until converged do the loop.
		while (keepGoing == true) {
			
			// Calculate the y value of each vertex.
			// vertex0
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(0.001)), 2);

			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Find the best (=lowest) and the worst (=highest) y value, and the respective t value.
			double bestY = vertex0Y;
			double bestt = vertex0t;
			double worstY = vertex1Y;
			double worstt = vertex1t;
			if (vertex1Y < bestY) {
				bestY = vertex1Y;
				bestt = vertex1t;
				worstY = vertex0Y;
				worstt = vertex0t;
			}

			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
			double reflectiont = bestt + alph * (bestt - worstt);
			
			// Calculate the y value of the reflection.
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Figure out which values to use as replacement for the values of the worst vertex.
				double replacet = 0;
				if (reflectionY < extensionY) {
					replacet = reflectiont;
				} else {
					replacet = extensiont;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0t = replacet;
				} else if (vertex1Y == worstY) {
					vertex1t = replacet;
				}

			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractiont = bestt + beta * (bestt - worstt);
				
				// Calculate the y value of the contraction.
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0t = contractiont;
					} else if (vertex1Y == worstY) {
						vertex1t = contractiont;
					}
				
				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
					// vertex0
					vertex0t = bestt + delt  * (vertex0t - bestt);
					
					// vertex1
					vertex1t = bestt + delt  * (vertex1t - bestt);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY) {
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0t - vertex1t) > 0.0001) {
				keepGoing = true;
			} 

		} // while (keepGoing == true) {

		// Memorize the result as the max_simulation_age
		double max_tree_duration = vertex0t;

		// Part 2: Select four time points and estimate clade age probabilities for these time points.
		// This is equivalent to part 4 of function run_fitted_cladeage_star.
		
		// Set the number of replicates for both ndr and epsilon variation.
		int replicates = 30;
		
		// Specify the values of ndr to be used for simulations.
		double[] ndrs = new double[replicates];
		if (ndr_min == ndr_max) {
			for (int x = 0; x < replicates; x++) {
				ndrs[x] = ndr_min; 
			}
		} else {
			for (int x = 0; x < replicates; x++) {
				ndrs[x] = ndr_min + (ndr_max-ndr_min)/(double) (2*replicates) + (ndr_max-ndr_min)*(x/(double) replicates); 
			}			
		}
		
		// Specify the values of epsilon to be used for simulations.
		double[] epsilons = new double[replicates];
		if (epsilon_min == epsilon_max) {
			for (int x = 0; x < replicates; x++) {
				epsilons[x] = epsilon_min; 
			}
		} else {
			for (int x = 0; x < replicates; x++) {
				epsilons[x] = epsilon_min + (epsilon_max-epsilon_min)/(double) (2*replicates) + (epsilon_max-epsilon_min)*(x/(double) replicates); 
			}			
		}
		
		// Initiate arrays xs and ys.
		double[] xs = new double[4];
		xs[0] = 0.0;
		xs[1] = 0.2*max_tree_duration;
		xs[2] = 0.4*max_tree_duration;
		xs[3] = 0.8*max_tree_duration;
		double[] ys = new double[4];
		ys[0] = (psi_min+psi_max)/2.0;
		
		// Start tree simulation loop for different values of x, ndr, and epsilon.
		for (int z = 1; z < 4; z++) {
			double x = xs[z];
			double y = 0;
			for (int zz = 0; zz < replicates; zz++) {
				double ndr = ndrs[zz];
				for (int zzz = 0; zzz < replicates; zzz++) {
					double epsilon = epsilons[zzz];
					
					// Calculate values of lambda and mu from the net diversification rate and the turnover.
					double mu = (ndr*epsilon)/(1.0 - epsilon);
					double lambda = ndr + mu;
					
					// Set success to false
					boolean success = false;
					
					// Start the tree generation loop with these settings for x, lambda, and mu.
					while (success == false) {
						
						// Initiate arrays for branch origin and branch termination.
						ArrayList<Double> branch_origin = new ArrayList<Double>();
						ArrayList<Double> branch_termination = new ArrayList<Double>();
						ArrayList<Double> new_branch_origin = new ArrayList<Double>();
						ArrayList<Double> new_branch_termination = new ArrayList<Double>();
						new_branch_origin.add(0.0);

						// Start the tree generation loop.
						while (new_branch_origin.size() > 0) {
							
							// For each new origin, add a new termination.
							for (int o = 0; o < new_branch_origin.size(); o++) {
								new_branch_termination.add(new_branch_origin.get(o) + (Math.log(Math.random())/(-(lambda+mu))));
							}

							// Add new origin and termination to the old collection.
							for (int o = 0; o < new_branch_origin.size(); o++) {
								branch_origin.add(new_branch_origin.get(o));
							}
							for (int t = 0; t < new_branch_termination.size(); t++) {
								branch_termination.add(new_branch_termination.get(t));
							}

							// Empty the new origin array.
							new_branch_origin = new ArrayList<Double>();

							// For each new termination, add it to the new origin array if it is < x and rand < lambda/(lambda+mu) - this represents a speciation event.
							for (int t = 0; t < new_branch_termination.size(); t++) {										
								if (new_branch_termination.get(t) < x) {
									if (Math.random() < lambda/(lambda+mu)) {
										new_branch_origin.add(new_branch_termination.get(t));
										new_branch_origin.add(new_branch_termination.get(t));
									}
								}
							}
							
							// Empty the new termination array.
							new_branch_termination = new ArrayList<Double>();
							
						} // while (new_branch_origin.size() > 0)
						
						// Calculate the sum of species durations and the number of extant taxa.
						double sum_of_species_durations = 0.0;
						int number_of_extant_taxa = 0;
						for (int o = branch_origin.size()-1; o >= 0; o--) {
							if (branch_termination.get(o) >= x) {
								branch_termination.set(o,x);
								number_of_extant_taxa += 1;
							}
							sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
						} // for (int o = branch_origin.size()-1; o >= 0; o--)

						// Calculate the probability density of first observing a fossil at this x.
						if (number_of_extant_taxa > 0) {
							if (psi_min == psi_max) {
								y += (psi_min * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa)/(double) Math.pow(replicates,2);
							} else {
								if (sum_of_species_durations == 0) {
									y += (((psi_max+psi_min)/2.0) * number_of_extant_taxa)/(double) Math.pow(replicates,2);
								} else {
									y += (((-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_max*sum_of_species_durations) * number_of_extant_taxa * (1+psi_max*sum_of_species_durations))-(-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa * (1+psi_min*sum_of_species_durations)))/(psi_max-psi_min))/(double) Math.pow(replicates,2);
								}
							}
							success = true;
						} // if (number_of_extant_taxa > 0)
						
					} // while (success == false)
										
				} // for (int zzz = 0; zzz < replicates; zzz++)
			} // for (int zz = 0; zz < replicates; zz++)
			
			ys[z] = y;
			
		} // for (int z = 1; z < 4; z++)

		// Part 3: Analytical distribution fitting.
		// This is equivalent to part 5 in function run_fitted_cladeage_star(...).
		
		// Now that probability densities have been calculated, the CladeAge distribution is analytically calculated to fit these.
        // S is the only parameter of the CladeAge distribution which can not be solved analytically, it has to be solved numerically with the following
        // Nelder-Mead Downhill Simplex.
		
		// For easier readability, copy values for x and y.
        double x1 = xs[1];
        double y0 = ys[0];
        double y1 = ys[1];
        double y2 = ys[2];
        double y3 = ys[3];

		// Initiate the simplex, find 2 vertices.
        // vertex0
        double vertex0s = 1.0;
        
        // vertex1
        double vertex1s = 10.0;
        
        // Prepare for the Nelder-Mead loop.
        keepGoing = true;
		
        // Until converged do the loop.
		while (keepGoing == true) {

			// Calculate the y value of each vertex.
			// vertex0
			double lhs_denom1denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(2*x1+vertex0s)),2);
			double lhs_denom1nom = Math.log((y2/y0)*((2*x1+vertex0s)/vertex0s));
			double lhs_denom2denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(x1+vertex0s)),2);
			double lhs_denom2nom = Math.log((y1/y0)*((x1+vertex0s)/vertex0s));
			double lhs_nom1denom = 2*Math.log((x1+vertex0s)/vertex0s);
			double lhs_nom1nom = lhs_denom2nom;
			double lhs_nom2denom = 2*Math.log((2*x1+vertex0s)/vertex0s);
			double lhs_nom2_nom = lhs_denom1nom;
			double lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
			double rhs_denom1denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(4*x1+vertex0s)),2);
			double rhs_denom1nom = Math.log((y3/y0)*((4*x1+vertex0s)/vertex0s));
			double rhs_denom2denom = Math.pow((Math.log(vertex0s)),2) - Math.pow((Math.log(x1+vertex0s)),2);
			double rhs_denom2nom = Math.log((y1/y0)*((x1+vertex0s)/vertex0s));
			double rhs_nom1denom = 2*Math.log((x1+vertex0s)/vertex0s);
			double rhs_nom1nom = rhs_denom2nom;
			double rhs_nom2denom = 2*Math.log((4*x1+vertex0s)/vertex0s);
			double rhs_nom2_nom = rhs_denom1nom;
			double rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
			double vertex0Y = Math.pow((lhs-rhs),2);

			// vertex1
			lhs_denom1denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(2*x1+vertex1s)),2);
			lhs_denom1nom = Math.log((y2/y0)*((2*x1+vertex1s)/vertex1s));
			lhs_denom2denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(x1+vertex1s)),2);
			lhs_denom2nom = Math.log((y1/y0)*((x1+vertex1s)/vertex1s));
			lhs_nom1denom = 2*Math.log((x1+vertex1s)/vertex1s);
			lhs_nom1nom = lhs_denom2nom;
			lhs_nom2denom = 2*Math.log((2*x1+vertex1s)/vertex1s);
			lhs_nom2_nom = lhs_denom1nom;
			lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
			rhs_denom1denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(4*x1+vertex1s)),2);
			rhs_denom1nom = Math.log((y3/y0)*((4*x1+vertex1s)/vertex1s));
			rhs_denom2denom = Math.pow((Math.log(vertex1s)),2) - Math.pow((Math.log(x1+vertex1s)),2);
			rhs_denom2nom = Math.log((y1/y0)*((x1+vertex1s)/vertex1s));
			rhs_nom1denom = 2*Math.log((x1+vertex1s)/vertex1s);
			rhs_nom1nom = rhs_denom2nom;
			rhs_nom2denom = 2*Math.log((4*x1+vertex1s)/vertex1s);
			rhs_nom2_nom = rhs_denom1nom;
			rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
			double vertex1Y = Math.pow((lhs-rhs),2);

			// Find the best (=lowest) y value.
			double bestY = vertex0Y;
			double worstY = vertex1Y;
			double bests = vertex0s;
			double worsts = vertex1s;
			if (vertex1Y < vertex0Y) {
				bestY = vertex1Y;
				worstY = vertex0Y;
				bests = vertex1s;
				worsts = vertex0s;
			}
			
			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
			double reflections = bests + alph * (bests - worsts);
			
			// Calculate the y value of the reflection.
			lhs_denom1denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(2*x1+reflections)),2);
			lhs_denom1nom = Math.log((y2/y0)*((2*x1+reflections)/reflections));
			lhs_denom2denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(x1+reflections)),2);
			lhs_denom2nom = Math.log((y1/y0)*((x1+reflections)/reflections));
			lhs_nom1denom = 2*Math.log((x1+reflections)/reflections);
			lhs_nom1nom = lhs_denom2nom;
			lhs_nom2denom = 2*Math.log((2*x1+reflections)/reflections);
			lhs_nom2_nom = lhs_denom1nom;
			lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
			rhs_denom1denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(4*x1+reflections)),2);
			rhs_denom1nom = Math.log((y3/y0)*((4*x1+reflections)/reflections));
			rhs_denom2denom = Math.pow((Math.log(reflections)),2) - Math.pow((Math.log(x1+reflections)),2);
			rhs_denom2nom = Math.log((y1/y0)*((x1+reflections)/reflections));
			rhs_nom1denom = 2*Math.log((x1+reflections)/reflections);
			rhs_nom1nom = rhs_denom2nom;
			rhs_nom2denom = 2*Math.log((4*x1+reflections)/reflections);
			rhs_nom2_nom = rhs_denom1nom;
			rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
			double reflectionY = Math.pow((lhs-rhs),2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
				double extensions = bests + gamm * (bests - worsts);
				
				// Calculate the y value of the extension.
				lhs_denom1denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(2*x1+extensions)),2);
				lhs_denom1nom = Math.log((y2/y0)*((2*x1+extensions)/extensions));
				lhs_denom2denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(x1+extensions)),2);
				lhs_denom2nom = Math.log((y1/y0)*((x1+extensions)/extensions));
				lhs_nom1denom = 2*Math.log((x1+extensions)/extensions);
				lhs_nom1nom = lhs_denom2nom;
				lhs_nom2denom = 2*Math.log((2*x1+extensions)/extensions);
				lhs_nom2_nom = lhs_denom1nom;
				lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
				rhs_denom1denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(4*x1+extensions)),2);
				rhs_denom1nom = Math.log((y3/y0)*((4*x1+extensions)/extensions));
				rhs_denom2denom = Math.pow((Math.log(extensions)),2) - Math.pow((Math.log(x1+extensions)),2);
				rhs_denom2nom = Math.log((y1/y0)*((x1+extensions)/extensions));
				rhs_nom1denom = 2*Math.log((x1+extensions)/extensions);
				rhs_nom1nom = rhs_denom2nom;
				rhs_nom2denom = 2*Math.log((4*x1+extensions)/extensions);
				rhs_nom2_nom = rhs_denom1nom;
				rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
				double extensionY = Math.pow((lhs-rhs),2);

				// Figure out which values to use as replacement for the values of the worst vertex.
				double replaces = 0;
				if (reflectionY < extensionY) {
					replaces = reflections;
				} else {
					replaces = extensions;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0s = replaces;
				} else {
					vertex1s = replaces;
				}
				
			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractions = bests + beta * (bests - worsts);
				
				// Calculate the y value of the contraction.
				lhs_denom1denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(2*x1+contractions)),2);
				lhs_denom1nom = Math.log((y2/y0)*((2*x1+contractions)/contractions));
				lhs_denom2denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(x1+contractions)),2);
				lhs_denom2nom = Math.log((y1/y0)*((x1+contractions)/contractions));
				lhs_nom1denom = 2*Math.log((x1+contractions)/contractions);
				lhs_nom1nom = lhs_denom2nom;
				lhs_nom2denom = 2*Math.log((2*x1+contractions)/contractions);
				lhs_nom2_nom = lhs_denom1nom;
				lhs = ((lhs_denom1denom/lhs_denom1nom)-(lhs_denom2denom/lhs_denom2nom))/((lhs_nom1denom/lhs_nom1nom)-(lhs_nom2denom/lhs_nom2_nom));
				rhs_denom1denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(4*x1+contractions)),2);
				rhs_denom1nom = Math.log((y3/y0)*((4*x1+contractions)/contractions));
				rhs_denom2denom = Math.pow((Math.log(contractions)),2) - Math.pow((Math.log(x1+contractions)),2);
				rhs_denom2nom = Math.log((y1/y0)*((x1+contractions)/contractions));
				rhs_nom1denom = 2*Math.log((x1+contractions)/contractions);
				rhs_nom1nom = rhs_denom2nom;
				rhs_nom2denom = 2*Math.log((4*x1+contractions)/contractions);
				rhs_nom2_nom = rhs_denom1nom;
				rhs = ((rhs_denom1denom/rhs_denom1nom)-(rhs_denom2denom/rhs_denom2nom))/((rhs_nom1denom/rhs_nom1nom)-(rhs_nom2denom/rhs_nom2_nom));
				double contractionY = Math.pow((lhs-rhs),2);

				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0s = contractions;
					} else {
						vertex1s = contractions;
					}

				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
					// vertex0
					vertex0s = bests + delt * (vertex0s - bests);
					
					// vertex1
					vertex1s = bests + delt * (vertex1s - bests);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY)
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0s - vertex1s) > 0.0001) {
				keepGoing = true;
			}
			 
		} // while (keepGoing == true)

		double fittedCladeAgeS = vertex0s;
		double m_denom1denom = Math.pow((Math.log(fittedCladeAgeS)),2) - Math.pow((Math.log(2*x1+fittedCladeAgeS)),2);
		double m_denom1nom = Math.log((y2/y0)*((2*x1+fittedCladeAgeS)/fittedCladeAgeS));
		double m_denom2denom = Math.pow((Math.log(fittedCladeAgeS)),2) - Math.pow((Math.log(x1+fittedCladeAgeS)),2);
		double m_denom2nom = Math.log((y1/y0)*((x1+fittedCladeAgeS)/fittedCladeAgeS));
		double m_nom1denom = 2*Math.log((x1+fittedCladeAgeS)/fittedCladeAgeS);
		double m_nom1nom = m_denom2nom;
		double m_nom2denom = 2*Math.log((2*x1+fittedCladeAgeS)/fittedCladeAgeS);
		double m_nom2_nom = m_denom1nom;
		double fittedCladeAgeM = ((m_denom1denom/m_denom1nom)-(m_denom2denom/m_denom2nom))/((m_nom1denom/m_nom1nom)-(m_nom2denom/m_nom2_nom));
		double fittedCladeAgeW = (Math.pow((Math.log(fittedCladeAgeS)-fittedCladeAgeM),2) - Math.pow((Math.log(x1+fittedCladeAgeS)-fittedCladeAgeM),2)) / Math.log((y1/y0)*((x1+fittedCladeAgeS)/fittedCladeAgeS));
		double fittedCladeAgeC = y0 * fittedCladeAgeS * Math.exp(Math.pow((Math.log(fittedCladeAgeS)-fittedCladeAgeM),2)/fittedCladeAgeW);
		// Calculate the root mean square deviation (in order to be fair with other distributions, do that based on the 101 original ages. This can be removed if the fast distribution calculation is required.
		double sd = 0;
		for (int a = 0; a < ages.length; a++) {
			sd += Math.pow((fix_probabilities[a]-((fittedCladeAgeC/(ages[a]-first_occurrence_age_min+fittedCladeAgeS)) * Math.exp( -Math.pow((Math.log(ages[a]-first_occurrence_age_min+fittedCladeAgeS)-fittedCladeAgeM),2)/fittedCladeAgeW ))),2);
		}
		double msd = sd/(double) ages.length;
		double fittedCladeAgeRmsd = Math.sqrt(msd);
		
		// Variable fittedCladeAgeC, as determined above has been adjusted to provide good fit with the empirical probability values, 
		// however, it may not be the right scale factor for a distribution, which must sum to 1. Therefore, the integral under
		// the distribution is calculated between fossil_age_min and Inf, and a scaling factor @fittedCladeAgeCcorr is determined so that
		// the distribution sums to 1.
		double erf_argument = -(Math.log(fittedCladeAgeS)-fittedCladeAgeM) / Math.sqrt(fittedCladeAgeW);
		double erf_t = 1/(1.0+0.5*Math.abs(erf_argument));
		double erf_polynomial = -Math.pow(erf_argument,2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t,2) + 0.09678418*Math.pow(erf_t,3) - 0.18628806*Math.pow(erf_t,4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
		double erf_tau = erf_t * Math.exp(erf_polynomial);
		double erf = 0;
		if (erf_argument >= 0) {
			erf = 1 - erf_tau;
		} else {
			erf = erf_tau - 1;
		}
		double fittedCladeAgeCcorr = 1/(0.5 * Math.sqrt(Math.PI) * Math.sqrt(fittedCladeAgeW) * (1+erf));
		
		// Return the fitted CladeAge distribution.
		return new FittedCladeAgeDistribution(first_occurrence_age_min, first_occurrence_age_max, fittedCladeAgeCcorr, fittedCladeAgeS, fittedCladeAgeM, fittedCladeAgeW, fittedCladeAgeRmsd);
		
	} // public FittedCladeAgeDistribution run_fitted_cladeage_rapid(...)

	public ContinuousDistribution run_standard(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, String distribution_type, JProgressBar dpb) {
		
		// Reset arrays.
		ages = new double[number_of_ages + 1];
		fix_probabilities = new double[number_of_ages + 1];
		int_probabilities = new double[number_of_ages + 1];
		double[] raw_probabilities = new double[number_of_ages + 1];

		// Part 1: Calculate the maximum age, and 100 time points between fossil age and maximum age.
		
		// Roughly calculate a maximum tree duration so that the probability for the maximum tree duration is 0.001 * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		double ndr_mean = (ndr_max+ndr_min)/2.0;
		double psi_mean = (psi_max+psi_min)/2.0;
		
		// Initiate the simplex, find 2 vertices.
		// vertex0
		// Math.log(ndr_mean/psi_mean)/ndr_mean is the peak of p(t), as determined through its first derivative.
		double vertex0t = 0;
		if (Math.log(ndr_mean/psi_mean)/ndr_mean > 0) {
			vertex0t = Math.log(ndr_mean/psi_mean)/ndr_mean;
		}
		
		// vertex1
		double vertex1t = vertex0t + 10;
		
		// Prepare for the Nelder-Mead loop.
		boolean keepGoing = true;
		
		// Until converged do the loop.
		while (keepGoing == true) {
			
			// Calculate the y value of each vertex.
			// vertex0
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(0.001)), 2);

			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Find the best (=lowest) and the worst (=highest) y value.
			double bestY = vertex0Y;
			double bestt = vertex0t;
			double worstY = vertex1Y;
			double worstt = vertex1t;
			if (vertex1Y < bestY) {
				bestY = vertex1Y;
				bestt = vertex1t;
				worstY = vertex0Y;
				worstt = vertex0t;
			}

			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alph).
			double reflectiont = bestt + alph * (bestt - worstt);
			
			// Calculate the y value of the reflection.
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(0.001)), 2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamm).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Figure out which values to use as replacement for the values of the worst vertex.
				double replacet = 0;
				if (reflectionY < extensionY) {
					replacet = reflectiont;
				} else {
					replacet = extensiont;
				}
				
				// Replace the parameter values of the worst vertex with the replacement values.
				if (vertex0Y == worstY) {
					vertex0t = replacet;
				} else if (vertex1Y == worstY) {
					vertex1t = replacet;
				}

			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
			} else {
				
				// Calculate the contraction.
				double contractiont = bestt + beta * (bestt - worstt);
				
				// Calculate the y value of the contraction.
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(0.001)), 2);
				
				// Consider two subcases of case iii:
				// iiia) The contraction is better than the worst vertex.
				// iiib) The contraction is not better than the worst vertex.
				
				// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
				if (contractionY < worstY) {
					
					// Replace the parameter values of the worst vertex with the contraction values.
					if (vertex0Y == worstY) {
						vertex0t = contractiont;
					} else if (vertex1Y == worstY) {
						vertex1t = contractiont;
					}
				
				// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that delt  = 0.5, the default setting)
					// vertex0
					vertex0t = bestt + delt  * (vertex0t - bestt);
					
					// vertex1
					vertex1t = bestt + delt  * (vertex1t - bestt);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY) {
			
			// Stop the loop when all parameter values are identical in the first 3 decimals.
			keepGoing = false;
			if (Math.abs(vertex0t - vertex1t) > 0.0001) {
				keepGoing = true;
			} 

		} // while (keepGoing == true) {

		// Memorize the result as the max_simulation_age.
		double max_simulation_age = vertex0t + first_occurrence_age_min;
		
		// Determine the ages.
		ages[0] = first_occurrence_age_min;
		for (int i = 0; i < number_of_ages; i++) {
			ages[i+1] = first_occurrence_age_min + ((i+1)/(double) number_of_ages)*(max_simulation_age-first_occurrence_age_min);
		}

		// Part 2: Estimate clade age probability densities.
		
		// Prepare the arrays raw_probabilities and successful simulations (those in which at least one taxon survives), and fill them with 0.
		int[] successful_simulations = new int[number_of_ages];
		
		while (successful_simulations[-1] < 10000) {

			if (cancel) {
				return null;
			}

			// Update the progress bar.
			if (dpb != null) {
				dpb.setValue(successful_simulations[-1]);
			}

			// Draw values for parameters ndr (net diversification rate, lambda - mu) and epsilon (turnover rate, mu/lambda) from uniform distributions, and calculate lambda and mu from them.
			double ndr = ndr_min + Math.random() * (ndr_max-ndr_min);
			double epsilon = epsilon_min + Math.random() * (epsilon_max-epsilon_min);
			double mu = (ndr*epsilon)/(1.0 - epsilon);
			double lambda = ndr + mu;
			
			// Determine the maximum tree duration.
			double max_tree_duration = ages[-1] - first_occurrence_age_min;
			
			// Initiate arrays for branch origin and branch termination.
			ArrayList<Double> branch_origin = new ArrayList<Double>();
			ArrayList<Double> branch_termination = new ArrayList<Double>();
			ArrayList<Double> new_branch_origin = new ArrayList<Double>();
			ArrayList<Double> new_branch_termination = new ArrayList<Double>();
			new_branch_origin.add(0.0);

			// Start the tree generation loop.
			while (new_branch_origin.size() > 0) {

				// For each new origin, add a new termination.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					new_branch_termination.add(new_branch_origin.get(o) + (Math.log(Math.random())/(-(lambda+mu))));
				}

				// Add new origin and termination to the old collection.
				for (int o = 0; o < new_branch_origin.size(); o++) {
					branch_origin.add(new_branch_origin.get(o));
				}
				for (int t = 0; t < new_branch_termination.size(); t++) {
					branch_termination.add(new_branch_termination.get(t));
				}
				
				// Empty the new origin array.
				new_branch_origin = new ArrayList<Double>();
				
				// For each new termination, add it to the new origin array if it is < max_tree_duration and rand < lambda/(lambda+mu) - this represents a speciation event.
				for (int t = 0; t < new_branch_termination.size(); t++) {										
					if (new_branch_termination.get(t) < max_tree_duration) {
						if (Math.random() < lambda/(lambda+mu)) {
							new_branch_origin.add(new_branch_termination.get(t));
							new_branch_origin.add(new_branch_termination.get(t));
						}
					}
				}
				
				// Empty the new termination array.
				new_branch_termination = new ArrayList<Double>();
				
			} // while (new_branch_origin.size() > 0)
			
			// Analyse the obtained set of branches.
			for (int i = ages.length-1; i >= 0; i--) {
				
				// For each age, determine the tree length.
				double tree_duration = ages[i] - first_occurrence_age_min;
				
				// Trim the tree so that no branches are older than ages[i] - first_occurrence_age. While trimming, make sure the tree contains at least one extant species.
				double sum_of_species_durations = 0.0;
				int number_of_extant_taxa = 0;
				for (int o = branch_origin.size()-1; o >= 0; o--) {
					boolean remove_this_branch = false;
					if (branch_origin.get(o) > tree_duration) {
						remove_this_branch = true;
					} else {
						if (branch_termination.get(o) >= tree_duration) {
							branch_termination.set(o,tree_duration);
							number_of_extant_taxa += 1;
						}
						sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
					}
					if (remove_this_branch == true) {
						branch_origin.remove(o);
						branch_termination.remove(o);
					}
				} // for (int o = branch_origin.size()-1; o >= 0; o--)
					
				// If at least one species is extant at this age, get the sum of lineage durations.
				if (number_of_extant_taxa > 0) {
					if (psi_min == psi_max) {
						raw_probabilities[i] += psi_min * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa;
					} else {
						if (sum_of_species_durations == 0) {
							raw_probabilities[i] += ((psi_max+psi_min)/2.0) * number_of_extant_taxa;
						} else {
							raw_probabilities[i] += ((-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_max*sum_of_species_durations) * number_of_extant_taxa * (1+psi_max*sum_of_species_durations))-(-(1/Math.pow(sum_of_species_durations,2)) * Math.exp(-psi_min*sum_of_species_durations) * number_of_extant_taxa * (1+psi_min*sum_of_species_durations)))/(psi_max-psi_min);
						}
					}
					successful_simulations[i] += 1;
				} // if (number_of_extant_taxa > 0)
				
			} // for (int i = 0; i < ages.length; i++)
			
		} // while (successful_simulations[-1] < 10000)
		
		if (cancel) {
			return null;
		}

		// Calculate the fix_probabilities (that's when first_occurrence_age_min == first_occurrence_age_max, and int_probabilities which are
		// different from fix_probabilities when this is not the case.		
		for (int i = 0; i < raw_probabilities.length; i++) {
			fix_probabilities[i] = raw_probabilities[i]/(double) successful_simulations[i];
		}
		
		// If first_occurrence_age_min == first_occurrence_age_max, fix_probabilities and int_probabilities are equal.
		if (first_occurrence_age_min == first_occurrence_age_max) {
			
			for (int i = 0; i < fix_probabilities.length; i++) {
				int_probabilities[i] = fix_probabilities[i];
			}
			
		// If not, int_probabilities have to be calculated from fix_probabilities.
		} else {
			
			int_probabilities[0] = 0.0;
			for (int a = 1; a < ages.length; a++) {
				double probability_mass = 0;
				// If the difference between first_occurrence_age_max and first_occurrence_age_min is smaller than half the difference
				// between two adjacent ages, the probability mass for a given age is simply the difference between first_occurrence_age_max
				// and first_occurrence_age_min times the probability of that age.
				if (first_occurrence_age_max-first_occurrence_age_min <= (ages[1]-ages[0])/2.0) {
					probability_mass += fix_probabilities[a] * (first_occurrence_age_max-first_occurrence_age_min);
				// If the difference between first_occurrence_age_max and first_occurrence_age_min is more than half the difference between
				// two adjacent ages (this means the integral extends over more than one segment), three steps are necessary.
				} else {
					// i.) Half of the last segment is added in all cases.
					probability_mass += fix_probabilities[a] * (ages[1]-ages[0])/2.0;
					
					// ii.) The first, or part of the first segment is added if the integral extends beyond or into it.
					if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[0]) {
						// Add half the first segment.
						probability_mass += fix_probabilities[0] * (ages[1]-ages[0])/2.0;
					} else if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[0] + (ages[1]-ages[0])/2.0) {
						// Add a proportion of the first segment.
						probability_mass += fix_probabilities[0] * ((ages[0]+(ages[1]-ages[0])/2.0) - (ages[a]-(first_occurrence_age_max-first_occurrence_age_min)));
					}
					
					// iii.) For each segment between the first and the last, the whole segment, or part of that segment are added if the integral
					// extends beyond or into it.
					for (int aa = 1; aa < a; a++) {
						if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) <= ages[aa] - (ages[1]-ages[0])/2.0) {
							// Add the full segment.
							probability_mass += fix_probabilities[aa] * (ages[1]-ages[0]);
						} else if (ages[a] - (first_occurrence_age_max-first_occurrence_age_min) < ages[aa] + (ages[1]-ages[0])/2.0) {
							// Add a proportion of this segment.
							probability_mass += fix_probabilities[aa] * ((ages[aa]+(ages[1]-ages[0])/2.0) - (ages[a]-(first_occurrence_age_max-first_occurrence_age_min)));
						}
					}
				}
				// The probability mass is divided by the difference between first_occurrence_age_max first_occurrence_age_min.
				int_probabilities[a] = probability_mass/(first_occurrence_age_max-first_occurrence_age_min);
			} // for (int a = 1; a < ages.length; a++) {
			
		} // if (first_occurrence_age_min == first_occurrence_age_max)

		// Part 3: Distribution fitting.
		
		if (distribution_type == "Lognormal") {

			// Set the number of Nelder Mead Downhill Simplex repetitions.
			int nmRepetitions = 10;
			
			// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
			double[] cs = new double[nmRepetitions];
			double[] sigmas = new double[nmRepetitions];
			double[] mus = new double[nmRepetitions];
			double[] ys = new double[nmRepetitions];
			
			for (int x = 0; x < nmRepetitions; x++) {
				
				if (cancel != false) {
					return null;
				}

				// Initiate the simplex, find 4 vertices.
				// vertex0
				double vertex0c____ = 0.5 + Math.random();
				double vertex0mu___ = 1 + Math.random()*2;
				double vertex0sigma = 1 + Math.random()*3;
				double vertex0Y = 0;
				
				// vertex1
				double vertex1c____ = 0.5 + Math.random();
				double vertex1mu___ = 1 + Math.random()*2;
				double vertex1sigma = 1 + Math.random()*3;
				double vertex1Y = 0;

				// vertex2
				double vertex2c____ = 0.5 + Math.random();
				double vertex2mu___ = 1 + Math.random()*2;
				double vertex2sigma = 1 + Math.random()*3;
				double vertex2Y = 0;

				// vertex3
				double vertex3c____ = 0.5 + Math.random();
				double vertex3mu___ = 1 + Math.random()*2;
				double vertex3sigma = 1 + Math.random()*3;
				double vertex3Y = 0;

				// Prepare for the Nelder-Mead loop.
				keepGoing = true;
				int stepCounter = 0;

				// Until converged do the loop.
				while (keepGoing == true) {

					if (cancel != false) {
						return null;
					}

					// Increment the stepCounter.
					stepCounter += 1;

					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-first_occurrence_age_min > 0) {
							temp = (1.0/((ages[i]-first_occurrence_age_min)*Math.sqrt(2*Math.PI*Math.pow(vertex0sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min)-vertex0mu___),2)/(2*Math.pow(vertex0sigma,2)) ));
						}
						vertex0Y += Math.pow((int_probabilities[i]-vertex0c____*temp),2);
					}

					// Calculate the y value of each vertex.
					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-first_occurrence_age_min > 0) {
							temp = (1.0/((ages[i]-first_occurrence_age_min)*Math.sqrt(2*Math.PI*Math.pow(vertex1sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min)-vertex1mu___),2)/(2*Math.pow(vertex1sigma,2)) ));
						}
						vertex1Y += Math.pow((int_probabilities[i]-vertex1c____*temp),2);
					}

					// Calculate the y value of each vertex.
					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-first_occurrence_age_min > 0) {
							temp = (1.0/((ages[i]-first_occurrence_age_min)*Math.sqrt(2*Math.PI*Math.pow(vertex2sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min)-vertex2mu___),2)/(2*Math.pow(vertex2sigma,2)) ));
						}
						vertex2Y += Math.pow((int_probabilities[i]-vertex2c____*temp),2);
					}

					// Calculate the y value of each vertex.
					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-first_occurrence_age_min > 0) {
							temp = (1.0/((ages[i]-first_occurrence_age_min)*Math.sqrt(2*Math.PI*Math.pow(vertex3sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min)-vertex3mu___),2)/(2*Math.pow(vertex3sigma,2)) ));
						}
						vertex3Y += Math.pow((int_probabilities[i]-vertex3c____*temp),2);
					}

					// Find the best (=lowest) y value.
					double bestY = vertex0Y;
					if (vertex1Y < bestY) {
						bestY = vertex1Y;
					}
					if (vertex2Y < bestY) {
						bestY = vertex2Y;
					}
					if (vertex3Y < bestY) {
						bestY = vertex3Y;
					}

					// Find the worst (=highest) y value.
					double worstY = vertex0Y;
					if (vertex1Y > worstY) {
						worstY = vertex1Y;
					}
					if (vertex2Y > worstY) {
						worstY = vertex2Y;
					}
					if (vertex3Y > worstY) {
						worstY = vertex3Y;
					}

					// Find the second-worst (=second-highest) y value.
					double secondWorstY = bestY;
					if (vertex0Y > secondWorstY && vertex0Y != worstY) {
						secondWorstY = vertex0Y;
					}
					if (vertex1Y > secondWorstY && vertex1Y != worstY) {
						secondWorstY = vertex1Y;
					}
					if (vertex2Y > secondWorstY && vertex2Y != worstY) {
						secondWorstY = vertex2Y;
					}
					if (vertex3Y > secondWorstY && vertex3Y != worstY) {
						secondWorstY = vertex3Y;
					}

					// Find the parameter values of the best vertex.
					double bestc____ = 0;
					double bestmu___ = 0;
					double bestsigma = 0;
					if (vertex0Y == bestY) {
						bestc____ = vertex0c____;
						bestmu___ = vertex0mu___;
						bestsigma = vertex0sigma;
					} else if (vertex1Y == bestY) {
						bestc____ = vertex1c____;
						bestmu___ = vertex1mu___;
						bestsigma = vertex1sigma;
					} else if (vertex2Y == bestY) {
						bestc____ = vertex2c____;
						bestmu___ = vertex2mu___;
						bestsigma = vertex2sigma;
					} else if (vertex3Y == bestY) {
						bestc____ = vertex3c____;
						bestmu___ = vertex3mu___;
						bestsigma = vertex3sigma;
					}

					// Find the parameter values of the worst vertex.
					double worstc____ = 0;
					double worstmu___ = 0;
					double worstsigma = 0;
					if (vertex0Y == worstY) {
						worstc____ = vertex0c____;
						worstmu___ = vertex0mu___;
						worstsigma = vertex0sigma;
					} else if (vertex1Y == worstY) {
						worstc____ = vertex1c____;
						worstmu___ = vertex1mu___;
						worstsigma = vertex1sigma;
					} else if (vertex2Y == worstY) {
						worstc____ = vertex2c____;
						worstmu___ = vertex2mu___;
						worstsigma = vertex2sigma;
					} else if (vertex3Y == worstY) {
						worstc____ = vertex3c____;
						worstmu___ = vertex3mu___;
						worstsigma = vertex3sigma;
					}

					// Calculate the sum of the parameters over all vertices.
					double sumc____ = vertex0c____ + vertex1c____ + vertex2c____ + vertex3c____;
					double summu___ = vertex0mu___ + vertex1mu___ + vertex2mu___ + vertex3mu___;
					double sumsigma = vertex0sigma + vertex1sigma + vertex2sigma + vertex3sigma;

					// Calculate the parameter values of the centroid.
					double centroidc____ = (sumc____ - worstc____)/3.0;
					double centroidmu___ = (summu___ - worstmu___)/3.0;
					double centroidsigma = (sumsigma - worstsigma)/3.0;

					// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
					double reflectionc____ = centroidc____ + alph * (centroidc____ - worstc____);
					double reflectionmu___ = centroidmu___ + alph * (centroidmu___ - worstmu___);
					if (reflectionmu___ <= 0) {
						reflectionmu___ = 1 + Math.random()*2;
					}
					double reflectionsigma = centroidsigma + alph * (centroidsigma - worstsigma);
					if (reflectionsigma <= 0) {
						reflectionsigma = 1 + Math.random()*3;
					}

					// Calculate the y value of the reflection.
					double reflectionY = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-first_occurrence_age_min > 0) {
							temp = (1.0/((ages[i]-first_occurrence_age_min)*Math.sqrt(2*Math.PI*Math.pow(reflectionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min)-reflectionmu___),2)/(2*Math.pow(reflectionsigma,2)) ));
						}
						reflectionY += Math.pow((int_probabilities[i]-reflectionc____*temp),2);
					}

					// Consider the three cases:
					// i.)   reflection is better than all vertices.
					// ii.)  reflection is better than the second-worst vertex.
					// iii.) reflection is worse than or equally good as the second-worst vertex.

					// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
					if (reflectionY < bestY) {

						// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
						double extensionc____ = centroidc____ + gamm * (centroidc____ - worstc____);
						double extensionmu___ = centroidmu___ + gamm * (centroidmu___ - worstmu___);
						if (extensionmu___ <= 0) {
							extensionmu___ = 1 + Math.random()*2;
						}
						double extensionsigma = centroidsigma + gamm * (centroidsigma - worstsigma);
						if (extensionsigma <= 0) {
							extensionsigma = 1 + Math.random()*3;
						}

						// Calculate the y value of the extension.
						double extensionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = 0;
							if (ages[i]-first_occurrence_age_min > 0) {
								temp = (1.0/((ages[i]-first_occurrence_age_min)*Math.sqrt(2*Math.PI*Math.pow(extensionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min)-extensionmu___),2)/(2*Math.pow(extensionsigma,2)) ));
							}
							extensionY += Math.pow((int_probabilities[i]-extensionc____*temp),2);
						}

						// Figure out which values to use as replacement for the values of the worst vertex.
						double replacec____ = 0;
						double replacemu___ = 0;
						double replacesigma = 0;
						if (reflectionY < extensionY) {
							replacec____ = reflectionc____;
							replacemu___ = reflectionmu___;
							replacesigma = reflectionsigma;
						} else {
							replacec____ = extensionc____;
							replacemu___ = extensionmu___;
							replacesigma = extensionsigma;
						}

						// Replace the parameter values of the worst vertex with the replacement values.
						if (vertex0Y == worstY) {
							vertex0c____ = replacec____;
							vertex0mu___ = replacemu___;
							vertex0sigma = replacesigma;
						} else if (vertex1Y == worstY) {
							vertex1c____ = replacec____;
							vertex1mu___ = replacemu___;
							vertex1sigma = replacesigma;
						} else if (vertex2Y == worstY) {
							vertex2c____ = replacec____;
							vertex2mu___ = replacemu___;
							vertex2sigma = replacesigma;
						} else if (vertex3Y == worstY) {
							vertex3c____ = replacec____;
							vertex3mu___ = replacemu___;
							vertex3sigma = replacesigma;
						}

						// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
					} else if (reflectionY < secondWorstY) {

						if (vertex0Y == worstY) {
							vertex0c____ = reflectionc____;
							vertex0mu___ = reflectionmu___;
							vertex0sigma = reflectionsigma;
						} else if  (vertex1Y == worstY) {
							vertex1c____ = reflectionc____;
							vertex1mu___ = reflectionmu___;
							vertex1sigma = reflectionsigma;
						} else if (vertex2Y == worstY) {
							vertex2c____ = reflectionc____;
							vertex2mu___ = reflectionmu___;
							vertex2sigma = reflectionsigma;
						} else if (vertex3Y == worstY) {
							vertex3c____ = reflectionc____;
							vertex3mu___ = reflectionmu___;
							vertex3sigma = reflectionsigma;
						}

						// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
					} else {

						// Calculate the contraction.
						double contractionc____ = centroidc____ + beta * (centroidc____ - worstc____);
						double contractionmu___ = centroidmu___ + beta * (centroidmu___ - worstmu___);
						double contractionsigma = centroidsigma + beta * (centroidsigma - worstsigma);

						// Calculate the y value of the contraction.
						double contractionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = 0;
							if (ages[i]-first_occurrence_age_min > 0) {
								temp = (1.0/((ages[i]-first_occurrence_age_min)*Math.sqrt(2*Math.PI*Math.pow(contractionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min)-contractionmu___),2)/(2*Math.pow(contractionsigma,2)) ));
							}
							contractionY += Math.pow((int_probabilities[i]-contractionc____*temp),2);
						}

						// Consider two subcases of case iii:
						// iiia) The contraction is better than the worst vertex.
						// iiib) The contraction is not better than the worst vertex.

						// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
						if (contractionY < worstY) {

							// Replace the parameter values of the worst vertex with the contraction values.
							if (vertex0Y == worstY) {
								vertex0c____ = contractionc____;
								vertex0mu___ = contractionmu___;
								vertex0sigma = contractionsigma;
							} else if (vertex1Y == worstY) {
								vertex1c____ = contractionc____;
								vertex1mu___ = contractionmu___;
								vertex1sigma = contractionsigma;
							} else if (vertex2Y == worstY) {
								vertex2c____ = contractionc____;
								vertex2mu___ = contractionmu___;
								vertex2sigma = contractionsigma;
							} else if (vertex3Y == worstY) {
								vertex3c____ = contractionc____;
								vertex3mu___ = contractionmu___;
								vertex3sigma = contractionsigma;
							}


						// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
						} else {

							// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that delt  = 0.5, the default setting)
							// vertex0
							vertex0c____ = bestc____ + delt  * (vertex0c____ - bestc____);
							vertex0mu___ = bestmu___ + delt  * (vertex0mu___ - bestmu___);
							vertex0sigma = bestsigma + delt  * (vertex0sigma - bestsigma);

							// vertex1
							vertex1c____ = bestc____ + delt  * (vertex1c____ - bestc____);
							vertex1mu___ = bestmu___ + delt  * (vertex1mu___ - bestmu___);
							vertex1sigma = bestsigma + delt  * (vertex1sigma - bestsigma);

							// vertex0
							vertex2c____ = bestc____ + delt  * (vertex2c____ - bestc____);
							vertex2mu___ = bestmu___ + delt  * (vertex2mu___ - bestmu___);
							vertex2sigma = bestsigma + delt  * (vertex2sigma - bestsigma);

							// vertex0
							vertex3c____ = bestc____ + delt  * (vertex3c____ - bestc____);
							vertex3mu___ = bestmu___ + delt  * (vertex3mu___ - bestmu___);
							vertex3sigma = bestsigma + delt  * (vertex3sigma - bestsigma);

						} // if (contractionY < worstY)

					} // if (reflectionY < bestY)

					// Stop the loop when all parameter values are identical in the first 10 decimals.
					keepGoing = false;
					if (Math.abs(vertex0c____- vertex1c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___- vertex1mu___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma- vertex1sigma) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____- vertex2c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___- vertex2mu___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma- vertex2sigma) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____- vertex3c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___- vertex3mu___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma- vertex3sigma) > 0.0000001) {
						keepGoing = true;
					}
					if (stepCounter > 1000) {
						keepGoing = false;
					}

				} // while (keepGoing == true) 

				// Report all parameters and y.
				cs[x] = vertex0c____;
				mus[x] = vertex0mu___;
				sigmas[x] = vertex0sigma;
				ys[x] = vertex0Y;

			} // for (int x = 0; x < nmRepetitions; x++)

			// Find the best result among the nmReplicates replicates.
			int index = 0;
			for (int x = 1; x < nmRepetitions; x++) {
				if (ys[x] < ys[index]) {
					index = x;
				}
			}

			double logConstant = cs[index];
			double logMean = mus[index];
			double logStdev = sigmas[index];
			double logRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));

			// Fill variables distribution_type, distribution_parameters, and distribution_rmsd
			distribution_parameters[0] = logMean;
			distribution_parameters[1] = logStdev;
			distribution_rmsd = logRmsd;

			// XXX this is just here for tests.
			System.out.println("Distribution type: " + distribution_type);
			System.out.println("Mean (log): " + distribution_parameters[0]);
			System.out.println("Stdev (log): " + distribution_parameters[1]);
			System.out.println("RMSD: " + distribution_rmsd);

			// XXX todo: Find a standard way to use the normaliser!
			normaliser = logConstant;			
			
			// XXX todo: The below distribution must be implemented, remove the below.
			// return new XXX(distribution_parameters[0], distribution_parameters[1], distribution_parameters[2], distribution_parameters[3], distribution_parameters[4], distribution_parameters[5]);
			return null;
	
		} else if (distribution_type == "Gamma") {

			// Set the number of Nelder Mead Downhill Simplex repetitions.
			int nmRepetitions = 10;

			// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
			double[] cs = new double[nmRepetitions];
			double[] ks = new double[nmRepetitions];
			double[] thetas = new double[nmRepetitions];
			double[] gammaKs = new double[nmRepetitions];
			double[] ys = new double[nmRepetitions];
			
			for (int x = 0; x < nmRepetitions; x++) {
				
				if (cancel != false) {
					return null;
				}
					
				// Initiate the simplex, find 4 vertices.
				// vertex0
				double vertex0c____ = 0.5 + Math.random();
				double vertex0k____ = 1 + Math.random()*3;
				double vertex0theta = (Math.random()*90 + 10)/vertex0k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				double k = vertex0k____;
				double lanczos_z = k-1;
				double lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				double lanczos_t = lanczos_z + lanczos_g + 0.5;
				double gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex0gammaK = gammaK;
				double vertex0Y = 0;

				// vertex1
				double vertex1c____ = 0.5 + Math.random();
				double vertex1k____ = 1 + Math.random()*3;
				double vertex1theta = (Math.random()*90 + 10)/vertex1k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				k = vertex1k____;
				lanczos_z = k-1;
				lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				lanczos_t = lanczos_z + lanczos_g + 0.5;
				gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex1gammaK = gammaK;
				double vertex1Y = 0;

				// vertex2
				double vertex2c____ = 0.5 + Math.random();
				double vertex2k____ = 1 + Math.random()*3;
				double vertex2theta = (Math.random()*90 + 10)/vertex2k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				k = vertex2k____;
				lanczos_z = k-1;
				lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				lanczos_t = lanczos_z + lanczos_g + 0.5;
				gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex2gammaK = gammaK;
				double vertex2Y = 0;

				// vertex3
				double vertex3c____ = 0.5 + Math.random();
				double vertex3k____ = 1 + Math.random()*3;
				double vertex3theta = (Math.random()*90 + 10)/vertex3k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				k = vertex3k____;
				lanczos_z = k-1;
				lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				lanczos_t = lanczos_z + lanczos_g + 0.5;
				gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex3gammaK = gammaK;
				double vertex3Y = 0;

				// Prepare for the Nelder-Mead loop.
				keepGoing = true;
				int stepCounter = 0;

				// Until converged do the loop.
				while (keepGoing == true) {

					// Increment the stepCounter.
					stepCounter += 1;

					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex0theta,vertex0k____)))*(1/vertex0gammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(vertex0k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/vertex0theta));
						vertex0Y += Math.pow((int_probabilities[i]-vertex0c____*temp),2);
					}

					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex1theta,vertex1k____)))*(1/vertex1gammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(vertex1k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/vertex1theta));
						vertex1Y += Math.pow((int_probabilities[i]-vertex1c____*temp),2);
					}

					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex2theta,vertex2k____)))*(1/vertex2gammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(vertex2k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/vertex2theta));
						vertex2Y += Math.pow((int_probabilities[i]-vertex2c____*temp),2);
					}

					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex3theta,vertex3k____)))*(1/vertex3gammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(vertex3k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/vertex3theta));
						vertex3Y += Math.pow((int_probabilities[i]-vertex3c____*temp),2);
					}

					// Find the best (=lowest) y value.
					double bestY = vertex0Y;
					if (vertex1Y < bestY) {
						bestY = vertex1Y;
					}
					if (vertex2Y < bestY) {
						bestY = vertex2Y;
					}
					if (vertex3Y < bestY) {
						bestY = vertex3Y;
					}

					// Find the worst (=highest) y value.
					double worstY = vertex0Y;
					if (vertex1Y > worstY) {
						worstY = vertex1Y;
					}
					if (vertex2Y > worstY) {
						worstY = vertex2Y;
					}
					if (vertex3Y > worstY) {
						worstY = vertex3Y;
					}

					// Find the second-worst (=second-highest) y value.
					double secondWorstY = bestY;
					if (vertex0Y > secondWorstY && vertex0Y != worstY) {
						secondWorstY = vertex0Y;
					}
					if (vertex1Y > secondWorstY && vertex1Y != worstY) {
						secondWorstY = vertex1Y;
					}
					if (vertex2Y > secondWorstY && vertex2Y != worstY) {
						secondWorstY = vertex2Y;
					}
					if (vertex3Y > secondWorstY && vertex3Y != worstY) {
						secondWorstY = vertex3Y;
					}

					// Find the parameter values of the best vertex.

					double bestc____ = 0;
					double bestk____ = 0;
					double besttheta = 0;
					if (vertex0Y == bestY) {
						bestc____ = vertex0c____;
						bestk____ = vertex0k____;
						besttheta = vertex0theta;
					} else if (vertex1Y == bestY) {
						bestc____ = vertex1c____;
						bestk____ = vertex1k____;
						besttheta = vertex1theta;
					} else if (vertex2Y == bestY) {
						bestc____ = vertex2c____;
						bestk____ = vertex2k____;
						besttheta = vertex2theta;
					} else if (vertex3Y == bestY) {
						bestc____ = vertex3c____;
						bestk____ = vertex3k____;
						besttheta = vertex3theta;
					}

					// Find the parameter values of the worst vertex.
					double worstc____ = 0;
					double worstk____ = 0;
					double worsttheta = 0;
					if (vertex0Y == worstY) {
						worstc____ = vertex0c____;
						worstk____ = vertex0k____;
						worsttheta = vertex0theta;
					} else if (vertex1Y == worstY) {
						worstc____ = vertex1c____;
						worstk____ = vertex1k____;
						worsttheta = vertex1theta;
					} else if (vertex2Y == worstY) {
						worstc____ = vertex2c____;
						worstk____ = vertex2k____;
						worsttheta = vertex2theta;
					} else if (vertex3Y == worstY) {
						worstc____ = vertex3c____;
						worstk____ = vertex3k____;
						worsttheta = vertex3theta;
					}

					// Calculate the sum of the parameters over all vertices.
					double sumc____ = vertex0c____ + vertex1c____ + vertex2c____ + vertex3c____;
					double sumk____ = vertex0k____ + vertex1k____ + vertex2k____ + vertex3k____;
					double sumtheta = vertex0theta + vertex1theta + vertex2theta + vertex3theta;

					// Calculate the parameter values of the centroid.
					double centroidc____ = (sumc____ - worstc____)/3.0;
					double centroidk____ = (sumk____ - worstk____)/3.0;
					double centroidtheta = (sumtheta - worsttheta)/3.0;

					// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
					double reflectionc____ = centroidc____ + alph * (centroidc____ - worstc____);
					double reflectionk____ = centroidk____ + alph * (centroidk____ - worstk____);
					if (reflectionk____ < 0) {
						reflectionk____ = 1 + Math.random()*3;
					}
					double reflectiontheta = centroidtheta + alph * (centroidtheta - worsttheta);
					if (reflectiontheta < 0) {
						reflectiontheta = (Math.random()*90 + 10)/reflectionk____;
					}
					// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
					// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
					k = reflectionk____;
					lanczos_z = k-1;
					lanczos_x = lanczos_p[0];
					for (int i = 1; i <= lanczos_g+1; i++){
						lanczos_x += lanczos_p[i]/(lanczos_z+i);
					}
					lanczos_t = lanczos_z + lanczos_g + 0.5;
					gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
					double reflectiongammaK = gammaK;

					// Calculate the y value of the reflection.
					double reflectionY = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(reflectiontheta,reflectionk____)))*(1/reflectiongammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(reflectionk____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/reflectiontheta));
						reflectionY += Math.pow((int_probabilities[i]-reflectionc____*temp),2);
					}

					// Consider the three cases:
					// i.)   reflection is better than all vertices.
					// ii.)  reflection is better than the second-worst vertex.
					// iii.) reflection is worse than or equally good as the second-worst vertex.

					// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
					if (reflectionY < bestY) {

						// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
						double extensionc____ = centroidc____ + gamm * (centroidc____ - worstc____);
						double extensionk____ = centroidk____ + gamm * (centroidk____ - worstk____);
						if (extensionk____ < 0) {
							extensionk____ = 1 + Math.random()*3;
						}
						double extensiontheta = centroidtheta + gamm * (centroidtheta - worsttheta);
						if (extensiontheta < 0) {
							extensiontheta = (Math.random()*90 + 10)/extensionk____;
						}
						// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
						// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
						k = extensionk____;
						lanczos_z = k-1;
						lanczos_x = lanczos_p[0];
						for (int i = 1; i <= lanczos_g+1; i++){
							lanczos_x += lanczos_p[i]/(lanczos_z+i);
						}
						lanczos_t = lanczos_z + lanczos_g + 0.5;
						gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
						double extensiongammaK = gammaK;

						// Calculate the y value of the reflection.
						double extensionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp =  (1.0/(Math.pow(extensiontheta,extensionk____)))*(1/extensiongammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(extensionk____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/extensiontheta));
							extensionY += Math.pow((int_probabilities[i]-extensionc____*temp),2);
						}

						// Figure out which values to use as replacement for the values of the worst vertex.
						double replacec____ = 0;
						double replacek____ = 0;
						double replacetheta = 0;
						double replacegammaK = 0;
						if (reflectionY < extensionY) {
							replacec____ = reflectionc____;
							replacek____ = reflectionk____;
							replacetheta = reflectiontheta;
							replacegammaK = reflectiongammaK;
						} else {
							replacec____ = extensionc____;
							replacek____ = extensionk____;
							replacetheta = extensiontheta;
							replacegammaK = extensiongammaK;
						}

						// Replace the parameter values of the worst vertex with the replacement values.
						if (vertex0Y == worstY) {
							vertex0c____ = replacec____;
							vertex0k____ = replacek____;
							vertex0theta = replacetheta;
							vertex0gammaK = replacegammaK;
						} else if (vertex1Y == worstY) {
							vertex1c____ = replacec____;
							vertex1k____ = replacek____;
							vertex1theta = replacetheta;
							vertex1gammaK = replacegammaK;
						} else if (vertex2Y == worstY) {
							vertex2c____ = replacec____;
							vertex2k____ = replacek____;
							vertex2theta = replacetheta;
							vertex2gammaK = replacegammaK;
						} else if (vertex3Y == worstY) {
							vertex3c____ = replacec____;
							vertex3k____ = replacek____;
							vertex3theta = replacetheta;
							vertex3gammaK = replacegammaK;
						}

					// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
					} else if (reflectionY < secondWorstY) {

						if (vertex0Y == worstY) {
							vertex0c____ = reflectionc____;
							vertex0k____ = reflectionk____;
							vertex0theta = reflectiontheta;
							vertex0gammaK = reflectiongammaK;
						} else if  (vertex1Y == worstY) {
							vertex1c____ = reflectionc____;
							vertex1k____ = reflectionk____;
							vertex1theta = reflectiontheta;
							vertex1gammaK = reflectiongammaK;
						} else if (vertex2Y == worstY) {
							vertex2c____ = reflectionc____;
							vertex2k____ = reflectionk____;
							vertex2theta = reflectiontheta;
							vertex2gammaK = reflectiongammaK;
						} else if (vertex3Y == worstY) {
							vertex3c____ = reflectionc____;
							vertex3k____ = reflectionk____;
							vertex3theta = reflectiontheta;
							vertex3gammaK = reflectiongammaK;
						}

					// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
					} else {

						// Calculate the contraction.
						double contractionc____ = centroidc____ + beta * (centroidc____ - worstc____);
						double contractionk____ = centroidk____ + beta * (centroidk____ - worstk____);
						double contractiontheta = centroidtheta + beta * (centroidtheta - worsttheta);
						// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
						// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
						k = contractionk____;
						lanczos_z = k-1;
						lanczos_x = lanczos_p[0];
						for (int i = 1; i <= lanczos_g+1; i++){
							lanczos_x += lanczos_p[i]/(lanczos_z+i);
						}
						lanczos_t = lanczos_z + lanczos_g + 0.5;
						gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
						double contractiongammaK = gammaK;

						// Calculate the y value of the contraction.
						double contractionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp =  (1.0/(Math.pow(contractiontheta,contractionk____)))*(1/contractiongammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(contractionk____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/contractiontheta));
							contractionY += Math.pow((int_probabilities[i]-contractionc____*temp),2);
						}

						// Consider two subcases of case iii:
						// iiia) The contraction is better than the worst vertex.
						// iiib) The contraction is not better than the worst vertex.

						// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
						if (contractionY < worstY) {

							// Replace the parameter values of the worst vertex with the contraction values.
							if (vertex0Y == worstY) {
								vertex0c____ = contractionc____;
								vertex0k____ = contractionk____;
								vertex0theta = contractiontheta;
								vertex0gammaK = contractiongammaK;
							} else if (vertex1Y == worstY) {
								vertex1c____ = contractionc____;
								vertex1k____ = contractionk____;
								vertex1theta = contractiontheta;
								vertex1gammaK = contractiongammaK;
							} else if (vertex2Y == worstY) {
								vertex2c____ = contractionc____;
								vertex2k____ = contractionk____;
								vertex2theta = contractiontheta;
								vertex2gammaK = contractiongammaK;
							} else if (vertex3Y == worstY) {
								vertex3c____ = contractionc____;
								vertex3k____ = contractionk____;
								vertex3theta = contractiontheta;
								vertex3gammaK = contractiongammaK;
							}

						// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
						} else {

							// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
							// vertex0
							vertex0c____ = bestc____ + delt  * (vertex0c____ - bestc____);
							vertex0k____ = bestk____ + delt  * (vertex0k____ - bestk____);
							vertex0theta = besttheta + delt  * (vertex0theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex0k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex0gammaK = gammaK;

							// vertex1
							vertex1c____ = bestc____ + delt  * (vertex1c____ - bestc____);
							vertex1k____ = bestk____ + delt  * (vertex1k____ - bestk____);
							vertex1theta = besttheta + delt  * (vertex1theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex1k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex1gammaK = gammaK;

							// vertex2
							vertex2c____ = bestc____ + delt  * (vertex2c____ - bestc____);
							vertex2k____ = bestk____ + delt  * (vertex2k____ - bestk____);
							vertex2theta = besttheta + delt  * (vertex2theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex2k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex2gammaK = gammaK;

							// vertex3
							vertex3c____ = bestc____ + delt  * (vertex3c____ - bestc____);
							vertex3k____ = bestk____ + delt  * (vertex3k____ - bestk____);
							vertex3theta = besttheta + delt  * (vertex3theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex3k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex3gammaK = gammaK;

						} // if (contractionY < worstY)

					} // if (reflectionY < bestY)

					// Stop the loop when all parameter values are identical in the first 10 decimals.
					keepGoing = false;
					if (Math.abs(vertex0c____- vertex1c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____- vertex1k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta- vertex1theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____- vertex2c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____- vertex2k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta- vertex2theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____- vertex3c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____- vertex3k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta- vertex3theta) > 0.0000001) {
						keepGoing = true;
					}
					if (stepCounter > 10000) {
						keepGoing = false;
					}

				} // while (keepGoing == true) 

				// Report all parameters and y.
				cs[x] = vertex0c____;
				ks[x] = vertex0k____;
				thetas[x] = vertex0theta;
				gammaKs[x] = vertex0gammaK;
				ys[x] = vertex0Y;
	                
			} // for (int x = 0; x < nmRepetitions; x++)
						
			// Find the best result among the nmReplicates replicates.
			int index = 0;
			for (int x = 1; x < nmRepetitions; x++) {
				if (ys[x] < ys[index]) {
					index = x;
				}
			}

			double gamConstant = cs[index];
			double gamShape = ks[index];
			double gamScale = thetas[index];
			double gamRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));

			// Fill variables distribution_type, distribution_parameters, and distribution_rmsd
			distribution_parameters[0] = gamShape;
			distribution_parameters[1] = gamScale;
			distribution_rmsd = gamRmsd;

			// XXX this is just here for tests.
			System.out.println("Distribution type: " + distribution_type);
			System.out.println("Shape: " + distribution_parameters[0]);
			System.out.println("Scale: " + distribution_parameters[1]);
			System.out.println("RMSD: " + distribution_rmsd);

			// XXX todo: Find a standard way to use the normaliser!
			normaliser = gamConstant;			
			
			return new GammaDistributionImpl(distribution_parameters[0], distribution_parameters[1]);						
	
		} else if (distribution_type == "Exponential") {

			// Set the number of Nelder Mead Downhill Simplex repetitions.
			int nmRepetitions = 10;

			// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
			double[] cs = new double[nmRepetitions];
			double[] means = new double[nmRepetitions];
			double[] ys = new double[nmRepetitions];
			
			for (int x = 0; x < nmRepetitions; x++) {

				if (cancel != false) {
					return null;
				}
					
					// Initiate the simplex, find 3 vertices.
					// vertex0
					double vertex0c___ = 0.5 + Math.random();
					double vertex0mean = 10 + Math.random()*50;
					double vertex0Y = 0;
				
					// vertex1
					double vertex1c___ = 0.5 + Math.random();
					double vertex1mean = 10 + Math.random()*50;
					double vertex1Y = 0;
				
					// vertex2
					double vertex2c___ = 0.5 + Math.random();
					double vertex2mean = 10 + Math.random()*50;
					double vertex2Y = 0;
					
	                // Prepare for the Nelder-Mead loop.
	                keepGoing = true;
	                int stepCounter = 0;

	                // Until converged do the loop.
					while (keepGoing == true) {

						// Increment the stepCounter.
						stepCounter += 1;
						
						// Calculate the y value of each vertex.
						// vertex0
						vertex0Y = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = (1/vertex0mean) * Math.exp(-(1/vertex0mean)*(ages[i]-first_occurrence_age_min));
							vertex0Y += Math.pow(int_probabilities[i]-vertex0c___*temp, 2);
						}

						// vertex1
						vertex1Y = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = (1/vertex1mean) * Math.exp(-(1/vertex1mean)*(ages[i]-first_occurrence_age_min));
							vertex1Y += Math.pow(int_probabilities[i]-vertex1c___*temp, 2);
						}

						// vertex2
						vertex2Y = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = (1/vertex2mean) * Math.exp(-(1/vertex2mean)*(ages[i]-first_occurrence_age_min));
							vertex2Y += Math.pow(int_probabilities[i]-vertex2c___*temp, 2);
						}
						
						// Find the best (=lowest) y value.
						double bestY = vertex0Y;
						if (vertex1Y < bestY) {
							bestY = vertex1Y;
						}
						if (vertex2Y < bestY) {
							bestY = vertex2Y;
						}
						
						// Find the worst (=highest) y value.
						double worstY = vertex0Y;
						if (vertex1Y > worstY) {
							worstY = vertex1Y;
						}
						if (vertex2Y > worstY) {
							worstY = vertex2Y;
						}
						
						// Find the second-worst (=second-highest) y value.
						double secondWorstY = bestY;
						if (vertex0Y > secondWorstY && vertex0Y != worstY) {
							secondWorstY = vertex0Y;
						}
						if (vertex1Y > secondWorstY && vertex1Y != worstY) {
							secondWorstY = vertex1Y;
						}
						if (vertex2Y > secondWorstY && vertex2Y != worstY) {
							secondWorstY = vertex2Y;
						}
						
						// Find the parameter values of the best vertex.
						double bestc___ = 0;
						double bestmean = 0;
						if (vertex0Y == bestY) {
							bestc___ = vertex0c___;
							bestmean = vertex0mean;						
						} else if (vertex1Y == bestY) {
							bestc___ = vertex1c___;
							bestmean = vertex1mean;												
						} else if (vertex2Y == bestY) {
							bestc___ = vertex2c___;
							bestmean = vertex2mean;												
						}
						
						// Find the parameter values of the worst vertex.
						double worstc___ = 0;
						double worstmean = 0;
						if (vertex0Y == worstY) {
							worstc___ = vertex0c___;
							worstmean = vertex0mean;
						} else if (vertex1Y == worstY) {
							worstc___ = vertex1c___;
							worstmean = vertex1mean;
						} else if (vertex2Y == worstY) {
							worstc___ = vertex2c___;
							worstmean = vertex2mean;
						}
						
						// Calculate the sum of the parameters over all vertices.
						double sumc___ = vertex0c___ + vertex1c___ + vertex2c___;
						double summean = vertex0mean + vertex1mean + vertex2mean;

						// Calculate the parameter values of the centroid.
						double centroidc___ = (sumc___ - worstc___)/2.0;
						double centroidmean = (summean - worstmean)/2.0;
						
						// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
						double reflectionc___ = centroidc___ + alph * (centroidc___ - worstc___);
						double reflectionmean = centroidmean + alph * (centroidmean - worstmean);

						// Calculate the y value of the reflection.
						double reflectionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = (1/reflectionmean) * Math.exp(-(1/reflectionmean)*(ages[i]-first_occurrence_age_min));
							reflectionY += Math.pow(int_probabilities[i]-reflectionc___*temp, 2);
						}

						// Consider the three cases:
						// i.)   reflection is better than all vertices.
						// ii.)  reflection is better than the second-worst vertex.
						// iii.) reflection is worse than or equally good as the second-worst vertex.
						
						// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
						if (reflectionY < bestY) {
							
							// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
							double extensionc___ = centroidc___ + gamm * (centroidc___ - worstc___);
							double extensionmean = centroidmean + gamm * (centroidmean - worstmean);
							
							// Calculate the y value of the extension.
							double extensionY = 0;
							for (int i = 0; i < ages.length; i++) {
								double temp = (1/extensionmean) * Math.exp(-(1/extensionmean)*(ages[i]-first_occurrence_age_min));
								reflectionY += Math.pow(int_probabilities[i]-extensionc___*temp, 2);
							}
							
							// Figure out which values to use as replacement for the values of the worst vertex.
							double replacec___ = 0;
							double replacemean = 0;
							if (reflectionY < extensionY) {
								replacec___ = reflectionc___;
								replacemean = reflectionmean;
							} else {
								replacec___ = extensionc___;
								replacemean = extensionmean;
							}
							
							// Replace the parameter values of the worst vertex with the replacement values.
							if (vertex0Y == worstY) {
								vertex0c___ = replacec___;
								vertex0mean = replacemean;
							} else if (vertex1Y == worstY) {
								vertex1c___ = replacec___;
								vertex1mean = replacemean;
							} else if (vertex2Y == worstY) {
								vertex2c___ = replacec___;
								vertex2mean = replacemean;
							}
							
						// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
						} else if (reflectionY < secondWorstY) {
							
							if (vertex0Y == worstY) {
								vertex0c___ = reflectionc___;
								vertex0mean = reflectionmean;
							} else if  (vertex1Y == worstY) {
								vertex1c___ = reflectionc___;
								vertex1mean = reflectionmean;
							} else if (vertex2Y == worstY) {
								vertex2c___ = reflectionc___;
								vertex2mean = reflectionmean;
							}
							
						// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
						} else {
							
							// Calculate the contraction.
							double contractionc___ = centroidc___ + beta * (centroidc___ - worstc___);
							double contractionmean = centroidmean + beta * (centroidmean - worstmean);
							
							// Calculate the y value of the contraction.
							double contractionY = 0;
							for (int i = 0; i < ages.length; i++) {
								double temp = (1/contractionmean) * Math.exp(-(1/contractionmean)*(ages[i]-first_occurrence_age_min));
								contractionY += Math.pow(int_probabilities[i]-contractionc___*temp, 2);
							}
							
							// Consider two subcases of case iii:
							// iiia) The contraction is better than the worst vertex.
							// iiib) The contraction is not better than the worst vertex.
							
							// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
							if (contractionY < worstY) {
								
								// Replace the parameter values of the worst vertex with the contraction values.
								if (vertex0Y == worstY) {
									vertex0c___ = contractionc___;
									vertex0mean = contractionmean;
								} else if (vertex1Y == worstY) {
									vertex1c___ = contractionc___;
									vertex1mean = contractionmean;
								} else if (vertex2Y == worstY) {
									vertex2c___ = contractionc___;
									vertex2mean = contractionmean;
								}
								
							// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
							} else {
								
								// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
								// vertex0
								vertex0c___ = bestc___ + delt * (vertex0c___ - bestc___);
								vertex0mean = bestmean + delt * (vertex0mean - bestmean);
								
								// vertex1
								vertex1c___ = bestc___ + delt * (vertex1c___ - bestc___);
								vertex1mean = bestmean + delt * (vertex1mean - bestmean);
								
								// vertex2
								vertex2c___ = bestc___ + delt * (vertex2c___ - bestc___);
								vertex2mean = bestmean + delt * (vertex2mean - bestmean);
								
							} // if (contractionY < worstY)
							
						} // if (reflectionY < bestY)

						// Stop the loop when all parameter values are identical in the first 10 decimals.
						keepGoing = false;
						if (Math.abs(vertex0c___- vertex1c___) > 0.0000001) {
							keepGoing = true;
						} else if (Math.abs(vertex0mean- vertex1mean) > 0.0000001) {
							keepGoing = true;
						} else if (Math.abs(vertex0c___- vertex2c___) > 0.0000001) {
							keepGoing = true;
						} else if (Math.abs(vertex0mean- vertex2mean) > 0.0000001) {
							keepGoing = true;
						}
						if (stepCounter > 1000) {
							keepGoing = false;
						}
						
					} // while (keepGoing == true) 
					
					// Report all parameters and y.
					cs[x] = vertex0c___;
					means[x] = vertex0mean;
					ys[x] = vertex0Y;
								
				//} if (cancel2 == false) {

			} // for (int x = 0; x < nmRepetitions; x++)

			// Find the best result among the nmReplicates replicates.
			int index = 0;
			for (int x = 1; x < nmRepetitions; x++) {
				if (ys[x] < ys[index]) {
					index = x;
				}
			}

			double expConstant = cs[index];
			double expMean = means[index];
			double expRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));

			
			
			// Fill variables distribution_parameters and distribution_rmsd.
			distribution_parameters[0] = expMean;
			distribution_rmsd = expRmsd;
						
			// XXX this is just here for tests.
			System.out.println("Distribution type: " + distribution_type);
			System.out.println("Mean: " + distribution_parameters[0]);
			System.out.println("RMSD: " + distribution_parameters);

			// XXX todo: Find a standard way to use the normaliser!			
			normaliser = expConstant;

			// XXX todo: Check whether this is the right way to return the distribution.
			return new ExponentialDistributionImpl(distribution_parameters[0]);		
	
		} else if (distribution_type == "TruncatedLognormal") {

			// Set the number of Nelder Mead Downhill Simplex repetitions.
			int nmRepetitions = 3;

			// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
			double[] cs = new double[nmRepetitions];
			double[] shifts = new double[nmRepetitions];
			double[] sigmas = new double[nmRepetitions];
			double[] mus = new double[nmRepetitions];
			double[] ys = new double[nmRepetitions];

			for (int x = 0; x < nmRepetitions; x++) {

				if (cancel != false) {
					return null;
				}

				// Initiate the simplex, find 4 vertices.
				// vertex0
				double vertex0c____ = 0.5 + Math.random();
				double vertex0shift = 10 + Math.random()*10;
				double vertex0mu___ = 1 + Math.random()*2;
				double vertex0sigma = 1 + Math.random()*3;
				double vertex0Y = 0;

				// vertex1
				double vertex1c____ = 0.5 + Math.random();
				double vertex1shift = 10 + Math.random()*10;
				double vertex1mu___ = 1 + Math.random()*2;
				double vertex1sigma = 1 + Math.random()*3;
				double vertex1Y = 0;

				// vertex2
				double vertex2c____ = 0.5 + Math.random();
				double vertex2shift = 10 + Math.random()*10;
				double vertex2mu___ = 1 + Math.random()*2;
				double vertex2sigma = 1 + Math.random()*3;
				double vertex2Y = 0;

				// vertex3
				double vertex3c____ = 0.5 + Math.random();
				double vertex3shift = 10 + Math.random()*10;
				double vertex3mu___ = 1 + Math.random()*2;
				double vertex3sigma = 1 + Math.random()*3;
				double vertex3Y = 0;

				// vertex4
				double vertex4c____ = 0.5 + Math.random();
				double vertex4shift = 10 + Math.random()*10;
				double vertex4mu___ = 1 + Math.random()*2;
				double vertex4sigma = 1 + Math.random()*3;
				double vertex4Y = 0;

				// Prepare for the Nelder-Mead loop.
				keepGoing = true;
				int stepCounter = 0;

				// Until converged do the loop.
				while (keepGoing == true) {

					if (cancel != false) {
						return null;
					}

					// Increment the stepCounter.
					stepCounter += 1;

					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/((ages[i]-first_occurrence_age_min+vertex0shift)*Math.sqrt(2*Math.PI*Math.pow(vertex0sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex0shift)-vertex0mu___),2)/(2*Math.pow(vertex0sigma,2)) ));
						vertex0Y += Math.pow((fix_probabilities[i]-vertex0c____*temp),2);
					}

					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/((ages[i]-first_occurrence_age_min+vertex1shift)*Math.sqrt(2*Math.PI*Math.pow(vertex1sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex1shift)-vertex1mu___),2)/(2*Math.pow(vertex1sigma,2)) ));
						vertex1Y += Math.pow((fix_probabilities[i]-vertex1c____*temp),2);
					}

					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/((ages[i]-first_occurrence_age_min+vertex2shift)*Math.sqrt(2*Math.PI*Math.pow(vertex2sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex2shift)-vertex2mu___),2)/(2*Math.pow(vertex2sigma,2)) ));
						vertex2Y += Math.pow((fix_probabilities[i]-vertex2c____*temp),2);
					}

					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/((ages[i]-first_occurrence_age_min+vertex3shift)*Math.sqrt(2*Math.PI*Math.pow(vertex3sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex3shift)-vertex3mu___),2)/(2*Math.pow(vertex3sigma,2)) ));
						vertex3Y += Math.pow((fix_probabilities[i]-vertex3c____*temp),2);
					}

					// vertex4
					vertex4Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/((ages[i]-first_occurrence_age_min+vertex4shift)*Math.sqrt(2*Math.PI*Math.pow(vertex4sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+vertex4shift)-vertex4mu___),2)/(2*Math.pow(vertex4sigma,2)) ));
						vertex4Y += Math.pow((fix_probabilities[i]-vertex4c____*temp),2);
					}

					// Find the best (=lowest) y value.
					double bestY = vertex0Y;
					if (vertex1Y < bestY) {
						bestY = vertex1Y;
					}
					if (vertex2Y < bestY) {
						bestY = vertex2Y;
					}
					if (vertex3Y < bestY) {
						bestY = vertex3Y;
					}
					if (vertex4Y < bestY) {
						bestY = vertex4Y;
					}

					// Find the worst (=highest) y value.
					double worstY = vertex0Y;
					if (vertex1Y > worstY) {
						worstY = vertex1Y;
					}
					if (vertex2Y > worstY) {
						worstY = vertex2Y;
					}
					if (vertex3Y > worstY) {
						worstY = vertex3Y;
					}
					if (vertex4Y > worstY) {
						worstY = vertex4Y;
					}

					// Find the second-worst (=second-highest) y value.
					double secondWorstY = bestY;
					if (vertex0Y > secondWorstY && vertex0Y != worstY) {
						secondWorstY = vertex0Y;
					}
					if (vertex1Y > secondWorstY && vertex1Y != worstY) {
						secondWorstY = vertex1Y;
					}
					if (vertex2Y > secondWorstY && vertex2Y != worstY) {
						secondWorstY = vertex2Y;
					}
					if (vertex3Y > secondWorstY && vertex3Y != worstY) {
						secondWorstY = vertex3Y;
					}
					if (vertex4Y > secondWorstY && vertex4Y != worstY) {
						secondWorstY = vertex4Y;
					}

					// Find the parameter values of the best vertex.
					double bestc____ = 0;
					double bestshift = 0;
					double bestmu___ = 0;
					double bestsigma = 0;
					if (vertex0Y == bestY) {
						bestc____ = vertex0c____;
						bestshift = vertex0shift;
						bestmu___ = vertex0mu___;
						bestsigma = vertex0sigma;
					} else if (vertex1Y == bestY) {
						bestc____ = vertex1c____;
						bestshift = vertex1shift;
						bestmu___ = vertex1mu___;
						bestsigma = vertex1sigma;
					} else if (vertex2Y == bestY) {
						bestc____ = vertex2c____;
						bestshift = vertex2shift;
						bestmu___ = vertex2mu___;
						bestsigma = vertex2sigma;
					} else if (vertex3Y == bestY) {
						bestc____ = vertex3c____;
						bestshift = vertex3shift;
						bestmu___ = vertex3mu___;
						bestsigma = vertex3sigma;
					} else if (vertex4Y == bestY) {
						bestc____ = vertex4c____;
						bestshift = vertex4shift;
						bestmu___ = vertex4mu___;
						bestsigma = vertex4sigma;
					}

					// Find the parameter values of the worst vertex.
					double worstc____ = 0;
					double worstshift = 0;
					double worstmu___ = 0;
					double worstsigma = 0;
					if (vertex0Y == worstY) {
						worstc____ = vertex0c____;
						worstshift = vertex0shift;
						worstmu___ = vertex0mu___;
						worstsigma = vertex0sigma;
					} else if (vertex1Y == worstY) {
						worstc____ = vertex1c____;
						worstshift = vertex1shift;
						worstmu___ = vertex1mu___;
						worstsigma = vertex1sigma;
					} else if (vertex2Y == worstY) {
						worstc____ = vertex2c____;
						worstshift = vertex2shift;
						worstmu___ = vertex2mu___;
						worstsigma = vertex2sigma;
					} else if (vertex3Y == worstY) {
						worstc____ = vertex3c____;
						worstshift = vertex3shift;
						worstmu___ = vertex3mu___;
						worstsigma = vertex3sigma;
					} else if (vertex4Y == worstY) {
						worstc____ = vertex4c____;
						worstshift = vertex4shift;
						worstmu___ = vertex4mu___;
						worstsigma = vertex4sigma;
					}

					// Calculate the sum of the parameters over all vertices.
					double sumc____ = vertex0c____ + vertex1c____ + vertex2c____ + vertex3c____ + vertex4c____;
					double sumshift = vertex0shift + vertex1shift + vertex2shift + vertex3shift + vertex4shift;
					double summu___ = vertex0mu___ + vertex1mu___ + vertex2mu___ + vertex3mu___ + vertex4mu___;
					double sumsigma = vertex0sigma + vertex1sigma + vertex2sigma + vertex3sigma + vertex4sigma;

					// Calculate the parameter values of the centroid.
					double centroidc____ = (sumc____ - worstc____)/4.0;
					double centroidshift = (sumshift - worstshift)/4.0;
					double centroidmu___ = (summu___ - worstmu___)/4.0;
					double centroidsigma = (sumsigma - worstsigma)/4.0;

					// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
					double reflectionc____ = centroidc____ + alph * (centroidc____ - worstc____);
					double reflectionshift = centroidshift + alph * (centroidshift - worstshift);
					if (reflectionshift <= 0) {
						reflectionshift = 10 + Math.random()*10;
					}
					double reflectionmu___ = centroidmu___ + alph * (centroidmu___ - worstmu___);
					if (reflectionmu___ <= 0) {
						reflectionmu___ = 1 + Math.random()*2;
					}
					double reflectionsigma = centroidsigma + alph * (centroidsigma - worstsigma);
					if (reflectionsigma <= 0) {
						reflectionsigma = 1 + Math.random()*3;
					}
					
					// Calculate the y value of the reflection.
					double reflectionY = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/((ages[i]-first_occurrence_age_min+reflectionshift)*Math.sqrt(2*Math.PI*Math.pow(reflectionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+reflectionshift)-reflectionmu___),2)/(2*Math.pow(reflectionsigma,2)) ));
						reflectionY += Math.pow((fix_probabilities[i]-reflectionc____*temp),2);
					}

					// Consider the three cases:
					// i.)   reflection is better than all vertices.
					// ii.)  reflection is better than the second-worst vertex.
					// iii.) reflection is worse than or equally good as the second-worst vertex.

					// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
					if (reflectionY < bestY) {

						// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
						double extensionc____ = centroidc____ + gamm * (centroidc____ - worstc____);
						double extensionshift = centroidshift + gamm * (centroidshift - worstshift);
						if (extensionshift <= 0) {
							extensionshift = 10 + Math.random()*10;
						}						
						double extensionmu___ = centroidmu___ + gamm * (centroidmu___ - worstmu___);
						if (extensionmu___ <= 0) {
							extensionmu___ = 1 + Math.random()*2;
						}
						double extensionsigma = centroidsigma + gamm * (centroidsigma - worstsigma);
						if (extensionsigma <= 0) {
							extensionsigma = 1 + Math.random()*3;
						}

						// Calculate the y value of the extension.
						double extensionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = (1.0/((ages[i]-first_occurrence_age_min+extensionshift)*Math.sqrt(2*Math.PI*Math.pow(extensionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+extensionshift)-extensionmu___),2)/(2*Math.pow(extensionsigma,2)) ));
							extensionY += Math.pow((int_probabilities[i]-extensionc____*temp),2);
						}
						
						// Figure out which values to use as replacement for the values of the worst vertex.
						double replacec____ = 0;
						double replaceshift = 0;
						double replacemu___ = 0;
						double replacesigma = 0;
						if (reflectionY < extensionY) {
							replacec____ = reflectionc____;
							replaceshift = reflectionshift;
							replacemu___ = reflectionmu___;
							replacesigma = reflectionsigma;
						} else {
							replacec____ = extensionc____;
							replaceshift = extensionshift;
							replacemu___ = extensionmu___;
							replacesigma = extensionsigma;
						}

						// Replace the parameter values of the worst vertex with the replacement values.
						if (vertex0Y == worstY) {
							vertex0c____ = replacec____;
							vertex0shift = replaceshift;
							vertex0mu___ = replacemu___;
							vertex0sigma = replacesigma;
						} else if (vertex1Y == worstY) {
							vertex1c____ = replacec____;
							vertex1shift = replaceshift;
							vertex1mu___ = replacemu___;
							vertex1sigma = replacesigma;
						} else if (vertex2Y == worstY) {
							vertex2c____ = replacec____;
							vertex2shift = replaceshift;
							vertex2mu___ = replacemu___;
							vertex2sigma = replacesigma;
						} else if (vertex3Y == worstY) {
							vertex3c____ = replacec____;
							vertex3shift = replaceshift;
							vertex3mu___ = replacemu___;
							vertex3sigma = replacesigma;
						} else if (vertex4Y == worstY) {
							vertex4c____ = replacec____;
							vertex4shift = replaceshift;
							vertex4mu___ = replacemu___;
							vertex4sigma = replacesigma;
						}
						
					// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
					} else if (reflectionY < secondWorstY) {

						if (vertex0Y == worstY) {
							vertex0c____ = reflectionc____;
							vertex0shift = reflectionshift;
							vertex0mu___ = reflectionmu___;
							vertex0sigma = reflectionsigma;
						} else if  (vertex1Y == worstY) {
							vertex1c____ = reflectionc____;
							vertex1shift = reflectionshift;
							vertex1mu___ = reflectionmu___;
							vertex1sigma = reflectionsigma;
						} else if (vertex2Y == worstY) {
							vertex2c____ = reflectionc____;
							vertex2shift = reflectionshift;
							vertex2mu___ = reflectionmu___;
							vertex2sigma = reflectionsigma;
						} else if (vertex3Y == worstY) {
							vertex3c____ = reflectionc____;
							vertex3shift = reflectionshift;
							vertex3mu___ = reflectionmu___;
							vertex3sigma = reflectionsigma;
						} else if (vertex4Y == worstY) {
							vertex4c____ = reflectionc____;
							vertex4shift = reflectionshift;
							vertex4mu___ = reflectionmu___;
							vertex4sigma = reflectionsigma;
						}

					// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor beta (=ro in the Wikipedia example)).
					} else {
						
						// Calculate the contraction.
						double contractionc____ = centroidc____ + beta * (centroidc____ - worstc____);
						double contractionshift = centroidshift + beta * (centroidshift - worstshift);
						double contractionmu___ = centroidmu___ + beta * (centroidmu___ - worstmu___);
						double contractionsigma = centroidsigma + beta * (centroidsigma - worstsigma);
						
						// Calculate the y value of the contraction.
						double contractionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp = (1.0/((ages[i]-first_occurrence_age_min+contractionshift)*Math.sqrt(2*Math.PI*Math.pow(contractionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-first_occurrence_age_min+contractionshift)-contractionmu___),2)/(2*Math.pow(contractionsigma,2)) ));
							contractionY += Math.pow((int_probabilities[i]-contractionc____*temp),2);
						}
						
						// Consider two subcases of case iii:
						// iiia) The contraction is better than the worst vertex.
						// iiib) The contraction is not better than the worst vertex.

						// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
						if (contractionY < worstY) {

							// Replace the parameter values of the worst vertex with the contraction values.
							if (vertex0Y == worstY) {
								vertex0c____ = contractionc____;
								vertex0shift = contractionshift;
								vertex0mu___ = contractionmu___;
								vertex0sigma = contractionsigma;
							} else if (vertex1Y == worstY) {
								vertex1c____ = contractionc____;
								vertex1shift = contractionshift;
								vertex1mu___ = contractionmu___;
								vertex1sigma = contractionsigma;
							} else if (vertex2Y == worstY) {
								vertex2c____ = contractionc____;
								vertex2shift = contractionshift;
								vertex2mu___ = contractionmu___;
								vertex2sigma = contractionsigma;
							} else if (vertex3Y == worstY) {
								vertex3c____ = contractionc____;
								vertex3shift = contractionshift;
								vertex3mu___ = contractionmu___;
								vertex3sigma = contractionsigma;
							} else if (vertex4Y == worstY) {
								vertex4c____ = contractionc____;
								vertex4shift = contractionshift;
								vertex4mu___ = contractionmu___;
								vertex4sigma = contractionsigma;
							}

						// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
						} else {

							// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that delt  = 0.5, the default setting)
							// vertex0
							vertex0c____ = bestc____ + delt * (vertex0c____ - bestc____);
							vertex0shift = bestshift + delt * (vertex0shift - bestshift); 
							vertex0mu___ = bestmu___ + delt * (vertex0mu___ - bestmu___);
							vertex0sigma = bestsigma + delt * (vertex0sigma - bestsigma);

							// vertex1
							vertex1c____ = bestc____ + delt * (vertex1c____ - bestc____);
							vertex1shift = bestshift + delt * (vertex1shift - bestshift);
							vertex1mu___ = bestmu___ + delt * (vertex1mu___ - bestmu___);
							vertex1sigma = bestsigma + delt * (vertex1sigma - bestsigma);

							// vertex2
							vertex2c____ = bestc____ + delt * (vertex2c____ - bestc____);
							vertex2shift = bestshift + delt * (vertex2shift - bestshift);
							vertex2mu___ = bestmu___ + delt * (vertex2mu___ - bestmu___);
							vertex2sigma = bestsigma + delt * (vertex2sigma - bestsigma);

							// vertex3
							vertex3c____ = bestc____ + delt * (vertex3c____ - bestc____);
							vertex3shift = bestshift + delt * (vertex3shift - bestshift);
							vertex3mu___ = bestmu___ + delt * (vertex3mu___ - bestmu___);
							vertex3sigma = bestsigma + delt * (vertex3sigma - bestsigma);
		
							// vertex4
							vertex4c____ = bestc____ + delt * (vertex4c____ - bestc____);
							vertex4shift = bestshift + delt * (vertex4shift - bestshift);
							vertex4mu___ = bestmu___ + delt * (vertex4mu___ - bestmu___);
							vertex4sigma = bestsigma + delt * (vertex4sigma - bestsigma);

						} // if (contractionY < worstY) {
												
					} // if (reflectionY < bestY) {
					
					// Stop the loop when all parameter values are identical in the first 10 decimals.
					keepGoing = false;
					if (Math.abs(vertex0c____ - vertex1c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex1shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___ - vertex1mu___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma - vertex1sigma) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____ - vertex2c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex2shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___ - vertex2mu___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma - vertex2sigma) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____ - vertex3c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex3shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___ - vertex3mu___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma - vertex3sigma) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____ - vertex4c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex4shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___ - vertex4mu___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma - vertex4sigma) > 0.0000001) {
						keepGoing = true;
					}
					if (stepCounter > 100000) {
						keepGoing = false;
					}
				
				}  // while (keepGoing == true)
				
				// Report all parameters and y.
				cs[x] = vertex0c____;
				shifts[x] = vertex0shift;
				mus[x] = vertex0mu___;
				sigmas[x] = vertex0sigma;
				ys[x] = vertex0Y;

			} // for (int x = 0; x < nmRepetitions; x++)
			
			// Find the best result among the nmReplicates replicates.
			int index = 0;
			for (int x = 1; x < nmRepetitions; x++) {
				if (ys[x] < ys[index]) {
					index = x;
				}
			}

			double truncLogConstant = cs[index];
			double truncLogShift = shifts[index];
			double truncLogMean = mus[index];
			double truncLogStdev = sigmas[index];
			double truncLogRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));

			// Fill variables distribution_type, distribution_parameters, and distribution_rmsd
			distribution_parameters[0] = truncLogShift;
			distribution_parameters[1] = truncLogMean;
			distribution_parameters[2] = truncLogStdev;
			distribution_rmsd = truncLogRmsd;

			// XXX this is just here for tests.
			System.out.println("Distribution type: " + distribution_type);
			System.out.println("Shift: " + distribution_parameters[0]);
			System.out.println("Mean (log): " + distribution_parameters[1]);
			System.out.println("Stdev (log): " + distribution_parameters[2]);
			System.out.println("RMSD: " + distribution_rmsd);

			// XXX todo: Find a standard way to use the normaliser!
			normaliser = truncLogConstant;

			// XXX todo: The below distribution must be implemented, remove the below.
			// return new XXX(distribution_parameters[0], distribution_parameters[1], distribution_parameters[2], distribution_parameters[3], distribution_parameters[4], distribution_parameters[5]);
			return null;
	
		} else if (distribution_type == "TruncatedGamma") {

			// Set the number of Nelder Mead Downhill Simplex repetitions.
			int nmRepetitions = 10;

			// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
			double[] cs = new double[nmRepetitions];
			double[] shifts = new double[nmRepetitions];
			double[] ks = new double[nmRepetitions];
			double[] thetas = new double[nmRepetitions];
			double[] gammaKs = new double[nmRepetitions];
			double[] ys = new double[nmRepetitions];
			
			for (int x = 0; x < nmRepetitions; x++) {

				if (cancel != false) {
					return null;
				}
	
				// Initiate the simplex, find 5 vertices.
				// vertex0
				double vertex0c____ = 0.5 + Math.random();
				double vertex0k____ = 1 + Math.random()*3;
                double vertex0shift = 10 + Math.random()*10;
                double vertex0theta = (Math.random()*90 + 10)/vertex0k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				double k = vertex0k____;
				double lanczos_z = k-1;
				double lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				double lanczos_t = lanczos_z + lanczos_g + 0.5;
				double gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex0gammaK = gammaK;
				double vertex0Y = 0;

				// vertex1
				double vertex1c____ = 0.5 + Math.random();
				double vertex1k____ = 1 + Math.random()*3;
                double vertex1shift = 10 + Math.random()*10;
                double vertex1theta = (Math.random()*90 + 10)/vertex1k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				k = vertex1k____;
				lanczos_z = k-1;
				lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				lanczos_t = lanczos_z + lanczos_g + 0.5;
				gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex1gammaK = gammaK;
				double vertex1Y = 0;

				// vertex2
				double vertex2c____ = 0.5 + Math.random();
				double vertex2k____ = 1 + Math.random()*3;
                double vertex2shift = 10 + Math.random()*10;
                double vertex2theta = (Math.random()*90 + 10)/vertex2k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				k = vertex2k____;
				lanczos_z = k-1;
				lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				lanczos_t = lanczos_z + lanczos_g + 0.5;
				gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex2gammaK = gammaK;
				double vertex2Y = 0;

				// vertex3
				double vertex3c____ = 0.5 + Math.random();
				double vertex3k____ = 1 + Math.random()*3;
                double vertex3shift = 10 + Math.random()*10;
                double vertex3theta = (Math.random()*90 + 10)/vertex3k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				k = vertex3k____;
				lanczos_z = k-1;
				lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				lanczos_t = lanczos_z + lanczos_g + 0.5;
				gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex3gammaK = gammaK;
				double vertex3Y = 0;

				// vertex4
				double vertex4c____ = 0.5 + Math.random();
				double vertex4k____ = 1 + Math.random()*3;
                double vertex4shift = 10 + Math.random()*10;
                double vertex4theta = (Math.random()*90 + 10)/vertex4k____;
				// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
				// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
				k = vertex4k____;
				lanczos_z = k-1;
				lanczos_x = lanczos_p[0];
				for (int i = 1; i <= lanczos_g+1; i++){
					lanczos_x += lanczos_p[i]/(lanczos_z+i);
				}
				lanczos_t = lanczos_z + lanczos_g + 0.5;
				gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
				double vertex4gammaK = gammaK;
				double vertex4Y = 0;

				// Prepare for the Nelder-Mead loop.
				keepGoing = true;
				int stepCounter = 0;

				// Until converged do the loop.
				while (keepGoing == true) {

					// Increment the stepCounter.
					stepCounter += 1;

					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/(Math.pow(vertex0theta,vertex0k____)))*(1/vertex0gammaK)*(Math.pow((ages[i]-first_occurrence_age_min+vertex0shift),(vertex0k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min+vertex0shift)/vertex0theta));
						vertex0Y += Math.pow((fix_probabilities[i]-vertex0c____*temp),2);
					}

					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/(Math.pow(vertex1theta,vertex1k____)))*(1/vertex1gammaK)*(Math.pow((ages[i]-first_occurrence_age_min+vertex1shift),(vertex1k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min+vertex1shift)/vertex1theta));
						vertex1Y += Math.pow((fix_probabilities[i]-vertex1c____*temp),2);
					}

					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/(Math.pow(vertex2theta,vertex2k____)))*(1/vertex2gammaK)*(Math.pow((ages[i]-first_occurrence_age_min+vertex2shift),(vertex2k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min+vertex2shift)/vertex2theta));
						vertex2Y += Math.pow((fix_probabilities[i]-vertex2c____*temp),2);
					}

					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/(Math.pow(vertex3theta,vertex3k____)))*(1/vertex3gammaK)*(Math.pow((ages[i]-first_occurrence_age_min+vertex3shift),(vertex3k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min+vertex3shift)/vertex3theta));
						vertex3Y += Math.pow((fix_probabilities[i]-vertex3c____*temp),2);
					}

					// vertex4
					vertex4Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1.0/(Math.pow(vertex4theta,vertex4k____)))*(1/vertex4gammaK)*(Math.pow((ages[i]-first_occurrence_age_min+vertex4shift),(vertex4k____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min+vertex4shift)/vertex4theta));
						vertex4Y += Math.pow((fix_probabilities[i]-vertex4c____*temp),2);
					}

					// Find the best (=lowest) y value.
					double bestY = vertex0Y;
					if (vertex1Y < bestY) {
						bestY = vertex1Y;
					}
					if (vertex2Y < bestY) {
						bestY = vertex2Y;
					}
					if (vertex3Y < bestY) {
						bestY = vertex3Y;
					}
					if (vertex4Y < bestY) {
						bestY = vertex4Y;
					}
					
					// Find the worst (=highest) y value.
					double worstY = vertex0Y;
					if (vertex1Y > worstY) {
						worstY = vertex1Y;
					}
					if (vertex2Y > worstY) {
						worstY = vertex2Y;
					}
					if (vertex3Y > worstY) {
						worstY = vertex3Y;
					}
					if (vertex4Y > worstY) {
						worstY = vertex4Y;
					}

					// Find the second-worst (=second-highest) y value.
					double secondWorstY = bestY;
					if (vertex0Y > secondWorstY && vertex0Y != worstY) {
						secondWorstY = vertex0Y;
					}
					if (vertex1Y > secondWorstY && vertex1Y != worstY) {
						secondWorstY = vertex1Y;
					}
					if (vertex2Y > secondWorstY && vertex2Y != worstY) {
						secondWorstY = vertex2Y;
					}
					if (vertex3Y > secondWorstY && vertex3Y != worstY) {
						secondWorstY = vertex3Y;
					}
					if (vertex4Y > secondWorstY && vertex4Y != worstY) {
						secondWorstY = vertex4Y;
					}

					// Find the parameter values of the best vertex.
					double bestc____ = 0;
					double bestk____ = 0;
					double bestshift = 0;
					double besttheta = 0;
					if (vertex0Y == bestY) {
						bestc____ = vertex0c____;
						bestk____ = vertex0k____;
						bestshift = vertex0shift;
						besttheta = vertex0theta;
					} else if (vertex1Y == bestY) {
						bestc____ = vertex1c____;
						bestk____ = vertex1k____;
						bestshift = vertex1shift;
						besttheta = vertex1theta;
					} else if (vertex2Y == bestY) {
						bestc____ = vertex2c____;
						bestk____ = vertex2k____;
						bestshift = vertex2shift;
						besttheta = vertex2theta;
					} else if (vertex3Y == bestY) {
						bestc____ = vertex3c____;
						bestk____ = vertex3k____;
						bestshift = vertex3shift;
						besttheta = vertex3theta;
					} else if (vertex4Y == bestY) {
						bestc____ = vertex4c____;
						bestk____ = vertex4k____;
						bestshift = vertex4shift;
						besttheta = vertex4theta;
					}

					// Find the parameter values of the worst vertex.
					double worstc____ = 0;
					double worstk____ = 0;
					double worstshift = 0;
					double worsttheta = 0;
					if (vertex0Y == worstY) {
						worstc____ = vertex0c____;
						worstk____ = vertex0k____;
						worstshift = vertex0shift;
						worsttheta = vertex0theta;
					} else if (vertex1Y == worstY) {
						worstc____ = vertex1c____;
						worstk____ = vertex1k____;
						worstshift = vertex1shift;
						worsttheta = vertex1theta;
					} else if (vertex2Y == worstY) {
						worstc____ = vertex2c____;
						worstk____ = vertex2k____;
						worstshift = vertex2shift;
						worsttheta = vertex2theta;
					} else if (vertex3Y == worstY) {
						worstc____ = vertex3c____;
						worstk____ = vertex3k____;
						worstshift = vertex3shift;
						worsttheta = vertex3theta;
					} else if (vertex4Y == worstY) {
						worstc____ = vertex4c____;
						worstk____ = vertex4k____;
						worstshift = vertex4shift;
						worsttheta = vertex4theta;
					}

					// Calculate the sum of the parameters over all vertices.
					double sumc____ = vertex0c____ + vertex1c____ + vertex2c____ + vertex3c____ + vertex4c____;
					double sumk____ = vertex0k____ + vertex1k____ + vertex2k____ + vertex3k____ + vertex4k____;
					double sumshift = vertex0shift + vertex1shift + vertex2shift + vertex3shift + vertex4shift;
					double sumtheta = vertex0theta + vertex1theta + vertex2theta + vertex3theta + vertex4theta;

					// Calculate the parameter values of the centroid.
					double centroidc____ = (sumc____ - worstc____)/3.0;
					double centroidk____ = (sumk____ - worstk____)/3.0;
					double centroidshift = (sumshift - worstshift)/3.0;
					double centroidtheta = (sumtheta - worsttheta)/3.0;

					// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
					double reflectionc____ = centroidc____ + alph * (centroidc____ - worstc____);
					double reflectionk____ = centroidk____ + alph * (centroidk____ - worstk____);
					if (reflectionk____ < 0) {
						reflectionk____ = 1 + Math.random()*3;
					}
					double reflectionshift = centroidshift + alph * (centroidshift - worstshift);
					if (reflectionshift < 0) {
						reflectionshift = 10 + Math.random()*10;
					}
					double reflectiontheta = centroidtheta + alph * (centroidtheta - worsttheta);
					if (reflectiontheta < 0) {
						reflectiontheta = (Math.random()*90 + 10)/reflectionk____;
					}
					// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
					// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
					k = reflectionk____;
					lanczos_z = k-1;
					lanczos_x = lanczos_p[0];
					for (int i = 1; i <= lanczos_g+1; i++){
						lanczos_x += lanczos_p[i]/(lanczos_z+i);
					}
					lanczos_t = lanczos_z + lanczos_g + 0.5;
					gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
					double reflectiongammaK = gammaK;

					// Calculate the y value of the reflection.
					double reflectionY = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(reflectiontheta,reflectionk____)))*(1/reflectiongammaK)*(Math.pow((ages[i]-first_occurrence_age_min),(reflectionk____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min)/reflectiontheta));
						reflectionY += Math.pow((fix_probabilities[i]-reflectionc____*temp),2);
					}

					// Consider the three cases:
					// i.)   reflection is better than all vertices.
					// ii.)  reflection is better than the second-worst vertex.
					// iii.) reflection is worse than or equally good as the second-worst vertex.

					// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
					if (reflectionY < bestY) {

						// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
						double extensionc____ = centroidc____ + gamm * (centroidc____ - worstc____);
						double extensionk____ = centroidk____ + gamm * (centroidk____ - worstk____);
						if (extensionk____ < 0) {
							extensionk____ = 1 + Math.random()*3;
						}
						double extensionshift = centroidshift + alph * (centroidshift - worstshift);
						if (extensionshift < 0) {
							extensionshift = 10 + Math.random()*10;
						}
						double extensiontheta = centroidtheta + gamm * (centroidtheta - worsttheta);
						if (extensiontheta < 0) {
							extensiontheta = (Math.random()*90 + 10)/extensionk____;
						}
						// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
						// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
						k = extensionk____;
						lanczos_z = k-1;
						lanczos_x = lanczos_p[0];
						for (int i = 1; i <= lanczos_g+1; i++){
							lanczos_x += lanczos_p[i]/(lanczos_z+i);
						}
						lanczos_t = lanczos_z + lanczos_g + 0.5;
						gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
						double extensiongammaK = gammaK;

						// Calculate the y value of the extension.
						double extensionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp =  (1.0/(Math.pow(extensiontheta,extensionk____)))*(1/extensiongammaK)*(Math.pow((ages[i]-first_occurrence_age_min+extensionshift),(extensionk____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min+extensionshift)/extensiontheta));
							extensionY += Math.pow((fix_probabilities[i]-extensionc____*temp),2);
						}

						// Figure out which values to use as replacement for the values of the worst vertex.
						double replacec____ = 0;
						double replacek____ = 0;
						double replaceshift = 0;
						double replacetheta = 0;
						double replacegammaK = 0;
						if (reflectionY < extensionY) {
							replacec____ = reflectionc____;
							replacek____ = reflectionk____;
							replaceshift = reflectionshift;
							replacetheta = reflectiontheta;
							replacegammaK = reflectiongammaK;
						} else {
							replacec____ = extensionc____;
							replacek____ = extensionk____;
							replaceshift = extensionshift;
							replacetheta = extensiontheta;
							replacegammaK = extensiongammaK;
						}

						// Replace the parameter values of the worst vertex with the replacement values.
						if (vertex0Y == worstY) {
							vertex0c____ = replacec____;
							vertex0k____ = replacek____;
							vertex0shift = replaceshift;
							vertex0theta = replacetheta;
							vertex0gammaK = replacegammaK;
						} else if (vertex1Y == worstY) {
							vertex1c____ = replacec____;
							vertex1k____ = replacek____;
							vertex1shift = replaceshift;
							vertex1theta = replacetheta;
							vertex1gammaK = replacegammaK;
						} else if (vertex2Y == worstY) {
							vertex2c____ = replacec____;
							vertex2k____ = replacek____;
							vertex2shift = replaceshift;
							vertex2theta = replacetheta;
							vertex2gammaK = replacegammaK;
						} else if (vertex3Y == worstY) {
							vertex3c____ = replacec____;
							vertex3k____ = replacek____;
							vertex3shift = replaceshift;
							vertex3theta = replacetheta;
							vertex3gammaK = replacegammaK;
						} else if (vertex4Y == worstY) {
							vertex4c____ = replacec____;
							vertex4k____ = replacek____;
							vertex4shift = replaceshift;
							vertex4theta = replacetheta;
							vertex4gammaK = replacegammaK;
						}
						
					// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
					} else if (reflectionY < secondWorstY) {

						if (vertex0Y == worstY) {
							vertex0c____ = reflectionc____;
							vertex0k____ = reflectionk____;
							vertex0shift = reflectionshift;
							vertex0theta = reflectiontheta;
							vertex0gammaK = reflectiongammaK;
						} else if  (vertex1Y == worstY) {
							vertex1c____ = reflectionc____;
							vertex1k____ = reflectionk____;
							vertex1shift = reflectionshift;
							vertex1theta = reflectiontheta;
							vertex1gammaK = reflectiongammaK;
						} else if (vertex2Y == worstY) {
							vertex2c____ = reflectionc____;
							vertex2k____ = reflectionk____;
							vertex2shift = reflectionshift;
							vertex2theta = reflectiontheta;
							vertex2gammaK = reflectiongammaK;
						} else if (vertex3Y == worstY) {
							vertex3c____ = reflectionc____;
							vertex3k____ = reflectionk____;
							vertex3shift = reflectionshift;
							vertex3theta = reflectiontheta;
							vertex3gammaK = reflectiongammaK;
						} else if (vertex4Y == worstY) {
							vertex4c____ = reflectionc____;
							vertex4k____ = reflectionk____;
							vertex4shift = reflectionshift;
							vertex4theta = reflectiontheta;
							vertex4gammaK = reflectiongammaK;
						}
						
					// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
					} else {

						// Calculate the contraction.
						double contractionc____ = centroidc____ + beta * (centroidc____ - worstc____);
						double contractionk____ = centroidk____ + beta * (centroidk____ - worstk____);
						double contractionshift = centroidshift + beta * (centroidshift - worstshift);
						double contractiontheta = centroidtheta + beta * (centroidtheta - worsttheta);
						// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
						// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
						k = contractionk____;
						lanczos_z = k-1;
						lanczos_x = lanczos_p[0];
						for (int i = 1; i <= lanczos_g+1; i++){
							lanczos_x += lanczos_p[i]/(lanczos_z+i);
						}
						lanczos_t = lanczos_z + lanczos_g + 0.5;
						gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
						double contractiongammaK = gammaK;
		
						// Calculate the y value of the contraction.
						double contractionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double temp =  (1.0/(Math.pow(contractiontheta,contractionk____)))*(1/contractiongammaK)*(Math.pow((ages[i]-first_occurrence_age_min+contractionshift),(contractionk____-1)))*(Math.exp(-(ages[i]-first_occurrence_age_min+contractionshift)/contractiontheta));
							contractionY += Math.pow((fix_probabilities[i]-contractionc____*temp),2);
						}

						// Consider two subcases of case iii:
						// iiia) The contraction is better than the worst vertex.
						// iiib) The contraction is not better than the worst vertex.

						// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
						if (contractionY < worstY) {

							// Replace the parameter values of the worst vertex with the contraction values.
							if (vertex0Y == worstY) {
								vertex0c____ = contractionc____;
								vertex0k____ = contractionk____;
								vertex0shift = contractionshift;
								vertex0theta = contractiontheta;
								vertex0gammaK = contractiongammaK;
							} else if (vertex1Y == worstY) {
								vertex1c____ = contractionc____;
								vertex1k____ = contractionk____;
								vertex1shift = contractionshift;
								vertex1theta = contractiontheta;
								vertex1gammaK = contractiongammaK;
							} else if (vertex2Y == worstY) {
								vertex2c____ = contractionc____;
								vertex2k____ = contractionk____;
								vertex2shift = contractionshift;
								vertex2theta = contractiontheta;
								vertex2gammaK = contractiongammaK;
							} else if (vertex3Y == worstY) {
								vertex3c____ = contractionc____;
								vertex3k____ = contractionk____;
								vertex3shift = contractionshift;
								vertex3theta = contractiontheta;
								vertex3gammaK = contractiongammaK;
							} else if (vertex4Y == worstY) {
								vertex4c____ = contractionc____;
								vertex4k____ = contractionk____;
								vertex4shift = contractionshift;
								vertex4theta = contractiontheta;
								vertex4gammaK = contractiongammaK;
							}
							
						// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
						} else {

							// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
							// vertex0
							vertex0c____ = bestc____ + delt * (vertex0c____ - bestc____);
							vertex0k____ = bestk____ + delt * (vertex0k____ - bestk____);
							vertex0shift = bestshift + delt * (vertex0shift - bestshift);
							vertex0theta = besttheta + delt * (vertex0theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex0k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex0gammaK = gammaK;

							// vertex1
							vertex1c____ = bestc____ + delt * (vertex1c____ - bestc____);
							vertex1k____ = bestk____ + delt * (vertex1k____ - bestk____);
							vertex1shift = bestshift + delt * (vertex1shift - bestshift);
							vertex1theta = besttheta + delt * (vertex1theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex1k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex1gammaK = gammaK;

							// vertex2
							vertex2c____ = bestc____ + delt * (vertex2c____ - bestc____);
							vertex2k____ = bestk____ + delt * (vertex2k____ - bestk____);
							vertex2shift = bestshift + delt * (vertex2shift - bestshift);
							vertex2theta = besttheta + delt * (vertex2theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex2k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex2gammaK = gammaK;

							// vertex3
							vertex3c____ = bestc____ + delt * (vertex3c____ - bestc____);
							vertex3k____ = bestk____ + delt * (vertex3k____ - bestk____);
							vertex3shift = bestshift + delt * (vertex3shift - bestshift);
							vertex3theta = besttheta + delt * (vertex3theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex3k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex3gammaK = gammaK;

							// vertex4
							vertex4c____ = bestc____ + delt * (vertex4c____ - bestc____);
							vertex4k____ = bestk____ + delt * (vertex4k____ - bestk____);
							vertex4shift = bestshift + delt * (vertex4shift - bestshift);
							vertex4theta = besttheta + delt * (vertex4theta - besttheta);
							// Wherever a new k is chosen, gamma(k) is also calculated to avoid time-consuming recalculation of it for every value of t.
							// This is done using the Lanczos approximation (http://en.wikipedia.org/wiki/Lanczos_approximation), and only positive values are considered.
							k = vertex4k____;
							lanczos_z = k-1;
							lanczos_x = lanczos_p[0];
							for (int i = 1; i <= lanczos_g+1; i++){
								lanczos_x += lanczos_p[i]/(lanczos_z+i);
							}
							lanczos_t = lanczos_z + lanczos_g + 0.5;
							gammaK = Math.sqrt(2*Math.PI) * Math.pow(lanczos_t,(lanczos_z+0.5)) * Math.exp(-lanczos_t) * lanczos_x;
							vertex4gammaK = gammaK;
							
						} // if (contractionY < worstY)
						
					} // if (reflectionY < bestY)
					
					// Stop the loop when all parameter values are identical in the first 10 decimals.
					keepGoing = false;
					if (Math.abs(vertex0c____ - vertex1c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex1k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex1shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex1theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____ - vertex2c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex2k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex2shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex2theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____ - vertex3c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex3k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex3shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex3theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____ - vertex4c____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex4k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0shift - vertex4shift) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex4theta) > 0.0000001) {
						keepGoing = true;
					}
					if (stepCounter > 100000) {
						keepGoing = false;
					}
					
				} // while (keepGoing == true)
			
				// Report all parameters and y.
				cs[x] = vertex0c____;
				ks[x] = vertex0k____;
				shifts[x] = vertex0shift;
				thetas[x] = vertex0theta;
				gammaKs[x] = vertex0gammaK;
				ys[x] = vertex0Y;
			
			} // for (int x = 0; x < nmRepetitions; x++)
				
			// Find the best result among the nmReplicates replicates.
			int index = 0;
			for (int x = 1; x < nmRepetitions; x++) {
				if (ys[x] < ys[index]) {
					index = x;
				}
			}

			double truncGamConstant = cs[index];
			double truncGamShift = shifts[index];
			double truncGamShape = ks[index];
			double truncGamScale = thetas[index];
			double truncGamRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));
				
			// Fill variables distribution_parameters and distribution_rmsd
			distribution_parameters[0] = truncGamShift;
			distribution_parameters[1] = truncGamShape;
			distribution_parameters[2] = truncGamScale;
			distribution_rmsd = truncGamRmsd;
				
			// XXX this is just here for tests.
			System.out.println("Distribution type: " + distribution_type);
			System.out.println("Shift: " + distribution_parameters[0]);
			System.out.println("Shape: " + distribution_parameters[1]);
			System.out.println("Scale: " + distribution_parameters[2]);
			System.out.println("RMSD: " + distribution_rmsd);

			// XXX todo: Find a standard way to use the normaliser!
			normaliser = truncGamConstant;
			
			// XXX todo: The below distribution must be implemented, remove the below.
			// return new XXX(distribution_parameters[0], distribution_parameters[1], distribution_parameters[2], distribution_parameters[3], distribution_parameters[4], distribution_parameters[5]);
			return null;
	
		} else if (distribution_type == "TruncatedNormal") {

			// Set the number of Nelder Mead Downhill Simplex repetitions.
			int nmRepetitions = 10;

			// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
			double[] cs = new double[nmRepetitions];
			double[] ms = new double[nmRepetitions];
			double[] ws = new double[nmRepetitions];
			double[] ys = new double[nmRepetitions];

			for (int x = 0; x < nmRepetitions; x++) {
				
				if (cancel != false) {
					return null;
				}

				// Initiate the simplex, find 4 vertices.
				// vertex0
				double vertex0c = 0.5 + Math.random();
				double vertex0m = 10 + Math.random()*20;
				double vertex0w = 10 + Math.random()*30;
				double vertex0Y = 0;

				// vertex1
				double vertex1c = 0.5 + Math.random();
				double vertex1m = 10 + Math.random()*20;
				double vertex1w = 10 + Math.random()*30;
				double vertex1Y = 0;

				// vertex2
				double vertex2c = 0.5 + Math.random();
				double vertex2m = 10 + Math.random()*20;
				double vertex2w = 10 + Math.random()*30;
				double vertex2Y = 0;

				// vertex3
				double vertex3c = 0.5 + Math.random();
				double vertex3m = 10 + Math.random()*20;
				double vertex3w = 10 + Math.random()*30;
				double vertex3Y = 0;

				// Prepare for the Nelder-Mead loop.
				keepGoing = true;
				int stepCounter = 0;

				// Until converged do the loop.
				while (keepGoing == true) {

					// Increment the stepCounter.
					stepCounter += 1;

					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						vertex0Y += Math.pow(fix_probabilities[i]-vertex0c*Math.exp(-Math.pow((ages[i]-first_occurrence_age_min)-vertex0m,2)/vertex0w),2);
					}

					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						vertex1Y += Math.pow(fix_probabilities[i]-vertex1c*Math.exp(-Math.pow((ages[i]-first_occurrence_age_min)-vertex1m,2)/vertex1w),2);
					}
					
					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						vertex2Y += Math.pow(fix_probabilities[i]-vertex2c*Math.exp(-Math.pow((ages[i]-first_occurrence_age_min)-vertex2m,2)/vertex2w),2);
					}
					
					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						vertex3Y += Math.pow(fix_probabilities[i]-vertex3c*Math.exp(-Math.pow((ages[i]-first_occurrence_age_min)-vertex3m,2)/vertex3w),2);
					}
					
					// Find the best (=lowest) y value.
					double bestY = vertex0Y;
					if (vertex1Y < bestY) {
						bestY = vertex1Y;
					}
					if (vertex2Y < bestY) {
						bestY = vertex2Y;
					}
					if (vertex3Y < bestY) {
						bestY = vertex3Y;
					}

					// Find the worst (=highest) y value.
					double worstY = vertex0Y;
					if (vertex1Y > worstY) {
						worstY = vertex1Y;
					}
					if (vertex2Y > worstY) {
						worstY = vertex2Y;
					}
					if (vertex3Y > worstY) {
						worstY = vertex3Y;
					}

					// Find the second-worst (=second-highest) y value.
					double secondWorstY = bestY;
					if (vertex0Y > secondWorstY && vertex0Y != worstY) {
						secondWorstY = vertex0Y;
					}
					if (vertex1Y > secondWorstY && vertex1Y != worstY) {
						secondWorstY = vertex1Y;
					}
					if (vertex2Y > secondWorstY && vertex2Y != worstY) {
						secondWorstY = vertex2Y;
					}
					if (vertex3Y > secondWorstY && vertex3Y != worstY) {
						secondWorstY = vertex3Y;
					}
					
					// Find the parameter values of the best vertex.
					double bestc = 0;
					double bestm = 0;
					double bestw = 0;
					if (vertex0Y == bestY) {
						bestc = vertex0c;
						bestm = vertex0m;
						bestw = vertex0w;
					} else if (vertex1Y == bestY) {
						bestc = vertex1c;
						bestm = vertex1m;
						bestw = vertex1w;
					} else if (vertex2Y == bestY) {
						bestc = vertex2c;
						bestm = vertex2m;
						bestw = vertex2w;
					} else if (vertex3Y == bestY) {
						bestc = vertex3c;
						bestm = vertex3m;
						bestw = vertex3w;
					}
					
					// Find the parameter values of the worst vertex.
					double worstc = 0;
					double worstm = 0;
					double worstw = 0;
					if (vertex0Y == worstY) {
						worstc = vertex0c;
						worstm = vertex0m;
						worstw = vertex0w;
					} else if (vertex1Y == worstY) {
						worstc = vertex1c;
						worstm = vertex1m;
						worstw = vertex1w;
					} else if (vertex2Y == worstY) {
						worstc = vertex2c;
						worstm = vertex2m;
						worstw = vertex2w;
					} else if (vertex3Y == worstY) {
						worstc = vertex3c;
						worstm = vertex3m;
						worstw = vertex3w;
					}

					// Calculate the sum of the parameters over all vertices.
					double sumc = vertex0c + vertex1c + vertex2c + vertex3c;
					double summ = vertex0m + vertex1m + vertex2m + vertex3m;
					double sumw = vertex0w + vertex1w + vertex2w + vertex3w;

					// Calculate the parameter values of the centroid.
					double centroidc = (sumc - worstc)/3.0;
					double centroidm = (summ - worstm)/3.0;
					double centroidw = (sumw - worstw)/3.0;

                    // Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
                    double reflectionc = centroidc + alph * (centroidc - worstc);
                    double reflectionm = centroidm + alph * (centroidm - worstm);
                    double reflectionw = centroidw + alph * (centroidw - worstw);

					// Calculate the y value of the reflection.
					// vertex0
					double reflectionY = 0;
					for (int i = 0; i < ages.length; i++) {
						reflectionY += Math.pow(fix_probabilities[i]-reflectionc*Math.exp(-Math.pow((ages[i]-first_occurrence_age_min)-reflectionm,2)/reflectionw),2);
					}

					// Consider the three cases:
					// i.)   reflection is better than all vertices.
					// ii.)  reflection is better than the second-worst vertex.
					// iii.) reflection is worse than or equally good as the second-worst vertex.

					// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
					if (reflectionY < bestY) {

                        // Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
						double extensionc = centroidc + gamm * (centroidc - worstc);
                        double extensionm = centroidm + gamm * (centroidm - worstm);
                        double extensionw = centroidw + gamm * (centroidw - worstw);

    					// Calculate the y value of the extension.
    					// vertex0
    					double extensionY = 0;
    					for (int i = 0; i < ages.length; i++) {
    						extensionY += Math.pow(fix_probabilities[i]-extensionc*Math.exp(-Math.pow((ages[i]-first_occurrence_age_min)-extensionm,2)/extensionw),2);
    					}
						
						// Figure out which values to use as replacement for the values of the worst vertex.
						double replacec = 0;
						double replacem = 0;
						double replacew = 0;
						if (reflectionY < extensionY) {
							replacec = reflectionc;
							replacem = reflectionm;
							replacew = reflectionw;
						} else {
							replacec = extensionc;
							replacem = extensionm;
							replacew = extensionw;
						}

						// Replace the parameter values of the worst vertex with the replacement values.
						if (vertex0Y == worstY) {
							vertex0c = replacec;
							vertex0m = replacem;
							vertex0w = replacew;
						} else if (vertex1Y == worstY) {
							vertex1c = replacec;
							vertex1m = replacem;
							vertex1w = replacew;
						} else if (vertex2Y == worstY) {
							vertex2c = replacec;
							vertex2m = replacem;
							vertex2w = replacew;
						} else if (vertex3Y == worstY) {
							vertex3c = replacec;
							vertex3m = replacem;
							vertex3w = replacew;
						}

					// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
					} else if (reflectionY < secondWorstY) {

						if (vertex0Y == worstY) {
							vertex0c = reflectionc;
							vertex0m = reflectionm;
							vertex0w = reflectionw;
						} else if  (vertex1Y == worstY) {
							vertex1c = reflectionc;
							vertex1m = reflectionm;
							vertex1w = reflectionw;
						} else if (vertex2Y == worstY) {
							vertex2c = reflectionc;
							vertex2m = reflectionm;
							vertex2w = reflectionw;
						} else if (vertex3Y == worstY) {
							vertex3c = reflectionc;
							vertex3m = reflectionm;
							vertex3w = reflectionw;
						}
												
					// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
					} else {

                        // Calculate the contraction.
                        double contractionc = centroidc + beta * (centroidc - worstc);
                        double contractionm = centroidm + beta * (centroidm - worstm);
                        double contractionw = centroidw + beta * (centroidw - worstw);

    					// Calculate the y value of the contraction.
    					// vertex0
    					double contractionY = 0;
    					for (int i = 0; i < ages.length; i++) {
    						contractionY += Math.pow(fix_probabilities[i]-contractionc*Math.exp(-Math.pow((ages[i]-first_occurrence_age_min)-contractionm,2)/contractionw),2);
    					}

						// Consider two subcases of case iii:
						// iiia) The contraction is better than the worst vertex.
						// iiib) The contraction is not better than the worst vertex.

						// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
						if (contractionY < worstY) {

							// Replace the parameter values of the worst vertex with the contraction values.
							if (vertex0Y == worstY) {
								vertex0c = contractionc;
								vertex0m = contractionm;
								vertex0w = contractionw;
							} else if (vertex1Y == worstY) {
								vertex1c = contractionc;
								vertex1m = contractionm;
								vertex1w = contractionw;
							} else if (vertex2Y == worstY) {
								vertex2c = contractionc;
								vertex2m = contractionm;
								vertex2w = contractionw;
							} else if (vertex3Y == worstY) {
								vertex3c = contractionc;
								vertex3m = contractionm;
								vertex3w = contractionw;
							}

						// Case iiib): If the contraction is not better than the worst vertex, bring all vertices closer to the best vertex.
						} else {

                            // Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
                            // vertex0
                            vertex0c = bestc + delt * (vertex0c - bestc);
                            vertex0m = bestm + delt * (vertex0m - bestm);
                            vertex0w = bestw + delt * (vertex0w - bestw);
		
                            // vertex1
                            vertex1c = bestc + delt * (vertex1c - bestc);
                            vertex1m = bestm + delt * (vertex1m - bestm);
                            vertex1w = bestw + delt * (vertex1w - bestw);
							
                            // vertex2
                            vertex2c = bestc + delt * (vertex2c - bestc);
                            vertex2m = bestm + delt * (vertex2m - bestm);
                            vertex2w = bestw + delt * (vertex2w - bestw);
							
                            // vertex3
                            vertex3c = bestc + delt * (vertex3c - bestc);
                            vertex3m = bestm + delt * (vertex3m - bestm);
                            vertex3w = bestw + delt * (vertex3w - bestw);
							
						} // if (contractionY < worstY)
    					
					} // if (reflectionY < bestY)

					// Stop the loop when all parameter values are identical in the first 10 decimals.
					keepGoing = false;
					if (Math.abs(vertex0c - vertex1c) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0m - vertex1m) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0w- vertex1w) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c - vertex2c) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0m - vertex2m) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0w - vertex2w) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c - vertex3c) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0m - vertex3m) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0w - vertex3w) > 0.0000001) {
						keepGoing = true;
					}
					if (stepCounter > 100000) {
						keepGoing = false;
					}
					
				} // while (keepGoing == true)
													
				// Report all parameters and y.
				cs[x] = vertex0c;
				ms[x] = vertex0m;
				ws[x] = vertex0w;
				ys[x] = vertex0Y;

			} // for (int x = 0; x < nmRepetitions; x++)

			// Find the best result among the nmReplicates replicates.
			int index = 0;
			for (int x = 1; x < nmRepetitions; x++) {
				if (ys[x] < ys[index]) {
					index = x;
				}
			}

			double normConstant = cs[index];
			double normMean = ms[index];
			double normStdev = Math.sqrt(ws[index]/2.0);
			double gamRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));

			// Fill variables distribution_parameters and distribution_rmsd
			distribution_parameters[0] = normMean;
			distribution_parameters[1] = normStdev;
			distribution_rmsd = gamRmsd;

			// XXX this is just here for tests.
			System.out.println("Distribution type: " + distribution_type);
			System.out.println("Mean: " + distribution_parameters[0]);
			System.out.println("StDev: " + distribution_parameters[1]);
			System.out.println("RMSD: " + distribution_rmsd);

			// XXX todo: Find a standard way to use the normaliser!
			normaliser = normConstant;
			
			// XXX todo: The below distribution must be implemented, remove the below.
			// return new XXX(distribution_parameters[0], distribution_parameters[1], distribution_parameters[2], distribution_parameters[3], distribution_parameters[4], distribution_parameters[5]);
			return null;
	
		} // if (distribution_type == "Lognormal")

		// XXX todo: The below distribution must be implemented, remove the below.
		// return new XXX(distribution_parameters[0], distribution_parameters[1], distribution_parameters[2], distribution_parameters[3], distribution_parameters[4], distribution_parameters[5]);
		return null;
				
	}

	public static void main(String[] args) {

		int z = 0;
		for (int x = 0; x < 1000000; x += 1) {
			if (z == 1000) {
				System.out.println(System.currentTimeMillis());
				z = 0;
			} else {
				z += 1;
			}
			CladeAgeProbabilities cladeAgeProbabilities = new CladeAgeProbabilities();
			cladeAgeProbabilities.run_fitted_cladeage_rapid(10.0,10.0,0.04,0.04,0.2,0.2,0.02,0.02);
		}

  }

	public double getNormaliser() {
		return normaliser;
	}
}
