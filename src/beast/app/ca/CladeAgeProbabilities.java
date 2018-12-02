package beast.app.ca;

import java.util.ArrayList;

import javax.swing.JProgressBar;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import beast.math.distributions.LogNormalDistributionModel;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.NormalDistributionImpl;

// Import CladeAge distributions.
import beast.math.distributions.CladeAgeDistribution;

public class CladeAgeProbabilities {

	// Initialize objects that are to become variables of class instances.
	private int number_of_ages = 100;
	private double[] ages = new double[number_of_ages+1];
	private double[] ages_younger = new double[number_of_ages+1];
	private double[] ages_older = new double[number_of_ages+1];
	private double[] fix_probabilities = new double[number_of_ages+1];
	private double[] int_probabilities = new double[number_of_ages+1];
	private double[] int_probabilities_younger = new double[number_of_ages+1];
	private double[] int_probabilities_older = new double[number_of_ages+1];
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

	public CladeAgeDistribution run_standard_cladeage(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, JProgressBar dpb) throws Exception {

		// Make sure cancel is set to false.
		cancel = false;

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		int_probabilities = new double[number_of_ages + 1];
		double[] raw_probabilities = new double[number_of_ages + 1];
		
		// Part 1: Calculate the maximum age, and 100 time points between fossil age and maximum age.
		
		// Roughly calculate a maximum tree duration so that the probability for the maximum tree duration is prob_ratio (e.g. 0.001) * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		double ndr_mean = (ndr_max+ndr_min)/2.0;
		double psi_mean = (psi_max+psi_min)/2.0;
		double prob_ratio = 0.001;

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
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(prob_ratio)), 2);

			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
			
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
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
				
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
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
				
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
		double max_simulation_age = vertex0t + first_occurrence_age_max;

		// Determine the ages.
		ages[0] = first_occurrence_age_min;
		for (int i = 0; i < number_of_ages; i++) {
			ages[i+1] = first_occurrence_age_min + ((i+1)/(double) number_of_ages)*(max_simulation_age-first_occurrence_age_min);
		}

		// Part 2: Estimate clade age probability densities
		
		// Prepare the arrays raw_probabilities and successful simulations (those in which at least one taxon survives), and fill them with 0.
		int[] successful_simulations = new int[number_of_ages + 1];
		
		while (successful_simulations[successful_simulations.length-1] < 10000) {
			if (cancel) {
				return null;
			}
			// Increment the progress indicator.
			if (dpb != null) {
				dpb.setValue(successful_simulations[successful_simulations.length-1]/100);
			}

			// Draw values for parameters ndr (net diversification rate, lambda - mu) and epsilon (turnover rate, mu/lambda) from uniform distributions, and calculate lambda and mu from them.
			double ndr = ndr_min + Math.random() * (ndr_max-ndr_min);
			double epsilon = epsilon_min + Math.random() * (epsilon_max-epsilon_min);
			double mu = (ndr*epsilon)/(1.0 - epsilon);
			double lambda = ndr + mu;

			// Draw a duration t for this simulation.
			double first_occurrence_age = first_occurrence_age_min + Math.random()*(first_occurrence_age_max-first_occurrence_age_min);
			
			// Determine the maximum tree duration.
			double max_tree_duration = ages[ages.length-1] - first_occurrence_age;
			
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
							sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
						}
						if (remove_this_branch == true) {
							branch_origin.remove(o);
							branch_termination.remove(o);
						}
					} // for (int o = branch_origin.size()-1; o >= 0; o--)
					
					// If at least one species is extant at this age, get the sum of lineage durations.
					if (number_of_extant_taxa > 0) {
						// Draw a value for the sampling rate psi.
						for (int pp = 0; pp < 10; pp++) {
							double psi = psi_min + Math.random()*(psi_max-psi_min);
							raw_probabilities[i] += (psi*Math.exp(-psi*sum_of_species_durations)*number_of_extant_taxa)/10.0;
						}
						successful_simulations[i] += 1;
					} // if (extant_at_this_age == true)

				} else {
					successful_simulations[i] += 1;
				} // if (tree_duration > 0)
				
			} // for (int i = 0; i < ages.length; i++)
			
		} // while (successful_simulations[successful_simulations.length-1] < 10000)
		
		for (int i = 0; i < raw_probabilities.length; i ++) {
			int_probabilities[i] = raw_probabilities[i]/(double) successful_simulations[i];
		}
		
		// Return the CladeAge distribution.
		return new CladeAgeDistribution(ages, int_probabilities, true);

	} // public void run_standard_cladeage(...)

	public CladeAgeDistribution run_duo_cladeage(double younger_first_occurrence_age_min, double younger_first_occurrence_age_max, double younger_first_occurrence_weight, double older_first_occurrence_age_min, double older_first_occurrence_age_max, double older_first_occurrence_weight, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, JProgressBar dpb) throws Exception {

		// Make sure cancel is set to false.
		cancel = false;

		// Determine the weight proportions for the younger and older fossils.
		double first_occurrence_weight_sum = younger_first_occurrence_weight + older_first_occurrence_weight;
		double younger_first_occurrence_weight_proportion = younger_first_occurrence_weight/first_occurrence_weight_sum;
		double older_first_occurrence_weight_proportion = older_first_occurrence_weight/first_occurrence_weight_sum;

		// Check that the younger first occurrence is actually younger than the older one.
		if (older_first_occurrence_age_min < younger_first_occurrence_age_min) {
			throw new Exception("The minimum age for the older first occurrence is younger than the minimum age for the younger first occurrence!");
		}

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		int_probabilities = new double[number_of_ages + 1];
		ages_younger = new double[number_of_ages + 1];
		int_probabilities_younger = new double[number_of_ages + 1];
		ages_older = new double[number_of_ages + 1];
		int_probabilities_older = new double[number_of_ages + 1];

		// Part 1: Calculate the maximum age, and 100 time points between fossil age and maximum age.
		
		// Roughly calculate a maximum tree duration so that the probability for the maximum tree duration is prob_ratio (e.g. 0.001) * the probability of age=0.
		// The approximation s(t) = exp(ndr*t) is used. This leads to
		// p(t) = psi * exp(ndr*t) * exp(-(psi/ndr) * (exp(ndr*t) - 1)) = 0.1 * p(0) = 0.1 * psi.
		
		double ndr_mean = (ndr_max+ndr_min)/2.0;
		double psi_mean = (psi_max+psi_min)/2.0;
		double prob_ratio = 0.001;
		
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
			double vertex0Y = Math.pow((Math.pow(ndr_mean,2)*vertex0t - psi_mean*(Math.exp(ndr_mean*vertex0t) - 1) - ndr_mean*Math.log(prob_ratio)), 2);

			// vertex1
			double vertex1Y = Math.pow((Math.pow(ndr_mean,2)*vertex1t - psi_mean*(Math.exp(ndr_mean*vertex1t) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
			
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
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
				double extensiont = bestt + gamm * (bestt - worstt);
				
				// Calculate the y value of the extension.
				double extensionY = Math.pow((Math.pow(ndr_mean,2)*extensiont - psi_mean*(Math.exp(ndr_mean*extensiont) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
				
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
				double contractionY = Math.pow((Math.pow(ndr_mean,2)*contractiont - psi_mean*(Math.exp(ndr_mean*contractiont) - 1) - ndr_mean*Math.log(prob_ratio)), 2);
				
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

		// Memorize the result as the max_simulation_ages for both the younger and the older fossil.
		double max_simulation_age_younger = vertex0t + younger_first_occurrence_age_max;
		double max_simulation_age_older = vertex0t + older_first_occurrence_age_max;

		// Determine the ages used for the younger and older fossils.
		ages_younger[0] = younger_first_occurrence_age_min;
		for (int i = 0; i < number_of_ages; i++) {
			ages_younger[i+1] = younger_first_occurrence_age_min + ((i+1)/(double) number_of_ages)*(max_simulation_age_younger-younger_first_occurrence_age_min);
		}
		ages_older[0] = older_first_occurrence_age_min;
		for (int i = 0; i < number_of_ages; i++) {
			ages_older[i+1] = older_first_occurrence_age_min + ((i+1)/(double) number_of_ages)*(max_simulation_age_older-older_first_occurrence_age_min);
		}

		// Part 2: Estimate clade age probability densities.
		
		// Repeat for both the younger and older fossil.
		String[] fossil_types = new String[] { "younger", "older" };
		for (String fossil_type : fossil_types) {

			// Reset the array for raw probabilities.
			double[] raw_probabilities = new double[number_of_ages + 1];

			// Prepare the arrays raw_probabilities and successful simulations (those in which at least one taxon survives), and fill them with 0.
			int[] successful_simulations = new int[number_of_ages + 1];
		
			while (successful_simulations[successful_simulations.length-1] < 10000) {
				if (cancel) {
					return null;
				}
				// Increment the progress indicator.
				if (dpb != null) {
					dpb.setValue(successful_simulations[successful_simulations.length-1]/100);
				}

				// Draw values for parameters ndr (net diversification rate, lambda - mu) and epsilon (turnover rate, mu/lambda) from uniform distributions, and calculate lambda and mu from them.
				double ndr = ndr_min + Math.random() * (ndr_max-ndr_min);
				double epsilon = epsilon_min + Math.random() * (epsilon_max-epsilon_min);
				double mu = (ndr*epsilon)/(1.0 - epsilon);
				double lambda = ndr + mu;

				// Draw a duration t for this simulation.
				double first_occurrence_age = 0;
				if (fossil_type == "younger") {
					first_occurrence_age = younger_first_occurrence_age_min + Math.random()*(younger_first_occurrence_age_max-younger_first_occurrence_age_min);
				} else if (fossil_type == "older") {
					first_occurrence_age = older_first_occurrence_age_min + Math.random()*(older_first_occurrence_age_max-older_first_occurrence_age_min);
				} else {
					throw new Exception("Expected fossil type to be either 'younger' or 'older'!");
				}
			
				// Determine the maximum tree duration.
				double max_tree_duration = 0;
				if (fossil_type == "younger") {
					max_tree_duration = ages_younger[ages_younger.length-1] - first_occurrence_age;
				} else if (fossil_type == "older") {
					max_tree_duration = ages_older[ages_older.length-1] - first_occurrence_age;
				} else {
					throw new Exception("Expected fossil type to be either 'younger' or 'older'!");
				}
			
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
				int ages_length = 0;
				if (fossil_type == "younger") {
					ages_length = ages_younger.length;
				} else if (fossil_type == "older") {
					ages_length = ages_older.length;
				} else {
					throw new Exception("Expected fossil type to be either 'younger' or 'older'!");
				}
				for (int i = ages_length-1; i >= 0; i--) {
				
					// For each age, determine the tree duration.
					double tree_duration = 0;
					if (fossil_type == "younger") {
						tree_duration = ages_younger[i] - first_occurrence_age;
					} else if (fossil_type == "older") {
						tree_duration = ages_older[i] - first_occurrence_age;
					} else {
						throw new Exception("Expected fossil type to be either 'younger' or 'older'!");
					}
				
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
								sum_of_species_durations += branch_termination.get(o)-branch_origin.get(o);
							}
							if (remove_this_branch == true) {
								branch_origin.remove(o);
								branch_termination.remove(o);
							}
						} // for (int o = branch_origin.size()-1; o >= 0; o--)
					
						// If at least one species is extant at this age, get the sum of lineage durations.
						if (number_of_extant_taxa > 0) {
							// Draw a value for the sampling rate psi.
							for (int pp = 0; pp < 10; pp++) {
								double psi = psi_min + Math.random()*(psi_max-psi_min);
								raw_probabilities[i] += (psi*Math.exp(-psi*sum_of_species_durations)*number_of_extant_taxa)/10.0;
							}
							successful_simulations[i] += 1;
						} // if (extant_at_this_age == true)

					} else {
						successful_simulations[i] += 1;
					} // if (tree_duration > 0)
				
				} // for (int i = 0; i < ages.length; i++)
			
			} // while (successful_simulations[successful_simulations.length-1] < 10000)
		
			if (fossil_type == "younger") {
				for (int i = 0; i < raw_probabilities.length; i ++) {
					int_probabilities_younger[i] = raw_probabilities[i]/(double) successful_simulations[i];
				}
			} else if (fossil_type == "older") {
				for (int i = 0; i < raw_probabilities.length; i ++) {
					int_probabilities_older[i] = raw_probabilities[i]/(double) successful_simulations[i];
				}
			} else {
				throw new Exception("Expected fossil type to be either 'younger' or 'older'!");
			}

		} // for (String fossil_type : fossil_types)

		// Determine the ages used for the joint distribution.
		double min_simulation_age = younger_first_occurrence_age_min;
		if (older_first_occurrence_age_min < younger_first_occurrence_age_min) {
			min_simulation_age = older_first_occurrence_age_min;
		}
		double max_simulation_age = max_simulation_age_younger;
		if (max_simulation_age_older > max_simulation_age_younger) {
			max_simulation_age = max_simulation_age_older;
		}
		ages[0] = min_simulation_age;
		for (int i = 0; i < number_of_ages; i++) {
			ages[i+1] = min_simulation_age + ((i+1)/(double) number_of_ages)*(max_simulation_age-min_simulation_age);
		}

		// Determine the added weighted probabilities of both distributions.
		for (int i = 0; i < ages.length; i ++) {
			double int_probability_younger = -1;
			double int_probability_older = -1;

			// Check if ages[i] is to the left of the leftmost ages of the distribution for the younger fossil  (this should actually never happen).
			if (ages[i] < ages_younger[0]) {
				int_probability_younger = 0;
			}

			// Check if ages[i] is to the left of the leftmost ages of the distribution for the older fossil.
			if (ages[i] < ages_older[0]) {
				int_probability_older = 0;
			}

			// Check if ages[i] is to the right of the rightmost ages of the distribution for the younger fossil.
			if (ages[i] > ages_younger[ages_younger.length-1]) {

				// Assume an exponential distribution defined by the last two points.
				// From the next two equations, exp_lambda and exp_scaler can be calculated.
				// - int_probabilities_younger[ages_younger.length-2] = exp_lambda * exp_scaler
				//   int_probabilities_younger[ages_younger.length-1] = exp_lambda * exp_scaler * Math.exp(-exp_lambda*delta_age)
				// - exp_scaler = int_probabilities_younger[ages_younger.length-2] / exp_lambda
				//   exp_scaler = int_probabilities_younger[ages_younger.length-1] / (exp_lambda * Math.exp(-exp_lambda*delta_age))
				// int_probabilities_younger[ages_younger.length-2] / exp_lambda = int_probabilities_younger[ages_younger.length-1] / (exp_lambda * Math.exp(-exp_lambda*delta_age))
				// int_probabilities_younger[ages_younger.length-2] * (exp_lambda * Math.exp(-exp_lambda*delta_age)) = int_probabilities_younger[ages_younger.length-1] * exp_lambda
				// (exp_lambda * Math.exp(-exp_lambda*delta_age)) / exp_lambda = int_probabilities_younger[ages_younger.length-1] / int_probabilities_younger[ages_younger.length-2]
				// Math.exp(-exp_lambda*delta_age) = 
				// -exp_lambda*delta_age = Math.log(int_probabilities_younger[ages_younger.length-1] / int_probabilities_younger[ages_younger.length-2])
				// -exp_lambda = Math.log(int_probabilities_younger[ages_younger.length-1] / int_probabilities_younger[ages_younger.length-2]) / delta_age
				// exp_lambda = -Math.log(int_probabilities_younger[ages_younger.length-1] / int_probabilities_younger[ages_younger.length-2]) / delta_age
				double delta_age1 = ages_younger[ages_younger.length-1] - ages_younger[ages_younger.length-2];
				double exp_lambda1 = Math.log(int_probabilities_younger[ages_younger.length-2] / int_probabilities_younger[ages_younger.length-1]) / delta_age1;
				double exp_scaler1 = int_probabilities_younger[ages_younger.length-2] / exp_lambda1;
				double int_probability_younger1 = exp_lambda1 *  exp_scaler1 * Math.exp(-exp_lambda1*(ages[i]-ages_younger[ages_younger.length-2]));
				double delta_age2 = ages_younger[ages_younger.length-1] - ages_younger[ages_younger.length-3];
				double exp_lambda2 = Math.log(int_probabilities_younger[ages_younger.length-3] / int_probabilities_younger[ages_younger.length-1]) / delta_age2;
				double exp_scaler2 = int_probabilities_younger[ages_younger.length-3] / exp_lambda2;
				double int_probability_younger2 = exp_lambda2 *  exp_scaler2 * Math.exp(-exp_lambda2*(ages[i]-ages_younger[ages_younger.length-3]));
				double delta_age3 = ages_younger[ages_younger.length-1] - ages_younger[ages_younger.length-4];
				double exp_lambda3 = Math.log(int_probabilities_younger[ages_younger.length-4] / int_probabilities_younger[ages_younger.length-1]) / delta_age3;
				double exp_scaler3 = int_probabilities_younger[ages_younger.length-4] / exp_lambda3;
				double int_probability_younger3 = exp_lambda3 *  exp_scaler3 * Math.exp(-exp_lambda3*(ages[i]-ages_younger[ages_younger.length-4]));
				int_probability_younger = (int_probability_younger1 + int_probability_younger2 + int_probability_younger3)/3.0;
			}

			// Check if ages[i] is to the right of the rightmost ages of the distribution for the older fossil (this should actually never happen).
			if (ages[i] > ages_older[ages_older.length-1]) {
				double delta_age1 = ages_older[ages_older.length-1] - ages_older[ages_older.length-2];
				double exp_lambda1 = Math.log(int_probabilities_older[ages_older.length-2] / int_probabilities_older[ages_older.length-1]) / delta_age1;
				double exp_scaler1 = int_probabilities_older[ages_older.length-2] / exp_lambda1;
				double int_probability_older1 = exp_lambda1 *  exp_scaler1 * Math.exp(-exp_lambda1*(ages[i]-ages_older[ages_older.length-2]));
				double delta_age2 = ages_older[ages_older.length-1] - ages_older[ages_older.length-3];
				double exp_lambda2 = Math.log(int_probabilities_older[ages_older.length-3] / int_probabilities_older[ages_older.length-1]) / delta_age2;
				double exp_scaler2 = int_probabilities_older[ages_older.length-3] / exp_lambda2;
				double int_probability_older2 = exp_lambda2 *  exp_scaler2 * Math.exp(-exp_lambda2*(ages[i]-ages_older[ages_older.length-3]));
				double delta_age3 = ages_older[ages_older.length-1] - ages_older[ages_older.length-4];
				double exp_lambda3 = Math.log(int_probabilities_older[ages_older.length-4] / int_probabilities_older[ages_older.length-1]) / delta_age3;
				double exp_scaler3 = int_probabilities_older[ages_older.length-4] / exp_lambda3;
				double int_probability_older3 = exp_lambda3 *  exp_scaler3 * Math.exp(-exp_lambda3*(ages[i]-ages_older[ages_older.length-4]));
				int_probability_older = (int_probability_older1 + int_probability_older2 + int_probability_older3)/3.0;
			}

			// If ages[i] is between two ages of the distribution of the younger fossil, calculate the probability from the probabilities of those two ages.
			if (int_probability_younger == -1) {
				double ages_younger_next_to_left = 0;
				double ages_younger_next_to_right = 0;
				double int_probability_younger_next_to_left = 0;
				double int_probability_younger_next_to_right = 0;
				for (int j = 0; j < ages_younger.length; j ++) {
					if (ages_younger[j] <= ages[i]) {
						ages_younger_next_to_left = ages_younger[j];
						int_probability_younger_next_to_left = int_probabilities_younger[j];
					}
				}

				for (int j = ages_younger.length-1; j >= 0; j --) {
					if (ages_younger[j] >= ages[i]) {
						ages_younger_next_to_right = ages_younger[j];
						int_probability_younger_next_to_right = int_probabilities_younger[j];
					}
				}
				if (ages_younger_next_to_left == ages_younger_next_to_right) {
					int_probability_younger = int_probability_younger_next_to_left;
				} else if (int_probability_younger_next_to_left == int_probability_younger_next_to_right) {
					int_probability_younger = int_probability_younger_next_to_left;
				} else {
					int_probability_younger = int_probability_younger_next_to_left + (ages[i]-ages_younger_next_to_left)/(ages_younger_next_to_right-ages_younger_next_to_left) * (int_probability_younger_next_to_right-int_probability_younger_next_to_left);
				}
			}

			// If ages[i] is between two ages of the distribution of the older fossil, calculate the probability from the probabilities of those two ages.
			if (int_probability_older == -1) {
				double ages_older_next_to_left = 0;
				double ages_older_next_to_right = 0;
				double int_probability_older_next_to_left = 0;
				double int_probability_older_next_to_right = 0;
				for (int j = 0; j < ages_older.length; j ++) {
					if (ages_older[j] <= ages[i]) {
						ages_older_next_to_left = ages_older[j];
						int_probability_older_next_to_left = int_probabilities_older[j];
					}
				}

				for (int j = ages_older.length-1; j >= 0; j --) {
					if (ages_older[j] >= ages[i]) {
						ages_older_next_to_right = ages_older[j];
						int_probability_older_next_to_right = int_probabilities_older[j];
					}
				}
				if (ages_older_next_to_left == ages_older_next_to_right) {
					int_probability_older = int_probability_older_next_to_left;
				} else if (int_probability_older_next_to_left == int_probability_older_next_to_right) {
					int_probability_older = int_probability_older_next_to_left;
				} else {
					int_probability_older = int_probability_older_next_to_left + (ages[i]-ages_older_next_to_left)/(ages_older_next_to_right-ages_older_next_to_left) * (int_probability_older_next_to_right-int_probability_older_next_to_left);
				}
			}

			// Keep the added weighted probability for this age.
			int_probabilities[i] = int_probability_younger * younger_first_occurrence_weight_proportion + int_probability_older * older_first_occurrence_weight_proportion;
		}


		System.out.println("  Clade age probability densities based on younger potential first occurrence:");
		for (int i = 0; i < ages_younger.length; i ++) {
			System.out.println("    " + ages_younger[i] + "\t" + int_probabilities_younger[i]);
		}
		System.out.println("  Clade age probability densities based on older potential first occurrence:");
		for (int i = 0; i < ages_older.length; i ++) {
			System.out.println("    " + ages_older[i] + "\t" + int_probabilities_older[i]);
		}
		System.out.println("  Clade age probability densitites based on both potential first occurrences:");
		for (int i = 0; i < ages.length; i ++) {
			System.out.println("    " + ages[i] + "\t" + int_probabilities[i]);
		}
		System.out.println("");

		// Return the CladeAge distribution.
		return new CladeAgeDistribution(ages, int_probabilities, true);

	} // public void run_duo_cladeage(...)

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
		}

	}

	public double getNormaliser() {
		return normaliser;
	}
}
