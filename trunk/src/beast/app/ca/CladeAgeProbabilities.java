package beast.app.ca;

import java.util.ArrayList;

import javax.swing.JProgressBar;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.math.distributions.ExpGamma;
import beast.math.distributions.LogNormalImpl;



public class CladeAgeProbabilities {

	// Initialize objects that are to become variables of class instances.
	private int number_of_ages = 100;
	private double[] ages = new double[number_of_ages+1];
	private double[] probabilities = new double[number_of_ages+1];
	private double[] fix_probabilities = new double[number_of_ages+1];
	private double[] int_probabilities = new double[number_of_ages+1];
	private double[] approx_ages = new double[number_of_ages+1];
	private double[] approx_probabilities = new double[number_of_ages+1];
	private double[] raw_probabilities = new double[number_of_ages+1];
	private int[] successful_simulations = new int[number_of_ages+1];
	private String approx_distribution_type;
	private double[] approx_distribution_parameters;
	private double approx_distribution_rmsd;
	private double offset = 0;
	private boolean cancel1 = false;
	private double mean_psi = 0;

	// Initialize constants for Lanczos approximation of the gamma function.
	private double lanczos_g = 7;
	private double[] lanczos_p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,771.32342877765313, -176.61502916214059, 12.507343278686905,-0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

	// Initialize constants for the Nelder-Mead Downhill Simplex (according to Wikipedia example).
	private double alph = 1.0;
	private double gamm = 2.0;
	private double beta = -0.5;
	private double delt  = 0.5;
	
	// normaliser for simulated probablities
	private double normaliser = 1.0;
	
	final static int DELTA_PROGRESS = 10;
	final static int nmRepetitions = 10;

	public double getOffset() {
		return offset;
	}

	public void setOffset(double offset) {
		this.offset = offset;
	}

	public double[] getAges() {
		return ages;
	}

	public void setAges(double[] ages) {
		this.ages = ages;
	}

	public double[] getProbabilities() {
		return probabilities;
	}

	public void setProbabilities(double[] probabilities) {
		this.probabilities = probabilities;
	}

	public double[] getApprox_ages() {
		return approx_ages;
	}

	public double[] getApprox_probabilities() {
		return approx_probabilities;
	}

	public String getApprox_distribution_type() {
		return approx_distribution_type;
	}

	public double[] getApprox_distribution_parameters() {
		return approx_distribution_parameters;
	}

	public double getApprox_distribution_rmsd() {
		return approx_distribution_rmsd;
	}

	// public void setSimulation_progress_indicator(???)
	// public ??? getSimulation_progress_indicator()

	// public void setSimulation_progress_indicator_text(???)
	// public ??? getSimulation_progress_indicator_text()

	// public void setSimulation_progress_indicator_box(???)
	// public ??? getSimulation_progress_indicator_box()

	public void setCancel1() {
		cancel1 = true;
	}
	public boolean getCancel1() {
		return cancel1;
	}

	// public void setCancel2(???)
	// public ??? getCancel2()

	public void bd_simulate(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, double sampling_gap_min, double sampling_gap_max, int bd_sample_size, int max_tree_size, int psi_sample_size, JProgressBar dpb) {
		cancel1 = false;

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		probabilities = new double[number_of_ages + 1];
		approx_ages = new double[number_of_ages + 1];
		approx_probabilities = new double[number_of_ages + 1];
		raw_probabilities = new double[number_of_ages + 1];
		
		
		// Memorize the offset and the mean sampling rate.
		offset = first_occurrence_age_min;
		mean_psi = (psi_min+psi_max)/2.0;

		// Roughly calculate a maximum age.
		double max_simulation_age = Math.log(9.210340371976182*((ndr_max+ndr_min)/(psi_max+psi_min)) + 1.0) * (2.0/(ndr_max+ndr_min)) + first_occurrence_age_max + sampling_gap_max;

		// Determine the ages.
		for (int i = 0; i < number_of_ages; i++) {
			ages[i] = max_simulation_age - (i/ (double) number_of_ages)*(max_simulation_age-first_occurrence_age_min);
		}
		ages[number_of_ages] = first_occurrence_age_min;

		while (successful_simulations[0] < bd_sample_size) {
			if (cancel1) {
				return;
			}
			// Update the progress bar.
			if (dpb != null) {
				dpb.setValue(successful_simulations[0]);
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
			double max_tree_duration = ages[0] - first_occurrence_age;
			
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
			
			// Set extant_at_this_age to false before looping through all ages. This can be done because we start with the oldest age, and as soon as we
			// find an age for which extant_at_this_age is true, this will also be true for all younger ages.
			for (int i = 0; i < ages.length; i++) {
				
				// For each age, determine the tree length.
				double tree_duration = ages[i] - first_occurrence_age;
				
				// Trim the tree so that no branches are older than @ages[i] - first_occurrence_age. While trimming, make sure the tree contains at least one extant species.
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
						for (int pp = 0; pp < psi_sample_size; pp++) {
							double psi = psi_min + Math.random()*(psi_max-psi_min);
							if (tree_duration >= sampling_gap) {
								raw_probabilities[i] += (psi*Math.exp(-psi*sum_of_species_durations)*number_of_extant_taxa)/(double) psi_sample_size;
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
			probabilities[i] = raw_probabilities[i]/(double) successful_simulations[i];
		}

	} // public void bd_simulate(...)

	public void bd_simulate_new(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, JProgressBar dpb) {
		cancel1 = false;

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		probabilities = new double[number_of_ages + 1];
		approx_ages = new double[number_of_ages + 1];
		approx_probabilities = new double[number_of_ages + 1];
		raw_probabilities = new double[number_of_ages + 1];
		
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
			double bestY = 0;
			double bestt = 0;
			double worstY = 0;
			double worstt = 0;
			if (vertex0Y < vertex1Y) {
				bestY = vertex0Y;
				bestt = vertex0t;
				worstY = vertex1Y;
				worstt = vertex1t;
			} else {
				bestY = vertex1Y;
				bestt = vertex1t;
				worstY = vertex0Y;
				worstt = vertex0t;
			}
			
			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
			double reflectiont = bestt + alph * (bestt - worstt);
			
			// Calculate the y value of the reflection.
			double reflectionY = Math.pow((Math.pow(ndr_mean,2)*reflectiont - psi_mean*(Math.exp(ndr_mean*reflectiont) - 1) - ndr_mean*Math.log(0.001)),2);
			
			// Consider the three cases:
			// i.)   reflection is better than all vertices.
			// ii.)  reflection is better than the second-worst vertex.
			// iii.) reflection is worse than or equally good as the second-worst vertex.
			
			// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
			if (reflectionY < bestY) {
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
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
				
			// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
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
					
				// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
				} else {
					
					// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
					// vertex0
					vertex0t = bestt + delt  * (vertex0t - bestt);
					
					// vertex1
					vertex1t = bestt + delt  * (vertex1t - bestt);
					
				} // if (contractionY < worstY)
				
			} // if (reflectionY < bestY)
			
			// Stop the loop when all parameter values are identical in the first 5 decimals.
			keepGoing = false;
			if (Math.abs(vertex0t-vertex1t) > 0.00001) {
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

		while (successful_simulations[-1] < 10000) {
			if (cancel1) {
				return;
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
			
			for (int i = ages.length-1; i >= 0; i--) {
				
				// For each age, determine the tree length.
				double tree_duration = ages[i] - first_occurrence_age_min;
				
				// Trim the tree so that no branches are older than @ages[i] - first_occurrence_age. While trimming, make sure the tree contains at least one extant species.
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
		
		// Calculate the fix_probabilities (that's when @first_occurrence_age_min == @first_occurrence_age_max, and int_probabilities which are
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

	} // public void bd_simulate_new(...)

	public ContinuousDistribution fitCladeAge(double first_occurrence_age_min, double first_occurrence_age_max, double psi_min, double psi_max, JProgressBar progress) {
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] ss = new double[nmRepetitions];
		double[] ms = new double[nmRepetitions];
		double[] ws = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];

		for (int x = 0; x < nmRepetitions; x++) {
			
			if (progress != null) {
				progress.setValue(progress.getValue() + DELTA_PROGRESS);
			}
			if (cancel1 != false) {
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
			boolean keepGoing = true;
			int stepCounter = 0;

			// Until converged do the loop.
			while (keepGoing == true) {
				
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

				// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
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
					
					// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
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
                        vertex0w = reflections;
                        vertex0m = reflectionm;
                        vertex0w = reflectionw;
                        vertex0c = reflectionc;
					} else if  (vertex1Y == worstY) {
                        vertex1w = reflections;
                        vertex1m = reflectionm;
                        vertex1w = reflectionw;
                        vertex1c = reflectionc;
					} else if (vertex2Y == worstY) {
                        vertex2w = reflections;
                        vertex2m = reflectionm;
                        vertex2w = reflectionw;
                        vertex2c = reflectionc;
					} else if (vertex3Y == worstY) {
                        vertex3w = reflections;
                        vertex3m = reflectionm;
                        vertex3w = reflectionw;
                        vertex3c = reflectionc;
					}
                	
				// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
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
						
						// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
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
				// if (cancel2 == true) {
				// 	keepGoing = false;
				// }

			} // while (keepGoing == true)

			// Report all parameters and y.
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

        double cladeAgeS = ss[index];
        double cladeAgeM = ms[index];
        double cladeAgeW = ws[index];
        double cladeAgeC = ((psi_min+psi_max)/2.0)/((1.0/(ages[0]-first_occurrence_age_min+cladeAgeS))*Math.exp( -Math.pow((Math.log(ages[0]-first_occurrence_age_min+cladeAgeS)-cladeAgeM),2)/cladeAgeW));
        double cladeAgeRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));

		// Fill arrays approx_ages and approx_probabilities.
		approx_ages = new double[1001];
		approx_probabilities = new double[1001];
		for (int a = 0; a < 1001; a++) {
			approx_ages[a] = first_occurrence_age_min + (a/1000.0)*(ages[-1]-first_occurrence_age_min);
		}
		for (int a = 0; a < 1001; a++) {
			if (first_occurrence_age_min == first_occurrence_age_max) {
				approx_probabilities[a] = (cladeAgeC/(approx_ages[a]-first_occurrence_age_min+cladeAgeS)) * Math.exp( -Math.pow((Math.log(approx_ages[a]-first_occurrence_age_min+cladeAgeS)-cladeAgeM),2)/cladeAgeW );
			} else {
				// Calculate the integral of the CladeAge distribution pdf for an upper bound, which is just the same as this age, a.
				// The CladeAge distribution pdf is c/(x+s) * exp(-(log(x+s)-m)**2/w), and the integral of this pdf is
				// -1/2 * c * sqrt(pi) * sqrt(w) * erf(-(log(x+s)-m)/w)
				// whereby erf is the error function (see http://en.wikipedia.org/wiki/Error_function)
				double x_upper = approx_ages[a]-first_occurrence_age_min;
				double erf_argument = -(Math.log(x_upper+cladeAgeS)-cladeAgeM) / Math.sqrt(cladeAgeW);
				// Calculate the error function using the Wikipedia approximation (which is from Numerical Recipes in Fortran 77).
				double erf_t = 1/(1.0+0.5*Math.abs(erf_argument));
				double erf_polynomial = -Math.pow(erf_argument,2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t,2) + 0.09678418*Math.pow(erf_t,3) - 0.18628806*Math.pow(erf_t,4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
				double erf_tau = erf_t * Math.exp(erf_polynomial);
				double erf = 0;
				if (erf_argument >= 0) {
					erf = 1 - erf_tau;
				} else {
					erf = erf_tau - 1;
				}
				double int_upper = -0.5 * cladeAgeC * Math.sqrt(Math.PI) * Math.sqrt(cladeAgeW) * erf;
				
				// Calculate the integral of the probabilities for a lower bound which is either first_occurrence_age_min, or this age, a,
				// minus the difference between @first_occurrence_age_max and first_occurrence_age_min, depending on which of the two is greater.
				double x_lower = 0;
				if (approx_ages[a]-first_occurrence_age_max > 0) {
					x_lower = approx_ages[a]-first_occurrence_age_max;
				}
				erf_argument = -(Math.log(x_lower+cladeAgeS)-cladeAgeM) / Math.sqrt(cladeAgeW);
				// Calculate the error function using the Wikipedia approximation (which is from Numerical Recipes in Fortran 77).
				erf_t = 1/(1.0+0.5*Math.abs(erf_argument));
				erf_polynomial = -Math.pow(erf_argument,2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t,2) + 0.09678418*Math.pow(erf_t,3) - 0.18628806*Math.pow(erf_t,4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
				erf_tau = erf_t * Math.exp(erf_polynomial);
				erf = 0;
				if (erf_argument >= 0) {
					erf = 1 - erf_tau;
				} else {
					erf = erf_tau - 1;
				}
				double int_lower = -0.5 * cladeAgeC * Math.sqrt(Math.PI) * Math.sqrt(cladeAgeW) * erf;
				
				// The approximated probability for this age is the integral of the upper bound minus the integral of the lower bound,
				// divided by the difference between first_occurrence_age_max and first_occurrence_age_min.
				approx_probabilities[a] = (int_upper-int_lower)/(first_occurrence_age_max-first_occurrence_age_min);
			}
		}
		
		// Variable cladeAgeC, as determined above is adjusted to provide good fit with the empirical probability values, however, it may
		// not be the right scale factor for a distribution, which must sum to 1. Therefore, the integral under the distribution is calculated
		// between fossil_age_min and Inf, and a scaling factor cladeAgeCcorr is determined so that the distribution sums to 1.
		double erf_argument = -(Math.log(cladeAgeS)-cladeAgeM) / Math.sqrt(cladeAgeW);
		double erf_t = 1/(1.0+0.5*Math.abs(erf_argument));
		double erf_polynomial = -Math.pow(erf_argument,2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t,2) + 0.09678418*Math.pow(erf_t,3) - 0.18628806*Math.pow(erf_t,4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
		double erf_tau = erf_t * Math.exp(erf_polynomial);
		double erf = 0;
		if (erf_argument >= 0) {
			erf = 1 - erf_tau;
		} else {
			erf = erf_tau - 1;
		}
		double cladeAgeCcorr = 1/(0.5 * Math.sqrt(Math.PI) * Math.sqrt(cladeAgeW) * (1+erf));

		// Fill variables approx_distribution_type, approx_distribution_parameters, and approx_distribution_rmsd
		String approx_distribution_type = "CladeAge";
		double[] approx_distribution_parameters = {first_occurrence_age_min,first_occurrence_age_max,cladeAgeCcorr,cladeAgeS,cladeAgeM,cladeAgeW};
		approx_distribution_rmsd = cladeAgeRmsd;
		
		// XXX this is just here for tests.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Offset: " + approx_distribution_parameters[0]);
		System.out.println("Cutoff: " + approx_distribution_parameters[1]);
		System.out.println("C: " + approx_distribution_parameters[2]);
		System.out.println("S: " + approx_distribution_parameters[3]);
		System.out.println("M: " + approx_distribution_parameters[4]);
		System.out.println("W: " + approx_distribution_parameters[5]);
		System.out.println("RMSD: " + approx_distribution_rmsd);

		// XXX todo: The below must be changed to return a CladeAge distribution
		return new LogNormalImpl(approx_distribution_parameters[0], approx_distribution_parameters[1]);
	}
	
	public void fitAnalyticalCladeAge(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max) {

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
				
				// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
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
		
		// Set the number of replicates for both ndr and epsilon variation.
		// XXX This should be tweaked depending on the precision and speed we need.
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
			double bestY = 0;
			double worstY = 0;
			double bests = 0;
			double worsts = 0;
			if (vertex1Y < vertex0Y) {
				bestY = vertex1Y;
				worstY = vertex0Y;
				bests = vertex1s;
				worsts = vertex0s;
			} else {
				bestY = vertex0Y;
				worstY = vertex1Y;
				bests = vertex0s;
				worsts = vertex1s;
			}
			
			// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
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
				
				// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
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

				// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
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

		double cladeAgeS = vertex0s;
		double m_denom1denom = Math.pow((Math.log(cladeAgeS)),2) - Math.pow((Math.log(2*x1+cladeAgeS)),2);
		double m_denom1nom = Math.log((y2/y0)*((2*x1+cladeAgeS)/cladeAgeS));
		double m_denom2denom = Math.pow((Math.log(cladeAgeS)),2) - Math.pow((Math.log(x1+cladeAgeS)),2);
		double m_denom2nom = Math.log((y1/y0)*((x1+cladeAgeS)/cladeAgeS));
		double m_nom1denom = 2*Math.log((x1+cladeAgeS)/cladeAgeS);
		double m_nom1nom = m_denom2nom;
		double m_nom2denom = 2*Math.log((2*x1+cladeAgeS)/cladeAgeS);
		double m_nom2_nom = m_denom1nom;
		double cladeAgeM = ((m_denom1denom/m_denom1nom)-(m_denom2denom/m_denom2nom))/((m_nom1denom/m_nom1nom)-(m_nom2denom/m_nom2_nom));
		double cladeAgeW = (Math.pow((Math.log(cladeAgeS)-cladeAgeM),2) - Math.pow((Math.log(x1+cladeAgeS)-cladeAgeM),2)) / Math.log((y1/y0)*((x1+cladeAgeS)/cladeAgeS));
		double cladeAgeC = y0 * cladeAgeS * Math.exp(Math.pow((Math.log(cladeAgeS)-cladeAgeM),2)/cladeAgeW);
		// Calculate the root mean square deviation (in order to be fair with other distributions, do that based on the 101 original ages.
		//double sd = 0;
		//for (int a = 0; a < ages.length; a++) {
		//	sd += Math.pow((probabilities[a]-((cladeAgeC/(ages[a]-first_occurrence_age_min+cladeAgeS)) * Math.exp( -Math.pow((Math.log(ages[a]-first_occurrence_age_min+cladeAgeS)-cladeAgeM),2)/cladeAgeW ))),2);
		//}
		//double msd = sd/(double) ages.length;
		//double cladeAgeRmsd = Math.sqrt(msd);
		
		// Variable @cladeAgeC, as determined above has been adjusted to provide good fit with the empirical probability values, 
		// however, it may not be the right scale factor for a distribution, which must sum to 1. Therefore, the integral under
		// the distribution is calculated between fossil_age_min and Inf, and a scaling factor @cladeAgeCcorr is determined so that
		// the distribution sums to 1.
		double erf_argument = -(Math.log(cladeAgeS)-cladeAgeM) / Math.sqrt(cladeAgeW);
		double erf_t = 1/(1.0+0.5*Math.abs(erf_argument));
		double erf_polynomial = -Math.pow(erf_argument,2) - 1.26551223 + 1.00002368*erf_t + 0.37409196*Math.pow(erf_t,2) + 0.09678418*Math.pow(erf_t,3) - 0.18628806*Math.pow(erf_t,4) + 0.27886807*Math.pow(erf_t,5) - 1.13520398*Math.pow(erf_t,6) + 1.48851587*Math.pow(erf_t,7) - 0.82215223*Math.pow(erf_t,8) + 0.17087277*Math.pow(erf_t,9);
		double erf_tau = erf_t * Math.exp(erf_polynomial);
		double erf = 0;
		if (erf_argument >= 0) {
			erf = 1 - erf_tau;
		} else {
			erf = erf_tau - 1;
		}
		double cladeAgeCcorr = 1/(0.5 * Math.sqrt(Math.PI) * Math.sqrt(cladeAgeW) * (1+erf));
		
		// Fill variables approx_distribution_type, approx_distribution_parameters, and approx_distribution_rmsd
		String approx_distribution_type = "CladeAge";
		double[] approx_distribution_parameters = {first_occurrence_age_min,first_occurrence_age_max,cladeAgeCcorr,cladeAgeS,cladeAgeM,cladeAgeW};
		// approx_distribution_rmsd = cladeAgeRmsd;
		
		// XXX this is just here for tests.
//		System.out.println("Distribution type: " + approx_distribution_type);
//		System.out.println("Offset: " + approx_distribution_parameters[0]);
//		System.out.println("Cutoff: " + approx_distribution_parameters[1]);
//		System.out.println("C: " + approx_distribution_parameters[2]);
//		System.out.println("S: " + approx_distribution_parameters[3]);
//		System.out.println("M: " + approx_distribution_parameters[4]);
//		System.out.println("W: " + approx_distribution_parameters[5]);
		//System.out.println("RMSD: " + approx_distribution_rmsd);

		// XXX todo: return a distribution.
		
	} // public void fitAnalyticalCladeAge(...)
	
	public ContinuousDistribution fitExponential(JProgressBar progress) {
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] cs = new double[nmRepetitions];
		double[] means = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			if (progress != null) {
				progress.setValue(progress.getValue() + DELTA_PROGRESS);
			}
			if (cancel1 != false) {
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
                boolean keepGoing = true;
                int stepCounter = 0;

                // Until converged do the loop.
				while (keepGoing == true) {
					// Increment the stepCounter.
					stepCounter += 1;
					
					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1/vertex0mean) * Math.exp(-(1/vertex0mean)*(ages[i]-offset));
						vertex0Y += Math.pow(probabilities[i]-vertex0c___*temp, 2);
					}

					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1/vertex1mean) * Math.exp(-(1/vertex1mean)*(ages[i]-offset));
						vertex1Y += Math.pow(probabilities[i]-vertex1c___*temp, 2);
					}

					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = (1/vertex2mean) * Math.exp(-(1/vertex2mean)*(ages[i]-offset));
						vertex2Y += Math.pow(probabilities[i]-vertex2c___*temp, 2);
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
						double temp = (1/reflectionmean) * Math.exp(-(1/reflectionmean)*(ages[i]-offset));
						reflectionY += Math.pow(probabilities[i]-reflectionc___*temp, 2);
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
							double temp = (1/extensionmean) * Math.exp(-(1/extensionmean)*(ages[i]-offset));
							reflectionY += Math.pow(probabilities[i]-extensionc___*temp, 2);
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
							double temp = (1/contractionmean) * Math.exp(-(1/contractionmean)*(ages[i]-offset));
							contractionY += Math.pow(probabilities[i]-contractionc___*temp, 2);
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
							
						// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
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
					// if (cancel2 == true) {
					// 	keepGoing = false;
					// }
					
				} // while (keepGoing == true) 
				
				// Report all parameters and y.
				cs[x] = vertex0c___;
				means[x] = vertex0mean;
				ys[x] = vertex0Y;
							
			//} if (cancel2 == false) {

		} // for (int x = 0; x < nmRepetitions; x++)
		
		// if (cancel2 == false) {
			
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
			
			// Fill arrays approx_ages and approx_probabilities.
			approx_ages = new double[1001];
			approx_probabilities = new double[1001];
			for (int x = 0; x < 1001; x++) {
				approx_ages[x] = offset + (x/1000.0)*(ages[0]-offset);
				approx_probabilities[x] = expConstant * (1/expMean) * Math.exp(-(1/expMean)*(approx_ages[x]-offset));
			}
			
			// Fill variables approx_distribution_type, approx_distribution_parameters, and approx_distribution_rmsd
			String approx_distribution_type = "Exponential";
			double[] approx_distribution_parameters = {expMean};
			approx_distribution_rmsd = expRmsd;
		
		// }
		
		// XXX this is just here for tests.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Mean: " + approx_distribution_parameters[0]);
		System.out.println("RMSD: " + approx_distribution_rmsd);

		
		normaliser = expConstant;

		return new ExponentialDistributionImpl(approx_distribution_parameters[0]);		
	
	}
	
	public ContinuousDistribution fitLognormal(JProgressBar progress) {
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] cs = new double[nmRepetitions];
		double[] sigmas = new double[nmRepetitions];
		double[] mus = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			
			if (progress != null) {
				progress.setValue(progress.getValue() + DELTA_PROGRESS);
			}
			if (cancel1 != false) {
				return null;
			}

				// Initiate the simplex, find 4 verteces.
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
				boolean keepGoing = true;
				int stepCounter = 0;
				
				// Until converged do the loop.
				while (keepGoing == true) {
					
					// Increment the stepCounter.
					stepCounter += 1;
					
					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-offset > 0) {
							temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex0sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex0mu___),2)/(2*Math.pow(vertex0sigma,2)) ));
						}
						vertex0Y += Math.pow((probabilities[i]-vertex0c____*temp),2);
					}

					// Calculate the y value of each vertex.
					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-offset > 0) {
							temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex1sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex1mu___),2)/(2*Math.pow(vertex1sigma,2)) ));
						}
						vertex1Y += Math.pow((probabilities[i]-vertex1c____*temp),2);
					}

					// Calculate the y value of each vertex.
					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-offset > 0) {
							temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex2sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex2mu___),2)/(2*Math.pow(vertex2sigma,2)) ));
						}
						vertex2Y += Math.pow((probabilities[i]-vertex2c____*temp),2);
					}

					// Calculate the y value of each vertex.
					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp = 0;
						if (ages[i]-offset > 0) {
							temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex3sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex3mu___),2)/(2*Math.pow(vertex3sigma,2)) ));
						}
						vertex3Y += Math.pow((probabilities[i]-vertex3c____*temp),2);
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
						if (ages[i]-offset > 0) {
							temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(reflectionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-reflectionmu___),2)/(2*Math.pow(reflectionsigma,2)) ));
						}
						reflectionY += Math.pow((probabilities[i]-reflectionc____*temp),2);
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
							if (ages[i]-offset > 0) {
								temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(extensionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-extensionmu___),2)/(2*Math.pow(extensionsigma,2)) ));
							}
							extensionY += Math.pow((probabilities[i]-extensionc____*temp),2);
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
							if (ages[i]-offset > 0) {
								temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(contractionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-contractionmu___),2)/(2*Math.pow(contractionsigma,2)) ));
							}
							contractionY += Math.pow((probabilities[i]-contractionc____*temp),2);
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

							
						// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
						} else {
							
							// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
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
					// if (cancel2 == true) {
					// 	keepGoing = false;
					// }

					
				} // while (keepGoing == true) 
				
				// Report all parameters and y.
				cs[x] = vertex0c____;
                mus[x] = vertex0mu___;
                sigmas[x] = vertex0sigma;
                ys[x] = vertex0Y;
			
			
			//} if (cancel2 == false)
			
		} // for (int x = 0; x < nmRepetitions; x++)
		
		// if (cancel2 == false) {
		
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

			// Fill arrays approx_ages and approx_probabilities.
			approx_ages = new double[1001];
			approx_probabilities = new double[1001];
			for (int x = 0; x < 1001; x++) {
				approx_ages[x] = offset + (x/1000.0)*(ages[0]-offset);
				if (approx_ages[x]-offset <= 0) {
					approx_probabilities[x] = 0.0;
				} else {
					approx_probabilities[x] = logConstant * (1.0/((approx_ages[x]-offset)*Math.sqrt(2*Math.PI*Math.pow(logStdev,2))))*(Math.exp( -Math.pow((Math.log(approx_ages[x]-offset)-logMean),2)/(2*Math.pow(logStdev,2)) ));
				}
			}
			
			// Fill variables approx_distribution_type, approx_distribution_parameters, and approx_distribution_rmsd
			String approx_distribution_type = "Lognormal";
			double[] approx_distribution_parameters = {logMean,logStdev};
			approx_distribution_rmsd = logRmsd;
			
		// }

		// XXX this is just here for tests.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Mean (log): " + approx_distribution_parameters[0]);
		System.out.println("Stdev (log): " + approx_distribution_parameters[1]);
		System.out.println("RMSD: " + approx_distribution_rmsd);

		normaliser = logConstant;

		return new LogNormalImpl(approx_distribution_parameters[0], approx_distribution_parameters[1]);
	}

	public ContinuousDistribution fitGamma(JProgressBar progress) {
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] cs = new double[nmRepetitions];
		double[] ks = new double[nmRepetitions];
		double[] thetas = new double[nmRepetitions];
		double[] gammaKs = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			
			if (progress != null) {
				progress.setValue(progress.getValue() + DELTA_PROGRESS);
			}
			if (cancel1 != false) {
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
				boolean keepGoing = true;
				int stepCounter = 0;

				// Until converged do the loop.
				while (keepGoing == true) {
					
					// Increment the stepCounter.
					stepCounter += 1;

					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex0theta,vertex0k____)))*(1/vertex0gammaK)*(Math.pow((ages[i]-offset),(vertex0k____-1)))*(Math.exp(-(ages[i]-offset)/vertex0theta));
						vertex0Y += Math.pow((probabilities[i]-vertex0c____*temp),2);
					}
					
					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex1theta,vertex1k____)))*(1/vertex1gammaK)*(Math.pow((ages[i]-offset),(vertex1k____-1)))*(Math.exp(-(ages[i]-offset)/vertex1theta));
						vertex1Y += Math.pow((probabilities[i]-vertex1c____*temp),2);
					}
					
					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex2theta,vertex2k____)))*(1/vertex2gammaK)*(Math.pow((ages[i]-offset),(vertex2k____-1)))*(Math.exp(-(ages[i]-offset)/vertex2theta));
						vertex2Y += Math.pow((probabilities[i]-vertex2c____*temp),2);
					}

					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double temp =  (1.0/(Math.pow(vertex3theta,vertex3k____)))*(1/vertex3gammaK)*(Math.pow((ages[i]-offset),(vertex3k____-1)))*(Math.exp(-(ages[i]-offset)/vertex3theta));
						vertex3Y += Math.pow((probabilities[i]-vertex3c____*temp),2);
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
						double temp =  (1.0/(Math.pow(reflectiontheta,reflectionk____)))*(1/reflectiongammaK)*(Math.pow((ages[i]-offset),(reflectionk____-1)))*(Math.exp(-(ages[i]-offset)/reflectiontheta));
						reflectionY += Math.pow((probabilities[i]-reflectionc____*temp),2);
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
							double temp =  (1.0/(Math.pow(extensiontheta,extensionk____)))*(1/extensiongammaK)*(Math.pow((ages[i]-offset),(extensionk____-1)))*(Math.exp(-(ages[i]-offset)/extensiontheta));
							extensionY += Math.pow((probabilities[i]-extensionc____*temp),2);
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
							double temp =  (1.0/(Math.pow(contractiontheta,contractionk____)))*(1/contractiongammaK)*(Math.pow((ages[i]-offset),(contractionk____-1)))*(Math.exp(-(ages[i]-offset)/contractiontheta));
							contractionY += Math.pow((probabilities[i]-contractionc____*temp),2);
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
							
						// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
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
					if (stepCounter > 1000) {
						keepGoing = false;
					}
					
				} // while (keepGoing == true) 

				// Report all parameters and y.
                cs[x] = vertex0c____;
                ks[x] = vertex0k____;
                thetas[x] = vertex0theta;
                gammaKs[x] = vertex0gammaK;
                ys[x] = vertex0Y;
				
			//} if (cancel2 == false)
                
		} // for (int x = 0; x < nmRepetitions; x++)
		
		// if (cancel2 == false) {
		
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
			double gamGammaK = gammaKs[index];
			double gamRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));
		
			// Fill arrays approx_ages and approx_probabilities.
			approx_ages = new double[1001];
			approx_probabilities = new double[1001];
			for (int x = 0; x < 1001; x++) {
				approx_ages[x] = offset + (x/1000.0)*(ages[0]-offset);
				approx_probabilities[x] = gamConstant * (1.0/(Math.pow(gamScale,gamShape)))*(1/gamGammaK)*(Math.pow((approx_ages[x]-offset),(gamShape-1)))*(Math.exp(-(approx_ages[x]-offset)/gamScale));
			}

			// Fill variables approx_distribution_type, approx_distribution_parameters, and approx_distribution_rmsd
			String approx_distribution_type = "Gamma";
			double[] approx_distribution_parameters = {gamShape,gamScale};
			approx_distribution_rmsd = gamRmsd;
			
		// }

		// XXX this is just here for tests.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Shape: " + approx_distribution_parameters[0]);
		System.out.println("Scale: " + approx_distribution_parameters[1]);
		System.out.println("RMSD: " + approx_distribution_rmsd);

		normaliser = gamConstant;

		return new GammaDistributionImpl(approx_distribution_parameters[0], approx_distribution_parameters[1]);			
	}
	
	public ContinuousDistribution fitExpGamma(JProgressBar progress) {
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] means = new double[nmRepetitions];
		double[] c2s = new double[nmRepetitions];
		double[] ks = new double[nmRepetitions];
		double[] thetas = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			
			if (progress != null) {
				progress.setValue(progress.getValue() + DELTA_PROGRESS);
			}
			if (cancel1 != false) {
				return null;
			}
			
				// Initiate the simplex, find 5 vertices.
				// vertex0
				double vertex0mean_ = 10 + Math.random()*50;
				double vertex0c2___ = 0.5*(0.5 + Math.random());
				double vertex0k____ = 1.5 + Math.random();
				double vertex0theta = (10 + Math.random()*50)/2.0;
				double vertex0Y = 0;
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

		
				// vertex1
				double vertex1mean_ = 10 + Math.random()*50;
				double vertex1c2___ = 0.5*(0.5 + Math.random());
				double vertex1k____ = 1.5 + Math.random();
				double vertex1theta = (10 + Math.random()*50)/2.0;
				double vertex1Y = 0;
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
		
				// vertex2
				double vertex2mean_ = 10 + Math.random()*50;
				double vertex2c2___ = 0.5*(0.5 + Math.random());
				double vertex2k____ = 1.5 + Math.random();
				double vertex2theta = (10 + Math.random()*50)/2.0;
				double vertex2Y = 0;
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
		
				// vertex3
				double vertex3mean_ = 10 + Math.random()*50;
				double vertex3c2___ = 0.5*(0.5 + Math.random());
				double vertex3k____ = 1.5 + Math.random();
				double vertex3theta = (10 + Math.random()*50)/2.0;
				double vertex3Y = 0;
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
						
				// vertex4
				double vertex4mean_ = 10 + Math.random()*50;
				double vertex4c2___ = 0.5*(0.5 + Math.random());
				double vertex4k____ = 1.5 + Math.random();
				double vertex4theta = (10 + Math.random()*50)/2.0;
				double vertex4Y = 0;
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
						
				// Prepare for the Nelder-Mead loop.
				boolean keepGoing = true;
				int stepCounter = 0;
		
				// Until converged do the loop.
				while (keepGoing == true) {
					
					// Increment the stepCounter.
					stepCounter += 1;

					// Calculate the y value of each vertex.
					// vertex0
					vertex0Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex0mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex0theta,vertex0k____)))*(1/vertex0gammaK)*(Math.pow((ages[i]-offset),(vertex0k____-1)))*(Math.exp(-(ages[i]-offset)/vertex0theta));
						vertex0Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex0c2___*gamma_part)),2);
					}
					
					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex1mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex1theta,vertex1k____)))*(1/vertex1gammaK)*(Math.pow((ages[i]-offset),(vertex1k____-1)))*(Math.exp(-(ages[i]-offset)/vertex1theta));
						vertex1Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex1c2___*gamma_part)),2);
					}
					
					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex2mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex2theta,vertex2k____)))*(1/vertex2gammaK)*(Math.pow((ages[i]-offset),(vertex2k____-1)))*(Math.exp(-(ages[i]-offset)/vertex2theta));
						vertex2Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex2c2___*gamma_part)),2);
					}
					
					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex3mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex3theta,vertex3k____)))*(1/vertex3gammaK)*(Math.pow((ages[i]-offset),(vertex3k____-1)))*(Math.exp(-(ages[i]-offset)/vertex3theta));
						vertex3Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex3c2___*gamma_part)),2);
					}
															
					// vertex4
					vertex4Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex4mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex4theta,vertex4k____)))*(1/vertex4gammaK)*(Math.pow((ages[i]-offset),(vertex4k____-1)))*(Math.exp(-(ages[i]-offset)/vertex4theta));
						vertex4Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex4c2___*gamma_part)),2);
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
					double bestmean_ = 0;
					double bestk____ = 0;
					double bestc2___ = 0;
					double besttheta = 0;
					if (vertex0Y == bestY) {
		                bestmean_ = vertex0mean_;
		                bestc2___ = vertex0c2___;
		                bestk____ = vertex0k____;
		                besttheta = vertex0theta;
					} else if (vertex1Y == bestY) {
		                bestmean_ = vertex1mean_;
		                bestc2___ = vertex1c2___;
		                bestk____ = vertex1k____;
		                besttheta = vertex1theta;
					} else if (vertex2Y == bestY) {
		                bestmean_ = vertex2mean_;
		                bestc2___ = vertex2c2___;
		                bestk____ = vertex2k____;
		                besttheta = vertex2theta;
					} else if (vertex3Y == bestY) {
		                bestmean_ = vertex3mean_;
		                bestc2___ = vertex3c2___;
		                bestk____ = vertex3k____;
		                besttheta = vertex3theta;
					} else if (vertex4Y == bestY) {
		                bestmean_ = vertex4mean_;
		                bestc2___ = vertex4c2___;
		                bestk____ = vertex4k____;
		                besttheta = vertex4theta;
					}
		
					// Find the parameter values of the worst vertex.
		            double worstmean_ = 0;
		            double worstc2___ = 0;
		            double worstk____ = 0;
		            double worsttheta = 0;
					if (vertex0Y == worstY) {
		                worstmean_ = vertex0mean_;
		                worstc2___ = vertex0c2___;
		                worstk____ = vertex0k____;
		                worsttheta = vertex0theta;
					} else if (vertex1Y == worstY) {
		                worstmean_ = vertex1mean_;
		                worstc2___ = vertex1c2___;
		                worstk____ = vertex1k____;
		                worsttheta = vertex1theta;
					} else if (vertex2Y == worstY) {
		                worstmean_ = vertex2mean_;
		                worstc2___ = vertex2c2___;
		                worstk____ = vertex2k____;
		                worsttheta = vertex2theta;
					} else if (vertex3Y == worstY) {
		                worstmean_ = vertex3mean_;
		                worstc2___ = vertex3c2___;
		                worstk____ = vertex3k____;
		                worsttheta = vertex3theta;
					} else if (vertex4Y == worstY) {
		                worstmean_ = vertex4mean_;
		                worstc2___ = vertex4c2___;
		                worstk____ = vertex4k____;
		                worsttheta = vertex4theta;
					}
		
					// Calculate the sum of the parameters over all vertices.
					double summean_ = vertex0mean_ + vertex1mean_ + vertex2mean_ + vertex3mean_ + vertex4mean_;
					double sumc2___ = vertex0c2___ + vertex1c2___ + vertex2c2___ + vertex3c2___ + vertex4c2___;
					double sumk____ = vertex0k____ + vertex1k____ + vertex2k____ + vertex3k____ + vertex4k____;
					double sumtheta = vertex0theta + vertex1theta + vertex2theta + vertex3theta + vertex4theta;
					
					// Calculate the parameter values of the centroid.
					double centroidmean_ = (summean_ - worstmean_)/4.0;
					double centroidc2___ = (sumc2___ - worstc2___)/4.0;
					double centroidk____ = (sumk____ - worstk____)/4.0;
					double centroidtheta = (sumtheta - worsttheta)/4.0;
					
					// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
					double reflectionmean_ = centroidmean_ + alph * (centroidmean_ - worstmean_);
					double reflectionc2___ = centroidc2___ + alph * (centroidc2___ - worstc2___);
					if (reflectionc2___ <= 0) {
						reflectionc2___ = 0.5*(0.5 + Math.random());
					}
					double reflectionk____ = centroidk____ + alph * (centroidk____ - worstk____);
					if (reflectionk____ < 0) {
						reflectionk____ = 1.5 + Math.random();
					}
					double reflectiontheta = centroidtheta + alph * (centroidtheta - worsttheta);
					if (reflectiontheta < 0) {
						reflectiontheta = (10 + Math.random()*50)/2.0;
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
						double exp_part = Math.exp(-(1/reflectionmean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(reflectiontheta,reflectionk____)))*(1/reflectiongammaK)*(Math.pow((ages[i]-offset),(reflectionk____-1)))*(Math.exp(-(ages[i]-offset)/reflectiontheta));
						reflectionY += Math.pow((probabilities[i]-(mean_psi*exp_part+reflectionc2___*gamma_part)),2);
					}
					
					// Consider the three cases:
					// i.)   reflection is better than all vertices.
					// ii.)  reflection is better than the second-worst vertex.
					// iii.) reflection is worse than or equally good as the second-worst vertex.
		
					// Case i): If reflection is better than all vertices, replace the worst vertex with an extension of the reflection (by factor gamm).
					if (reflectionY < bestY) {
						
						// Calculate the extension of the worst vertex at the centroid (with extension coefficient gamma).
						double extensionmean_ = centroidmean_ + gamm * (centroidmean_ - worstmean_);
						double extensionc2___ = centroidc2___ + gamm * (centroidc2___ - worstc2___);
						if (extensionc2___ <= 0) {
							extensionc2___ = 0.5*(0.5 + Math.random());
						}
						double extensionk____ = centroidk____ + gamm * (centroidk____ - worstk____);
						if (extensionk____ < 0) {
							extensionk____ = 1.5 + Math.random();
						}
						double extensiontheta = centroidtheta + gamm * (centroidtheta - worsttheta);
						if (extensiontheta < 0) {
							extensiontheta = (10 + Math.random()*50)/2.0;
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
							double exp_part = Math.exp(-(1/extensionmean_)*(ages[i]-offset));
							double gamma_part =  (1.0/(Math.pow(extensiontheta,extensionk____)))*(1/extensiongammaK)*(Math.pow((ages[i]-offset),(extensionk____-1)))*(Math.exp(-(ages[i]-offset)/extensiontheta));
							extensionY += Math.pow((probabilities[i]-(mean_psi*exp_part+extensionc2___*gamma_part)),2);
						}
						
						// Figure out which values to use as replacement for the values of the worst vertex.
		                double replacemean_ = 0;
		                double replacec2___ = 0;
		                double replacek____ = 0;
		                double replacetheta = 0;
						if (reflectionY < extensionY) {
		                    replacemean_ = reflectionmean_;
		                    replacec2___ = reflectionc2___;
		                    replacek____ = reflectionk____;
		                    replacetheta = reflectiontheta;
						} else {
		                    replacemean_ = extensionmean_;
		                    replacec2___ = extensionc2___;
		                    replacek____ = extensionk____;
		                    replacetheta = extensiontheta;
						}
						
						// Replace the parameter values of the worst vertex with the replacement values.
						if (vertex0Y == worstY) {
		                    vertex0mean_ = replacemean_;
		                    vertex0c2___ = replacec2___;
		                    vertex0k____ = replacek____;
		                    vertex0theta = replacetheta;
						} else if (vertex1Y == worstY) {
		                    vertex1mean_ = replacemean_;
		                    vertex1c2___ = replacec2___;
		                    vertex1k____ = replacek____;
		                    vertex1theta = replacetheta;
						} else if (vertex2Y == worstY) {
		                    vertex2mean_ = replacemean_;
		                    vertex2c2___ = replacec2___;
		                    vertex2k____ = replacek____;
		                    vertex2theta = replacetheta;
						} else if (vertex3Y == worstY) {
		                    vertex3mean_ = replacemean_;
		                    vertex3c2___ = replacec2___;
		                    vertex3k____ = replacek____;
		                    vertex3theta = replacetheta;
						} else if (vertex4Y == worstY) {
		                    vertex4mean_ = replacemean_;
		                    vertex4c2___ = replacec2___;
		                    vertex4k____ = replacek____;
		                    vertex4theta = replacetheta;
						}
		
					// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
					} else if (reflectionY < secondWorstY) {
		
						if (vertex0Y == worstY) {
		                    vertex0mean_ = reflectionmean_;
		                    vertex0c2___ = reflectionc2___;
		                    vertex0k____ = reflectionk____;
		                    vertex0theta = reflectiontheta;
						} else if  (vertex1Y == worstY) {
		                    vertex1mean_ = reflectionmean_;
		                    vertex1c2___ = reflectionc2___;
		                    vertex1k____ = reflectionk____;
		                    vertex1theta = reflectiontheta;
						} else if (vertex2Y == worstY) {
		                    vertex2mean_ = reflectionmean_;
		                    vertex2c2___ = reflectionc2___;
		                    vertex2k____ = reflectionk____;
		                    vertex2theta = reflectiontheta;
						} else if (vertex3Y == worstY) {
		                    vertex3mean_ = reflectionmean_;
		                    vertex3c2___ = reflectionc2___;
		                    vertex3k____ = reflectionk____;
		                    vertex3theta = reflectiontheta;
						} else if (vertex4Y == worstY) {
		                    vertex4mean_ = reflectionmean_;
		                    vertex4c2___ = reflectionc2___;
		                    vertex4k____ = reflectionk____;
		                    vertex4theta = reflectiontheta;
						}
		
					// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
					} else {
		
						// Calculate the contraction.
						double contractionmean_ = centroidmean_ + beta * (centroidmean_ - worstmean_);
						double contractionc2___ = centroidc2___ + beta * (centroidc2___ - worstc2___);
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
							double exp_part = Math.exp(-(1/contractionmean_)*(ages[i]-offset));
							double gamma_part = (1.0/(Math.pow(contractiontheta,contractionk____)))*(1/contractiongammaK)*(Math.pow((ages[i]-offset),(contractionk____-1)))*(Math.exp(-(ages[i]-offset)/contractiontheta));
							contractionY += Math.pow((probabilities[i]-(mean_psi*exp_part+contractionc2___*gamma_part)),2);
						}
		
						// Consider two subcases of case iii:
						// iiia) The contraction is better than the worst vertex.
						// iiib) The contraction is not better than the worst vertex.
						
						// Case iiia): If the contraction is better than the worst vertex, replace the worst vertex with the contraction.
						if (contractionY < worstY) {
							
							// Replace the parameter values of the worst vertex with the contraction values.
							if (vertex0Y == worstY) {
		                        vertex0mean_ = contractionmean_;
		                        vertex0c2___ = contractionc2___;
		                        vertex0k____ = contractionc2___;
		                        vertex0theta = contractiontheta;
							} else if (vertex1Y == worstY) {
		                        vertex1mean_ = contractionmean_;
		                        vertex1c2___ = contractionc2___;
		                        vertex1k____ = contractionc2___;
		                        vertex1theta = contractiontheta;
							} else if (vertex2Y == worstY) {
		                        vertex2mean_ = contractionmean_;
		                        vertex2c2___ = contractionc2___;
		                        vertex2k____ = contractionc2___;
		                        vertex2theta = contractiontheta;
							} else if (vertex3Y == worstY) {
		                        vertex3mean_ = contractionmean_;
		                        vertex3c2___ = contractionc2___;
		                        vertex3k____ = contractionc2___;
		                        vertex3theta = contractiontheta;
							} else if (vertex4Y == worstY) {
		                        vertex4mean_ = contractionmean_;
		                        vertex4c2___ = contractionc2___;
		                        vertex4k____ = contractionc2___;
		                        vertex4theta = contractiontheta;
							}
		
						// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
						} else {
		
							// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
							// vertex0
		                    vertex0mean_ = bestmean_ + delt  * (vertex0mean_ - bestmean_);
		                    vertex0c2___ = bestc2___ + delt  * (vertex0c2___ - bestc2___);
		                    vertex0k____ = bestk____ + delt  * (vertex0k____ - bestk____);
		                    vertex0theta = besttheta + delt  * (vertex0theta - besttheta);
		
							// vertex1
		                    vertex1mean_ = bestmean_ + delt  * (vertex1mean_ - bestmean_);
		                    vertex1c2___ = bestc2___ + delt  * (vertex1c2___ - bestc2___);
		                    vertex1k____ = bestk____ + delt  * (vertex1k____ - bestk____);
		                    vertex1theta = besttheta + delt  * (vertex1theta - besttheta);
		
							// vertex2
		                    vertex2mean_ = bestmean_ + delt  * (vertex2mean_ - bestmean_);
		                    vertex2c2___ = bestc2___ + delt  * (vertex2c2___ - bestc2___);
		                    vertex2k____ = bestk____ + delt  * (vertex2k____ - bestk____);
		                    vertex2theta = besttheta + delt  * (vertex2theta - besttheta);
		
							// vertex3
		                    vertex3mean_ = bestmean_ + delt  * (vertex3mean_ - bestmean_);
		                    vertex3c2___ = bestc2___ + delt  * (vertex3c2___ - bestc2___);
		                    vertex3k____ = bestk____ + delt  * (vertex3k____ - bestk____);
		                    vertex3theta = besttheta + delt  * (vertex3theta - besttheta);
																	
							// vertex4
		                    vertex4mean_ = bestmean_ + delt  * (vertex4mean_ - bestmean_);
		                    vertex4c2___ = bestc2___ + delt  * (vertex4c2___ - bestc2___);
		                    vertex4k____ = bestk____ + delt  * (vertex4k____ - bestk____);
		                    vertex4theta = besttheta + delt  * (vertex4theta - besttheta);
																	
						} // if (contractionY < worstY)	
											
					} // if (reflectionY < bestY)
		
					// Stop the loop when all parameter values are identical in the first 10 decimals.
					keepGoing = false;
					if (Math.abs(vertex0mean_ - vertex1mean_) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c2___ - vertex1c2___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex1k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex1theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mean_ - vertex2mean_) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c2___ - vertex2c2___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex2k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex2theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mean_ - vertex3mean_) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c2___ - vertex3c2___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex3k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex3theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mean_ - vertex4mean_) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c2___ - vertex4c2___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0k____ - vertex4k____) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex4theta) > 0.0000001) {
						keepGoing = true;
					}
					if (stepCounter > 1000) {
						keepGoing = false;
					}
					// if (cancel2 == true) {
					// 	keepGoing = false;
					// }
		
				} // while (keepGoing == true) 
		
				// Report all parameters and y.
				means[x] = vertex0mean_;
				c2s[x] = vertex0c2___;
				ks[x] = vertex0k____;
				thetas[x] = vertex0theta;
				ys[x] = vertex0Y;
				
			//} if (cancel2 == false)
			
		} // for (int x = 0; x < nmRepetitions; x++)
				
		// if (cancel2 == false) {
		
			// Find the best result among the nmReplicates replicates.
			int index = 0;
			for (int x = 1; x < nmRepetitions; x++) {
				if (ys[x] < ys[index]) {
					index = x;
				}
			}
			
			double expGamMean = means[index];
			double expGamConstant1 = expGamMean * mean_psi;
			double expGamConstant2 = c2s[index];
			double expGamShape = ks[index];
			double expGamScale = thetas[index];
			double expGamRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1));

			// Fill arrays approx_ages and approx_probabilities.
			approx_ages = new double[1001];
			approx_probabilities = new double[1001];
			for (int x = 0; x < 1001; x++) {
				approx_ages[x] = offset + (x/1000.0)*(ages[0]-offset);				
				double exp_part = Math.exp(-(1/expGamMean)*(approx_ages[x]-offset));
				double gamma_part = (1.0/(Math.pow(expGamScale,2.0)))*((approx_ages[x]-offset))*(Math.exp(-(approx_ages[x]-offset)/expGamScale));
				approx_probabilities[x] = mean_psi*exp_part + expGamConstant2*gamma_part;
			}

			// Fill variables approx_distribution_type, approx_distribution_parameters, and approx_distribution_rmsd
			String approx_distribution_type = "ExpGamma";
			double[] approx_distribution_parameters = {expGamMean,(expGamConstant1/expGamConstant2),expGamScale};
			approx_distribution_rmsd = expGamRmsd;
			
		// }
		
		// XXX this is just here for tests.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Mean: " + approx_distribution_parameters[0]);
		System.out.println("ExpGamConstant1: " + expGamConstant1);
		System.out.println("ExpGamConstant2: " + expGamConstant2);
		System.out.println("Ratio: " + approx_distribution_parameters[1]);
		System.out.println("Shape: " + expGamShape);
		System.out.println("Scale: " + approx_distribution_parameters[2]);
		System.out.println("RMSD: " + approx_distribution_rmsd);
		System.out.println("Weight: " + expGamConstant2/(expGamConstant1+expGamConstant2));
		System.out.println("Normaliser: " + normaliser);

		normaliser = expGamConstant1 + expGamConstant2;
		
		return new ExpGamma(expGamConstant1/(expGamConstant1+expGamConstant2), 
				approx_distribution_parameters[0], 
				expGamShape, 
				approx_distribution_parameters[2]); 
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
			cladeAgeProbabilities.fitAnalyticalCladeAge(10.0,10.0,0.04,0.04,0.2,0.2,0.02,0.02);
		}

  }

	public double getNormaliser() {
		return normaliser;
	}
}
