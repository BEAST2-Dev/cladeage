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
			dpb.setValue(successful_simulations[0]);

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
			boolean extant_at_this_age = false;
			for (int i = 0; i < ages.length; i++) {
				
				// For each age, determine the tree length.
				double tree_duration = ages[i] - first_occurrence_age;
				
				// Trim the tree so that no branches are older than @ages[i] - first_occurrence_age. While trimming, make sure the tree contains at least one extant species.
				if (tree_duration >= 0) {
					for (int o = branch_origin.size()-1; o >= 0; o--) {
						boolean remove_this_branch = false;
						if (branch_origin.get(o) >= tree_duration) {
							remove_this_branch = true;
							extant_at_this_age = true;
						} else if (branch_termination.get(o) < sampling_gap) {
							remove_this_branch = true;
						} else {
							if (branch_termination.get(o) >= tree_duration) {
								branch_termination.set(o,tree_duration);
								extant_at_this_age = true;
							}
							if (branch_origin.get(o) < sampling_gap) {
								branch_origin.set(o,sampling_gap); 
							}
						}
						if (remove_this_branch == true) {
							branch_origin.remove(o);
							branch_termination.remove(o);
						}
					} // for (int o = branch_origin.size()-1; o >= 0; o--)
					
					// If at least one species is extant at this age, get the sum of lineage durations, excluding the sampling gap.
					if (extant_at_this_age == true) {
						double sum_of_species_duration = 0;
						for (int o = 0; o < branch_origin.size(); o++) {
							if (branch_termination.get(o) > sampling_gap) {
								if (branch_origin.get(o) > sampling_gap) {
									sum_of_species_duration += branch_termination.get(o) - branch_origin.get(o);
								} else {
									sum_of_species_duration += branch_termination.get(o) - sampling_gap;
								}
							}
						}
						// Draw a value for the sampling rate psi.
						for (int pp = 0; pp < psi_sample_size; pp++) {
							double psi = psi_min + Math.random()*(psi_max-psi_min);
							if (tree_duration >= sampling_gap) {
								raw_probabilities[i] += (psi*Math.exp(-psi*sum_of_species_duration))/(double) psi_sample_size;
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

		/// XXX
		System.out.println("From CladeAgeProbabilities:");
		for (int ccc = 0; ccc < probabilities.length; ccc++) {
			System.out.println("Age: " + ages[ccc] + "\tProbability: " + probabilities[ccc]);
		}

		
	} // public void bd_simulate(...)

	public ContinuousDistribution fitExponential() {
		int nmRepetitions = 10;
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] cs = new double[nmRepetitions];
		double[] means = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			// if (cancel2 == false) {
				
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

// XXX
System.out.println("Normaliser: " + normaliser);

		return new ExponentialDistributionImpl(approx_distribution_parameters[0]);		
	
	}
	
	public ContinuousDistribution fitLognormal() {
		int nmRepetitions = 10;
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] cs = new double[nmRepetitions];
		double[] sigmas = new double[nmRepetitions];
		double[] mus = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			
			// if (cancel2 == false) {

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

	public ContinuousDistribution fitGamma() {
		int nmRepetitions = 10;
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] cs = new double[nmRepetitions];
		double[] ks = new double[nmRepetitions];
		double[] thetas = new double[nmRepetitions];
		double[] gammaKs = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			
			// if (cancel2 == false) {
				
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
	
	public ContinuousDistribution fitExpGamma() {
		int nmRepetitions = 10;
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] means = new double[nmRepetitions];
		double[] c2s = new double[nmRepetitions];
		double[] thetas = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			
			// if (cancel2 == false) {
			
				// Initiate the simplex, find 4 vertices.
				// vertex0
				double vertex0mean_ = 10 + Math.random()*50;
				double vertex0c2___ = 0.5*(0.5 + Math.random());
				double vertex0theta = (10 + Math.random()*50)/2.0;
				double vertex0Y = 0;
		
				// vertex1
				double vertex1mean_ = 10 + Math.random()*50;
				double vertex1c2___ = 0.5*(0.5 + Math.random());
				double vertex1theta = (10 + Math.random()*50)/2.0;
				double vertex1Y = 0;
		
				// vertex2
				double vertex2mean_ = 10 + Math.random()*50;
				double vertex2c2___ = 0.5*(0.5 + Math.random());
				double vertex2theta = (10 + Math.random()*50)/2.0;
				double vertex2Y = 0;
		
				// vertex3
				double vertex3mean_ = 10 + Math.random()*50;
				double vertex3c2___ = 0.5*(0.5 + Math.random());
				double vertex3theta = (10 + Math.random()*50)/2.0;
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
						double exp_part = Math.exp(-(1/vertex0mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex0theta,2.0)))*((ages[i]-offset))*(Math.exp(-(ages[i]-offset)/vertex0theta));
						vertex0Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex0c2___*gamma_part)),2);
					}
					
					// Calculate the y value of each vertex.
					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex1mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex1theta,2.0)))*((ages[i]-offset))*(Math.exp(-(ages[i]-offset)/vertex1theta));
						vertex1Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex1c2___*gamma_part)),2);
					}
					
					// Calculate the y value of each vertex.
					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex2mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex2theta,2.0)))*((ages[i]-offset))*(Math.exp(-(ages[i]-offset)/vertex2theta));
						vertex2Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex2c2___*gamma_part)),2);
					}
					
					// Calculate the y value of each vertex.
					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/vertex3mean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(vertex3theta,2.0)))*((ages[i]-offset))*(Math.exp(-(ages[i]-offset)/vertex3theta));
						vertex3Y += Math.pow((probabilities[i]-(mean_psi*exp_part+vertex3c2___*gamma_part)),2);
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
					double bestmean_ = 0;
					double bestc2___ = 0;
					double besttheta = 0;
					if (vertex0Y == bestY) {
		                bestmean_ = vertex0mean_;
		                bestc2___ = vertex0c2___;
		                besttheta = vertex0theta;
					} else if (vertex1Y == bestY) {
		                bestmean_ = vertex1mean_;
		                bestc2___ = vertex1c2___;
		                besttheta = vertex1theta;
					} else if (vertex2Y == bestY) {
		                bestmean_ = vertex2mean_;
		                bestc2___ = vertex2c2___;
		                besttheta = vertex2theta;
					} else if (vertex3Y == bestY) {
		                bestmean_ = vertex3mean_;
		                bestc2___ = vertex3c2___;
		                besttheta = vertex3theta;
					}
		
					// Find the parameter values of the worst vertex.
		            double worstmean_ = 0;
		            double worstc2___ = 0;
		            double worsttheta = 0;
					if (vertex0Y == worstY) {
		                worstmean_ = vertex0mean_;
		                worstc2___ = vertex0c2___;
		                worsttheta = vertex0theta;
					} else if (vertex1Y == worstY) {
		                worstmean_ = vertex1mean_;
		                worstc2___ = vertex1c2___;
		                worsttheta = vertex1theta;
					} else if (vertex2Y == worstY) {
		                worstmean_ = vertex2mean_;
		                worstc2___ = vertex2c2___;
		                worsttheta = vertex2theta;
					} else if (vertex3Y == worstY) {
		                worstmean_ = vertex3mean_;
		                worstc2___ = vertex3c2___;
		                worsttheta = vertex3theta;
					}
		
					// Calculate the sum of the parameters over all vertices.
					double summean_ = vertex0mean_ + vertex1mean_ + vertex2mean_ + vertex3mean_;
					double sumc2___ = vertex0c2___ + vertex1c2___ + vertex2c2___ + vertex3c2___;
					double sumtheta = vertex0theta + vertex1theta + vertex2theta + vertex3theta;
					
					// Calculate the parameter values of the centroid.
					double centroidmean_ = (summean_ - worstmean_)/3.0;
					double centroidc2___ = (sumc2___ - worstc2___)/3.0;
					double centroidtheta = (sumtheta - worsttheta)/3.0;
					
					// Calculate the reflection of the worst vertex at the centroid (with reflection coefficient alpha).
					double reflectionmean_ = centroidmean_ + alph * (centroidmean_ - worstmean_);
					double reflectionc2___ = centroidc2___ + alph * (centroidc2___ - worstc2___);
					if (reflectionc2___ <= 0) {
						reflectionc2___ = 0.5*(0.5 + Math.random());
					}
					double reflectiontheta = centroidtheta + alph * (centroidtheta - worsttheta);
					if (reflectiontheta < 0) {
						reflectiontheta = (10 + Math.random()*50)/2.0;
					}
		
					// Calculate the y value of the reflection.
					double reflectionY = 0;
					for (int i = 0; i < ages.length; i++) {
						double exp_part = Math.exp(-(1/reflectionmean_)*(ages[i]-offset));
						double gamma_part = (1.0/(Math.pow(reflectiontheta,2.0)))*((ages[i]-offset))*(Math.exp(-(ages[i]-offset)/reflectiontheta));
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
						double extensiontheta = centroidtheta + gamm * (centroidtheta - worsttheta);
						if (extensiontheta < 0) {
							extensiontheta = (10 + Math.random()*50)/2.0;
						}
		
						// Calculate the y value of the extension.
						double extensionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double exp_part = Math.exp(-(1/extensionmean_)*(ages[i]-offset));
							double gamma_part = (1.0/(Math.pow(extensiontheta,2.0)))*((ages[i]-offset))*(Math.exp(-(ages[i]-offset)/extensiontheta));
							extensionY += Math.pow((probabilities[i]-(mean_psi*exp_part+extensionc2___*gamma_part)),2);
						}
						
						// Figure out which values to use as replacement for the values of the worst vertex.
		                double replacemean_ = 0;
		                double replacec2___ = 0;
		                double replacetheta = 0;
						if (reflectionY < extensionY) {
		                    replacemean_ = reflectionmean_;
		                    replacec2___ = reflectionc2___;
		                    replacetheta = reflectiontheta;
						} else {
		                    replacemean_ = extensionmean_;
		                    replacec2___ = extensionc2___;
		                    replacetheta = extensiontheta;
						}
						
						// Replace the parameter values of the worst vertex with the replacement values.
						if (vertex0Y == worstY) {
		                    vertex0mean_ = replacemean_;
		                    vertex0c2___ = replacec2___;
		                    vertex0theta = replacetheta;
						} else if (vertex1Y == worstY) {
		                    vertex1mean_ = replacemean_;
		                    vertex1c2___ = replacec2___;
		                    vertex1theta = replacetheta;
						} else if (vertex2Y == worstY) {
		                    vertex2mean_ = replacemean_;
		                    vertex2c2___ = replacec2___;
		                    vertex2theta = replacetheta;
						} else if (vertex3Y == worstY) {
		                    vertex3mean_ = replacemean_;
		                    vertex3c2___ = replacec2___;
		                    vertex3theta = replacetheta;
						}
		
					// Case ii): If the reflection is better than the second worst vertex, replace the worst vertex with the reflection.
					} else if (reflectionY < secondWorstY) {
		
						if (vertex0Y == worstY) {
		                    vertex0mean_ = reflectionmean_;
		                    vertex0c2___ = reflectionc2___;
		                    vertex0theta = reflectiontheta;
						} else if  (vertex1Y == worstY) {
		                    vertex1mean_ = reflectionmean_;
		                    vertex1c2___ = reflectionc2___;
		                    vertex1theta = reflectiontheta;
						} else if (vertex2Y == worstY) {
		                    vertex2mean_ = reflectionmean_;
		                    vertex2c2___ = reflectionc2___;
		                    vertex2theta = reflectiontheta;
						} else if (vertex3Y == worstY) {
		                    vertex3mean_ = reflectionmean_;
		                    vertex3c2___ = reflectionc2___;
		                    vertex3theta = reflectiontheta;
						}
		
					// Case iii): If the reflection is worse than or equally good as the second-worst vertex, calculate the contraction (by factor @@beta (=ro in the Wikipedia example)).
					} else {
		
						// Calculate the contraction.
						double contractionmean_ = centroidmean_ + beta * (centroidmean_ - worstmean_);
						double contractionc2___ = centroidc2___ + beta * (centroidc2___ - worstc2___);
						double contractiontheta = centroidtheta + beta * (centroidtheta - worsttheta);
		
						// Calculate the y value of the contraction.
						double contractionY = 0;
						for (int i = 0; i < ages.length; i++) {
							double exp_part = Math.exp(-(1/contractionmean_)*(ages[i]-offset));
							double gamma_part = (1.0/(Math.pow(contractiontheta,2.0)))*((ages[i]-offset))*(Math.exp(-(ages[i]-offset)/contractiontheta));
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
		                        vertex0theta = contractiontheta;
							} else if (vertex1Y == worstY) {
		                        vertex1mean_ = contractionmean_;
		                        vertex1c2___ = contractionc2___;
		                        vertex1theta = contractiontheta;
							} else if (vertex2Y == worstY) {
		                        vertex2mean_ = contractionmean_;
		                        vertex2c2___ = contractionc2___;
		                        vertex2theta = contractiontheta;
							} else if (vertex3Y == worstY) {
		                        vertex3mean_ = contractionmean_;
		                        vertex3c2___ = contractionc2___;
		                        vertex3theta = contractiontheta;
							}
		
						// Case iiib): If the contraction is not better than the worst vertex, bring all verteces closer to the best vertex.
						} else {
		
							// Replace each vertex by a new vertex that is halfway between the original vertex and the best vertex (given that @@delt  = 0.5, the default setting)
							// vertex0
		                    vertex0mean_ = bestmean_ + delt  * (vertex0mean_ - bestmean_);
		                    vertex0c2___ = bestc2___ + delt  * (vertex0c2___ - bestc2___);
		                    vertex0theta = besttheta + delt  * (vertex0theta - besttheta);
		
							// vertex1
		                    vertex1mean_ = bestmean_ + delt  * (vertex1mean_ - bestmean_);
		                    vertex1c2___ = bestc2___ + delt  * (vertex1c2___ - bestc2___);
		                    vertex1theta = besttheta + delt  * (vertex1theta - besttheta);
		
							// vertex2
		                    vertex2mean_ = bestmean_ + delt  * (vertex2mean_ - bestmean_);
		                    vertex2c2___ = bestc2___ + delt  * (vertex2c2___ - bestc2___);
		                    vertex2theta = besttheta + delt  * (vertex2theta - besttheta);
		
							// vertex3
		                    vertex3mean_ = bestmean_ + delt  * (vertex3mean_ - bestmean_);
		                    vertex3c2___ = bestc2___ + delt  * (vertex3c2___ - bestc2___);
		                    vertex3theta = besttheta + delt  * (vertex3theta - besttheta);
																	
						} // if (contractionY < worstY)	
											
					} // if (reflectionY < bestY)
		
					// Stop the loop when all parameter values are identical in the first 10 decimals.
					keepGoing = false;
					if (Math.abs(vertex0mean_ - vertex1mean_) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c2___ - vertex1c2___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex1theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mean_ - vertex2mean_) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c2___ - vertex2c2___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex2theta) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mean_ - vertex3mean_) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c2___ - vertex3c2___) > 0.0000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0theta - vertex3theta) > 0.0000001) {
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
// XXX			
for (int ccc = 0; ccc < approx_probabilities.length; ccc++) {
	System.out.println("Approx ages: " + approx_ages[ccc] + "\tApprox probabilities: " + approx_probabilities[ccc]);
}

			// Fill variables approx_distribution_type, approx_distribution_parameters, and approx_distribution_rmsd
			String approx_distribution_type = "ExpGamma";
			double[] approx_distribution_parameters = {expGamMean,(expGamConstant1/expGamConstant2),expGamScale};
			approx_distribution_rmsd = expGamRmsd;
			
		// }
		
		// XXX this is just here for tests.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Mean: " + approx_distribution_parameters[0]);
		System.out.println("Ratio: " + approx_distribution_parameters[1]);
		System.out.println("Scale: " + approx_distribution_parameters[2]);
		System.out.println("RMSD: " + approx_distribution_rmsd);
		
		normaliser = expGamConstant1 + expGamConstant2;

		return new ExpGamma(expGamConstant1/(expGamConstant1+expGamConstant2), 
				approx_distribution_parameters[0], 
				approx_distribution_parameters[1], 
				approx_distribution_parameters[2]); 
	}
	
	public static void main(String[] args) {

		System.out.println(System.currentTimeMillis());
		CladeAgeProbabilities cladeAgeProbabilities = new CladeAgeProbabilities();
		cladeAgeProbabilities.bd_simulate(10.0,10.0,0.01,0.01,0.1,0.1,0.01,0.01,0,0,1000,100000,10, new JProgressBar());
		System.out.println(System.currentTimeMillis());
		cladeAgeProbabilities.fitExponential();
		System.out.println(System.currentTimeMillis());

  }

	public double getNormaliser() {
		return normaliser;
	}
}
