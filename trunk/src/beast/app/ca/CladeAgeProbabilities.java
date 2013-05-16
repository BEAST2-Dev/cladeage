package beast.app.ca;

import java.util.ArrayList;

import javax.swing.JProgressBar;

public class CladeAgeProbabilities {

	// Initialize objects that are to become variables of class instances.
	private int number_of_ages = 100;
	private double[] ages = new double[number_of_ages+1];
	private double[] probabilities = new double[number_of_ages+1];
	private double[] approx_ages = new double[number_of_ages+1];
	private double[] approx_probabilities = new double[number_of_ages+1];
	private boolean[] one_of_the_trees_too_large = new boolean[number_of_ages + 1];
	private String approx_distribution_type;
	private double[] approx_distribution_parameters;
	private double approx_distribution_rmsd;
	private int number_of_ages_with_trees_that_are_too_large;
	private double offset = 0;
	private boolean cancel1 = false;

	// Initialize constants for Lanczos approximation of the gamma function.
	private double lanczos_g = 7;
	private double[] lanczos_p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,771.32342877765313, -176.61502916214059, 12.507343278686905,-0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

	// Initialize constants for the Nelder-Mead Downhill Simplex (according to Wikipedia example).
	private double alph = 1.0;
	private double gamm = 2.0;
	private double beta = -0.5;
	private double delt  = 0.5;

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

	public int getNumber_of_ages_with_trees_that_are_too_large() {
		return number_of_ages_with_trees_that_are_too_large;
	}

	public void bd_simulate(double first_occurrence_age_min, double first_occurrence_age_max, double ndr_min, double ndr_max, double epsilon_min, double epsilon_max, double psi_min, double psi_max, double sampling_gap_min, double sampling_gap_max, int bd_sample_size, int max_tree_size, int psi_sample_size, JProgressBar dpb) {
		cancel1 = false;

		// Unhide the progress indicator.
		// simulation_progress_indicator.setHidden(false);

		// Reset arrays.
		ages = new double[number_of_ages + 1];
		probabilities = new double[number_of_ages + 1];
		approx_ages = new double[number_of_ages + 1];
		approx_probabilities = new double[number_of_ages + 1];
		one_of_the_trees_too_large = new boolean[number_of_ages + 1];

		// Memorize the offset.
		offset = first_occurrence_age_min;

		// Roughly calculate a maximum age.
		double max_simulation_age = Math.log(9.210340371976182*((ndr_max+ndr_min)/(psi_max+psi_min)) + 1.0) * (2.0/(ndr_max+ndr_min)) + first_occurrence_age_max + sampling_gap_max;

		// Determine the ages.
		for (int i = 0; i < number_of_ages; i++) {
			ages[i] = max_simulation_age - (i/ (double) number_of_ages)*(max_simulation_age-first_occurrence_age_min);
		}
		ages[number_of_ages] = first_occurrence_age_min;

		// Initiate percentage and old_percentage for the progress indicator.
		// double multi_percentage = 0.0;
		// double old_multi_percentage = 0.0;

		for (int bd = 0; bd < bd_sample_size; bd++) {

			// Feedback.
			System.out.println(bd + 1);

			if (cancel1) {
				return;
			}
			dpb.setValue(bd);

				// Increment the progress indicator and update the text.
				// lines 88..102 of Probabilities.rb

				// For each age, perform birth death simulations.
				for (int i = 0; i < number_of_ages+1; i++) {

					if (one_of_the_trees_too_large[i] == false) {
						double age = ages[i];
						boolean successful_bd_simulation = false;
						double result_sum = 0;
						boolean keepGoing = true;

						while (keepGoing == true) {

							// Draw values for parameters ndr (net diversification rate, lambda - mu) and epsilon (turnover rate, mu/lambda) from uniform distributions, and calculate lambda and mu from them.
							double ndr = ndr_min + Math.random() * (ndr_max-ndr_min);
							double epsilon = epsilon_min + Math.random() * (epsilon_max-epsilon_min);
							double mu = (ndr*epsilon)/(1.0 - epsilon);
							double lambda = ndr + mu;

							// Draw a duration t for this simulation.
							double tree_duration = age - first_occurrence_age_max + Math.random() * (first_occurrence_age_max-first_occurrence_age_min);

							// Draw a sampling gap.
							double sampling_gap = sampling_gap_min + Math.random() * (sampling_gap_max-sampling_gap_min);

							// Run birth-death simulation with those parameter values, only if t > 0, else result = 0.
							if (tree_duration - sampling_gap >= 0) {
								
								// Initiate arrays for branch origin and branch termination.
								ArrayList<Double> branch_origin = new ArrayList<Double>();
								ArrayList<Double> branch_termination = new ArrayList<Double>();
								ArrayList<Double> new_branch_origin = new ArrayList<Double>();
								ArrayList<Double> new_branch_termination = new ArrayList<Double>();
//								double[] branch_origin = new double[max_tree_size + 1];
//								int branch_origin_counter = 0;
//								double[] branch_termination = new double[max_tree_size + 1];
//								int branch_termination_counter = 0;
//								double[] new_branch_origin = new double[max_tree_size + 1];
//								new_branch_origin[0] = tree_duration;
								new_branch_origin.add(tree_duration);
//								int new_branch_origin_counter = 1;
//								double[] new_branch_termination = new double[max_tree_size + 1];
//								int new_branch_termination_counter = 0;

								// Start the tree generation loop.
//								while (new_branch_origin_counter > 0 && branch_origin_counter < max_tree_size) {
								while (new_branch_origin.size() > 0 && branch_origin.size() < max_tree_size) {

									// For each new origin, add a new termination.
//									for (int o = 0; o < new_branch_origin_counter; o++) {
									for (int o = 0; o < new_branch_origin.size(); o++) {
//										new_branch_termination[new_branch_termination_counter] = new_branch_origin[o] - (Math.log(Math.random())/(-(lambda+mu)));
										new_branch_termination.add(new_branch_origin.get(o) - (Math.log(Math.random())/(-(lambda+mu))));
//										new_branch_termination_counter += 1;
									}

									// Add new origin and termination to the old collection.
//									for (int o = 0; o < new_branch_origin_counter; o++) {
									for (int o = 0; o < new_branch_origin.size(); o++) {
//										branch_origin[branch_origin_counter] = new_branch_origin[o];
										branch_origin.add(new_branch_origin.get(o));
//										branch_origin_counter += 1;
									}
//									for (int t = 0; t < new_branch_termination_counter; t++) {
									for (int t = 0; t < new_branch_termination.size(); t++) {
//										branch_termination[branch_termination_counter] = new_branch_termination[t];
										branch_termination.add(new_branch_termination.get(t));
//										branch_termination_counter += 1;
									}

									// Empty the new origin array.
//									new_branch_origin = new double[max_tree_size + 1];
									new_branch_origin = new ArrayList<Double>();
//									new_branch_origin_counter = 0;

									// For each new termination, add it to the new origin array if it is > 0 and rand < lambda/(lambda+mu) - this represents a speciation event.
//									for (int t = 0; t < new_branch_termination_counter; t++) {
									for (int t = 0; t < new_branch_termination.size(); t++) {										
//										if (new_branch_termination[t] > 0) {
										if (new_branch_termination.get(t) > 0) {
											if (Math.random() < lambda/(lambda+mu)) {
//												new_branch_origin[new_branch_origin_counter] = new_branch_termination[t];
												new_branch_origin.add(new_branch_termination.get(t));
//												new_branch_origin_counter += 1;
//												new_branch_origin[new_branch_origin_counter] = new_branch_termination[t];
												new_branch_origin.add(new_branch_termination.get(t));
//												new_branch_origin_counter += 1;
											}
										}
									}

									// Empty the new termination array.
//									new_branch_termination = new double[max_tree_size + 1];
									new_branch_termination = new ArrayList<Double>();
//									new_branch_termination_counter = 0;

									// Set one_of_the_trees_too_large to true if this tree has become too large.
//									if (branch_origin_counter >= max_tree_size) {
									if (branch_origin.size() >= max_tree_size) {
										one_of_the_trees_too_large[i] = true;
									}

								} // while (new_branch_origin_counter > 0 && branch_origin_counter < max_tree_size)

								// Count the number of extant species and set branch terminations that are less than 0 back to 0.
								int extant_taxa = 0;
								if (one_of_the_trees_too_large[i] == false) {
//									for (int b = 0; b < branch_termination_counter; b++) {
									for (int b = 0; b < branch_termination.size(); b++) {
//										if (branch_termination[b] < 0) {
										if (branch_termination.get(b) < 0) {
//											branch_termination[b] = 0.0;
											branch_termination.set(b,0.0);
											extant_taxa += 1;
										}
									}
								}

								// Set branch origins and terminations that are more than tree_duration-sampling_gap back to tree_duration-sampling_gap.
								if (one_of_the_trees_too_large[i] == false) {
//									for (int b = 0; b < branch_origin_counter; b++) {
									for (int b = 0; b < branch_origin.size(); b++) {
//										if (branch_origin[b] > tree_duration-sampling_gap) {
										if (branch_origin.get(b) > tree_duration-sampling_gap) {
//											branch_origin[b] = tree_duration-sampling_gap;
											branch_origin.set(b,tree_duration-sampling_gap);
										}
//										if (branch_termination[b] > tree_duration-sampling_gap) {
										if (branch_termination.get(b) > tree_duration-sampling_gap) {
//											branch_termination[b] = tree_duration-sampling_gap;
											branch_termination.set(b,tree_duration-sampling_gap);
										}
									}
								}

								// If the tree contains at least one extant species, calculate the sum of species durations.
								if (one_of_the_trees_too_large[i] == false) {
									if (extant_taxa > 0) {
										double sum_of_species_duration = 0;
//										for (int b = 0; b < branch_origin_counter; b++) {
										for (int b = 0; b < branch_origin.size(); b++) {
//											sum_of_species_duration += branch_origin[b] - branch_termination[b];
											sum_of_species_duration += branch_origin.get(b) - branch_termination.get(b);
										}

										// Once a suitable tree is found, calculate the result for multiple different values of psi.
										for (int p = 0; p < psi_sample_size; p++) {
											double psi = psi_min + Math.random() * (psi_max-psi_min);
											result_sum += psi * Math.exp(-psi*sum_of_species_duration);
										}

										// Mark the birth-death simulation as successful.
										successful_bd_simulation = true;
									}
								}

							} else if (tree_duration - sampling_gap < 0) { 

								// Mark the birth-death simulation as successful.
								successful_bd_simulation = true;

							} // if tree_duration >= 0 ... else if (tree_duration - sampling_gap < 0) 

							// Check whether the loop should go on.
							if (successful_bd_simulation == true || one_of_the_trees_too_large[i] == true) {
								keepGoing = false;
							}

						} // while (keepGoing == true)

						// Scale the result_sum and add it to the probabilities array.
						if (one_of_the_trees_too_large[i] == false) {
							probabilities[i] += result_sum/ (double) (bd_sample_size*psi_sample_size);
						}

					} // if (one_of_the_trees_too_large[i] == false)

				} // for (int i = 0; i < number_of_ages+1; i++)

			// } // if (cancel1 == false) {

		} // for (int b = 0; b < bd_sample_size; b++)
		
		// if (cancel1 == false) {
			int number_of_ages_with_trees_that_are_too_large = 0;
			for (int i = 0; i < number_of_ages+1; i++) {
				if (one_of_the_trees_too_large[i] == true) {
					number_of_ages_with_trees_that_are_too_large += 1;
					probabilities[i] = -1;
				}
			}
		//}

	} // public void bd_simulate(...)

	public double fitExponential() {
		int nmRepetitions = 10;
		
		// Prepare arrays for parameters that are to be filled with each Nelder-Mead Downhill Simplex run.
		double[] cs = new double[nmRepetitions];
		double[] means = new double[nmRepetitions];
		double[] ys = new double[nmRepetitions];
		
		for (int x = 0; x < nmRepetitions; x++) {
			// if (cancel2 == false) {
				
				// Initiate the simplex, find 3 verteces.
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
						if (one_of_the_trees_too_large[i] == false) {
							double temp = (1/vertex0mean) * Math.exp(-(1/vertex0mean)*(ages[i]-offset));
							vertex0Y += Math.pow(probabilities[i]-vertex0c___*temp, 2);
						}
					}

					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						if (one_of_the_trees_too_large[i] == false) {
							double temp = (1/vertex1mean) * Math.exp(-(1/vertex1mean)*(ages[i]-offset));
							vertex1Y += Math.pow(probabilities[i]-vertex1c___*temp, 2);
						}
					}

					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						if (one_of_the_trees_too_large[i] == false) {
							double temp = (1/vertex2mean) * Math.exp(-(1/vertex2mean)*(ages[i]-offset));
							vertex2Y += Math.pow(probabilities[i]-vertex2c___*temp, 2);
						}
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
						if (one_of_the_trees_too_large[i] == false) {
							double temp = (1/reflectionmean) * Math.exp(-(1/reflectionmean)*(ages[i]-offset));
							reflectionY += Math.pow(probabilities[i]-reflectionc___*temp, 2);
						}
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
							if (one_of_the_trees_too_large[i] == false) {
								double temp = (1/extensionmean) * Math.exp(-(1/extensionmean)*(ages[i]-offset));
								reflectionY += Math.pow(probabilities[i]-extensionc___*temp, 2);
							}
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
							if (one_of_the_trees_too_large[i] == false) {
								double temp = (1/contractionmean) * Math.exp(-(1/contractionmean)*(ages[i]-offset));
								contractionY += Math.pow(probabilities[i]-contractionc___*temp, 2);
							}
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
					if (Math.abs(vertex0c___- vertex1c___) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mean- vertex1mean) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c___- vertex2c___) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mean- vertex2mean) > 0.000000001) {
						keepGoing = true;
					}
					if (stepCounter > 100000) {
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
			double expRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1-number_of_ages_with_trees_that_are_too_large));
			
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
			double approx_distribution_rmsd = expRmsd;
		
		// }
		
		// XXX needs to be changed so that an object of class ParameterDistribution is returned.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Mean: " + approx_distribution_parameters[0]);
		System.out.println("RMSD: " + approx_distribution_rmsd);
		return approx_distribution_parameters[0];
	}
	
	public void fitLognormal() {
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
						if (one_of_the_trees_too_large[i] == false) {
							double temp = 0;
							if (ages[i]-offset > 0) {
								temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex0sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex0mu___),2)/(2*Math.pow(vertex0sigma,2)) ));
							}
							vertex0Y += Math.pow((probabilities[i]-vertex0c____*temp),2);
						}
					}

					// Calculate the y value of each vertex.
					// vertex1
					vertex1Y = 0;
					for (int i = 0; i < ages.length; i++) {
						if (one_of_the_trees_too_large[i] == false) {
							double temp = 0;
							if (ages[i]-offset > 0) {
								temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex1sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex1mu___),2)/(2*Math.pow(vertex1sigma,2)) ));
							}
							vertex1Y += Math.pow((probabilities[i]-vertex1c____*temp),2);
						}
					}

					// Calculate the y value of each vertex.
					// vertex2
					vertex2Y = 0;
					for (int i = 0; i < ages.length; i++) {
						if (one_of_the_trees_too_large[i] == false) {
							double temp = 0;
							if (ages[i]-offset > 0) {
								temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex2sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex2mu___),2)/(2*Math.pow(vertex2sigma,2)) ));
							}
							vertex2Y += Math.pow((probabilities[i]-vertex2c____*temp),2);
						}
					}

					// Calculate the y value of each vertex.
					// vertex3
					vertex3Y = 0;
					for (int i = 0; i < ages.length; i++) {
						if (one_of_the_trees_too_large[i] == false) {
							double temp = 0;
							if (ages[i]-offset > 0) {
								temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(vertex3sigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-vertex3mu___),2)/(2*Math.pow(vertex3sigma,2)) ));
							}
							vertex3Y += Math.pow((probabilities[i]-vertex3c____*temp),2);
						}
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
						if (one_of_the_trees_too_large[i] == false) {
							double temp = 0;
							if (ages[i]-offset > 0) {
								temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(reflectionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-reflectionmu___),2)/(2*Math.pow(reflectionsigma,2)) ));
							}
							reflectionY += Math.pow((probabilities[i]-reflectionc____*temp),2);
						}
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
							if (one_of_the_trees_too_large[i] == false) {
								double temp = 0;
								if (ages[i]-offset > 0) {
									temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(extensionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-extensionmu___),2)/(2*Math.pow(extensionsigma,2)) ));
								}
								extensionY += Math.pow((probabilities[i]-extensionc____*temp),2);
							}
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
							if (one_of_the_trees_too_large[i] == false) {
								double temp = 0;
								if (ages[i]-offset > 0) {
									temp = (1.0/((ages[i]-offset)*Math.sqrt(2*Math.PI*Math.pow(contractionsigma,2))))*(Math.exp( -Math.pow((Math.log(ages[i]-offset)-contractionmu___),2)/(2*Math.pow(contractionsigma,2)) ));
								}
								contractionY += Math.pow((probabilities[i]-contractionc____*temp),2);
							}
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
					if (Math.abs(vertex0c____- vertex1c____) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___- vertex1mu___) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma- vertex1sigma) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____- vertex2c____) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___- vertex2mu___) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma- vertex2sigma) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0c____- vertex3c____) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0mu___- vertex3mu___) > 0.000000001) {
						keepGoing = true;
					} else if (Math.abs(vertex0sigma- vertex3sigma) > 0.000000001) {
						keepGoing = true;
					}
					if (stepCounter > 100000) {
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
			double logRmsd = Math.sqrt(ys[index]/(double) (number_of_ages+1-number_of_ages_with_trees_that_are_too_large));

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
			double approx_distribution_rmsd = logRmsd;
			
		// }

		// XXX needs to be changed so that an object of class ParameterDistribution is returned.
		System.out.println("Distribution type: " + approx_distribution_type);
		System.out.println("Mean (log): " + approx_distribution_parameters[0]);
		System.out.println("Stdev (log): " + approx_distribution_parameters[1]);
		System.out.println("RMSD: " + approx_distribution_rmsd);

	}
	
	public static void main(String[] args) {

		CladeAgeProbabilities cladeAgeProbabilities = new CladeAgeProbabilities();
		cladeAgeProbabilities.bd_simulate(10.0,13.0,0.04,0.08,0.001,0.4,0.008,0.018,0,2,100,100000,10, new JProgressBar());
		cladeAgeProbabilities.fitLognormal();

  }
}
