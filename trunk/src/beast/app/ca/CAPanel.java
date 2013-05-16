package beast.app.ca;

import javax.swing.JPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;

import javax.swing.JLabel;
import java.awt.Insets;
import javax.swing.JTextField;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

import java.awt.Color;
import javax.swing.border.LineBorder;
import java.awt.Component;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;

import javax.swing.Box;
import javax.swing.border.TitledBorder;
import javax.swing.border.BevelBorder;

import beast.app.beauti.BeautiPanel;
import beast.math.distributions.Exponential;
import beast.math.distributions.ParametricDistribution;

public class CAPanel extends JPanel {
	private static final long serialVersionUID = 1L;

	static final String CA_ICON = "beast/app/ca/icons/cladeage_256x256px.png";

    // GUI components
	private JTextField textField_maxOccuranceAge;
	private JTextField textField_maxDivRate;
	private JTextField textField_maxTurnoverRate;
	private JTextField textField_maxSamplingRate;
	private JTextField textField_maxSamplingGap;
	private JTextField textField_minOccuranceAge;
	private JTextField textField_minDivRate;
	private JTextField textField_minTurnoverRate;
	private JTextField textField_minSamplingRate;
	private JTextField textField_minSamplginGap;
	private JTextField textField_NumberOfTreeSimulations;
	private JTextField textField_MaxNrOfBranches;
	private JTextField textField_SamplingReplicatesPerTree;
	JButton btnFindApproximation;
	JComboBox comboBox;
	JPanel panel_1;
	
	CladeAgeProbabilities probs = null;
	
	private double minOccuranceAge = 0;
	private double minDivRate = 0.01;
	private double minTurnoverRate = 0.1;
	private double minSamplingRate = 0.01;
	private double minSamplingGap = 0;

	private double maxOccuranceAge = minOccuranceAge;
	private double maxDivRate = minDivRate;
	private double maxTurnoverRate = minTurnoverRate;
	private double maxSamplingRate = minSamplingRate;
	private double maxSamplingGap = minSamplingGap;
	private int NumberOfTreeSimulations = 100;
	private int MaxNrOfBranches = 100000;
	private int SamplingReplicatesPerTree = 10;

	// used for graphing
	double [] ages;
	double [] probabilities;
	
    ParametricDistribution m_distr = null;


	
	public CAPanel() {
		GridBagLayout gridBagLayout = new GridBagLayout();
		gridBagLayout.columnWidths = new int[]{525, 375, 100};
		gridBagLayout.rowHeights = new int[]{150, 500, 20};
		gridBagLayout.columnWeights = new double[]{1.0, 1.0, 1.0};
		gridBagLayout.rowWeights = new double[]{1.0, 1.0, 1.0};
		setLayout(gridBagLayout);
		
		JPanel panel = new JPanel();
		panel.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Model parameters", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel = new GridBagConstraints();
		gbc_panel.fill = GridBagConstraints.BOTH;
		gbc_panel.gridx = 0;
		gbc_panel.gridy = 0;
		add(panel, gbc_panel);
		GridBagLayout gbl_panel = new GridBagLayout();
		gbl_panel.columnWidths = new int[]{0, 0, 0, 0, 0, 0};
		gbl_panel.rowHeights = new int[]{0, 0, 0, 0, 0, 0, 0};
		gbl_panel.columnWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		gbl_panel.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		panel.setLayout(gbl_panel);

		JLabel lblNewLabel_2 = new JLabel("Minimum");
		GridBagConstraints gbc_lblNewLabel_2 = new GridBagConstraints();
		gbc_lblNewLabel_2.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel_2.gridx = 1;
		gbc_lblNewLabel_2.gridy = 0;
		panel.add(lblNewLabel_2, gbc_lblNewLabel_2);
		
		JLabel lblNewLabel_3 = new JLabel("Maximum");
		GridBagConstraints gbc_lblNewLabel_3 = new GridBagConstraints();
		gbc_lblNewLabel_3.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel_3.gridx = 3;
		gbc_lblNewLabel_3.gridy = 0;
		panel.add(lblNewLabel_3, gbc_lblNewLabel_3);
		
		
		JLabel lblBirt = new JLabel("First occurance age:");
		GridBagConstraints gbc_lblBirt = new GridBagConstraints();
		gbc_lblBirt.anchor = GridBagConstraints.EAST;
		gbc_lblBirt.insets = new Insets(0, 0, 5, 5);
		gbc_lblBirt.gridx = 0;
		gbc_lblBirt.gridy = 1;
		panel.add(lblBirt, gbc_lblBirt);
		
		textField_minOccuranceAge = new JTextField();
		textField_minOccuranceAge.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_minOccuranceAge.setColumns(10);
		GridBagConstraints gbc_textField_7 = new GridBagConstraints();
		gbc_textField_7.insets = new Insets(0, 0, 5, 5);
		gbc_textField_7.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_7.gridx = 1;
		gbc_textField_7.gridy = 1;
		panel.add(textField_minOccuranceAge, gbc_textField_7);
		
		JLabel label = new JLabel("-");
		GridBagConstraints gbc_label = new GridBagConstraints();
		gbc_label.insets = new Insets(0, 0, 5, 5);
		gbc_label.anchor = GridBagConstraints.EAST;
		gbc_label.gridx = 2;
		gbc_label.gridy = 1;
		panel.add(label, gbc_label);
		
		textField_maxOccuranceAge = new JTextField();
		textField_maxOccuranceAge.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_maxOccuranceAge.setColumns(10);
		GridBagConstraints gbc_textField_1_1 = new GridBagConstraints();
		gbc_textField_1_1.insets = new Insets(0, 0, 5, 0);
		gbc_textField_1_1.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_1_1.gridx = 3;
		gbc_textField_1_1.gridy = 1;
		panel.add(textField_maxOccuranceAge, gbc_textField_1_1);
		
		JLabel lblNetDiversificationRate = new JLabel("<html>Net&nbsp;diversification&nbsp;rate&nbsp;&lambda;&minus;&mu;:</html>");
		GridBagConstraints gbc_lblNetDiversificationRate = new GridBagConstraints();
		gbc_lblNetDiversificationRate.anchor = GridBagConstraints.EAST;
		gbc_lblNetDiversificationRate.insets = new Insets(0, 0, 5, 5);
		gbc_lblNetDiversificationRate.gridx = 0;
		gbc_lblNetDiversificationRate.gridy = 2;
		panel.add(lblNetDiversificationRate, gbc_lblNetDiversificationRate);
		
		textField_minDivRate = new JTextField();
		textField_minDivRate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_minDivRate.setColumns(10);
		GridBagConstraints gbc_textField_8 = new GridBagConstraints();
		gbc_textField_8.insets = new Insets(0, 0, 5, 5);
		gbc_textField_8.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_8.gridx = 1;
		gbc_textField_8.gridy = 2;
		panel.add(textField_minDivRate, gbc_textField_8);
		
		JLabel label_1 = new JLabel("-");
		GridBagConstraints gbc_label_1 = new GridBagConstraints();
		gbc_label_1.insets = new Insets(0, 0, 5, 5);
		gbc_label_1.gridx = 2;
		gbc_label_1.gridy = 2;
		panel.add(label_1, gbc_label_1);
		
		textField_maxDivRate = new JTextField();
		textField_maxDivRate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_maxDivRate.setColumns(10);
		GridBagConstraints gbc_textField_2 = new GridBagConstraints();
		gbc_textField_2.insets = new Insets(0, 0, 5, 0);
		gbc_textField_2.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_2.gridx = 3;
		gbc_textField_2.gridy = 2;
		panel.add(textField_maxDivRate, gbc_textField_2);
		
		JLabel lblRurnoverRateDb = new JLabel("<html>Turnover&nbsp;rate&nbsp;&mu;/&lambda;:</html>");
		GridBagConstraints gbc_lblRurnoverRateDb = new GridBagConstraints();
		gbc_lblRurnoverRateDb.anchor = GridBagConstraints.EAST;
		gbc_lblRurnoverRateDb.insets = new Insets(0, 0, 5, 5);
		gbc_lblRurnoverRateDb.gridx = 0;
		gbc_lblRurnoverRateDb.gridy = 3;
		panel.add(lblRurnoverRateDb, gbc_lblRurnoverRateDb);
		
		textField_minTurnoverRate = new JTextField();
		textField_minTurnoverRate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_minTurnoverRate.setColumns(10);
		GridBagConstraints gbc_textField_9 = new GridBagConstraints();
		gbc_textField_9.insets = new Insets(0, 0, 5, 5);
		gbc_textField_9.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_9.gridx = 1;
		gbc_textField_9.gridy = 3;
		panel.add(textField_minTurnoverRate, gbc_textField_9);
		
		JLabel label_2 = new JLabel("-");
		GridBagConstraints gbc_label_2 = new GridBagConstraints();
		gbc_label_2.insets = new Insets(0, 0, 5, 5);
		gbc_label_2.gridx = 2;
		gbc_label_2.gridy = 3;
		panel.add(label_2, gbc_label_2);
		
		textField_maxTurnoverRate = new JTextField();
		textField_maxTurnoverRate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_maxTurnoverRate.setColumns(10);
		GridBagConstraints gbc_textField_3 = new GridBagConstraints();
		gbc_textField_3.insets = new Insets(0, 0, 5, 0);
		gbc_textField_3.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_3.gridx = 3;
		gbc_textField_3.gridy = 3;
		panel.add(textField_maxTurnoverRate, gbc_textField_3);
		
		JLabel lblSamplingRate = new JLabel("<html>Sampling&nbsp;rate&nbsp;&psi;:</html>");
		GridBagConstraints gbc_lblSamplingRate = new GridBagConstraints();
		gbc_lblSamplingRate.anchor = GridBagConstraints.EAST;
		gbc_lblSamplingRate.insets = new Insets(0, 0, 5, 5);
		gbc_lblSamplingRate.gridx = 0;
		gbc_lblSamplingRate.gridy = 4;
		panel.add(lblSamplingRate, gbc_lblSamplingRate);
		
		textField_minSamplingRate = new JTextField();
		textField_minSamplingRate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_minSamplingRate.setColumns(10);
		GridBagConstraints gbc_textField_10 = new GridBagConstraints();
		gbc_textField_10.insets = new Insets(0, 0, 5, 5);
		gbc_textField_10.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_10.gridx = 1;
		gbc_textField_10.gridy = 4;
		panel.add(textField_minSamplingRate, gbc_textField_10);
		
		JLabel label_3 = new JLabel("-");
		GridBagConstraints gbc_label_3 = new GridBagConstraints();
		gbc_label_3.insets = new Insets(0, 0, 5, 5);
		gbc_label_3.gridx = 2;
		gbc_label_3.gridy = 4;
		panel.add(label_3, gbc_label_3);
		
		textField_maxSamplingRate = new JTextField();
		textField_maxSamplingRate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_maxSamplingRate.setColumns(10);
		GridBagConstraints gbc_textField_4 = new GridBagConstraints();
		gbc_textField_4.insets = new Insets(0, 0, 5, 0);
		gbc_textField_4.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_4.gridx = 3;
		gbc_textField_4.gridy = 4;
		panel.add(textField_maxSamplingRate, gbc_textField_4);
		
		JLabel lblNewLabel = new JLabel("Sampling gap:");
		GridBagConstraints gbc_lblNewLabel = new GridBagConstraints();
		gbc_lblNewLabel.insets = new Insets(0, 0, 0, 5);
		gbc_lblNewLabel.anchor = GridBagConstraints.EAST;
		gbc_lblNewLabel.gridx = 0;
		gbc_lblNewLabel.gridy = 5;
		panel.add(lblNewLabel, gbc_lblNewLabel);
		
		textField_minSamplginGap = new JTextField();
		textField_minSamplginGap.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_minSamplginGap.setColumns(10);
		GridBagConstraints gbc_textField_11 = new GridBagConstraints();
		gbc_textField_11.insets = new Insets(0, 0, 0, 5);
		gbc_textField_11.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_11.gridx = 1;
		gbc_textField_11.gridy = 5;
		panel.add(textField_minSamplginGap, gbc_textField_11);
		
		JLabel label_4 = new JLabel("-");
		GridBagConstraints gbc_label_4 = new GridBagConstraints();
		gbc_label_4.insets = new Insets(0, 0, 0, 5);
		gbc_label_4.gridx = 2;
		gbc_label_4.gridy = 5;
		panel.add(label_4, gbc_label_4);
		
		textField_maxSamplingGap = new JTextField();
		textField_maxSamplingGap.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_maxSamplingGap.setColumns(10);
		GridBagConstraints gbc_textField_5 = new GridBagConstraints();
		gbc_textField_5.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_5.gridx = 3;
		gbc_textField_5.gridy = 5;
		panel.add(textField_maxSamplingGap, gbc_textField_5);

		// help buttons for panel1
		JButton btnHelpButtonFirstOccurance = new JButton("?");
		GridBagConstraints gbc_btnNewButton = new GridBagConstraints();
		gbc_btnNewButton.gridx = 5;
		gbc_btnNewButton.gridy = 1;
		panel.add(btnHelpButtonFirstOccurance, gbc_btnNewButton);
		btnHelpButtonFirstOccurance.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(null, "First occurance is ...");
			}
		});

		JButton btnHelpButton2= new JButton("?");
		GridBagConstraints gbc_btnNewButton2 = new GridBagConstraints();
		gbc_btnNewButton2.gridx = 5;
		gbc_btnNewButton2.gridy = 2;
		panel.add(btnHelpButton2, gbc_btnNewButton2);
		btnHelpButton2.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(null, "First occurance is ...");
			}
		});

		JButton btnHelpButton3= new JButton("?");
		GridBagConstraints gbc_btnNewButton3 = new GridBagConstraints();
		gbc_btnNewButton3.gridx = 5;
		gbc_btnNewButton3.gridy = 3;
		panel.add(btnHelpButton3, gbc_btnNewButton3);
		btnHelpButton3.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(null, "First occurance is ...");
			}
		});

		JButton btnHelpButton4= new JButton("?");
		GridBagConstraints gbc_btnNewButton4 = new GridBagConstraints();
		gbc_btnNewButton4.gridx = 5;
		gbc_btnNewButton4.gridy = 4;
		panel.add(btnHelpButton4, gbc_btnNewButton4);
		btnHelpButton4.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(null, "First occurance is ...");
			}
		});
		
		JButton btnHelpButton5= new JButton("?");
		GridBagConstraints gbc_btnNewButton5 = new GridBagConstraints();
		gbc_btnNewButton5.gridx = 5;
		gbc_btnNewButton5.gridy = 5;
		panel.add(btnHelpButton5, gbc_btnNewButton5);
		btnHelpButton5.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(null, "First occurance is ...");
			}
		});
		
		JPanel panel2 = new JPanel();
		panel2.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Simulation settings", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel2 = new GridBagConstraints();
		gbc_panel2.fill = GridBagConstraints.BOTH;
		gbc_panel2.gridx = 1;
		gbc_panel2.gridy = 0;
		add(panel2, gbc_panel2);
		GridBagLayout gbl_panel2 = new GridBagLayout();
		gbl_panel2.columnWidths = new int[]{0, 0};
		gbl_panel2.rowHeights = new int[]{0, 0, 0, 0, 0, 0, 0};
		gbl_panel2.columnWeights = new double[]{0, 0};
		gbl_panel2.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		panel2.setLayout(gbl_panel2);

		
		JLabel lblIcon = new MyJLabel(BeautiPanel.getIcon(CA_ICON));
		lblIcon.setMinimumSize(new Dimension(160,160));
		lblIcon.setPreferredSize(new Dimension(128,128));
		GridBagConstraints gbc_lblIcon = new GridBagConstraints();
		gbc_lblIcon.insets = new Insets(5, 0, 5, 5);
		gbc_lblIcon.anchor = GridBagConstraints.EAST;
		gbc_lblIcon.gridx = 2;
		gbc_lblIcon.gridy = 0;
		gbc_lblIcon.gridwidth = 2;
		add(lblIcon, gbc_lblIcon);
		
		textField_NumberOfTreeSimulations = new JTextField();
		textField_NumberOfTreeSimulations.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		GridBagConstraints gbc_textField = new GridBagConstraints();
		gbc_textField.insets = new Insets(0, 0, 5, 0);
		gbc_textField.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField.gridx = 1;
		gbc_textField.gridy = 1;
		panel2.add(textField_NumberOfTreeSimulations, gbc_textField);
		textField_NumberOfTreeSimulations.setColumns(10);
		
		textField_MaxNrOfBranches = new JTextField();
		textField_MaxNrOfBranches.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_MaxNrOfBranches.setColumns(10);
		GridBagConstraints gbc_textField_12 = new GridBagConstraints();
		gbc_textField_12.insets = new Insets(0, 0, 5, 0);
		gbc_textField_12.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_12.gridx = 1;
		gbc_textField_12.gridy = 2;
		panel2.add(textField_MaxNrOfBranches, gbc_textField_12);
		
		JLabel lblSamplingReplicatesPer = new JLabel("Sampling replicates per tree:");
		GridBagConstraints gbc_lblSamplingReplicatesPer = new GridBagConstraints();
		gbc_lblSamplingReplicatesPer.insets = new Insets(0, 0, 5, 5);
		gbc_lblSamplingReplicatesPer.anchor = GridBagConstraints.EAST;
		gbc_lblSamplingReplicatesPer.gridx = 0;
		gbc_lblSamplingReplicatesPer.gridy = 3;
		panel2.add(lblSamplingReplicatesPer, gbc_lblSamplingReplicatesPer);
		
		textField_SamplingReplicatesPerTree = new JTextField();
		textField_SamplingReplicatesPerTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		textField_SamplingReplicatesPerTree.setColumns(10);
		GridBagConstraints gbc_textField_13 = new GridBagConstraints();
		gbc_textField_13.insets = new Insets(0, 0, 5, 0);
		gbc_textField_13.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField_13.gridx = 1;
		gbc_textField_13.gridy = 3;
		panel2.add(textField_SamplingReplicatesPerTree, gbc_textField_13);
		
		JButton btnCalculate = new JButton("Estimate clade age probabilities");
		btnCalculate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				guiToData();
				probs = new CladeAgeProbabilities();
				
				Frame parentFrame = Frame.getFrames()[0];
			    final JDialog dlg = new JDialog(parentFrame, "Progress Dialog", true);
			    final JProgressBar dpb = new JProgressBar(0, NumberOfTreeSimulations);
			    final JButton cancelButton = new JButton("Cancel");
			    cancelButton.addActionListener(new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						try {
							System.err.println("Trying to stop");
							probs.setCancel1();
							dlg.setVisible(false);
						} catch (Exception ex) {
							ex.printStackTrace();
						}
						
					}
				});
			    dlg.add(BorderLayout.CENTER, dpb);
			    dlg.add(BorderLayout.NORTH, new JLabel("Progress..."));
			    dlg.add(BorderLayout.EAST, cancelButton);
			    dlg.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
			    dlg.setSize(300, 75);
			    dlg.setLocationRelativeTo(parentFrame);

			    final Thread t = new Thread(new Runnable() {
			      public void run() {
					probs.bd_simulate(
							minOccuranceAge, maxOccuranceAge,
							minDivRate, maxDivRate,
							minTurnoverRate, maxTurnoverRate,
							minSamplingRate, maxSamplingRate,
							minSamplingGap, maxSamplingGap,
							NumberOfTreeSimulations, MaxNrOfBranches, SamplingReplicatesPerTree, dpb);
			        dlg.setVisible(false);
			      }
			    });
			    t.start();
		        dlg.setVisible(true);

			    while (dlg.isVisible()) {
			    	try {
			    		Thread.sleep(500);
			    	} catch (Exception ex) {
			    		
			    	}
			    }
			    if (probs.getCancel1()) {
			    	return;
			    }
		        
				ages = probs.getAges();
				probabilities =  probs.getProbabilities();
				// normalize
				double sum = probabilities[0] * (ages[1] - ages[0]);;
				for (int i = 1; i < probabilities.length-1; i++) {
					sum += probabilities[i] * (ages[i-1] - ages[i+1])/2.0;
				}
				sum += probabilities[probabilities.length-1] * (ages[probabilities.length-2] - ages[probabilities.length-1]);
				
				for (int i = 0; i < probabilities.length; i++) {
					probabilities[i] /= sum;
				}
				

				btnFindApproximation.setEnabled(true);
				comboBox.setEnabled(true);
				panel_1.repaint();
			}
		});
		GridBagConstraints gbc_btnCalculate = new GridBagConstraints();
		gbc_btnCalculate.gridwidth = 2;
		gbc_btnCalculate.insets = new Insets(0, 0, 5, 0);
		gbc_btnCalculate.gridx = 0;
		gbc_btnCalculate.gridy = 4;
		panel2.add(btnCalculate, gbc_btnCalculate);
		
		JLabel lblMaximumNumberOf = new JLabel("Maximum number of branches:");
		GridBagConstraints gbc_lblMaximumNumberOf = new GridBagConstraints();
		gbc_lblMaximumNumberOf.anchor = GridBagConstraints.EAST;
		gbc_lblMaximumNumberOf.insets = new Insets(0, 0, 5, 5);
		gbc_lblMaximumNumberOf.gridx = 0;
		gbc_lblMaximumNumberOf.gridy = 2;
		panel2.add(lblMaximumNumberOf, gbc_lblMaximumNumberOf);
		
		JLabel lblNewLabel_1 = new JLabel("Number of tree simulations:");
		GridBagConstraints gbc_lblNewLabel_1 = new GridBagConstraints();
		gbc_lblNewLabel_1.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel_1.anchor = GridBagConstraints.EAST;
		gbc_lblNewLabel_1.gridx = 0;
		gbc_lblNewLabel_1.gridy = 1;
		panel2.add(lblNewLabel_1, gbc_lblNewLabel_1);
		
		Component verticalGlue = Box.createVerticalGlue();
		GridBagConstraints gbc_verticalGlue = new GridBagConstraints();
		gbc_verticalGlue.gridwidth = 2;
		gbc_verticalGlue.insets = new Insets(0, 0, 5, 5);
		gbc_verticalGlue.gridx = 0;
		gbc_verticalGlue.gridy = 5;
		panel2.add(verticalGlue, gbc_verticalGlue);

		JPanel panel3 = new JPanel();
		panel3.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Clade age probabilities", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel3 = new GridBagConstraints();
		gbc_panel3.anchor = GridBagConstraints.WEST;
		gbc_panel3.fill = GridBagConstraints.BOTH;
		gbc_panel3.gridx = 0;
		gbc_panel3.gridy = 1;
		gbc_panel3.gridwidth = 3;
		add(panel3, gbc_panel3);
		GridBagLayout gbl_panel3 = new GridBagLayout();
		gbl_panel3.columnWidths = new int[]{0};
		gbl_panel3.rowHeights = new int[]{0, 103};
		gbl_panel3.columnWeights = new double[]{1.0};
		gbl_panel3.rowWeights = new double[]{Double.MIN_VALUE, 1.0};
		panel3.setLayout(gbl_panel3);
		
		panel_1 = new JPanel() {
			protected void paintComponent(java.awt.Graphics g) {
					g.setColor(new Color(241, 241, 241));
				g.fillRect(0, 0, getWidth(), getHeight());
				if (ages != null) {
		            final int width = getWidth();
		            final int height = getHeight();
		            final int graphoffset = 30;
		            final int labeloffset = 0;
		            int nGraphWidth = width - graphoffset * 2;
		            int nGraphHeight = height - graphoffset * 2 - labeloffset;
		            g.setColor(Color.WHITE);
		            g.fillRect(graphoffset, graphoffset, nGraphWidth, nGraphHeight);
		            g.setColor(Color.BLACK);
		            g.drawRect(graphoffset, graphoffset, nGraphWidth, nGraphHeight);
		            try {
		            	if (m_distr != null) {
		            		m_distr.initAndValidate();
		            	}
		            } catch (Exception e1) {
		                // ignore
		            }
		            int nPoints = ages.length;
		            int[] xPoints = new int[nPoints];
		            int[] yPoints = new int[nPoints];
		            int[] yPoints2 = new int[nPoints];
		            double[] fyPoints = new double[nPoints];
		            double[] fyPoints2 = new double[nPoints];
		            Font font = g.getFont();
		            double fMinValue = 0.1;
		            double fMaxValue = 1;
	                fMinValue = ages[ages.length-1];
	                fMaxValue = ages[0];
		            double fXRange = fMaxValue - fMinValue;
		            // adjust fYMax so that the ticks come out right
		            double fX0 = fMinValue;
		            int k = 0;
		            double f = fXRange;
		            double f2 = fX0;
		            while (f > 10) {
		                f /= 10;
		                f2 /= 10;
		                k++;
		            }
		            while (f < 1 && f > 0) {
		                f *= 10;
		                f2 *= 10;
		                k--;
		            }
		            f = Math.ceil(f);
		            f2 = Math.floor(f2);
//					final int NR_OF_TICKS_X = NR_OF_TICKS[(int) f];
		            for (int i = 0; i < k; i++) {
		                f *= 10;
		                f2 *= 10;
		            }
		            for (int i = k; i < 0; i++) {
		                f /= 10;
		                f2 /= 10;
		            }
		            //double fAdjXRange = f;

		            fXRange = fXRange + fMinValue - f2;
		            fXRange = adjust(fXRange);
		            final int NR_OF_TICKS_X = m_nTicks;

		            fMinValue = f2; //fXRange = fAdjXRange;

		            double fYMax = 0;
		            for (int i = 0; i < nPoints; i++) {
		                xPoints[i] = graphoffset + nGraphWidth * i / nPoints;
		                fyPoints[i] = probabilities[nPoints - i - 1];
		                if (m_distr != null) {
		                    try {
		                        fyPoints2[i] = m_distr.density(fMinValue + (fXRange * i) / nPoints);
		                    } catch (Exception e) {
		                        fyPoints2[i] = 0;
		                    }
		                }
		                if (Double.isInfinite(fyPoints[i]) || Double.isNaN(fyPoints[i])) {
		                    fyPoints[i] = 0;
		                }
		                //fyPoints[i] = Math.exp(m_distr.logDensity(fMinValue + (fXRange * i)/nPoints));
		                fYMax = Math.max(fYMax, fyPoints[i]);
		            }

		            fYMax = adjust(fYMax);
		            final int NR_OF_TICKS_Y = m_nTicks;

		            
		            for (int i = 0; i < nPoints; i++) {
		                yPoints[i] = 1 + (int) (graphoffset + nGraphHeight - nGraphHeight * fyPoints[i] / fYMax);
		                yPoints2[i] = 1 + (int) (graphoffset + nGraphHeight - nGraphHeight * fyPoints2[i] / fYMax);
		                g.drawLine(xPoints[i]+2, yPoints[i], xPoints[i]-2, yPoints[i]);
		                g.drawLine(xPoints[i], yPoints[i]-2, xPoints[i], yPoints[i]+2);
		            }
		            if (m_distr != null) {
		            	g.drawPolyline(xPoints, yPoints2, nPoints);
		            }
		            

		            // draw ticks on edge
		            Font smallFont = new Font(font.getName(), font.getStyle(), 8);
		            g.setFont(smallFont);
		            for (int i = 0; i <= NR_OF_TICKS_X; i++) {
		                int x = graphoffset + i * nGraphWidth / NR_OF_TICKS_X;
		                g.drawLine(x, graphoffset + nGraphHeight, x, graphoffset + nGraphHeight + 5);
		                g.drawString(format(fMinValue + fXRange * i / NR_OF_TICKS_X), x + 2, graphoffset + nGraphHeight + 5 + 4);
		            }
		            for (int i = 0; i <= NR_OF_TICKS_Y; i++) {
		                int y = graphoffset + nGraphHeight - i * nGraphHeight / NR_OF_TICKS_Y;
		                g.drawLine(graphoffset - 5, y, graphoffset, y);
		                g.drawString(format(fYMax * i / NR_OF_TICKS_Y), 0, y - 2);
		            }
					
				}
			};
		    /**
		     * maps most significant digit to nr of ticks on graph *
		     */
		    final int[] NR_OF_TICKS = new int[]{5, 10, 8, 6, 8, 10, 6, 7, 8, 9, 10};
			int m_nTicks = 5;
			
	        private double adjust(double fYMax) {
	            // adjust fYMax so that the ticks come out right
	            int k = 0;
	            double fY = fYMax;
	            while (fY > 10) {
	                fY /= 10;
	                k++;
	            }
	            while (fY < 1 && fY > 0) {
	                fY *= 10;
	                k--;
	            }
	            fY = Math.ceil(fY);
	            m_nTicks = NR_OF_TICKS[(int) fY];
	            for (int i = 0; i < k; i++) {
	                fY *= 10;
	            }
	            for (int i = k; i < 0; i++) {
	                fY /= 10;
	            }
	            return fY;
	        }
	        
	        private String format(double value) {
	            StringWriter writer = new StringWriter();
	            PrintWriter pw = new PrintWriter(writer);
	            pw.printf("%.3g", value);
	            pw.flush();
	            return writer.toString();
	        }

		};
		panel_1.setBorder(new LineBorder(new Color(0, 0, 0)));
		panel_1.setBackground(Color.gray);
		panel_1.setPreferredSize(new Dimension(1024,400));
		GridBagConstraints gbc_panel_1 = new GridBagConstraints();
		gbc_panel_1.fill = GridBagConstraints.BOTH;
		gbc_panel_1.gridx = 0;
		gbc_panel_1.gridy = 1;
		panel3.add(panel_1, gbc_panel_1);

		JPanel panel4 = new JPanel();
		panel4.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Approximation", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(51, 51, 51)));
		GridBagConstraints gbc_panel4 = new GridBagConstraints();
		gbc_panel4.fill = GridBagConstraints.BOTH;
		gbc_panel4.gridx = 0;
		gbc_panel4.gridy = 2;
		gbc_panel4.gridwidth = 3;
		add(panel4, gbc_panel4);
		GridBagLayout gbl_panel4 = new GridBagLayout();
		gbl_panel4.columnWidths = new int[]{0, 0, 0, 0};
		gbl_panel4.rowHeights = new int[]{0, 0, 0};
		gbl_panel4.columnWeights = new double[]{0.0, 0.0, 0.0, 0.0};
		gbl_panel4.rowWeights = new double[]{Double.MIN_VALUE, 0.0, 0.0};
		panel4.setLayout(gbl_panel4);
		
		JLabel lblDistributionType = new JLabel("Distribution type:");
		GridBagConstraints gbc_lblDistributionType = new GridBagConstraints();
		gbc_lblDistributionType.insets = new Insets(0, 0, 0, 5);
		gbc_lblDistributionType.anchor = GridBagConstraints.EAST;
		gbc_lblDistributionType.gridx = 0;
		gbc_lblDistributionType.gridy = 2;
		panel4.add(lblDistributionType, gbc_lblDistributionType);
		
		comboBox = new JComboBox(new String[]{"Best fit","Exponential","Gamma","Log Normal"});
		GridBagConstraints gbc_comboBox = new GridBagConstraints();
		gbc_comboBox.insets = new Insets(0, 0, 0, 5);
		gbc_comboBox.fill = GridBagConstraints.HORIZONTAL;
		gbc_comboBox.gridx = 1;
		gbc_comboBox.gridy = 2;
		panel4.add(comboBox, gbc_comboBox);
		
		btnFindApproximation = new JButton("Find approximation");
		btnFindApproximation.setEnabled(false);
		comboBox.setEnabled(false);
		btnFindApproximation.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				String type = (String) comboBox.getSelectedItem();
				if (type.equals("Best fit")) {
					JOptionPane.showMessageDialog(getParent(), "Not implemented yet");
					return;
				}
				if (type.equals("Exponential")) {
					double mean = probs.fitExponential();
					m_distr = new Exponential();
					try {
						m_distr.initByName("mean", mean + "");
					} catch (Exception e2) {
						e2.printStackTrace();
					}
					panel_1.repaint();
					return;
				}
				if (type.equals("Gamma")) {
					JOptionPane.showMessageDialog(getParent(), "Not implemented yet");
					return;
				}
				if (type.equals("Log Normal")) {
					JOptionPane.showMessageDialog(getParent(), "Not implemented yet");
					return;
				}
				
				
				
			}
		});
		GridBagConstraints gbc_btnFindApproximation = new GridBagConstraints();
		gbc_btnFindApproximation.insets = new Insets(0, 0, 0, 5);
		gbc_btnFindApproximation.gridx = 2;
		gbc_btnFindApproximation.gridy = 2;
		panel4.add(btnFindApproximation, gbc_btnFindApproximation);
		
		Component horizontalGlue = Box.createHorizontalGlue();
		GridBagConstraints gbc_horizontalGlue = new GridBagConstraints();
		gbc_horizontalGlue.gridx = 3;
		gbc_horizontalGlue.gridy = 2;
		panel4.add(horizontalGlue, gbc_horizontalGlue);
		
		dataToGUI();
	}

	void dataToGUI() {
		textField_minOccuranceAge.setText(minOccuranceAge + "");
		textField_minDivRate.setText(minDivRate + "");
		textField_minTurnoverRate.setText(minTurnoverRate + "");
		textField_minSamplingRate.setText(minSamplingRate + "");
		textField_minSamplginGap.setText(minSamplingGap + "");
		textField_NumberOfTreeSimulations.setText(NumberOfTreeSimulations + "");
		textField_MaxNrOfBranches.setText(MaxNrOfBranches + "");
		textField_SamplingReplicatesPerTree.setText(SamplingReplicatesPerTree + "");
		if (maxOccuranceAge != minOccuranceAge) {
			textField_maxOccuranceAge.setText(maxOccuranceAge + "");
		} else {
			textField_maxOccuranceAge.setText("As minimum");
		}
		if (maxDivRate != minDivRate) {
			textField_maxDivRate.setText(maxDivRate + "");
		} else {
			textField_maxDivRate.setText("As minimum");
		}
		if (maxTurnoverRate != minTurnoverRate) {
			textField_maxTurnoverRate.setText(maxTurnoverRate + "");
		} else {
			textField_maxTurnoverRate.setText("As minimum");
		}
		if (maxSamplingRate != minSamplingRate) {
			textField_maxSamplingRate.setText(maxSamplingRate + "");
		} else {
			textField_maxSamplingRate.setText("As minimum");
		}
		if (maxSamplingGap != minSamplingGap) {
			textField_maxSamplingGap.setText(maxSamplingGap + "");
		} else {
			textField_maxSamplingGap.setText("As minimum");
		}
	}
	
	void guiToData() {
		maxOccuranceAge = parseDouble(textField_maxOccuranceAge.getText());
		maxDivRate = parseDouble(textField_maxDivRate.getText());
		maxTurnoverRate = parseDouble(textField_maxTurnoverRate.getText());
		maxSamplingRate = parseDouble(textField_maxSamplingRate.getText());
		maxSamplingGap = parseDouble(textField_maxSamplingGap.getText());
	
		minOccuranceAge = parseDouble(textField_minOccuranceAge.getText());
		if (minOccuranceAge < 0) {
			minOccuranceAge = 0;
		}
		if (Double.isInfinite(maxOccuranceAge) || maxOccuranceAge < minOccuranceAge) {
			maxOccuranceAge = minOccuranceAge;
		}
		minDivRate = parseDouble(textField_minDivRate.getText());
		if (minDivRate < 0) {
			minDivRate = 0;
		}
		if (Double.isInfinite(maxDivRate) || maxDivRate < minDivRate) {
			maxDivRate = minDivRate;
		}
		minTurnoverRate = parseDouble(textField_minTurnoverRate.getText());
		if (minTurnoverRate < 0) {
			minTurnoverRate = 0;
		}
		if (Double.isInfinite(maxTurnoverRate) || maxTurnoverRate < minTurnoverRate) {
			maxTurnoverRate = minTurnoverRate;
		}
		minSamplingRate = parseDouble(textField_minSamplingRate.getText());
		if (minSamplingRate < 1e-10) {
			minSamplingRate = 1e-10;
		}
		if (Double.isInfinite(maxSamplingRate) || maxSamplingRate < minSamplingRate) {
			maxSamplingRate = minSamplingRate;
		}
		minSamplingGap = parseDouble(textField_minSamplginGap.getText());
		if (minSamplingGap < 0) {
			minSamplingGap = 0;
		}
		if (Double.isInfinite(maxSamplingGap) || maxSamplingGap < minSamplingGap) {
			maxSamplingGap = minSamplingGap;
		}
		NumberOfTreeSimulations = parseInt(textField_NumberOfTreeSimulations.getText());
		MaxNrOfBranches = parseInt(textField_MaxNrOfBranches.getText());
		SamplingReplicatesPerTree = parseInt(textField_SamplingReplicatesPerTree.getText());

		// something changed, so the probabilities are not valid any more
		btnFindApproximation.setEnabled(false);
		comboBox.setEnabled(false);
		m_distr = null;
		dataToGUI();
	}

	private int parseInt(String text) {
		try {
			return Integer.parseInt(text);
		} catch (Exception e) {
			
		}
		return 0;
	}

	private double parseDouble(String text) {
		try {
			return Double.parseDouble(text);
		} catch (Exception e) {
			
		}
		return 0;
	}
	
	class MyJLabel extends JLabel {
		private static final long serialVersionUID = 1L;
		ImageIcon imageIcon;
	    public MyJLabel(ImageIcon icon)
	    {
	        super();
	        this.imageIcon = icon;
	    }
	    @Override
	    public void paintComponent(Graphics g)
	    {
	        super.paintComponent(g);
	        g.drawImage(imageIcon.getImage(),0,0,getWidth(),getHeight(),this);
	    }
	}

	public static void main(String[] args) {
		JFrame frame = new JFrame();
		frame.setSize(1024, 728);
        ImageIcon icon = BeautiPanel.getIcon(CA_ICON);
        if (icon != null) {
            frame.setIconImage(icon.getImage());
        }
		CAPanel pane = new CAPanel();
		frame.getContentPane().add(pane);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);	
		frame.setVisible(true);
	}
	
}