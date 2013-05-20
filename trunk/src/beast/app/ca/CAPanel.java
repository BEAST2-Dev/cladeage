package beast.app.ca;



import jam.framework.DocumentFrame;

import javax.swing.JPanel;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Cursor;
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
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JProgressBar;
import javax.swing.ToolTipManager;

import java.awt.Color;
import javax.swing.border.LineBorder;
import java.awt.Component;
import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Method;
import java.net.URL;

import javax.swing.Box;
import javax.swing.border.TitledBorder;
import javax.swing.border.BevelBorder;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.app.beauti.BeautiPanel;
import beast.app.draw.ModelBuilder;
import beast.app.draw.MyAction;
import beast.app.util.Utils;
import beast.math.distributions.ExpGamma;
import beast.math.distributions.LogNormalImpl;

public class CAPanel extends JPanel {
	private static final long serialVersionUID = 1L;

	static final String CA_ICON = "beast/app/ca/icons/cladeage_256x256px.png";

	final public static String OCCURRENCE_AGE_HELP = "<html>First occurrence age:<br/>"+
"<br/>"+
"The age of the oldest fossil of the clade. <br/>"+
"If this age is known exactly, it should be specified in the 'Minimum' <br/>"+
"field. If not, you should specify both a minimum and maximum age.<br/>" +
"<br/>" +
"Consider a fosil described in the literature to have come from the oligocene,<br/>" +
"which extends from about 34 million to 23 million years before the present <br/>" +
"(33.9±0.1 to 23.03±0.05 Ma). The uncertainties in the boundaries are small <br/>" +
"enough to be ignored, so the minimum is 23.03 and maximum 33.9 if you want <br/>" +
"to express ages in millions of years. <br/>" +
"</html>";

	final public static String DIV_RATE_HELP = "<html>Net diversification rate:<br/>"+
"<br/>"+
"The net diversification rate is the difference between <br/>"+
"speciation and extinction rate. Estimates can be obtained <br/>"+
"e.g. with MEDUSA (Alfaro et al. 2009) or TreePar (Stadler 2011).<br/>"+
"-<br/>"+
"Alfaro et al. (2009) PNAS 106, 13410-13414 <br/>"+
"<a href='http://www.ncbi.nlm.nih.gov/pubmed/19633192'>http://www.ncbi.nlm.nih.gov/pubmed/19633192</a><br/>"+
"Stadler (2009) PNAS 108, 6187-6192, <br/>"+
"<a href='http://www.ncbi.nlm.nih.gov/pubmed/19631666'>http://www.ncbi.nlm.nih.gov/pubmed/19631666</a></html>";

	final public static String TURNOVER_RATE_HELP = "<html>Turnover rate:<br/>"+
"<br/>"+
"The turnover rate is the ratio of extinction and speciation <br/>"+
"rate. Estimates can be obtained e.g. with MEDUSA (Alfaro et al. <br/>"+
"2009) or TreePar (Stadler 2011).<br/>"+
"-<br/>"+
"Alfaro et al. (2009) PNAS 106, 13410-13414 <br/>"+
"<a href='http://www.ncbi.nlm.nih.gov/pubmed/19633192'>http://www.ncbi.nlm.nih.gov/pubmed/19633192</a><br/>"+
"Stadler (2009) PNAS 108, 6187-6192, <br/>"+
"<a href='http://www.ncbi.nlm.nih.gov/pubmed/19631666'>http://www.ncbi.nlm.nih.gov/pubmed/19631666</a></html>";

	final public static String SAMPLING_RATE_HELP = "<html>Sampling rate:<br/>"+
"<br/>"+
"The sampling rate (sometimes called 'preservation rate') <br/>"+
"includes all processes that result in the publication of a fossil, <br/>"+
"including fossilization, discovery, identification and description <br/>"+
"(Friedman & Brazeau 2011). Estimates have been reported for many <br/>"+
"taxonomic groups (e.g. Foote et al. 1999, Foote & Sepkoski 1999) <br/>"+
"and can be obtained by methods outlined in Foote (1997).<br/>"+
"-<br/>"+
"Foote (1997) Paleobiology 23, 278-300 Proc R Soc B 278, 432-439, <br/>"+
"<a href='http://www.jstor.org/stable/2401105'>http://www.jstor.org/stable/2401105</a><br/>"+
"Foote et al. (1999) Science 283, 1310-1314, <br/>"+
"<a href='http://www.ncbi.nlm.nih.gov/pubmed/10037598'>http://www.ncbi.nlm.nih.gov/pubmed/10037598</a><br/>"+
"Foote & Sepkoski (1999) Nature 398, 415-417, <br/>"+
"<a href='http://www.ncbi.nlm.nih.gov/pubmed/11536900'>http://www.ncbi.nlm.nih.gov/pubmed/11536900</a><br/>"+
"Friedman & Brazeau (2011), <br/>"+
"<a href='http://www.ncbi.nlm.nih.gov/pubmed/20739322‎'>http://www.ncbi.nlm.nih.gov/pubmed/20739322‎</a></html>";

	final public static String SAMPLING_GAP_HELP = "<html>Sampling gap:<br/>"+
"<br/>"+
"The sampling gap represents the time period after a clade's <br/>"+
"origin during which it could not have fossilized (possibly due to <br/>"+
"small population size or limited geographic distribution), or its <br/>"+
"earliest fossils could not be recognized as part of this clade (as <br/>"+
"no apomorphies may have evolved yet). Specifying a sampling gap <br/>"+
"is optional. A sampling gap of 0.0-2.0 Ma may be a <br/>"+
"reasonable assumption.</html>";

	final public static String NR_SIMULATIONS_HELP = "<html>Number of tree simulations:<br/>"+
"<br/>"+
"The number of times birth-death trees are <br/>"+
"simulated to estimate the total unobserved lineage duration of a <br/>"+
"clade given the net diversification rate and the turnover rate. The <br/>"+
"higher this number, the better the estimates of the unobserved history, <br/>"+
"and thus of clade age probabilities. The default of 1000 <br/>"+
"tree simulations usually works well. But, if the simulated distribution <br/>" +
"looks bumpy, increasing this may help. If the calculation time is too long <br/>" +
"decreasing this number may help.</html>";

	final public static String MAX_NR_TREES_HELP = "<html>Maximum number of trees:<br/>"+
"<br/>"+
"If the net diversification rate is high (> 0.5), <br/>"+
"the simulated trees could potentially become very large, which can <br/>"+
"substantially slow down the estimation of clade age probabilities. <br/>"+
"The maximum number of branches limits the size of simulated trees. <br/>"+
"If the limit should be hit for a given age in one of the simulations, <br/>"+
"this age will not be further analysed and it <br/>"+
"will not appear in the output. There is no reason to decrease <br/>"+
"this number, but in extreme cases you might want <br/>"+
"to increase it.</html>";

	final public static String REPS_PER_TREE_HELP = "<html>Sampling replicates per tree:<br/>"+
"<br/>"+
"If a range of sampling rates has been specified, <br/>"+
"this number specifies how often rates will be drawn at random <br/>"+
"from this range for each simulated tree to estimate the probabilities <br/>"+
"of unobserved lineage durations. The default of 10 sampling <br/>"+
"replicates per tree usually works well.<br/></html>";

	final public static String ABOUT_HELP = "<html>CladeAge:<br/><br/>" +
			"Copyright 2013<br/><br/>" +
			"Michael Matschiner<br/>michaelmatschiner@mac.com<br/>" +
			"<br/>and <br/><br/>" +
			"Remco Bouckaert<br/>remco@cs.auckland.ac.nz<br/>" +
			"</html>";
	
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
	JPanel panel;
	JPanel panel2;
	JPanel panel_1;
	GridBagConstraints gbc_panel2;
	
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
	private int NumberOfTreeSimulations = 1000;
	private int MaxNrOfBranches = 100000;
	private int SamplingReplicatesPerTree = 10;

	// used for graphing
	double [] ages;
	double [] probabilities;
	
	ContinuousDistribution m_distr = null;
	double m_rmsd = 0;


	GridBagLayout gridBagLayout = new GridBagLayout();
	
	public CAPanel() {
		GridBagLayout gridBagLayout = new GridBagLayout();
		//gridBagLayout.columnWidths = new int[]{525, 375, 100};
		gridBagLayout.columnWidths = new int[]{800, 100, 200};
		gridBagLayout.rowHeights = new int[]{180, 450};
		gridBagLayout.columnWeights = new double[]{1.0, 0.0, 1.0};
		gridBagLayout.rowWeights = new double[]{1.0, 0.0, 1.0};
		setLayout(gridBagLayout);
		
		panel = new JPanel();
		panel.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Model parameters", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel = new GridBagConstraints();
		gbc_panel.fill = GridBagConstraints.BOTH;
		gbc_panel.gridx = 0;
		gbc_panel.gridy = 0;
		add(panel, gbc_panel);
		GridBagLayout gbl_panel = new GridBagLayout();
		gbl_panel.columnWidths = new int[]{0, 0, 0, 0, 0, 0};
		gbl_panel.rowHeights = new int[]{0, 0, 0, 0, 0, 0, 0};
		gbl_panel.columnWeights = new double[]{0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
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
		JButton btnHelpButtonFirstOccurance = newHelpButton();
		GridBagConstraints gbc_btnNewButton = new GridBagConstraints();
		gbc_btnNewButton.gridx = 5;
		gbc_btnNewButton.gridy = 1;
		panel.add(btnHelpButtonFirstOccurance, gbc_btnNewButton);
		btnHelpButtonFirstOccurance.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(OCCURRENCE_AGE_HELP);
			}

		});

		JButton btnHelpButton2= newHelpButton();
		GridBagConstraints gbc_btnNewButton2 = new GridBagConstraints();
		gbc_btnNewButton2.gridx = 5;
		gbc_btnNewButton2.gridy = 2;
		panel.add(btnHelpButton2, gbc_btnNewButton2);
		btnHelpButton2.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(DIV_RATE_HELP);
			}
		});

		JButton btnHelpButton3= newHelpButton();
		GridBagConstraints gbc_btnNewButton3 = new GridBagConstraints();
		gbc_btnNewButton3.gridx = 5;
		gbc_btnNewButton3.gridy = 3;
		panel.add(btnHelpButton3, gbc_btnNewButton3);
		btnHelpButton3.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(TURNOVER_RATE_HELP);
			}
		});

		JButton btnHelpButton4= newHelpButton();
		GridBagConstraints gbc_btnNewButton4 = new GridBagConstraints();
		gbc_btnNewButton4.gridx = 5;
		gbc_btnNewButton4.gridy = 4;
		panel.add(btnHelpButton4, gbc_btnNewButton4);
		btnHelpButton4.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(SAMPLING_RATE_HELP);
			}
		});
		
		JButton btnHelpButton5= newHelpButton();
		GridBagConstraints gbc_btnNewButton5 = new GridBagConstraints();
		gbc_btnNewButton5.gridx = 5;
		gbc_btnNewButton5.gridy = 5;
		panel.add(btnHelpButton5, gbc_btnNewButton5);
		btnHelpButton5.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(SAMPLING_GAP_HELP);
			}
		});
		
		panel2 = new JPanel();
		panel2.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Simulation settings", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		gbc_panel2 = new GridBagConstraints();
		gbc_panel2.fill = GridBagConstraints.BOTH;
		gbc_panel2.gridx = 1;
		gbc_panel2.gridy = 0;
		//add(panel2, gbc_panel2);
		
		GridBagLayout gbl_panel2 = new GridBagLayout();
		gbl_panel2.columnWidths = new int[]{0, 0, 0};
		gbl_panel2.rowHeights = new int[]{0, 0, 0, 0, 0, 0, 0};
		gbl_panel2.columnWeights = new double[]{0, 1.0, 0};
		gbl_panel2.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		panel2.setLayout(gbl_panel2);
		
		JPanel panel2b = new JPanel();
		panel2b.setLayout(new BorderLayout());
		
		JLabel lblIcon = new MyJLabel(BeautiPanel.getIcon(CA_ICON));
		lblIcon.setMinimumSize(new Dimension(160,160));
		lblIcon.setPreferredSize(new Dimension(128,128));
		panel2b.add(lblIcon, BorderLayout.NORTH);
		
		JButton btnCalculate = new JButton("Run");
		btnCalculate.setMinimumSize(new Dimension(128, 20));
		btnCalculate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				setCursor(new Cursor(Cursor.WAIT_CURSOR));
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
		        
			    // XXX todo
				ages = probs.getAges();
				probabilities =  probs.getProbabilities();

				calcFit();
				panel_1.repaint();
				setCursor(new Cursor(Cursor.DEFAULT_CURSOR));				
			}
		});
//		GridBagConstraints gbc_btnCalculate = new GridBagConstraints();
//		gbc_btnCalculate.gridwidth = 2;
//		gbc_btnCalculate.insets = new Insets(0, 0, 5, 0);
//		gbc_btnCalculate.gridx = 0;
//		gbc_btnCalculate.gridy = 4;
		panel2b.add(btnCalculate, BorderLayout.SOUTH);
		
		GridBagConstraints gbc_lblIcon = new GridBagConstraints();
		gbc_lblIcon.insets = new Insets(5, 0, 5, 5);
		gbc_lblIcon.anchor = GridBagConstraints.EAST;
		gbc_lblIcon.gridx = 2;
		gbc_lblIcon.gridy = 0;
		gbc_lblIcon.gridwidth = 2;
		
		add(panel2b, gbc_lblIcon);
		
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

		// help buttons for panel2
		JButton btnHelpNrSimulateions = newHelpButton();
		GridBagConstraints gbc_btnNewButton12 = new GridBagConstraints();
		gbc_btnNewButton12.gridx = 3;
		gbc_btnNewButton12.gridy = 1;
		panel2.add(btnHelpNrSimulateions, gbc_btnNewButton12);
		btnHelpNrSimulateions.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(NR_SIMULATIONS_HELP);
			}
		});
		
		JButton btnHelpNrBranches = newHelpButton();
		GridBagConstraints gbc_btnNewButton13 = new GridBagConstraints();
		gbc_btnNewButton13.gridx = 3;
		gbc_btnNewButton13.gridy = 2;
		panel2.add(btnHelpNrBranches, gbc_btnNewButton13);
		btnHelpNrBranches.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(MAX_NR_TREES_HELP);
			}
		});

		JButton btnHelpReplicates = newHelpButton();
		GridBagConstraints gbc_btnNewButton14 = new GridBagConstraints();
		gbc_btnNewButton14.gridx = 3;
		gbc_btnNewButton14.gridy = 3;
		panel2.add(btnHelpReplicates, gbc_btnNewButton14);
		btnHelpReplicates.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				showHelp(REPS_PER_TREE_HELP);
			}
		});

		comboBox = new JComboBox(new String[]{"Best fit","Exponential","Gamma","Log Normal", "Exp Gamma"});
		GridBagConstraints gbc_comboBox = new GridBagConstraints();
		gbc_comboBox.insets = new Insets(0, 0, 0, 5);
		gbc_comboBox.fill = GridBagConstraints.HORIZONTAL;
		gbc_comboBox.gridx = 1;
		gbc_comboBox.gridy = 5;
		panel2.add(comboBox, gbc_comboBox);
		
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
//		            try {
//		            	if (m_distr != null) {
//		            		m_distr.initAndValidate();
//		            	}
//		            } catch (Exception e1) {
//		                // ignore
//		            }
		            int nPoints = ages.length;
		            int[] xPoints = new int[nPoints];
		            int[] yPoints = new int[nPoints];
		            int[] xPoints2 = new int[nPoints];
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
		                xPoints2[i] = graphoffset + nGraphWidth * i / nPoints;  
		                xPoints[i] = (int)(graphoffset + nGraphWidth * (ages[ages.length - 1 - i] - f2)/fXRange);  

// XXX
//for (int ccc = 0; ccc < probabilities.length; ccc++) {
//	System.out.println("Age: " + ages[ccc] + "\tProbability: " + probabilities[ccc]);
//}

// XXX
//System.out.println("Normaliser from CAPanel: " + probs.getNormaliser());

		                fyPoints[i] = probabilities[nPoints - i - 1]/probs.getNormaliser();

// XXX
//for (int ccc = 0; ccc < probabilities.length; ccc++) {
//	System.out.println("FYPoints: " + fyPoints[ccc]);
//}

		                if (m_distr != null) {
		                    try {
		                        fyPoints2[i] = m_distr.density(fMinValue + (fXRange * i) / nPoints);
		                        if (Double.isInfinite(fyPoints2[i]) || Double.isNaN(fyPoints2[i])) {
		                        	fyPoints2[i] = 0;
		                        }
		                    } catch (Exception e) {
		                        fyPoints2[i] = 0;
		                    }
		                }

 // XXX
//for (int ccc = 0; ccc < probabilities.length; ccc++) {
//	System.out.println("FYPoints2: " + fyPoints2[ccc]);
//}

		                
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
		            	g.drawPolyline(xPoints2, yPoints2, nPoints);
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
		            
		            // draw statistics
		            if (m_distr != null) {
			            int statoffsetx = getWidth() - 120;
			            int statoffsety = graphoffset + 10;
		            	String text = getFitParameters();
		            	String [] strs = text.split("\n");
		            	for (int i = 0; i < strs.length; i++) {
			                g.drawString(strs[i], statoffsetx, statoffsety + i * 10);
		            	}
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
	        

		};
		panel_1.setToolTipText("Click to show parameters of fit");
		panel_1.addMouseListener(new MouseListener() {
			
			@Override
			public void mouseReleased(MouseEvent e) {}
			
			@Override
			public void mousePressed(MouseEvent e) {}
			
			@Override
			public void mouseExited(MouseEvent e) {}
			
			@Override
			public void mouseEntered(MouseEvent e) {}
			
			@Override
			public void mouseClicked(MouseEvent e) {
				String text = getFitParameters();
				text = text.replaceAll("\n", "<br/>");
				showHelp("<html>Distribution fit:<br/>" + text + "</html>");
			}
		});
		
		panel_1.setBorder(new LineBorder(new Color(0, 0, 0)));
		panel_1.setBackground(Color.gray);
		panel_1.setPreferredSize(new Dimension(1024,400));
		GridBagConstraints gbc_panel_1 = new GridBagConstraints();
		gbc_panel_1.fill = GridBagConstraints.BOTH;
		gbc_panel_1.gridx = 0;
		gbc_panel_1.gridy = 1;
		panel3.add(panel_1, gbc_panel_1);

//		JPanel panel4 = new JPanel();
//		panel4.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Approximation", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(51, 51, 51)));
//		GridBagConstraints gbc_panel4 = new GridBagConstraints();
//		gbc_panel4.fill = GridBagConstraints.BOTH;
//		gbc_panel4.gridx = 0;
//		gbc_panel4.gridy = 2;
//		gbc_panel4.gridwidth = 3;
//		add(panel4, gbc_panel4);
//		GridBagLayout gbl_panel4 = new GridBagLayout();
//		gbl_panel4.columnWidths = new int[]{0, 0, 0, 0};
//		gbl_panel4.rowHeights = new int[]{0, 0, 0};
//		gbl_panel4.columnWeights = new double[]{0.0, 0.0, 0.0, 0.0};
//		gbl_panel4.rowWeights = new double[]{Double.MIN_VALUE, 0.0, 0.0};
//		panel4.setLayout(gbl_panel4);
//		
//		JLabel lblDistributionType = new JLabel("Distribution type:");
//		GridBagConstraints gbc_lblDistributionType = new GridBagConstraints();
//		gbc_lblDistributionType.insets = new Insets(0, 0, 0, 5);
//		gbc_lblDistributionType.anchor = GridBagConstraints.EAST;
//		gbc_lblDistributionType.gridx = 0;
//		gbc_lblDistributionType.gridy = 2;
//		panel4.add(lblDistributionType, gbc_lblDistributionType);
//		
//		
//		btnFindApproximation = new JButton("Find approximation");
//		btnFindApproximation.setEnabled(false);
//		btnFindApproximation.addActionListener(new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//				String type = (String) comboBox.getSelectedItem();
//				if (type.equals("Best fit")) {
//					m_distr = probs.fitExponential();
//					m_rmsd = probs.getApprox_distribution_rmsd();
////					normalise();
//					ContinuousDistribution tmp = probs.fitGamma();
//					if (probs.getApprox_distribution_rmsd() < m_rmsd) {
//						m_rmsd = probs.getApprox_distribution_rmsd();
//						m_distr = tmp;
////						normalise();
//					}
//					tmp = probs.fitLognormal();
//					if (probs.getApprox_distribution_rmsd() < m_rmsd) {
//						m_rmsd = probs.getApprox_distribution_rmsd();
//						m_distr = tmp;
////						normalise();
//					}
//					tmp = probs.fitExpGamma();
//					if (probs.getApprox_distribution_rmsd() < m_rmsd) {
//						m_rmsd = probs.getApprox_distribution_rmsd();
//						m_distr = tmp;
////						normalise();
//					}
//				}
//				if (type.equals("Exponential")) {
//					m_distr = probs.fitExponential();
//					m_rmsd = probs.getApprox_distribution_rmsd();
////					normalise();
//				}
//				if (type.equals("Gamma")) {
//					m_distr = probs.fitGamma();
//					m_rmsd = probs.getApprox_distribution_rmsd();
////					normalise();
//				}
//				if (type.equals("Log Normal")) {
//					m_distr = probs.fitLognormal();
//					m_rmsd = probs.getApprox_distribution_rmsd();
////					normalise();
//				}
//				if (type.equals("Exp Gamma")) {
//					m_distr = probs.fitExpGamma();
//					m_rmsd = probs.getApprox_distribution_rmsd();
////					normalise();
//				}
//				panel_1.repaint();
//			}
//
//		});
//		GridBagConstraints gbc_btnFindApproximation = new GridBagConstraints();
//		gbc_btnFindApproximation.insets = new Insets(0, 0, 0, 5);
//		gbc_btnFindApproximation.gridx = 2;
//		gbc_btnFindApproximation.gridy = 2;
//		panel4.add(btnFindApproximation, gbc_btnFindApproximation);
//		
//		Component horizontalGlue = Box.createHorizontalGlue();
//		GridBagConstraints gbc_horizontalGlue = new GridBagConstraints();
//		gbc_horizontalGlue.gridx = 3;
//		gbc_horizontalGlue.gridy = 2;
//		panel4.add(horizontalGlue, gbc_horizontalGlue);
		
		dataToGUI();

		lblBirt.setToolTipText(OCCURRENCE_AGE_HELP);
		lblNetDiversificationRate.setToolTipText(DIV_RATE_HELP);
		lblRurnoverRateDb.setToolTipText(TURNOVER_RATE_HELP);
		lblSamplingRate.setToolTipText(SAMPLING_RATE_HELP);
		lblNewLabel.setToolTipText(SAMPLING_GAP_HELP);
		lblNewLabel_1.setToolTipText(NR_SIMULATIONS_HELP);
		lblMaximumNumberOf.setToolTipText(MAX_NR_TREES_HELP);
		lblSamplingReplicatesPer.setToolTipText(REPS_PER_TREE_HELP);
		textField_maxOccuranceAge.setToolTipText(OCCURRENCE_AGE_HELP);
		textField_maxDivRate.setToolTipText(DIV_RATE_HELP);
		textField_maxTurnoverRate.setToolTipText(TURNOVER_RATE_HELP);
		textField_maxSamplingRate.setToolTipText(SAMPLING_RATE_HELP);
		textField_maxSamplingGap.setToolTipText(SAMPLING_GAP_HELP);
		textField_minOccuranceAge.setToolTipText(OCCURRENCE_AGE_HELP);
		textField_minDivRate.setToolTipText(DIV_RATE_HELP);
		textField_minTurnoverRate.setToolTipText(TURNOVER_RATE_HELP);
		textField_minSamplingRate.setToolTipText(SAMPLING_RATE_HELP);
		textField_minSamplginGap.setToolTipText(SAMPLING_GAP_HELP);
		textField_NumberOfTreeSimulations.setToolTipText(NR_SIMULATIONS_HELP);
		textField_MaxNrOfBranches.setToolTipText(MAX_NR_TREES_HELP);
		textField_SamplingReplicatesPerTree.setToolTipText(REPS_PER_TREE_HELP);
	}

	private JButton newHelpButton() {
		return new HelpButton("?",true);
	}

	
	public void calcFit() {
		String type = (String) comboBox.getSelectedItem();
		if (type.equals("Best fit")) {
			m_distr = probs.fitExponential();
			m_rmsd = probs.getApprox_distribution_rmsd();
//			normalise();
			ContinuousDistribution tmp = probs.fitGamma();
			if (probs.getApprox_distribution_rmsd() < m_rmsd) {
				m_rmsd = probs.getApprox_distribution_rmsd();
				m_distr = tmp;
//				normalise();
			}
			tmp = probs.fitLognormal();
			if (probs.getApprox_distribution_rmsd() < m_rmsd) {
				m_rmsd = probs.getApprox_distribution_rmsd();
				m_distr = tmp;
//				normalise();
			}
			tmp = probs.fitExpGamma();
			if (probs.getApprox_distribution_rmsd() < m_rmsd) {
				m_rmsd = probs.getApprox_distribution_rmsd();
				m_distr = tmp;
//				normalise();
			}
		}
		if (type.equals("Exponential")) {
			m_distr = probs.fitExponential();
			m_rmsd = probs.getApprox_distribution_rmsd();
//			normalise();
		}
		if (type.equals("Gamma")) {
			m_distr = probs.fitGamma();
			m_rmsd = probs.getApprox_distribution_rmsd();
//			normalise();
		}
		if (type.equals("Log Normal")) {
			m_distr = probs.fitLognormal();
			m_rmsd = probs.getApprox_distribution_rmsd();
//			normalise();
		}
		if (type.equals("Exp Gamma")) {
			m_distr = probs.fitExpGamma();
			m_rmsd = probs.getApprox_distribution_rmsd();
//			normalise();
		}
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

	static void showHelp(String text) {
//		JTextPane k = new JTextPane();
//		k.setContentType("text/html");
//		k.setText(text);
		
		   final JEditorPane pane = new JEditorPane();
		   pane.setContentType("text/html");
		   pane.setText(text);
		    pane.setEditable(false);
		    ToolTipManager.sharedInstance().registerComponent(pane);

//		    HyperlinkListener l = new HyperlinkListener() {
//		        @Override
//		        public void hyperlinkUpdate(HyperlinkEvent e) {
//		            if (HyperlinkEvent.EventType.ACTIVATED == e.getEventType()) {
//		                try {
//		                    pane.setPage(e.getURL());
//		                } catch (Exception e1) {
//		                    e1.printStackTrace();
//		                }
//		            }
//
//		        }
//
//		    };
//		    pane.addHyperlinkListener(l);
		    String title = text.substring(0, text.indexOf(":"));
		    title = title.replaceAll("<html>", "");
			//JOptionPane.showMessageDialog(null, pane, title, JOptionPane.PLAIN_MESSAGE);
			JOptionPane.showMessageDialog(null, pane, title, JOptionPane.PLAIN_MESSAGE, BeautiPanel.getIcon(CA_ICON));
		
	}

    private String getFitParameters() {
    	if (m_distr == null) {
    		return "No model fit yet";
    	}
    	String text = "";
        String distr = m_distr.getClass().getName();
        distr = distr.substring(distr.lastIndexOf('.') + 1);
        distr = distr.replace("Impl", "");
        distr = distr.replaceAll("([A-Z])", " $1").trim();
        text += distr +"\n\n";
        text += "RMDS: " + format(m_rmsd,5) +"\n";
        if (m_distr instanceof ExponentialDistributionImpl) {
        	double mean = ((ExponentialDistributionImpl) m_distr).getMean();
            text += "mean: " + format(mean,5) + "\n";
        }
        if (m_distr instanceof LogNormalImpl) {
        	double mean = ((LogNormalImpl) m_distr).getMean();
        	double sigma = ((LogNormalImpl) m_distr).getSigma();
            text += "mean: " + format(mean,5) + "\n";
            text += "sigma: " + format(sigma,5) + "\n";
        }
        if (m_distr instanceof GammaDistributionImpl) {
        	double alpha = ((GammaDistributionImpl) m_distr).getAlpha();
        	double beta = ((GammaDistributionImpl) m_distr).getBeta();
            text += "alpha: " + format(alpha,5) + "\n";
            text += "beta: " + format(beta,5) + "\n";
        }
        if (m_distr instanceof ExpGamma) {
        	double mean = ((ExpGamma) m_distr).getMean();
        	double alpha = ((ExpGamma) m_distr).getAlpha();
        	double beta = ((ExpGamma) m_distr).getBeta();
        	double weight = ((ExpGamma) m_distr).getWeight();
            text += "mean: " + format(mean,5) + "\n";
            text += "alpha: " + format(alpha,5) + "\n";
            text += "beta: " + format(beta,5) + "\n";
            text += "weight: " + format(weight,5) + "\n";
        }
        try {
	        text += "\nmedian: " + format(m_distr.inverseCumulativeProbability(0.5), 5) + "\n";
	        text += "\n95% HPD: " + format(m_distr.inverseCumulativeProbability(0.025), 5) + "\n" +
	        " to " + format(m_distr.inverseCumulativeProbability(0.975)) + "\n";
        } catch (Exception e) {
        	// ignore
        }
        return text;
	}

    private String format(double value) {
        return format(value,3);
    }

    private String format(double value, int digits) {
        StringWriter writer = new StringWriter();
        PrintWriter pw = new PrintWriter(writer);
        pw.printf("%." + digits +"g", value);
        pw.flush();
        return writer.toString();
    }

    public static void main(String[] args) {
		JFrame frame = new JFrame();
        Utils.loadUIManager();

        if (Utils.isMac()) {
            // set up application about-menu for Mac
            // Mac-only stuff
        	try {
            URL url = ClassLoader.getSystemResource(ModelBuilder.ICONPATH + "beauti.png");
            Icon icon = null;
            if (url != null) {
                icon = new ImageIcon(url);
            } else {
                System.err.println("Unable to find image: " + ModelBuilder.ICONPATH + "beauti.png");
            }
            jam.framework.Application application = new jam.framework.MultiDocApplication(null, "CladeAge", "about" , icon) {

                @Override
                protected JFrame getDefaultFrame() {
                    return null;
                }

                @Override
                public void doQuit() {
                    System.exit(0);
                }

                @Override
                public void doAbout() {
                    showHelp(ABOUT_HELP);
                }

				@Override
                public DocumentFrame doOpenFile(File file) {
                    return null;
                }

                @Override
                public DocumentFrame doNew() {
                    return null;
                }
            };
            jam.mac.Utils.macOSXRegistration(application);
        	} catch (Exception e) {
        		// ignore
        	}
            try {
            	Class<?> class_ = Class.forName("jam.maconly.OSXAdapter");
                Method method = class_.getMethod("enablePrefs", boolean.class);
                method.invoke(null, false);
            } catch (java.lang.Exception e) {
            	// ignore
            }
        }
        
		frame.setSize(1024, 728);
		

        ImageIcon icon = BeautiPanel.getIcon(CA_ICON);
        if (icon != null) {
            frame.setIconImage(icon.getImage());
        }
		CAPanel pane = new CAPanel();
        JMenuBar menuBar = pane.makeMenuBar();
        frame.setJMenuBar(menuBar);
		frame.getContentPane().add(pane);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);	
		frame.setVisible(true);
	}
    
    boolean bAdvancedFeatures = false;

	private JMenuBar makeMenuBar() {
        JMenuBar menuBar = new JMenuBar();
        JMenu fileMenu = new JMenu("File");
        fileMenu.setMnemonic('F');
        menuBar.add(fileMenu);
        fileMenu.add(new MyAction("New", "Start new Clade Age", "new", KeyEvent.VK_N) {
        	public void actionPerformed(ActionEvent ae) {
                main(new String[0]);
            }
        });

        if (!Utils.isMac()) {
            fileMenu.addSeparator();
            fileMenu.add(new MyAction("Close", "Close Window", "close", KeyEvent.VK_W) {
			    public void actionPerformed(ActionEvent ae) {
			        JMenuItem menuItem = (JMenuItem) ae.getSource();
			        JPopupMenu popupMenu = (JPopupMenu) menuItem.getParent();
			        Component invoker = popupMenu.getInvoker();
			        JComponent invokerAsJComponent = (JComponent) invoker;
			        Container topLevel = invokerAsJComponent.getTopLevelAncestor();
			        if (topLevel != null) {
			            ((JFrame) topLevel).dispose();
			        }
			    }
			});
            
            fileMenu.add(new MyAction("Exit", "Exit Program", "exit", KeyEvent.VK_F4) {
            	public void actionPerformed(ActionEvent ae) {
            		System.exit(0);
            	}
            });
        }
        
        JMenu modeMenu = new JMenu("Mode");
        menuBar.add(modeMenu);
        modeMenu.setMnemonic('M');

        JCheckBoxMenuItem advancedMode = new JCheckBoxMenuItem("Show advanced settings", bAdvancedFeatures);
        advancedMode.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
            	JCheckBoxMenuItem advancedMode = (JCheckBoxMenuItem) ae.getSource();
                bAdvancedFeatures = advancedMode.getState();
//                panel2.setVisible(bAdvancedFeatures);
                if (bAdvancedFeatures) {
                	gridBagLayout.columnWidths = new int[]{525, 375, 100};
            		add(panel2, gbc_panel2);
                	
                } else {
                    gridBagLayout.columnWidths = new int[]{800, 0, 100};
            		remove(panel2);
                	gridBagLayout.removeLayoutComponent(panel2);
//                    panel2.setMinimumSize(new Dimension(0,0));
//                    panel2.setSize(new Dimension(0,0));
//                    panel2.setMaximumSize(new Dimension(0,0));
//                    gridBagLayout.removeLayoutComponent(panel2);
                }
                
		        JMenuItem menuItem = (JMenuItem) ae.getSource();
		        JPopupMenu popupMenu = (JPopupMenu) menuItem.getParent();
		        Component invoker = popupMenu.getInvoker();
		        JComponent invokerAsJComponent = (JComponent) invoker;
		        Container topLevel = invokerAsJComponent.getTopLevelAncestor();
		        if (topLevel != null) {
		        	JFrame frame = (JFrame) topLevel;
		        	Dimension size = frame.getSize(); 
		            ((JFrame) topLevel).setSize(new Dimension(size.width, size.height - 1));
		            ((JFrame) topLevel).setSize(size);
		        }

                repaint();
            }
        });
        modeMenu.add(advancedMode);
        
        if (!Utils.isMac()) {
            JMenu helpMenu = new JMenu("Help");
            helpMenu.setMnemonic('H');
            menuBar.add(helpMenu);
            helpMenu.add(new MyAction("About", "Help about", "help", -1) {
		        public void actionPerformed(ActionEvent ae) {
		        	showHelp(ABOUT_HELP);
		        }
            });
        }

        
		return menuBar;
	}
	
}
