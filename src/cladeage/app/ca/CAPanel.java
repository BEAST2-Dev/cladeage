package cladeage.app.ca;





import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;

import beastfx.app.inputeditor.SmallButton;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.Cursor;
import javafx.scene.Node;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Control;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.image.ImageView;
import javafx.scene.layout.Background;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import cladeage.math.distributions.CladeAgeDistribution;
import cladeage.math.distributions.FossilCalibration.CladeAgeMethod;


public class CAPanel extends VBox {

	static final String CA_ICON = "cladeage/app/ca/icons/cladeage_256x256px.png";
	static final String CA_ICON2 = "cladeage/app/ca/icons/cladeage_128x128px.png";

	final public static String OCCURRENCE_AGE_HELP = "<html>First occurrence age:<br/>"+
		"<br/>"+
		"The age of the oldest fossil of the clade. <br/>"+
		"If this age is known exactly, it should be specified in the 'Minimum' <br/>"+
		"field. If not, you should specify both a minimum and maximum age.<br/>" +
		"<br/>" +
		"Consider a fosil described in the literature to have come from the oligocene,<br/>" +
		"which extends from about 34 million to 23 million years before the present <br/>" +
		"(33.9-0.1 to 23.03-0.05 Ma). The uncertainties in the boundaries are small <br/>" +
		"enough to be ignored, so the minimum is 23.03 and maximum 33.9 if you want <br/>" +
		"to express ages in millions of years. <br/>" +
		"</html>";

	final public static String YOUNGER_OCCURRENCE_AGE_HELP = "<html>Younger potential first occurrence age:<br/>"+
		"<br/>"+
		"In cases where it is uncertain which of two occurrences represents the first <br/>"+
		"occurrence of a clade, this is the age of the younger of the two occurrences. <br/>"+
		"</html>";

	final public static String WEIGHT_YOUNGER_OCCURRENCE_HELP = "<html>Weight of the younger potential first occurrence:<br/>"+
		"<br/>"+
		"In cases where it is uncertain which of two occurrences represents the first <br/>"+
		"occurrence of a clade, this is the relative probability for the younger of <br/>"+
		"the two occurrences. <br/>"+
		"</html>";

	final public static String WEIGHT_OLDER_OCCURRENCE_HELP = "<html>Weight of the older potential first occurrence:<br/>"+
		"<br/>"+
		"In cases where it is uncertain which of two occurrences represents the first <br/>"+
		"occurrence of a clade, this is the relative probability for the older of <br/>"+
		"the two occurrences. <br/>"+
		"</html>";

	final public static String OLDER_OCCURRENCE_AGE_HELP = "<html>Older potential first occurrence age:<br/>"+
		"<br/>"+
		"In cases where it is uncertain which of two occurrences represents the first <br/>"+
		"occurrence of a clade, this is the age of the older of the two occurrences. <br/>"+
		"</html>";

	final static String RATE_HELP = 
		"See e.g. Alfaro et al. (2009), Santini et al. (2009), Jetz et al. <br/>" +
		"(2012), and Stadler (2011) for rate estimates for vertebrates, <br/>" +
		"teleost fishes, birds, and mammals.<br/>" +
		"<br/>---<br/>" +
		"Alfaro ME, Santini F, Brock CD et al. (2009) Nine exceptional <br/>" +
		"radiations plus high turnover explain species diversity in jawed <br/>" +
		"vertebrates. Proceedings of the National Academy of Sciences USA, <br/>" +
		"106, 13410-13414.<br/>" +
		"<br/>" +
		"Jetz W, Thomas GH, Joy JB, Hartmann K, Mooers A (2012) The global <br/>" +
		"diversity of birds in space and time. Nature, 491, 444-448.<br/>" +
		"<br/>" +
		"Santini F, Harmon LJ, Carnevale G, Alfaro ME (2009) Did genome <br/>" +
		"duplication drive the origin of teleosts? A comparative study of <br/>" +
		"diversification in ray-finned fishes. BMC Evolutionary Biology, <br/>" +
		"9, 194.<br/>" +
		"<br/>" +
		"Stadler T (2011) Mammalian phylogeny reveals recent diversification <br/>" +
		"rate shifts. Proceedings of the National Academy of Sciences, <br/>" +
		"108, 6187-6192.<br/>";
	
	final public static String DIV_RATE_HELP = "<html>Net diversification rate:<br/>"+
		"<br/>"+
		"The net diversification rate is the difference between <br/>"+
		"speciation and extinction rate.<br/>" + RATE_HELP + "</html>";

	final public static String TURNOVER_RATE_HELP = "<html>Turnover rate:<br/>"+
		"<br/>"+
		"The turnover rate is the ratio of extinction and speciation rate.<br/>"+
		RATE_HELP + "</html>";

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
		"<a href='http://www.ncbi.nlm.nih.gov/pubmed/20739322'>http://www.ncbi.nlm.nih.gov/pubmed/20739322</a></html>";

	final public static String ABOUT_HELP = "<html>CladeAge:<br/><br/>" +
		"Copyright 2013<br/><br/>" +
		"Michael Matschiner<br/>michaelmatschiner@mac.com<br/>" +
		"<br/>and <br/><br/>" +
		"Remco Bouckaert<br/>remco@cs.auckland.ac.nz<br/>" +
		"</html>";
	
    // GUI components
	private TextField textField_maxOccuranceAge;
	private TextField textField_maxDivRate;
	private TextField textField_maxTurnoverRate;
	private TextField textField_maxSamplingRate;
	private TextField textField_minOccuranceAge;
	private TextField textField_minDivRate;
	private TextField textField_minTurnoverRate;
	private TextField textField_minSamplingRate;
	Button btnFindApproximation;
	Button btnCalculate;
	public void setCalculateButtonText(String text) {btnCalculate.setText(text);}
	ComboBox comboBox;
//	public JComboBox getComboBox() {return comboBox;}
	GridPane panel;
	//JPanel panel2;
	GridPane panel_1;
	
	CladeAgeProbabilities probs = null;
	
	private double minOccuranceAge = 0;
	private double minDivRate = 0.01;
	private double minTurnoverRate = 0.1;
	private double minSamplingRate = 0.01;

	private double maxOccuranceAge = minOccuranceAge;
	private double maxDivRate = minDivRate;
	private double maxTurnoverRate = minTurnoverRate;
	private double maxSamplingRate = minSamplingRate;

	private CladeAgeMethod cladeAgeMethod = CladeAgeMethod.standard;
	
	public double getMinOccuranceAge() { return	minOccuranceAge;}
	public double getMinDivRate() { return	minDivRate;}
	public double getMinTurnoverRate() { return	minTurnoverRate;}
	public double getMinSamplingRate() { return	minSamplingRate;}

	public double getMaxOccuranceAge() { return	maxOccuranceAge;}
	public double getMaxDivRate() { return	maxDivRate;}
	public double getMaxTurnoverRate() { return	maxTurnoverRate;}
	public double getMaxSamplingRate() { return	maxSamplingRate;}
	
	public void setMinOccuranceAge(double minOccuranceAge) {this.minOccuranceAge = minOccuranceAge;}
	public void setMinDivRate(double minDivRate) {this.minDivRate = minDivRate;}
	public void setMinTurnoverRate(double minTurnoverRate) {this.minTurnoverRate = minTurnoverRate;}
	public void setMinSamplingRate(double minSamplingRate) {this.minSamplingRate = minSamplingRate;}
	   	
	public void setMaxOccuranceAge(double maxOccuranceAge) {this.maxOccuranceAge = maxOccuranceAge;}
	public void setMaxDivRate(double maxDivRate) {this.maxDivRate = maxDivRate;}
	public void setMaxTurnoverRate(double maxTurnoverRate) {this.maxTurnoverRate = maxTurnoverRate;}
	public void setMaxSamplingRate(double maxSamplingRate) {this.maxSamplingRate = maxSamplingRate;}
	
	private LineChart chart;
	LineChart.Series series;
	
	public void setMethod(CladeAgeMethod method) {
		cladeAgeMethod = method;
		if (comboBox != null) {
			comboBox.setValue(method);
		}
	}
	public CladeAgeMethod getMethod() {
		if (comboBox != null) {
			cladeAgeMethod = (CladeAgeMethod) comboBox.getValue();
		}
		return cladeAgeMethod;
	}
	
	boolean processingDataToGui = false;

	// used for graphing
	double [] ages;
	double [] probabilities;
	
	ContinuousDistribution m_distr = null;
	double m_rmsd = 0;
	
	
	public final static int MODE_STAND_ALONE = 0;
	public final static int MODE_BEAUTI_TOP = 1;
	public final static int MODE_BEAUTI_BOTTOM = 2;
	int mode = MODE_STAND_ALONE;

	
	public CAPanel(int mode) {
		this.mode = mode;
//		switch (mode) {
//		case MODE_STAND_ALONE:
////			gridBagLayout.columnWidths = new int[]{800, 100, 200};
////			gridBagLayout.columnWidths = new int[]{450, 80, 200};
////			gridBagLayout.columnWidths = new int[]{350, 40, 150};
////			gridBagLayout.columnWidths = new int[]{0};
////			gridBagLayout.rowHeights = new int[]{180, 450};
////			gridBagLayout.columnWeights = new double[]{1.0, 0.0, 1.0};
////			gridBagLayout.rowWeights = new double[]{1.0, 0.0, 1.0};
//			break;
//		case MODE_BEAUTI_BOTTOM:
////			gridBagLayout.columnWidths = new int[]{0};
////			gridBagLayout.rowHeights = new int[]{80, 450};
////			gridBagLayout.columnWeights = new double[]{1.0, 0.0, 1.0};
////			gridBagLayout.rowWeights = new double[]{1.0, 0.0, 1.0};
//			break;
//		case MODE_BEAUTI_TOP:
////			gridBagLayout.columnWidths = new int[]{0};
////			gridBagLayout.rowHeights = new int[]{180};
////			gridBagLayout.columnWeights = new double[]{1.0};
////			gridBagLayout.rowWeights = new double[]{1.0};
//			break;			
//		}
////		setLayout(gridBagLayout);
		
		panel = new GridPane();
		panel.setHgap(5);
		panel.setVgap(5);
		panel.setPadding(new Insets(0,0,5,0));
		//panel.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Model parameters", TitledBorder.LEADING, TitledBorder.TOP, null, null));
//		//gbc_panel.fill = GridBagConstraints.BOTH;
//		//gbc_panel.gridx = 0;
//		//gbc_panel.gridy = 0;
		getChildren().add(panel); //gbc_panel);
//		gbl_panel.columnWidths = new int[]{0, 0, 0, 0, 0, 0};
//		gbl_panel.rowHeights = new int[]{0, 0, 0, 0, 0, 0, 0};
//		gbl_panel.columnWeights = new double[]{0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
//		gbl_panel.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//		panel.setLayout(gbl_panel);

		Label lblNewLabel_2 = new Label("Minimum");
//		GridBagConstraints //gbc_lblNewLabel_2 = new GridBagConstraints();
//		//gbc_lblNewLabel_2.insets = new Insets(0, 0, 5, 5);
//		//gbc_lblNewLabel_2.gridx = 1;
//		//gbc_lblNewLabel_2.gridy = 0;
		panel.add(lblNewLabel_2, 1, 0, 1, 1); //gbc_lblNewLabel_2);
		
		Label lblNewLabel_3 = new Label("Maximum");
		//GridBagConstraints //gbc_lblNewLabel_3 = new GridBagConstraints();
		//gbc_lblNewLabel_3.insets = new Insets(0, 0, 5, 5);
		//gbc_lblNewLabel_3.gridx = 3;
		//gbc_lblNewLabel_3.gridy = 0;
		panel.add(lblNewLabel_3, 3, 0, 3, 1); //gbc_lblNewLabel_3);
		
		
		
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_BOTTOM) {
			Label lblBirt = new Label("First occurance age:");
			lblBirt.setTooltip(new Tooltip(OCCURRENCE_AGE_HELP));
			//GridBagConstraints //gbc_lblBirt = new GridBagConstraints();
			//gbc_lblBirt.anchor = GridBagConstraints.EAST;
			//gbc_lblBirt.insets = new Insets(0, 0, 5, 5);
			//gbc_lblBirt.gridx = 0;
			//gbc_lblBirt.gridy = 1;
			panel.add(lblBirt, 0, 1, 1, 1); //gbc_lblBirt);

			textField_minOccuranceAge = newTextField();
			textField_minOccuranceAge.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_7 = new GridBagConstraints();
			//gbc_textField_7.insets = new Insets(0, 0, 5, 5);
			//gbc_textField_7.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_7.gridx = 1;
			//gbc_textField_7.gridy = 1;
			panel.add(textField_minOccuranceAge, 1, 1, 1, 1); //gbc_textField_7);
			
			Label label = new Label("-");
			//GridBagConstraints //gbc_label = new GridBagConstraints();
			//gbc_label.insets = new Insets(0, 0, 5, 5);
			//gbc_label.anchor = GridBagConstraints.EAST;
			//gbc_label.gridx = 2;
			//gbc_label.gridy = 1;
			panel.add(label, 2 ,1 ,1, 1); //gbc_label);
			
			textField_maxOccuranceAge = newTextField();
			textField_maxOccuranceAge.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_1_1 = new GridBagConstraints();
			//gbc_textField_1_1.insets = new Insets(0, 0, 5, 0);
			//gbc_textField_1_1.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_1_1.gridx = 3;
			//gbc_textField_1_1.gridy = 1;
			panel.add(textField_maxOccuranceAge, 3, 1, 1, 1); //gbc_textField_1_1);
		}
		
		// always create the combobox, even when it is not displayed
		CladeAgeMethod[] values = CladeAgeMethod.values();
		comboBox = new ComboBox();
		for (CladeAgeMethod value : values) {
			comboBox.getItems().add(value);
		}
		comboBox.getSelectionModel().select(cladeAgeMethod);

		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_TOP) {
					//new String[]{"standard CladeAge"});
			//GridBagConstraints //gbc_comboBox = new GridBagConstraints();
			//gbc_comboBox.insets = new Insets(0, 0, 0, 5);
			//gbc_comboBox.fill = GridBagConstraints.HORIZONTAL;
			//gbc_comboBox.gridx = 0;
			//gbc_comboBox.gridy = 0;
			panel.add(comboBox, 0, 0, 1, 1);//gbc_comboBox);
			comboBox.setOnAction(e->guiToData());

			Label lblNetDiversificationRate = new Label("Net diversification rate \u03BB-\u03BC:");
			//GridBagConstraints //gbc_lblNetDiversificationRate = new GridBagConstraints();
			//gbc_lblNetDiversificationRate.anchor = GridBagConstraints.EAST;
			//gbc_lblNetDiversificationRate.insets = new Insets(0, 0, 5, 5);
			//gbc_lblNetDiversificationRate.gridx = 0;
			//gbc_lblNetDiversificationRate.gridy = 2;
			panel.add(lblNetDiversificationRate, 0, 2, 1, 1); //gbc_lblNetDiversificationRate);
			
			textField_minDivRate = newTextField();
			textField_minDivRate.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_8 = new GridBagConstraints();
			//gbc_textField_8.insets = new Insets(0, 0, 5, 5);
			//gbc_textField_8.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_8.gridx = 1;
			//gbc_textField_8.gridy = 2;
			panel.add(textField_minDivRate, 1, 2, 1, 1); //gbc_textField_8);
			
			Label label_1 = new Label("-");
			//GridBagConstraints //gbc_label_1 = new GridBagConstraints();
			//gbc_label_1.insets = new Insets(0, 0, 5, 5);
			//gbc_label_1.gridx = 2;
			//gbc_label_1.gridy = 2;
			panel.add(label_1, 2, 2, 1, 1); //gbc_label_1);
			
			textField_maxDivRate = newTextField();
			textField_maxDivRate.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_2 = new GridBagConstraints();
			//gbc_textField_2.insets = new Insets(0, 0, 5, 0);
			//gbc_textField_2.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_2.gridx = 3;
			//gbc_textField_2.gridy = 2;
			panel.add(textField_maxDivRate, 3, 2, 1, 1); //gbc_textField_2);
			
			Label lblRurnoverRateDb = new Label("Turnover rate \u03BC/\u03BB:");
			//GridBagConstraints //gbc_lblRurnoverRateDb = new GridBagConstraints();
			//gbc_lblRurnoverRateDb.anchor = GridBagConstraints.EAST;
			//gbc_lblRurnoverRateDb.insets = new Insets(0, 0, 5, 5);
			//gbc_lblRurnoverRateDb.gridx = 0;
			//gbc_lblRurnoverRateDb.gridy = 3;
			panel.add(lblRurnoverRateDb, 0, 3, 1, 1); //gbc_lblRurnoverRateDb);
			
			textField_minTurnoverRate = newTextField();
			textField_minTurnoverRate.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_9 = new GridBagConstraints();
			//gbc_textField_9.insets = new Insets(0, 0, 5, 5);
			//gbc_textField_9.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_9.gridx = 1;
			//gbc_textField_9.gridy = 3;
			panel.add(textField_minTurnoverRate, 1, 3, 1, 1); //gbc_textField_9);
			
			Label label_2 = new Label("-");
			//GridBagConstraints //gbc_label_2 = new GridBagConstraints();
			//gbc_label_2.insets = new Insets(0, 0, 5, 5);
			//gbc_label_2.gridx = 2;
			//gbc_label_2.gridy = 3;
			panel.add(label_2, 2, 3, 1, 1); //gbc_label_2);
			
			textField_maxTurnoverRate = newTextField();
			textField_maxTurnoverRate.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_3 = new GridBagConstraints();
			//gbc_textField_3.insets = new Insets(0, 0, 5, 0);
			//gbc_textField_3.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_3.gridx = 3;
			//gbc_textField_3.gridy = 3;
			panel.add(textField_maxTurnoverRate, 3, 3, 1, 1); //gbc_textField_3);
			
			Label lblSamplingRate = new Label("Sampling rate \u03a8:");
			//GridBagConstraints //gbc_lblSamplingRate = new GridBagConstraints();
			//gbc_lblSamplingRate.anchor = GridBagConstraints.EAST;
			//gbc_lblSamplingRate.insets = new Insets(0, 0, 5, 5);
			//gbc_lblSamplingRate.gridx = 0;
			//gbc_lblSamplingRate.gridy = 4;
			panel.add(lblSamplingRate, 0, 4, 1, 1); //gbc_lblSamplingRate);
			
			textField_minSamplingRate = newTextField();
			textField_minSamplingRate.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_10 = new GridBagConstraints();
			//gbc_textField_10.insets = new Insets(0, 0, 5, 5);
			//gbc_textField_10.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_10.gridx = 1;
			//gbc_textField_10.gridy = 4;
			panel.add(textField_minSamplingRate, 1, 4, 1, 1); //gbc_textField_10);
			
			Label label_3 = new Label("-");
			//GridBagConstraints //gbc_label_3 = new GridBagConstraints();
			//gbc_label_3.insets = new Insets(0, 0, 5, 5);
			//gbc_label_3.gridx = 2;
			//gbc_label_3.gridy = 4;
			panel.add(label_3, 2, 4, 1, 1); //gbc_label_3);
			
			textField_maxSamplingRate = newTextField();
			textField_maxSamplingRate.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_4 = new GridBagConstraints();
			//gbc_textField_4.insets = new Insets(0, 0, 5, 0);
			//gbc_textField_4.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_4.gridx = 3;
			//gbc_textField_4.gridy = 4;
			panel.add(textField_maxSamplingRate, 3, 4, 1, 1); //gbc_textField_4);
			
			setToolTipText(lblNetDiversificationRate, DIV_RATE_HELP);
			setToolTipText(lblRurnoverRateDb, TURNOVER_RATE_HELP);
			setToolTipText(lblSamplingRate, SAMPLING_RATE_HELP);
		}
		
				
		
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_BOTTOM) {
			// Label lblNewLabel = new Label("Sampling gap:");
			// lblNewLabel.setTooltip(new Tooltip(SAMPLING_GAP_HELP);
			// GridBagConstraints //gbc_lblNewLabel = new GridBagConstraints();
			// //gbc_lblNewLabel.insets = new Insets(0, 0, 0, 5);
			// //gbc_lblNewLabel.anchor = GridBagConstraints.EAST;
			// //gbc_lblNewLabel.gridx = 0;
			// //gbc_lblNewLabel.gridy = 5;
			// panel.add(lblNewLabel, //gbc_lblNewLabel);

			// textField_minSamplingGap = newTextField();
			// textField_minSamplingGap.addActionListener(new ActionListener() {
			// 	public void actionPerformed(ActionEvent e) {
			// 	}
			// });
			// textField_minSamplingGap.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_11 = new GridBagConstraints();
			//gbc_textField_11.insets = new Insets(0, 0, 0, 5);
			//gbc_textField_11.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_11.gridx = 1;
			//gbc_textField_11.gridy = 5;
			// panel.add(textField_minSamplingGap, //gbc_textField_11);
			
			//Label label_4 = new Label("-");
			//GridBagConstraints //gbc_label_4 = new GridBagConstraints();
			//gbc_label_4.insets = new Insets(0, 0, 0, 5);
			//gbc_label_4.gridx = 2;
			//gbc_label_4.gridy = 5;
			//panel.add(label_4, 2, 5, 1, 1); //gbc_label_4);
			
			// textField_maxSamplingGap = newTextField();
			// textField_maxSamplingGap.addActionListener(new ActionListener() {
			// 	public void actionPerformed(ActionEvent e) {
			// 	}
			// });
			// textField_maxSamplingGap.setPrefColumnCount(10);
			//GridBagConstraints //gbc_textField_5 = new GridBagConstraints();
			//gbc_textField_5.fill = GridBagConstraints.HORIZONTAL;
			//gbc_textField_5.gridx = 3;
			//gbc_textField_5.gridy = 5;
			// panel.add(textField_maxSamplingGap, //gbc_textField_5);
	
			// help buttons for panel1
			Button btnHelpButtonFirstOccurance = newHelpButton();
			//GridBagConstraints //gbc_btnNewButton = new GridBagConstraints();
			//gbc_btnNewButton.gridx = 5;
			//gbc_btnNewButton.gridy = 1;
			panel.add(btnHelpButtonFirstOccurance, 5, 1, 1, 1); //gbc_btnNewButton);
			btnHelpButtonFirstOccurance.setOnAction(e ->
					showHelp(OCCURRENCE_AGE_HELP)
			);

		}
		
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_TOP) {
			Button btnHelpButton2= newHelpButton();
			//GridBagConstraints //gbc_btnNewButton2 = new GridBagConstraints();
			//gbc_btnNewButton2.gridx = 5;
			//gbc_btnNewButton2.gridy = 2;
			panel.add(btnHelpButton2, 5, 2, 1, 1); //gbc_btnNewButton2);
			btnHelpButton2.setOnAction(e ->
					showHelp(DIV_RATE_HELP)
			);

			Button btnHelpButton3= newHelpButton();
			//GridBagConstraints //gbc_btnNewButton3 = new GridBagConstraints();
			//gbc_btnNewButton3.gridx = 5;
			//gbc_btnNewButton3.gridy = 3;
			panel.add(btnHelpButton3, 5, 3, 1, 1); //gbc_btnNewButton3);
			btnHelpButton3.setOnAction(e ->
					showHelp(TURNOVER_RATE_HELP)
			);
	
			Button btnHelpButton4= newHelpButton();
			//GridBagConstraints //gbc_btnNewButton4 = new GridBagConstraints();
			//gbc_btnNewButton4.gridx = 5;
			//gbc_btnNewButton4.gridy = 4;
			panel.add(btnHelpButton4, 5, 4, 1, 1); //gbc_btnNewButton4);
			btnHelpButton4.setOnAction(e ->
					showHelp(SAMPLING_RATE_HELP)
			);
		}
		
		// if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_BOTTOM) {
		// 	Button btnHelpButton5= newHelpButton();
		// 	GridBagConstraints //gbc_btnNewButton5 = new GridBagConstraints();
		// 	//gbc_btnNewButton5.gridx = 5;
		// 	//gbc_btnNewButton5.gridy = 5;
		// 	panel.add(btnHelpButton5, //gbc_btnNewButton5);
		// 	btnHelpButton5.addActionListener(new ActionListener() {
		// 		@Override
		// 		public void actionPerformed(ActionEvent e) {
		// 			showHelp(SAMPLING_GAP_HELP);
		// 		}
		// 	});
		// }
				
		HBox panel2b = FXUtils.newHBox();
		// panel2b.setLayout(new BorderLayout());
		
		ImageView icon = FXUtils.getIcon("CA", CA_ICON2);
		Label lblIcon = new Label();
		lblIcon.setGraphic(icon);
		if (mode == MODE_BEAUTI_BOTTOM) {
			lblIcon.setMinSize(4,4);
			lblIcon.setPrefSize(84,84);
			lblIcon.setMaxSize(4,4);
			panel.add(lblIcon, 6, 2, 1, 1);
		} else {
			lblIcon.setMinSize(160,160);
			lblIcon.setPrefSize(160,160);
			panel2b.getChildren().add(lblIcon);
		}
		
		btnCalculate = new Button("Run");
		btnCalculate.setMinSize(128, 20);
		btnCalculate.setOnAction(e -> {
			cursorProperty().set(Cursor.WAIT);
				guiToData();
				dataToGUI();
				ages = null;
				probs = new CladeAgeProbabilities();
				
				Dialog dlg = new Dialog<>();
				dlg.getDialogPane().setHeaderText("Progress Dialog");
			    final ProgressBar dpb = new ProgressBar();
			    dlg.getDialogPane().setContent(dpb);
			    dlg.getDialogPane().getButtonTypes().add(ButtonType.CANCEL);
			    
			    
			    final boolean [] finished = new boolean[1];
			    final Thread t = new Thread(new Runnable() {
			      public void run() {
//					try {
//						String method = comboBox.getSelectedItem().toString();
//						if (method.equals("standard CladeAge")) {
//							probs.run_standard_cladeage(
//									minOccuranceAge, maxOccuranceAge,
//									minDivRate, maxDivRate,
//									minTurnoverRate, maxTurnoverRate,
//									minSamplingRate, maxSamplingRate,
//									minSamplingGap, maxSamplingGap, dpb);
//							m_rmsd = 0.0;
//						}
//						if (!method.equals("standard CladeAge")) {
//							m_rmsd = probs.getDistribution_rmsd();
//						}
//
//					} catch (Exception e) {
//						e.printStackTrace();
//						throw new RuntimeException(e);
//					}
					calcFit(dpb);
					finished[0] = true;
					Platform.runLater(() ->dlg.close());
			      }
			    });
			    t.start();
		        Optional optional = dlg.showAndWait();
		        if (!finished[0] && optional.toString().toLowerCase().contains("cancel")) {
					try {
						System.err.println("Trying to stop");
						probs.setCancel();
						dlg.close();
					} catch (Exception ex) {
						ex.printStackTrace();
					}
		        }

			    if (probs.getCancel()) {
					refresh(); //panel_1.repaint();
					cursorProperty().set(Cursor.DEFAULT);
			    	return;
			    }
		        
				ages = probs.getAges();
				probabilities =  probs.getInt_probabilities();

				refresh(); //panel_1.repaint();
				cursorProperty().set(Cursor.DEFAULT);
			}
		);
//		GridBagConstraints //gbc_btnCalculate = new GridBagConstraints();
//		//gbc_btnCalculate.gridwidth = 2;
//		//gbc_btnCalculate.insets = new Insets(0, 0, 5, 0);
//		//gbc_btnCalculate.gridx = 0;
//		//gbc_btnCalculate.gridy = 4;
		if (mode == MODE_BEAUTI_BOTTOM) {
			panel.add(btnCalculate, 8, 1, 1, 1);
		} else {
			panel2b.getChildren().add(btnCalculate);
		}
		
		//GridBagConstraints gbc_lblIcon = new GridBagConstraints();
		//gbc_lblIcon.insets = new Insets(5, 0, 5, 5);
		//gbc_lblIcon.anchor = GridBagConstraints.EAST;
		//gbc_lblIcon.gridx = 2;
		//gbc_lblIcon.gridy = 0;
		//gbc_lblIcon.gridwidth = 2;
		
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_BOTTOM) {
			getChildren().add(panel2b); //gbc_lblIcon);
		}

// The below is not needed anymore, as parameters numberOfTreeSimulations, maxNrOfBranches, and samplingReplicatesPer don't exist anymore.
//
//		textField_NumberOfTreeSimulations = newTextField();
//		textField_NumberOfTreeSimulations.addActionListener(new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//			}
//		});
//		GridBagConstraints //gbc_textField = new GridBagConstraints();
//		//gbc_textField.insets = new Insets(0, 0, 5, 0);
//		//gbc_textField.fill = GridBagConstraints.HORIZONTAL;
//		//gbc_textField.gridx = 1;
//		//gbc_textField.gridy = 1;
//		panel2.add(textField_NumberOfTreeSimulations, //gbc_textField);
//		textField_NumberOfTreeSimulations.setPrefColumnCount(10);
//		
//		textField_MaxNrOfBranches = newTextField();
//		textField_MaxNrOfBranches.addActionListener(new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//			}
//		});
//		textField_MaxNrOfBranches.setPrefColumnCount(10);
//		GridBagConstraints //gbc_textField_12 = new GridBagConstraints();
//		//gbc_textField_12.insets = new Insets(0, 0, 5, 0);
//		//gbc_textField_12.fill = GridBagConstraints.HORIZONTAL;
//		//gbc_textField_12.gridx = 1;
//		//gbc_textField_12.gridy = 2;
//		panel2.add(textField_MaxNrOfBranches, //gbc_textField_12);
//		
//		Label lblSamplingReplicatesPer = new Label("Sampling replicates per tree:");
//		GridBagConstraints //gbc_lblSamplingReplicatesPer = new GridBagConstraints();
//		//gbc_lblSamplingReplicatesPer.insets = new Insets(0, 0, 5, 5);
//		//gbc_lblSamplingReplicatesPer.anchor = GridBagConstraints.EAST;
//		//gbc_lblSamplingReplicatesPer.gridx = 0;
//		//gbc_lblSamplingReplicatesPer.gridy = 3;
//		panel2.add(lblSamplingReplicatesPer, //gbc_lblSamplingReplicatesPer);
//		
//		textField_SamplingReplicatesPerTree = newTextField();
//		textField_SamplingReplicatesPerTree.addActionListener(new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//			}
//		});
//		textField_SamplingReplicatesPerTree.setPrefColumnCount(10);
//		GridBagConstraints //gbc_textField_13 = new GridBagConstraints();
//		//gbc_textField_13.insets = new Insets(0, 0, 5, 0);
//		//gbc_textField_13.fill = GridBagConstraints.HORIZONTAL;
//		//gbc_textField_13.gridx = 1;
//		//gbc_textField_13.gridy = 3;
//		panel2.add(textField_SamplingReplicatesPerTree, //gbc_textField_13);
//				
//		Label lblMaximumNumberOf = new Label("Maximum number of branches:");
//		GridBagConstraints //gbc_lblMaximumNumberOf = new GridBagConstraints();
//		//gbc_lblMaximumNumberOf.anchor = GridBagConstraints.EAST;
//		//gbc_lblMaximumNumberOf.insets = new Insets(0, 0, 5, 5);
//		//gbc_lblMaximumNumberOf.gridx = 0;
//		//gbc_lblMaximumNumberOf.gridy = 2;
//		panel2.add(lblMaximumNumberOf, //gbc_lblMaximumNumberOf);
//		
//		Label lblNewLabel_1 = new Label("Number of tree simulations:");
//		GridBagConstraints //gbc_lblNewLabel_1 = new GridBagConstraints();
//		//gbc_lblNewLabel_1.insets = new Insets(0, 0, 5, 5);
//		//gbc_lblNewLabel_1.anchor = GridBagConstraints.EAST;
//		//gbc_lblNewLabel_1.gridx = 0;
//		//gbc_lblNewLabel_1.gridy = 1;
//		panel2.add(lblNewLabel_1, //gbc_lblNewLabel_1);
//		
//		Component verticalGlue = Box.createVerticalGlue();
//		GridBagConstraints //gbc_verticalGlue = new GridBagConstraints();
//		//gbc_verticalGlue.gridwidth = 2;
//		//gbc_verticalGlue.insets = new Insets(0, 0, 5, 5);
//		//gbc_verticalGlue.gridx = 0;
//		//gbc_verticalGlue.gridy = 5;
//		panel2.add(verticalGlue, //gbc_verticalGlue);
//
//		// help buttons for panel2
//		Button btnHelpNrSimulateions = newHelpButton();
//		GridBagConstraints //gbc_btnNewButton12 = new GridBagConstraints();
//		//gbc_btnNewButton12.gridx = 3;
//		//gbc_btnNewButton12.gridy = 1;
//		panel2.add(btnHelpNrSimulateions, //gbc_btnNewButton12);
//		btnHelpNrSimulateions.addActionListener(new ActionListener() {
//			@Override
//			public void actionPerformed(ActionEvent e) {
//				showHelp(NR_SIMULATIONS_HELP);
//			}
//		});
//		
//		Button btnHelpNrBranches = newHelpButton();
//		GridBagConstraints //gbc_btnNewButton13 = new GridBagConstraints();
//		//gbc_btnNewButton13.gridx = 3;
//		//gbc_btnNewButton13.gridy = 2;
//		panel2.add(btnHelpNrBranches, //gbc_btnNewButton13);
//		btnHelpNrBranches.addActionListener(new ActionListener() {
//			@Override
//			public void actionPerformed(ActionEvent e) {
//				showHelp(MAX_NR_TREES_HELP);
//			}
//		});
//
//		Button btnHelpReplicates = newHelpButton();
//		GridBagConstraints //gbc_btnNewButton14 = new GridBagConstraints();
//		//gbc_btnNewButton14.gridx = 3;
//		//gbc_btnNewButton14.gridy = 3;
//		panel2.add(btnHelpReplicates, //gbc_btnNewButton14);
//		btnHelpReplicates.addActionListener(new ActionListener() {
//			@Override
//			public void actionPerformed(ActionEvent e) {
//				showHelp(REPS_PER_TREE_HELP);
//			}
//		});
		
		GridPane panel3 = new GridPane();
		//panel3.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Clade age probabilities", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		//GridBagConstraints //gbc_panel3 = new GridBagConstraints();
		//gbc_panel3.anchor = GridBagConstraints.WEST;
		//gbc_panel3.fill = GridBagConstraints.BOTH;
		//gbc_panel3.gridx = 0;
		//gbc_panel3.gridy = 1;
		//gbc_panel3.gridwidth = 3;
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_BOTTOM) {
			getChildren().add(panel3); //gbc_panel3);
		}
//		GridBagLayout gbl_panel3 = new GridBagLayout();
//		gbl_panel3.columnWidths = new int[]{0};
//		gbl_panel3.rowHeights = new int[]{0, 103};
//		gbl_panel3.columnWeights = new double[]{1.0};
//		gbl_panel3.rowWeights = new double[]{Double.MIN_VALUE, 1.0};
//		panel3.setLayout(gbl_panel3);
		
		
		NumberAxis xAxis = new NumberAxis();
		xAxis.setForceZeroInRange(false);
        //xAxis.setLabel("x");                
        NumberAxis yAxis = new NumberAxis();        
        yAxis.setLabel("p(x)");
        chart = new LineChart<Number,Number>(xAxis,yAxis);
        //chart.setAnimated(false);
        chart.setLegendVisible(false);
        chart.setCreateSymbols(false);
        chart.getXAxis().setAutoRanging(true);
        chart.getYAxis().setAutoRanging(true);
        series = new LineChart.Series<>();
        if (ages != null) {
        	for (int i = 0; i < ages.length; i++) {
        		series.getData().add(new XYChart.Data<Number,Number>(0,0));
        	}
        }
        chart.getData().add(series);

		panel_1 = new GridPane();// {
		panel3.add(panel2b, 0, 2, 1, 1);
		if (mode == MODE_BEAUTI_BOTTOM) {
			panel.add(chart, 0, 2, 6, 1);
			chart.setOnMouseClicked(e->{
				String text = getFitParameters();
				text = text.replaceAll("\n", "<br/>");
				showHelp("<html>Distribution fit:<br/>" + text + "</html>");
			}
		);
		} else {
			panel_1.add(chart, 0, 1, 1, 1);
		}
//			
//			
//			protected void paintComponent(java.awt.Graphics g) {
//				g.setColor(new Color(241, 241, 241));
//				g.fillRect(0, 0, getWidth(), getHeight());
//				if (ages != null) {
//		            final int width = getWidth();
//		            final int height = getHeight();
//		            final int graphoffset = 30;
//		            final int labeloffset = 0;
//		            int nGraphWidth = width - graphoffset * 2;
//		            int nGraphHeight = height - graphoffset * 2 - labeloffset;
//		            g.setColor(Color.WHITE);
//		            g.fillRect(graphoffset, graphoffset, nGraphWidth, nGraphHeight);
//		            g.setColor(Color.BLACK);
//		            g.drawRect(graphoffset, graphoffset, nGraphWidth, nGraphHeight);
////		            try {
////		            	if (m_distr != null) {
////		            		m_distr.initAndValidate();
////		            	}
////		            } catch (Exception e1) {
////		                // ignore
////		            }
//		            int nPoints = ages.length;
//		            int nPoints2 = getWidth();
//		            int[] xPoints = new int[nPoints];
//		            int[] yPoints = new int[nPoints];
//		            int[] xPoints2 = new int[nPoints2];
//		            int[] yPoints2 = new int[nPoints2];
//		            double[] fyPoints = new double[nPoints];
//		            double[] fyPoints2 = new double[nPoints2];
//		            Font font = g.getFont();
//		            double fMinValue = 0.1;
//		            double fMaxValue = 1;
//	                fMinValue = Math.min(ages[0], ages[ages.length-1]);
//	                fMaxValue = Math.max(ages[0], ages[ages.length-1]);
//		            double fXRange = fMaxValue - fMinValue;
//		            // adjust fYMax so that the ticks come out right
//		            double fX0 = fMinValue;
//		            int k = 0;
//		            double f = fXRange;
//		            double f2 = fX0;
//		            while (f > 10) {
//		                f /= 10;
//		                f2 /= 10;
//		                k++;
//		            }
//		            while (f < 1 && f > 0) {
//		                f *= 10;
//		                f2 *= 10;
//		                k--;
//		            }
//		            f = Math.ceil(f);
//		            f2 = Math.floor(f2);
////					final int NR_OF_TICKS_X = NR_OF_TICKS[(int) f];
//		            for (int i = 0; i < k; i++) {
//		                f *= 10;
//		                f2 *= 10;
//		            }
//		            for (int i = k; i < 0; i++) {
//		                f /= 10;
//		                f2 /= 10;
//		            }
//		            //double fAdjXRange = f;
//
//		            fXRange = fXRange + fMinValue - f2;
//		            fXRange = adjust(fXRange);
//		            final int NR_OF_TICKS_X = m_nTicks;
//
//		            fMinValue = f2; //fXRange = fAdjXRange;
//
//		            double fYMax = 0;
//		            for (int i = 0; i < nPoints2; i++) {
//		                xPoints2[i] = graphoffset + nGraphWidth * i / nPoints2;  
//		                if (m_distr != null) {
//		                    try {
//		                    	if (m_distr instanceof CladeAgeDistribution) {
//			                        fyPoints2[i] = m_distr.density(fMinValue + (fXRange * i) / nPoints2);
//		                    	} else {
//		                    		fyPoints2[i] = m_distr.density(fMinValue + (fXRange * i) / nPoints2);// - minOccuranceAge);X
//		                    	}
//		                        if (Double.isInfinite(fyPoints2[i]) || Double.isNaN(fyPoints2[i])) {
//		                        	fyPoints2[i] = 0;
//		                        }
//		                    } catch (Exception e) {
//		                        fyPoints2[i] = 0;
//		                    }
//		                }
//		            }
//		            
//		            for (int i = 0; i < nPoints; i++) {
//		                xPoints[i] = (int)(graphoffset + nGraphWidth * (ages[ages.length - 1 - i] - f2)/fXRange);  
//
//		                fyPoints[i] = probabilities[nPoints - i - 1]/probs.getNormaliser();
//
//		                
//		                if (Double.isInfinite(fyPoints[i]) || Double.isNaN(fyPoints[i])) {
//		                    fyPoints[i] = 0;
//		                }
//		                //fyPoints[i] = Math.exp(m_distr.logDensity(fMinValue + (fXRange * i)/nPoints));
//		                fYMax = Math.max(fYMax, fyPoints[i]);
//		            }
//
//		            fYMax = adjust(fYMax);
//		            final int NR_OF_TICKS_Y = m_nTicks;
//
//		            
//		            for (int i = 0; i < nPoints; i++) {
//		                yPoints[i] = 1 + (int) (graphoffset + nGraphHeight - nGraphHeight * fyPoints[i] / fYMax);
//		                g.drawLine(xPoints[i]+2, yPoints[i], xPoints[i]-2, yPoints[i]);
//		                g.drawLine(xPoints[i], yPoints[i]-2, xPoints[i], yPoints[i]+2);
//		            }
//		            for (int i = 0; i < nPoints2; i++) {
//		                yPoints2[i] = 1 + (int) (graphoffset + nGraphHeight - nGraphHeight * fyPoints2[i] / fYMax);
//		            }
//		            if (m_distr != null) {
//		            	g.drawPolyline(xPoints2, yPoints2, nPoints2);
//		            }
//		            
//
//		            // draw ticks on edge
//		            Font smallFont = new Font(font.getName(), font.getStyle(), 8);
//		            g.setFont(smallFont);
//		            for (int i = 0; i <= NR_OF_TICKS_X; i++) {
//		                int x = graphoffset + i * nGraphWidth / NR_OF_TICKS_X;
//		                g.drawLine(x, graphoffset + nGraphHeight, x, graphoffset + nGraphHeight + 5);
//		                g.drawString(format(fMinValue + fXRange * i / NR_OF_TICKS_X), x + 2, graphoffset + nGraphHeight + 5 + 4);
//		            }
//		            for (int i = 0; i <= NR_OF_TICKS_Y; i++) {
//		                int y = graphoffset + nGraphHeight - i * nGraphHeight / NR_OF_TICKS_Y;
//		                g.drawLine(graphoffset - 5, y, graphoffset, y);
//		                g.drawString(format(fYMax * i / NR_OF_TICKS_Y), 0, y - 2);
//		            }
//		            
//		            // draw statistics
//		            if (m_distr != null) {
//			            int statoffsetx = getWidth() - 120;
//			            int statoffsety = graphoffset + 10;
////		            	String text = getFitParameters();
////		            	String [] strs = text.split("\n");
////		            	for (int i = 0; i < strs.length; i++) {
////			                g.drawString(strs[i], statoffsetx, statoffsety + i * 10);
////		            	}
//		            }
//				}
//			};
//
//			/**
//		     * maps most significant digit to nr of ticks on graph *
//		     */
//		    final int[] NR_OF_TICKS = new int[]{5, 10, 8, 6, 8, 10, 6, 7, 8, 9, 10};
//			int m_nTicks = 5;
//			
//	        private double adjust(double fYMax) {
//	            // adjust fYMax so that the ticks come out right
//	            int k = 0;
//	            double fY = fYMax;
//	            while (fY > 10) {
//	                fY /= 10;
//	                k++;
//	            }
//	            while (fY < 1 && fY > 0) {
//	                fY *= 10;
//	                k--;
//	            }
//	            fY = Math.ceil(fY);
//	            m_nTicks = NR_OF_TICKS[(int) fY];
//	            for (int i = 0; i < k; i++) {
//	                fY *= 10;
//	            }
//	            for (int i = k; i < 0; i++) {
//	                fY /= 10;
//	            }
//	            return fY;
//	        }
//	        
//
//		};
        for (Node n : panel_1.getChildren()) {
        	if (n instanceof Control) {
        		((Control)n).setTooltip(new Tooltip("Click to show parameters of fit"));
        	}
        }
		panel_1.setOnMouseClicked(e->{
				String text = getFitParameters();
				text = text.replaceAll("\n", "<br/>");
				showHelp("<html>Distribution fit:<br/>" + text + "</html>");
			}
		);
		
//		panel_1.setBorder(new LineBorder(new Color(0, 0, 0)));
		panel_1.setBackground(Background.fill(Color.GRAY));
		panel_1.setPrefSize(1024,400);
		//GridBagConstraints //gbc_panel_1 = new GridBagConstraints();
		//gbc_panel_1.fill = GridBagConstraints.BOTH;
		//gbc_panel_1.gridx = 0;
		//gbc_panel_1.gridy = 1;
		if (mode != MODE_BEAUTI_BOTTOM) {
			panel3.getChildren().add(panel_1); //gbc_panel_1);
		}
//		JPanel panel4 = new JPanel();
//		panel4.setBorder(new TitledBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null), "Approximation", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(51, 51, 51)));
//		GridBagConstraints //gbc_panel4 = new GridBagConstraints();
//		//gbc_panel4.fill = GridBagConstraints.BOTH;
//		//gbc_panel4.gridx = 0;
//		//gbc_panel4.gridy = 2;
//		//gbc_panel4.gridwidth = 3;
//		add(panel4, //gbc_panel4);
//		GridBagLayout gbl_panel4 = new GridBagLayout();
//		gbl_panel4.columnWidths = new int[]{0, 0, 0, 0};
//		gbl_panel4.rowHeights = new int[]{0, 0, 0};
//		gbl_panel4.columnWeights = new double[]{0.0, 0.0, 0.0, 0.0};
//		gbl_panel4.rowWeights = new double[]{Double.MIN_VALUE, 0.0, 0.0};
//		panel4.setLayout(gbl_panel4);
//		
//		Label lblDistributionType = new Label("Distribution type:");
//		GridBagConstraints //gbc_lblDistributionType = new GridBagConstraints();
//		//gbc_lblDistributionType.insets = new Insets(0, 0, 0, 5);
//		//gbc_lblDistributionType.anchor = GridBagConstraints.EAST;
//		//gbc_lblDistributionType.gridx = 0;
//		//gbc_lblDistributionType.gridy = 2;
//		panel4.add(lblDistributionType, //gbc_lblDistributionType);
//		
//		
//		btnFindApproximation = new Button("Find approximation");
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
//		GridBagConstraints //gbc_btnFindApproximation = new GridBagConstraints();
//		//gbc_btnFindApproximation.insets = new Insets(0, 0, 0, 5);
//		//gbc_btnFindApproximation.gridx = 2;
//		//gbc_btnFindApproximation.gridy = 2;
//		panel4.add(btnFindApproximation, //gbc_btnFindApproximation);
//		
//		Component horizontalGlue = Box.createHorizontalGlue();
//		GridBagConstraints //gbc_horizontalGlue = new GridBagConstraints();
//		//gbc_horizontalGlue.gridx = 3;
//		//gbc_horizontalGlue.gridy = 2;
//		panel4.add(horizontalGlue, //gbc_horizontalGlue);
		
		dataToGUI();

		setToolTipText(textField_maxOccuranceAge, OCCURRENCE_AGE_HELP);
		setToolTipText(textField_maxDivRate, DIV_RATE_HELP);
		setToolTipText(textField_maxTurnoverRate, TURNOVER_RATE_HELP);
		setToolTipText(textField_maxSamplingRate, SAMPLING_RATE_HELP);
		// setToolTipText(textField_maxSamplingGap, SAMPLING_GAP_HELP);
		setToolTipText(textField_minOccuranceAge, OCCURRENCE_AGE_HELP);
		setToolTipText(textField_minDivRate, DIV_RATE_HELP);
		setToolTipText(textField_minTurnoverRate, TURNOVER_RATE_HELP);
		setToolTipText(textField_minSamplingRate, SAMPLING_RATE_HELP);
		// setToolTipText(textField_minSamplingGap, SAMPLING_GAP_HELP);
	}

	
	
	protected void refresh() {
	if (ages != null) {
        final int width = 1024;
        final int height = 1024;
        final int graphoffset = 30;
        final int labeloffset = 0;
        int nGraphWidth = width - graphoffset * 2;
        int nGraphHeight = height - graphoffset * 2 - labeloffset;
        int nPoints = ages.length;
        int nPoints2 = width;
        int[] xPoints = new int[nPoints];
        int[] yPoints = new int[nPoints];
        int[] xPoints2 = new int[nPoints2];
        int[] yPoints2 = new int[nPoints2];
        double[] fyPoints = new double[nPoints];
        double[] fyPoints2 = new double[nPoints2];
        double fMinValue = 0.1;
        double fMaxValue = 1;
        fMinValue = Math.min(ages[0], ages[ages.length-1]);
        fMaxValue = Math.max(ages[0], ages[ages.length-1]);
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
//		final int NR_OF_TICKS_X = NR_OF_TICKS[(int) f];
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
//        fXRange = adjust(fXRange);
//        final int NR_OF_TICKS_X = m_nTicks;

        fMinValue = f2; //fXRange = fAdjXRange;

        double fYMax = 0;
        for (int i = 0; i < nPoints2; i++) {
            xPoints2[i] = graphoffset + nGraphWidth * i / nPoints2;  
            if (m_distr != null) {
                try {
                	if (m_distr instanceof CladeAgeDistribution) {
                        fyPoints2[i] = m_distr.density(fMinValue + (fXRange * i) / nPoints2);
                	} else {
                		fyPoints2[i] = m_distr.density(fMinValue + (fXRange * i) / nPoints2);// - minOccuranceAge);X
                	}
                    if (Double.isInfinite(fyPoints2[i]) || Double.isNaN(fyPoints2[i])) {
                    	fyPoints2[i] = 0;
                    }
                } catch (Exception e) {
                    fyPoints2[i] = 0;
                }
            }
        }
        
        for (int i = 0; i < nPoints; i++) {
            xPoints[i] = (int)(graphoffset + nGraphWidth * (ages[ages.length - 1 - i] - f2)/fXRange);  

            fyPoints[i] = probabilities[nPoints - i - 1]/probs.getNormaliser();

            
            if (Double.isInfinite(fyPoints[i]) || Double.isNaN(fyPoints[i])) {
                fyPoints[i] = 0;
            }
            //fyPoints[i] = Math.exp(m_distr.logDensity(fMinValue + (fXRange * i)/nPoints));
            fYMax = Math.max(fYMax, fyPoints[i]);
        }

//        fYMax = adjust(fYMax);
//        final int NR_OF_TICKS_Y = m_nTicks;

        series.getData().clear();
        for (int i = 0; i < ages.length; i++) {
        	series.getData().add(new XYChart.Data<Number,Number>(ages[ages.length - 1 - i], fyPoints[i]));
        }

        
//        for (int i = 0; i < nPoints; i++) {
//            yPoints[i] = 1 + (int) (graphoffset + nGraphHeight - nGraphHeight * fyPoints[i] / fYMax);
//            g.drawLine(xPoints[i]+2, yPoints[i], xPoints[i]-2, yPoints[i]);
//            g.drawLine(xPoints[i], yPoints[i]-2, xPoints[i], yPoints[i]+2);
//        }
//        for (int i = 0; i < nPoints2; i++) {
//            yPoints2[i] = 1 + (int) (graphoffset + nGraphHeight - nGraphHeight * fyPoints2[i] / fYMax);
//        }
//        if (m_distr != null) {
//        	g.drawPolyline(xPoints2, yPoints2, nPoints2);
//        }
//        
//
//        // draw ticks on edge
//        Font smallFont = new Font(font.getName(), font.getStyle(), 8);
//        g.setFont(smallFont);
//        for (int i = 0; i <= NR_OF_TICKS_X; i++) {
//            int x = graphoffset + i * nGraphWidth / NR_OF_TICKS_X;
//            g.drawLine(x, graphoffset + nGraphHeight, x, graphoffset + nGraphHeight + 5);
//            g.drawString(format(fMinValue + fXRange * i / NR_OF_TICKS_X), x + 2, graphoffset + nGraphHeight + 5 + 4);
//        }
//        for (int i = 0; i <= NR_OF_TICKS_Y; i++) {
//            int y = graphoffset + nGraphHeight - i * nGraphHeight / NR_OF_TICKS_Y;
//            g.drawLine(graphoffset - 5, y, graphoffset, y);
//            g.drawString(format(fYMax * i / NR_OF_TICKS_Y), 0, y - 2);
//        }
//        
//        // draw statistics
//        if (m_distr != null) {
//            int statoffsetx = getWidth() - 120;
//            int statoffsety = graphoffset + 10;
////        	String text = getFitParameters();
////        	String [] strs = text.split("\n");
////        	for (int i = 0; i < strs.length; i++) {
////                g.drawString(strs[i], statoffsetx, statoffsety + i * 10);
////        	}
//        }
	}
};

	
	private void setToolTipText(Control component, String text) {
		if (component != null) {
			text = normalise(text);
			component.setTooltip(new Tooltip(text));
		}
	}
	
	private Button newHelpButton() {
		return new SmallButton("?",true);
	}
	
	public void calcFit(ProgressBar progress) {
		String type = comboBox.getValue().toString();
		
		if (type.equals("standard CladeAge")) {
			try {
				m_distr = probs.run_standard_cladeage(minOccuranceAge,maxOccuranceAge,minDivRate,maxDivRate,minTurnoverRate,maxTurnoverRate,minSamplingRate,maxSamplingRate, progress);
				// m_distr = probs.run_standard_cladeage(minOccuranceAge,maxOccuranceAge,minDivRate,maxDivRate,minTurnoverRate,maxTurnoverRate,minSamplingRate,maxSamplingRate,minSamplingGap,maxSamplingGap, progress);
			} catch (Exception e) {
				return;
			}
			m_rmsd = 0;
			return;
		}
		
	}
	
	public void dataToGUI() {
		processingDataToGui = true;
			
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_BOTTOM) {
			textField_minOccuranceAge.setText(minOccuranceAge + "");
			// textField_minSamplingGap.setText(minSamplingGap + "");

			if (maxOccuranceAge != minOccuranceAge) {
				textField_maxOccuranceAge.setText(maxOccuranceAge + "");
			} else {
				textField_maxOccuranceAge.setText("As minimum");
			}
			// if (maxSamplingGap != minSamplingGap) {
			// 	textField_maxSamplingGap.setText(maxSamplingGap + "");
			// } else {
			// 	textField_maxSamplingGap.setText("As minimum");
			// }
		}
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_TOP) {
			textField_minDivRate.setText(minDivRate + "");
			textField_minTurnoverRate.setText(minTurnoverRate + "");
			textField_minSamplingRate.setText(minSamplingRate + "");
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
			comboBox.setValue(cladeAgeMethod);
		}

		processingDataToGui = false;
	}
	
	void guiToData() {
		
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_BOTTOM) {
			maxOccuranceAge = parseDouble(textField_maxOccuranceAge.getText());
			// maxSamplingGap = parseDouble(textField_maxSamplingGap.getText());
			minOccuranceAge = parseDouble(textField_minOccuranceAge.getText());
			if (minOccuranceAge < 0) {
				minOccuranceAge = 0;
			}
			if (Double.isInfinite(maxOccuranceAge) || maxOccuranceAge < minOccuranceAge) {
				maxOccuranceAge = minOccuranceAge;
			}
			// minSamplingGap = parseDouble(textField_minSamplingGap.getText());
			// if (minSamplingGap < 0) {
			// 	minSamplingGap = 0;
			// }
			// if (Double.isInfinite(maxSamplingGap) || maxSamplingGap < minSamplingGap) {
			// 	maxSamplingGap = minSamplingGap;
			// }
		}
		if (mode == MODE_STAND_ALONE || mode == MODE_BEAUTI_TOP) {
			maxDivRate = parseDouble(textField_maxDivRate.getText());
			maxTurnoverRate = parseDouble(textField_maxTurnoverRate.getText());
			maxSamplingRate = parseDouble(textField_maxSamplingRate.getText());
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
			cladeAgeMethod = (CladeAgeMethod) comboBox.getValue();
		}	

		// something changed, so the probabilities are not valid any more
		m_distr = null;
		
		for (CAPanelListener listener : listeners) {
			listener.update();
		}
		//dataToGUI();
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
	
//	class MyLabel extends Label {
//		private static final long serialVersionUID = 1L;
//		ImageIcon imageIcon;
//	    public MyLabel(ImageIcon icon)
//	    {
//	        super();
//	        this.imageIcon = icon;
//	    }
//	    @Override
//	    public void paintComponent(Graphics g)
//	    {
//	        super.paintComponent(g);
//	        g.drawImage(imageIcon.getImage(),0,0,getWidth(),getHeight(),this);
//	    }
//	}

	static void showHelp(String text) {
			text = normalise(text);
		
		    String title = text.substring(0, text.indexOf(":"));

		    TextArea textArea = new TextArea(text);
	    	textArea.setPrefRowCount(30);
	    	textArea.setPrefColumnCount(70);
	    	textArea.setEditable(true);
	    	
			HBox box = FXUtils.newHBox();
	    	box.getChildren().add(FXUtils.getIcon(CA_ICON2));
            box.getChildren().add(textArea);
            
            Alert.showMessageDialog(null, box, title, Alert.INFORMATION_MESSAGE);
	}

    private static String normalise(String text) {
		text = text.replace("<br>", "\n");
		text = text.replace("<br/>", "\n");
		text = text.replace("&nbsp;", " ");
		text = text.replace("<html>", "");
		text = text.replace("</html>", "");
		text = text.replaceAll("<a[^>]+>", "");
		text = text.replace("</a>", "\n");
		return text;
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
//        if (m_distr instanceof LogNormalImpl) {
//        	double mean = ((LogNormalImpl) m_distr).getMean();
//        	double sigma = ((LogNormalImpl) m_distr).getSigma();
//            text += "mean: " + format(mean,5) + "\n";
//            text += "sigma: " + format(sigma,5) + "\n";
//        }
//        if (m_distr instanceof GammaDistributionImpl) {
//        	double alpha = ((GammaDistributionImpl) m_distr).getAlpha();
//        	double beta = ((GammaDistributionImpl) m_distr).getBeta();
//            text += "alpha: " + format(alpha,5) + "\n";
//            text += "beta: " + format(beta,5) + "\n";
//        }
//        if (m_distr instanceof ExpGamma) {
//        	double mean = ((ExpGamma) m_distr).getMean();
//        	double alpha = ((ExpGamma) m_distr).getAlpha();
//        	double beta = ((ExpGamma) m_distr).getBeta();
//        	double weight = ((ExpGamma) m_distr).getWeight();
//            text += "mean: " + format(mean,5) + "\n";
//            text += "alpha: " + format(alpha,5) + "\n";
//            text += "beta: " + format(beta,5) + "\n";
//            text += "weight: " + format(weight,5) + "\n";
//        }
        try {
        	if (m_distr instanceof CladeAgeDistribution) {
    	        text += "\nmedian: " + format(m_distr.inverseCumulativeProbability(0.5), 5) + "\n";
    	        text += "\n95% HPD: " + "\n" + 
    	        		format(m_distr.inverseCumulativeProbability(0.025), 5) + 
    	        		" to " + format(m_distr.inverseCumulativeProbability(0.975)) + "\n";
        	} else {
		        text += "\nmedian: " + format(minOccuranceAge + m_distr.inverseCumulativeProbability(0.5), 5) + "\n";
		        text += "\n95% HPD: " + "\n" + 
		        		format(minOccuranceAge + m_distr.inverseCumulativeProbability(0.025), 5) + 
		        		" to " + format(minOccuranceAge + m_distr.inverseCumulativeProbability(0.975)) + "\n";
        	}
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

    
//    boolean bAdvancedFeatures = false;
//
//	private JMenuBar makeMenuBar() {
//        JMenuBar menuBar = new JMenuBar();
//        JMenu fileMenu = new JMenu("File");
//        fileMenu.setMnemonic('F');
//        menuBar.add(fileMenu);
//        fileMenu.add(new MyAction("New", "Start new Clade Age", "new", KeyEvent.VK_N) {
//        	public void actionPerformed(ActionEvent ae) {
//                main(new String[0]);
//            }
//        });
//
//        if (!Utils.isMac()) {
//            fileMenu.addSeparator();
//            fileMenu.add(new MyAction("Close", "Close Window", "close", KeyEvent.VK_W) {
//			    public void actionPerformed(ActionEvent ae) {
//			        JMenuItem menuItem = (JMenuItem) ae.getSource();
//			        JPopupMenu popupMenu = (JPopupMenu) menuItem.getParent();
//			        Component invoker = popupMenu.getInvoker();
//			        JComponent invokerAsJComponent = (JComponent) invoker;
//			        Container topLevel = invokerAsJComponent.getTopLevelAncestor();
//			        if (topLevel != null) {
//			            ((JFrame) topLevel).dispose();
//			        }
//			    }
//			});
//            
//            fileMenu.add(new MyAction("Exit", "Exit Program", "exit", KeyEvent.VK_F4) {
//            	public void actionPerformed(ActionEvent ae) {
//            		System.exit(0);
//            	}
//            });
//        }
//        
//        
//        if (!Utils.isMac()) {
//            JMenu helpMenu = new JMenu("Help");
//            helpMenu.setMnemonic('H');
//            menuBar.add(helpMenu);
//            helpMenu.add(new MyAction("About", "Help about", "help", -1) {
//		        public void actionPerformed(ActionEvent ae) {
//		        	showHelp(ABOUT_HELP);
//		        }
//            });
//        }
//
//        
//		return menuBar;
//	}

	List<CAPanelListener> listeners = new ArrayList<CAPanelListener>();
	public void addChangeListener(CAPanelListener o) {
		listeners.add(o);
	}

	private TextField newTextField() {
		TextField entry = new TextField();
		entry.setOnKeyReleased(e -> guiToData());
	    return entry;
	}
	
//    public static void main(String[] args) {
//		JFrame frame = new JFrame();
//        Utils.loadUIManager();
//
//        if (Utils.isMac()) {
//            // set up application about-menu for Mac
//            // Mac-only stuff
//        	try {
//            ImageIcon icon = Utils.getIcon(ModelBuilder.ICONPATH + "beauti.png");
//            if (icon == null) {
//                System.err.println("Unable to find image: " + ModelBuilder.ICONPATH + "beauti.png");
//            }
//            jam.framework.Application application = new jam.framework.MultiDocApplication(null, "CladeAge", "about" , icon) {
//
//                @Override
//                protected JFrame getDefaultFrame() {
//                    return null;
//                }
//
//                @Override
//                public void doQuit() {
//                    System.exit(0);
//                }
//
//                @Override
//                public void doAbout() {
//                    showHelp(ABOUT_HELP);
//                }
//
//				@Override
//                public DocumentFrame doOpenFile(File file) {
//                    return null;
//                }
//
//                @Override
//                public DocumentFrame doNew() {
//                    return null;
//                }
//            };
//            jam.mac.Utils.macOSXRegistration(application);
//        	} catch (Exception e) {
//        		// ignore
//        	}
//            try {
//            	Class<?> class_ = Class.forName("jam.maconly.OSXAdapter");
//                Method method = class_.getMethod("enablePrefs", boolean.class);
//                method.invoke(null, false);
//            } catch (java.lang.Exception e) {
//            	// ignore
//            }
//        }
//        
//		frame.setSize(1024, 728);
//		
//        ImageIcon icon = Utils.getIcon(CA_ICON);
//        if (icon != null) {
//            frame.setIconImage(icon.getImage());
//        }
//		CAPanel pane = new CAPanel(CAPanel.MODE_STAND_ALONE);
//		
//        JMenuBar menuBar = pane.makeMenuBar();
//        frame.setJMenuBar(menuBar);
//		frame.getContentPane().add(pane);
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);	
//		frame.setVisible(true);
//	} // main

}
