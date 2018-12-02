package beast.app.beauti;





import java.awt.Dimension;

import beast.app.beauti.BeautiDoc;
import beast.app.ca.CAPanel;
import beast.app.ca.CAPanelListener;
import beast.app.draw.BEASTObjectInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.FossilCalibration;
import beast.math.distributions.FossilCalibration.CladeAgeMethod;


public class FossilCalibrationInputEditor extends BEASTObjectInputEditor implements CAPanelListener {

	private static final long serialVersionUID = 1L;

	public FossilCalibrationInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return FossilCalibration.class;
	}

	CAPanel panel;
	FossilCalibration calibration;
	
	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
        m_input = input;
        m_beastObject = plugin;
        calibration = (FossilCalibration) m_input.get();
        panel = new CAPanel(CAPanel.MODE_BEAUTI_BOTTOM);
        panel.setMethod(calibration.cladeAgeMethodInput.get());
        panel.setMinOccuranceAge(calibration.minOccuranceAgeInput.get().getValue());
//        panel.setMinDivRate(calibration.minDivRateInput.get().getValue());
//        panel.setMinTurnoverRate(calibration.minTurnoverRateInput.get().getValue());
//        panel.setMinSamplingRate(calibration.minSamplingRateInput.get().getValue());

        panel.setMaxOccuranceAge(calibration.maxOccuranceAgeInput.get().getValue());
//        panel.setMaxDivRate(calibration.maxDivRateInput.get().getValue());
//        panel.setMaxTurnoverRate(calibration.maxTurnoverRateInput.get().getValue());
//        panel.setMaxSamplingRate(calibration.maxSamplingRateInput.get().getValue());
        //panel.setMethod(calibration.cladeAgeMethodInput.get());
        panel.dataToGUI();
        panel.setCalculateButtonText("Preview");

        panel.addChangeListener(this);
        Dimension size = new Dimension(panel.getPreferredSize().width - 40, panel.getPreferredSize().height);

        panel.setSize(size);
        panel.setPreferredSize(size);
        panel.setMinimumSize(size);
        panel.setMaximumSize(size);
        add (panel);
	}

	@Override
	public void update() {
		//setValue(calibration.m_offset, panel.getMinOccuranceAge());
		setValue(calibration.minOccuranceAgeInput, panel.getMinOccuranceAge());
		setValue(calibration.maxOccuranceAgeInput, panel.getMaxOccuranceAge());
//		setValue(calibration.minDivRateInput, panel.getMinDivRate());
//		setValue(calibration.maxDivRateInput, panel.getMaxDivRate());
//		setValue(calibration.minTurnoverRateInput, panel.getMinTurnoverRate());
//		setValue(calibration.maxTurnoverRateInput, panel.getMaxTurnoverRate());
//		setValue(calibration.minSamplingRateInput, panel.getMinSamplingRate());
//		setValue(calibration.maxSamplingRateInput, panel.getMaxSamplingRate());
//		setValue(calibration.cladeAgeMethodInput, panel.getMethod());

		panel.setMethod(calibration.cladeAgeMethodInput.get());
		panel.setMinDivRate(calibration.minDivRateInput.get().getValue());
		panel.setMaxDivRate(calibration.maxDivRateInput.get().getValue());
		panel.setMinTurnoverRate(calibration.minTurnoverRateInput.get().getValue());
		panel.setMaxTurnoverRate(calibration.maxTurnoverRateInput.get().getValue());
		panel.setMinSamplingRate(calibration.minSamplingRateInput.get().getValue());
		panel.setMaxSamplingRate(calibration.maxSamplingRateInput.get().getValue());
	
	}

	private void setValue(Input<RealParameter> input, double value) {
		try {
			input.get().valuesInput.setValue(value+"", calibration);
			input.get().setValue(value);
		} catch (Exception e) {
			
		}
	}
	private void setValue(Input<CladeAgeMethod> input, CladeAgeMethod value) {
		try {
			input.setValue(value, calibration);
		} catch (Exception e) {
			
		}
	}
	
	
	
} // class FossilCalibrationInputEditor
