package beast.app.beauti;



import beast.app.draw.PluginInputEditor;
import beast.core.Input;
import beast.core.Plugin;
import beast.math.distributions.FossilCalibration;

public class FossilCalibrationInputEditor extends PluginInputEditor {

	private static final long serialVersionUID = 1L;

	public FossilCalibrationInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return FossilCalibration.class;
	}

	@Override
	public void init(Input<?> input, Plugin plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
		//JOptionPane.showMessageDialog(null, "NOT IMPLEMENTED YET");
		super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
	}
	
	
	
} // class FossilCalibrationInputEditor
