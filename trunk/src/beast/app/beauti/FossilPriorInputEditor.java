package beast.app.beauti;

import beast.math.distributions.FossilPrior;

public class FossilPriorInputEditor extends MRCAPriorInputEditor {

	public FossilPriorInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return FossilPrior.class;
	}

}
