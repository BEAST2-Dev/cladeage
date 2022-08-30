package beast.app.beauti;



import java.util.ArrayList;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JOptionPane;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiSubTemplate;
import beast.base.parser.PartitionContext;
import beastfx.app.beauti.PriorListInputEditor;
import beastfx.app.inputeditor.TaxonSetDialog;
import beast.app.ca.CAPanel;
import beast.app.ca.CAPanelListener;
import beastfx.app.inputeditor.BEASTObjectPanel;
import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.CompoundDistribution;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Tree;
import beast.math.distributions.FossilCalibration;
import beast.math.distributions.FossilPrior;
import beast.base.inference.distribution.OneOnX;


public class FossilPriorListInputEditor extends PriorListInputEditor implements CAPanelListener {

		public FossilPriorListInputEditor(BeautiDoc doc) {
			super(doc);
		}

	    @Override
	    public Class<?> type() {
	        return List.class;
	    }

	    @Override
	    public Class<?> baseType() {
	        return FossilPrior.class;
	    }
	    
	    CAPanel panel;
	    FossilCalibration calibration;
	    
	    @Override
	    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
	    	super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
	    	List<?> list = (List) input.get();
	    	if (list.size() > 0) {
	    		calibration = ((FossilPrior) list.get(0)).calibrationDistr.get();
	    		panel = new CAPanel(CAPanel.MODE_BEAUTI_TOP);
	            panel.setMinDivRate(calibration.minDivRateInput.get().getValue());
	            panel.setMinTurnoverRate(calibration.minTurnoverRateInput.get().getValue());
	            panel.setMinSamplingRate(calibration.minSamplingRateInput.get().getValue());
	            panel.setMaxDivRate(calibration.maxDivRateInput.get().getValue());
	            panel.setMaxTurnoverRate(calibration.maxTurnoverRateInput.get().getValue());
	            panel.setMaxSamplingRate(calibration.maxSamplingRateInput.get().getValue());
	            panel.setMethod(calibration.cladeAgeMethodInput.get());
	            panel.dataToGUI();
	            panel.addChangeListener(this);
	            panel.setVisible(true);
	            getChildren().add(panel);
	    	}
	    }
	    
	    @Override
	    public List<BEASTInterface> pluginSelector(Input<?> input, BEASTInterface parent, List<String> sTabuList) {
	        FossilPrior prior = new FossilPrior();
	        try {
	        	
	            List<Tree> trees = new ArrayList<Tree>();
	            getDoc().scrubAll(true, false);
	            State state = (State) doc.pluginmap.get("state");
	            for (StateNode node : state.stateNodeInput.get()) {
	                if (node instanceof Tree) { // && ((Tree) node).m_initial.get() != null) {
	                    trees.add((Tree) node);
	                }
	            }
	            int iTree = 0;
	            if (trees.size() > 1) {
	                String[] sTreeIDs = new String[trees.size()];
	                for (int j = 0; j < sTreeIDs.length; j++) {
	                    sTreeIDs[j] = trees.get(j).getID();
	                }
	                iTree = JOptionPane.showOptionDialog(null, "Select a tree", "MRCA selector",
	                        JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null,
	                        sTreeIDs, trees.get(0));
	            }
	            if (iTree < 0) {
	                return null;
	            }
	            prior.treeInput.setValue(trees.get(iTree), prior);
	            TaxonSet taxonSet = new TaxonSet();

	            TaxonSetDialog dlg = new TaxonSetDialog(taxonSet, getTaxonCandidates(prior, doc), doc);
	            if (!dlg.showDialog() || dlg.taxonSet.getID() == null) {
	                return null;
	            }
	            taxonSet = dlg.taxonSet;
	            BEASTObjectPanel.addPluginToMap(taxonSet, doc);
	            prior.taxonsetInput.setValue(taxonSet, prior);
	            //prior.setID(taxonSet.getID() + ".fossilprior");
	            prior.setID(taxonSet.getID()+".fossilprior");
	            // this sets up the type
	            prior.distInput.setValue(new OneOnX(), prior);
	            // this removes the parametric distribution
	            prior.distInput.setValue(null, prior);

	            Logger logger = (Logger) doc.pluginmap.get("tracelog");
	            logger.loggersInput.setValue(prior, logger);
//	            CompoundFossilPrior compoundPrior = (CompoundFossilPrior) doc.pluginmap.get("fossilCalibrations");
//	            compoundPrior.distributionsInput.setValue(prior, compoundPrior);
	            CompoundDistribution compoundPrior = (CompoundDistribution) doc.pluginmap.get("fossilCalibrations");
	            compoundPrior.pDistributions.setValue(prior, compoundPrior);

//				FossilCalibration fossilCalibration = new FossilCalibration();
//				fossilCalibration.setInputValue("minOccuranceAge",new RealParameter("0.0"));
//				fossilCalibration.setInputValue("maxOccuranceAge",new RealParameter("0.0"));
//				fossilCalibration.setInputValue("minDivRate",     new RealParameter("0.01"));
//				fossilCalibration.setInputValue("maxDivRate",     new RealParameter("0.01"));
//				fossilCalibration.setInputValue("maxSamplingRate",new RealParameter("0.01"));
//				fossilCalibration.setInputValue("minSamplingRate",new RealParameter("0.01"));
//				fossilCalibration.setInputValue("minTurnoverRate",new RealParameter("0.001"));
//				fossilCalibration.setInputValue("maxTurnoverRate",new RealParameter("0.001"));

				List<BeautiSubTemplate> availablePlugins = doc.getInputEditorFactory().getAvailableTemplates(
						new Input<FossilCalibration>("fossil","",FossilCalibration.class), doc, null, doc);
				FossilCalibration fossilCalibration = (FossilCalibration) availablePlugins.get(0).createSubNet(new PartitionContext(), false);
	        	fossilCalibration.setID("FossilCalibration.0");
	        	doc.registerPlugin(fossilCalibration);
		        prior.calibrationDistr.setValue(fossilCalibration, prior);

	            
	        } catch (Exception e) {
	        	e.printStackTrace();
	            // TODO: handle exception
	        }
	        List<BEASTInterface> selectedPlugins = new ArrayList<BEASTInterface>();
	        selectedPlugins.add(prior);
	        g_collapsedIDs.add(prior.getID());	        
	        return selectedPlugins;
	    }

		@Override
		public void update() {
			//setValue(calibration.m_offset, panel.getMinOccuranceAge());
			setValue(calibration.minDivRateInput, panel.getMinDivRate());
			setValue(calibration.maxDivRateInput, panel.getMaxDivRate());
			setValue(calibration.minTurnoverRateInput, panel.getMinTurnoverRate());
			setValue(calibration.maxTurnoverRateInput, panel.getMaxTurnoverRate());
			setValue(calibration.minSamplingRateInput, panel.getMinSamplingRate());
			setValue(calibration.maxSamplingRateInput, panel.getMaxSamplingRate());
			try {
				calibration.cladeAgeMethodInput.setValue(panel.getMethod(), calibration);
			} catch (Exception e) {
				e.printStackTrace();
			}
	}

		private void setValue(Input<RealParameter> input, double value) {
			try {
				input.get().valuesInput.setValue(value+"", calibration);
				input.get().setValue(value);
			} catch (Exception e) {
				
			}
		}

}
