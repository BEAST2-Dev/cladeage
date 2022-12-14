package cladeage.app.beauti;



import java.util.ArrayList;
import java.util.List;


import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiSubTemplate;
import beast.base.parser.PartitionContext;
import beastfx.app.beauti.PriorListInputEditor;
import beastfx.app.inputeditor.TaxonSetDialog;
import beastfx.app.util.Alert;
import cladeage.app.ca.CAPanel;
import cladeage.app.ca.CAPanelListener;
import beastfx.app.inputeditor.BEASTObjectPanel;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.CompoundDistribution;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Tree;
import cladeage.math.distributions.FossilCalibration;
import cladeage.math.distributions.FossilPrior;
import javafx.scene.layout.Pane;
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
	    
	    CAPanel caPanel;
	    FossilCalibration calibration;
	    
	    @Override
	    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
	    	super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
	    	
	    	// add panel at end
	    	List<?> list = (List) input.get();
	    	if (list.size() > 0) {
	    		calibration = ((FossilPrior) list.get(0)).calibrationDistr.get();
	    		caPanel = new CAPanel(CAPanel.MODE_BEAUTI_TOP);
	    		caPanel.setMinDivRate(calibration.minDivRateInput.get().getValue());
	    		caPanel.setMinTurnoverRate(calibration.minTurnoverRateInput.get().getValue());
	    		caPanel.setMinSamplingRate(calibration.minSamplingRateInput.get().getValue());
	    		caPanel.setMaxDivRate(calibration.maxDivRateInput.get().getValue());
	    		caPanel.setMaxTurnoverRate(calibration.maxTurnoverRateInput.get().getValue());
	    		caPanel.setMaxSamplingRate(calibration.maxSamplingRateInput.get().getValue());
	    		caPanel.setMethod(calibration.cladeAgeMethodInput.get());
	            caPanel.dataToGUI();
	            caPanel.addChangeListener(this);
	            caPanel.setVisible(true);
	            caPanel.setPrefWidth(640);
	            Pane p = ((Pane)pane.getChildren().get(0));
	            int pos = p.getChildren().size() - 1;
	            p.getChildren().add(pos, caPanel);
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
	                
	                String treeID = (String) Alert.showInputDialog(null, "Select a tree", "MRCA selector", Alert.QUESTION_MESSAGE, null, sTreeIDs, trees.get(0));
	                iTree = 0;
	                while (iTree < sTreeIDs.length && !sTreeIDs[iTree].equals(treeID)) {
	                	iTree++;
	                }
	                if (iTree == sTreeIDs.length) {
	                	iTree = -1;
	                }
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
			setValue(calibration.minDivRateInput, caPanel.getMinDivRate());
			setValue(calibration.maxDivRateInput, caPanel.getMaxDivRate());
			setValue(calibration.minTurnoverRateInput, caPanel.getMinTurnoverRate());
			setValue(calibration.maxTurnoverRateInput, caPanel.getMaxTurnoverRate());
			setValue(calibration.minSamplingRateInput, caPanel.getMinSamplingRate());
			setValue(calibration.maxSamplingRateInput, caPanel.getMaxSamplingRate());
			try {
				calibration.cladeAgeMethodInput.setValue(caPanel.getMethod(), calibration);
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
