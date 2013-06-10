package beast.app.beauti;


import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;

import beast.app.draw.PluginPanel;
import beast.core.Input;
import beast.core.Logger;
import beast.core.Plugin;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;
import beast.math.distributions.FossilCalibration;
import beast.math.distributions.FossilPrior;
import beast.math.distributions.OneOnX;

public class FossilPriorListInputEditor extends PriorListInputEditor {
	    private static final long serialVersionUID = 1L;

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
	    
	    @Override
	    public List<Plugin> pluginSelector(Input<?> input, Plugin parent, List<String> sTabuList) {
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
	            prior.m_treeInput.setValue(trees.get(iTree), prior);
	            TaxonSet taxonSet = new TaxonSet();

	            TaxonSetDialog dlg = new TaxonSetDialog(taxonSet, getTaxonCandidates(prior), doc);
	            if (!dlg.showDialog() || dlg.taxonSet.getID() == null) {
	                return null;
	            }
	            taxonSet = dlg.taxonSet;
	            PluginPanel.addPluginToMap(taxonSet, doc);
	            prior.m_taxonset.setValue(taxonSet, prior);
	            //prior.setID(taxonSet.getID() + ".fossilprior");
	            prior.setID(taxonSet.getID()+".fossilprior");
	            // this sets up the type
	            prior.m_distInput.setValue(new OneOnX(), prior);
	            // this removes the parametric distribution
	            prior.m_distInput.setValue(null, prior);

	            Logger logger = (Logger) doc.pluginmap.get("tracelog");
	            logger.m_pLoggers.setValue(prior, logger);
//	            CompoundFossilPrior compoundPrior = (CompoundFossilPrior) doc.pluginmap.get("fossilCalibrations");
//	            compoundPrior.distributionsInput.setValue(prior, compoundPrior);
	            CompoundDistribution compoundPrior = (CompoundDistribution) doc.pluginmap.get("fossilCalibrations");
	            compoundPrior.pDistributions.setValue(prior, compoundPrior);

//				FossilCalibration fossilCallibration = new FossilCalibration();
//				fossilCallibration.setInputValue("minOccuranceAge",new RealParameter("0.0"));
//				fossilCallibration.setInputValue("maxOccuranceAge",new RealParameter("0.0"));
//				fossilCallibration.setInputValue("minDivRate",     new RealParameter("0.01"));
//				fossilCallibration.setInputValue("maxDivRate",     new RealParameter("0.01"));
//				fossilCallibration.setInputValue("maxSamplingRate",new RealParameter("0.01"));
//				fossilCallibration.setInputValue("minSamplingRate",new RealParameter("0.01"));
//				fossilCallibration.setInputValue("minTurnoverRate",new RealParameter("0.001"));
//				fossilCallibration.setInputValue("maxTurnoverRate",new RealParameter("0.001"));
//				fossilCallibration.setInputValue("maxSamplingGap", new RealParameter("0.0"));
//				fossilCallibration.setInputValue("minSamplingGap", new RealParameter("0.0"));

				List<BeautiSubTemplate> availablePlugins = doc.getInpuEditorFactory().getAvailableTemplates(
						new Input<FossilCalibration>("fossil","",FossilCalibration.class), doc, null, doc);
				FossilCalibration fossilCallibration = (FossilCalibration) availablePlugins.get(0).createSubNet(new PartitionContext(), false);
	        	fossilCallibration.setID("FossilCallibration.0");
	        	doc.registerPlugin(fossilCallibration);
		        prior.callibrationDistr.setValue(fossilCallibration, prior);

	            
	        } catch (Exception e) {
	        	e.printStackTrace();
	            // TODO: handle exception
	        }
	        List<Plugin> selectedPlugins = new ArrayList<Plugin>();
	        selectedPlugins.add(prior);
	        g_collapsedIDs.add(prior.getID());	        
	        return selectedPlugins;
	    }

}
