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
	        	FossilCalibration fossilCallibration = new FossilCalibration();
	        	fossilCallibration.setID("FossilCallibration.0");
	        	doc.registerPlugin(fossilCallibration);
		        prior.callibrationDistr.setValue(fossilCallibration, prior);

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
	            prior.setID(taxonSet.getID() + ".prior");
	            // this sets up the type
	            prior.m_distInput.setValue(new OneOnX(), prior);
	            // this removes the parametric distribution
	            prior.m_distInput.setValue(null, prior);

	            Logger logger = (Logger) doc.pluginmap.get("tracelog");
	            logger.m_pLoggers.setValue(prior, logger);
	        } catch (Exception e) {
	            // TODO: handle exception
	        }
	        List<Plugin> selectedPlugins = new ArrayList<Plugin>();
	        selectedPlugins.add(prior);
	        g_collapsedIDs.add(prior.getID());
	        return selectedPlugins;
	    }

}
