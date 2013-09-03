package beast.app.beauti;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;

import beast.app.beauti.BeautiDoc;
import beast.app.beauti.MRCAPriorInputEditor;
import beast.app.beauti.PriorInputEditor;
import beast.app.beauti.TaxonSetDialog;
import beast.core.Input;
import beast.core.BEASTObject;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.math.distributions.FossilPrior;
import beast.math.distributions.MRCAPrior;
import beast.math.distributions.OneOnX;


public class FossilPriorInputEditor extends MRCAPriorInputEditor {

	public FossilPriorInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return FossilPrior.class;
	}

	@Override
	public void init(Input<?> input, BEASTObject plugin, final int listItemNr, ExpandOption bExpandOption, boolean bAddButtons) {
		doc.beautiConfig.suppressPlugins.add("beast.math.distributions.FossilPrior.tree");
		doc.beautiConfig.suppressPlugins.add("beast.math.distributions.FossilPrior.distr");
		doc.beautiConfig.suppressPlugins.add("beast.math.distributions.FossilPrior.monophyletic");
		doc.beautiConfig.suppressPlugins.add("beast.math.distributions.FossilPrior.useOriginate");
		doc.beautiConfig.suppressPlugins.add("beast.math.distributions.FossilPrior.tipsonly");
		doc.beautiConfig.suppressPlugins.add("beast.math.distributions.FossilPrior.taxonset");
		doc.beautiConfig.suppressPlugins.add("beast.math.distributions.FossilCalibration.offset");

        m_bAddButtons = bAddButtons;
        m_input = input;
        m_plugin = plugin;
        this.itemNr= listItemNr;
		
        Box itemBox = Box.createHorizontalBox();

        MRCAPrior prior = (MRCAPrior) plugin;
        String sText = prior.taxonsetInput.get().getID();

        JButton taxonButton = new JButton(sText);
        taxonButton.setMinimumSize(PriorInputEditor.PREFERRED_SIZE);
        taxonButton.setPreferredSize(PriorInputEditor.PREFERRED_SIZE);
        itemBox.add(taxonButton);
        taxonButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JButton taxonButton = (JButton) e.getSource();
                List<?> list = (List<?>) m_input.get();
                MRCAPrior prior = (MRCAPrior) list.get(itemNr);
                try {
                    TaxonSet taxonset = prior.taxonsetInput.get();
                    Set<Taxon> candidates = getTaxonCandidates(prior);
                    TaxonSetDialog dlg = new TaxonSetDialog(taxonset, candidates, doc);
                    if (dlg.showDialog()) {
                        prior.setID(dlg.taxonSet.getID()+".fossilprior");
                        prior.taxonsetInput.setValue(dlg.taxonSet, prior);
                    }
                    } catch (Exception e1) {
                    // TODO Auto-generated catch block
                    e1.printStackTrace();
                }
                refreshPanel();
            }
        });


        if (prior.distInput.getType() == null) {
            try {
                prior.distInput.setValue(new OneOnX(), prior);
                prior.distInput.setValue(null, prior);
            } catch (Exception e) {
                // TODO: handle exception
            }

        }

        JCheckBox isMonophyleticdBox = new JCheckBox(doc.beautiConfig.getInputLabel(prior, prior.isMonophyleticInput.getName()));
        isMonophyleticdBox.setName(sText+".isMonophyletic");
        isMonophyleticdBox.setSelected(prior.isMonophyleticInput.get());
        isMonophyleticdBox.setToolTipText(prior.isMonophyleticInput.getTipText());
        isMonophyleticdBox.addActionListener(new MRCAPriorActionListener(prior));
        itemBox.add(isMonophyleticdBox);
        itemBox.add(Box.createGlue());

        add(itemBox);
	}
	
    Set<Taxon> getTaxonCandidates(FossilPrior prior) {
        Set<Taxon> candidates = new HashSet<Taxon>();
        for (String sTaxon : prior.treeInput.get().getTaxaNames()) {
            Taxon taxon = null;
            for (Taxon taxon2 : doc.taxaset) {
                if (taxon2.getID().equals(sTaxon)) {
                    taxon = taxon2;
                    break;
                }
            }
            if (taxon == null) {
                taxon = new Taxon();
                taxon.setID(sTaxon);
                doc.taxaset.add(taxon);
            }
            candidates.add(taxon);
        }
        return candidates;
    }

}
