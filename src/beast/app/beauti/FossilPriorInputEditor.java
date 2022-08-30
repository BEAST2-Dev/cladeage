package beast.app.beauti;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.MRCAPriorInputEditor;
import beastfx.app.beauti.PriorInputEditor;
import beastfx.app.inputeditor.TaxonSetDialog;
import beastfx.app.util.FXUtils;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.math.distributions.FossilPrior;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.inference.distribution.OneOnX;


public class FossilPriorInputEditor extends MRCAPriorInputEditor {

	public FossilPriorInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return FossilPrior.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, final int listItemNr, ExpandOption bExpandOption, boolean bAddButtons) {
		doc.beautiConfig.suppressBEASTObjects.add("beast.math.distributions.FossilPrior.tree");
		doc.beautiConfig.suppressBEASTObjects.add("beast.math.distributions.FossilPrior.distr");
		doc.beautiConfig.suppressBEASTObjects.add("beast.math.distributions.FossilPrior.monophyletic");
		doc.beautiConfig.suppressBEASTObjects.add("beast.math.distributions.FossilPrior.useOriginate");
		doc.beautiConfig.suppressBEASTObjects.add("beast.math.distributions.FossilPrior.tipsonly");
		doc.beautiConfig.suppressBEASTObjects.add("beast.math.distributions.FossilPrior.taxonset");
		doc.beautiConfig.suppressBEASTObjects.add("beast.math.distributions.FossilCalibration.offset");

        m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = plugin;
        this.itemNr= listItemNr;
		
        HBox itemBox = FXUtils.newHBox();

        MRCAPrior prior = (MRCAPrior) plugin;
        String sText = prior.taxonsetInput.get().getID();

        Button taxonButton = new Button(sText);
        taxonButton.setMinSize(PriorInputEditor.PREFERRED_SIZE.getWidth(), PriorInputEditor.PREFERRED_SIZE.getHeight());
        taxonButton.setPrefSize(PriorInputEditor.PREFERRED_SIZE.getWidth(), PriorInputEditor.PREFERRED_SIZE.getHeight());
        itemBox.getChildren().add(taxonButton);
        taxonButton.setOnAction(e-> {
                //JButton taxonButton = (JButton) e.getSource();
                List<?> list = (List<?>) m_input.get();
                //MRCAPrior prior = (MRCAPrior) list.get(itemNr);
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
            });


        if (prior.distInput.getType() == null) {
            try {
                prior.distInput.setValue(new OneOnX(), prior);
                prior.distInput.setValue(null, prior);
            } catch (Exception e) {
                // TODO: handle exception
            }

        }

        CheckBox isMonophyleticdBox = new CheckBox(doc.beautiConfig.getInputLabel(prior, prior.isMonophyleticInput.getName()));
        isMonophyleticdBox.setId(sText+".isMonophyletic");
        isMonophyleticdBox.setSelected(prior.isMonophyleticInput.get());
        isMonophyleticdBox.setTooltip(new Tooltip(prior.isMonophyleticInput.getTipText()));
        isMonophyleticdBox.setOnAction(e->{
            prior.isMonophyleticInput.setValue(isMonophyleticdBox.isSelected(), prior);
        });
        itemBox.getChildren().add(isMonophyleticdBox);
        //itemBox.getChildren().add(Box.createGlue());

        getChildren().add(itemBox);
	}
	
    Set<Taxon> getTaxonCandidates(FossilPrior prior) {
        Set<Taxon> candidates = new HashSet<Taxon>();
        for (String sTaxon : prior.treeInput.get().getTaxaNames()) {
            Taxon taxon = null;
            for (Taxon taxon2 : doc.taxaset.values()) {
                if (taxon2.getID().equals(sTaxon)) {
                    taxon = taxon2;
                    break;
                }
            }
            if (taxon == null) {
                taxon = new Taxon();
                taxon.setID(sTaxon);
                doc.taxaset.put(sTaxon, taxon);
            }
            candidates.add(taxon);
        }
        return candidates;
    }

}
