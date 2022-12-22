package cladeage.app.beauti;

import java.util.HashSet;
import java.util.List;
import java.util.Set;


import beastfx.app.inputeditor.BeautiDoc;
// import beastfx.app.inputeditor.Expandable;
import beastfx.app.inputeditor.ListInputEditor;
import beastfx.app.inputeditor.MRCAPriorInputEditor;
import beastfx.app.beauti.PriorInputEditor;
import beastfx.app.inputeditor.TaxonSetDialog;
import beastfx.app.util.FXUtils;
import javafx.geometry.Insets;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import cladeage.math.distributions.FossilCalibration;
import cladeage.math.distributions.FossilPrior;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.inference.distribution.OneOnX;


public class FossilPriorInputEditor extends MRCAPriorInputEditor { //implements Expandable {

	public FossilPriorInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return FossilPrior.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, final int listItemNr, ExpandOption bExpandOption, boolean bAddButtons) {
		String className = FossilPrior.class.getName();
		doc.beautiConfig.suppressBEASTObjects.add(className + ".tree");
		doc.beautiConfig.suppressBEASTObjects.add(className + ".distr");
		doc.beautiConfig.suppressBEASTObjects.add(className + ".monophyletic");
		doc.beautiConfig.suppressBEASTObjects.add(className + ".useOriginate");
		doc.beautiConfig.suppressBEASTObjects.add(className + ".tipsonly");
		doc.beautiConfig.suppressBEASTObjects.add(className + ".taxonset");

		className = FossilCalibration.class.getName();
		doc.beautiConfig.suppressBEASTObjects.add(className + ".offset");

        m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = plugin;
        this.itemNr= listItemNr;
        
		
        HBox itemBox = new HBox();
        itemBox.setSpacing(5);

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
        isMonophyleticdBox.setPadding(new Insets(5));
        isMonophyleticdBox.setOnAction(e->{
            prior.isMonophyleticInput.setValue(isMonophyleticdBox.isSelected(), prior);
        });
        itemBox.getChildren().add(isMonophyleticdBox);
        //itemBox.getChildren().add(Box.createGlue());

        getChildren().add(itemBox);
        setPrefHeight(20);
        setMaxHeight(20);
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

    VBox expandBox = null;
	public void setExpandBox(VBox expandBox) {
		this.expandBox = expandBox;
	}
	
	@Override
	public void refreshPanel() {
		if (expandBox != null) {
			ListInputEditor.updateExpandBox(doc, expandBox, m_beastObject, this);
		}
		super.refreshPanel();
	}	

}
