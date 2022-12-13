package cladeage.app.ca;

import beast.base.core.ProgramStatus;
import beastfx.app.beauti.ThemeProvider;
import javafx.application.Application;
import javafx.scene.Scene;
import javafx.stage.Stage;

public class CladeAgeApp extends Application {

	@Override
	public void start(Stage primaryStage) throws Exception {
		CAPanel panel = new CAPanel(CAPanel.MODE_STAND_ALONE);
        Scene scene = new Scene(panel);
        ThemeProvider.loadStyleSheet(scene);
        primaryStage.setScene(scene);
		primaryStage.show();
	}

	public static void main(String[] args) {
    	ProgramStatus.name = "CladeAge";
    	
		launch(CladeAgeApp.class, args);
	}
}
