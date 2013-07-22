package statalign.postprocess.plugins;

import statalign.base.State;
import statalign.postprocess.Postprocess;

import javax.swing.*;
import java.awt.*;

public class ParameterHistogram extends Postprocess {
    private HistogramGUI gui = new HistogramGUI();

    private StringBuilder RMeasurements = new StringBuilder();
    private StringBuilder lambdaMeasurements = new StringBuilder();
    private StringBuilder muMeasurements = new StringBuilder();

    @Override
    public String getTabName() {
        return "Parameter histograms";
    }

    @Override
    public Icon getIcon() {
        return new ImageIcon(ClassLoader.getSystemResource("icons/loglikelihood.gif"));
    }

    @Override
    public JPanel getJPanel() {
        return gui;
    }

    @Override
    public String getTip() {
        return "A histogram view of the model parameters used in the MCMC process.";
    }

    @Override
    public void setSampling(boolean enabled) {
        // DJ Ötzi - Live is Life
    }

    @Override
    public void newSample(State state, int no, int total) {
        RMeasurements.append(state.indelParams[0] + "\n");
        lambdaMeasurements.append(state.indelParams[1] + "\n");
        muMeasurements.append(state.indelParams[2] + "\n");
    }

    @Override
    public void afterLastSample() {
        System.out.println(lambdaMeasurements.toString());
    }

    private void writeAndGenerateHistogram(String name, StringBuilder measurements) {

    }

    public static class HistogramGUI extends JPanel
    {
        private boolean isHistogramReady = false;

        @Override
        public void paintComponent(Graphics graphics) {
            Graphics2D g = (Graphics2D) graphics;

            if (isHistogramReady) {
                // Show the histogram
            } else {
                Font font = new Font("Andale Mono", Font.BOLD, 36);
                g.setFont(font);

                g.drawString("Waiting for data...", 10, 70);
            }
        }
    }
}
