package statalign.postprocess.plugins;

import statalign.base.Mcmc;
import statalign.base.State;
import statalign.postprocess.Postprocess;

import javax.swing.*;
import java.awt.*;

public class DiagnosticsView extends Postprocess {
    DiagnosticsGUI gui = new DiagnosticsGUI();

    @Override
    public String getTabName() {
        return "Diagnostics";
    }

    @Override
    public Icon getIcon() {
        return new ImageIcon(ClassLoader.getSystemResource("icons/loglikelihood.gif"));
    }

    @Override
    public JPanel getJPanel() {
        return gui;
    }

    public void newSample(State state, int no, int total) {
        gui.update(mcmc, state);
        gui.repaint();
    }

    @Override
    public String getTip() {
        return "Diagnostics view for monitoring behaviour of MCMC.";
    }

    @Override
    public void setSampling(boolean enabled) {
        // Ignore this. We always want to be enabled!
        // Uncle Sam can't kill us!
    }

    public static class DiagnosticsGUI extends JPanel {
        private double logLike = Double.NEGATIVE_INFINITY;

        private int alignmentSamples = 0;
        private double alignmentRatio = Double.NaN;
        private int edgeSamples = 0;
        private double edgeRatio = Double.NaN;
        private int topologySamples = 0;
        private double topologyRatio = Double.NaN;
        private int indelSamples = 0;
        private double indelRatio = Double.NaN;
        private int substSamples = 0;
        private double substRatio = Double.NaN;

        @Override
        public void paintComponent(Graphics g) {
            Graphics2D graphics = (Graphics2D) g;
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, getWidth(), getHeight());

            final Font titleFont = new Font("Andale Mono", Font.BOLD, 36);
            final Font textFont = new Font("Andale Mono", Font.PLAIN, 20);

            graphics.setColor(Color.BLACK);
            graphics.setFont(titleFont);
            graphics.drawString("General", 10, 40);

            graphics.setFont(textFont);
            graphics.drawString("logLike: " + logLike, 10, 70);

            graphics.setFont(titleFont);
            graphics.drawString("MCMC Sample Ratios", 10, 140);

            graphics.setFont(textFont);
            drawInfoLine(graphics, "Alignment", 170, alignmentRatio, alignmentSamples);
            drawInfoLine(graphics, "Edge", 200, edgeRatio, edgeSamples);
            drawInfoLine(graphics, "Topology", 230, topologyRatio, topologySamples);
            drawInfoLine(graphics, "Indel", 260, indelRatio, indelSamples);
            drawInfoLine(graphics, "Subst", 290, substRatio, substSamples);
        }

        private void drawInfoLine(Graphics2D graphics, String text, int y, double ratio, int total) {
            graphics.drawString(text + ":", 10, y);
            graphics.drawString("" + ratio, 170, y);
            graphics.drawString("Total samples: " + total, 450, y);
        }

        public void update(Mcmc mcmc, State state) {
            logLike = state.logLike;

            alignmentSamples = mcmc.alignmentSampled;
            alignmentRatio = (double) mcmc.alignmentAccepted / mcmc.alignmentSampled;

            edgeSamples = mcmc.edgeSampled;
            edgeRatio = (double) mcmc.edgeAccepted / mcmc.edgeSampled;

            topologySamples = mcmc.topologySampled;
            topologyRatio = (double) mcmc.topologyAccepted / mcmc.topologySampled;

            indelSamples = mcmc.indelSampled;
            indelRatio = (double) mcmc.indelAccepted / mcmc.indelSampled;

            substSamples = mcmc.substSampled;
            substRatio = (double) mcmc.substAccepted / mcmc.substSampled;
        }
    }
}
