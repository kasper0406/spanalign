package statalign.postprocess.plugins;

import statalign.base.InputData;
import statalign.base.Mcmc;
import statalign.base.State;
import statalign.postprocess.Postprocess;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

public class DiagnosticsView extends Postprocess {
    DiagnosticsGUI gui = new DiagnosticsGUI();

    private InputData inputData;
    private AlignmentScorer scorer;

    private MpdAlignment mpdAlignment;

    public DiagnosticsView() {
        screenable = true;
        outputable = true;
        postprocessable = true;
        postprocessWrite = true;
    }

    @Override
    public String getTabName() {
        return "Diagnostics";
    }

    @Override
    public String[] getDependences() {
        return new String[] { "statalign.postprocess.plugins.CurrentAlignment", "statalign.postprocess.plugins.MpdAlignment"};
    }

    @Override
    public void refToDependences(Postprocess[] plugins) {
        // curAlig = (CurrentAlignment) plugins[0];
        mpdAlignment = (MpdAlignment) plugins[1];
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
    public String getFileExtension() {
        return "stats";
    }

    static Comparator<String[]> compStringArr = new Comparator<String[]>() {
        @Override
        public int compare(String[] a1, String[] a2) {
            return a1[0].compareTo(a2[0]);
        }};

    @Override
    public void beforeFirstSample(InputData inputData) {
        this.inputData = inputData;

        try {
            if (postprocessWrite)
                outputFile.write("Sample\tn\tk\tLogLike\tMPD.bali\n"); // \tComponents\tAvg.comp. size\tMax comp.size\tMin comp. size\tAvg. edgelen (leaf)\tMin edgelen (leaf)\t#edges<0.1 (leaf)");
        } catch (IOException e) {
            // Ignore..
        }
    }

    @Override
    public void newSample(State state, int no, int total) {
        if (show) {
            gui.update(mcmc, state);
            gui.repaint();
        }

        if (postprocessWrite) {
            final int n = inputData.seqs.size();
            final int k = inputData.pars.componentSize;
            final double logLike = state.logLike;
            final double mpdBali = mpdAlignment.getScore();

            try {
                outputFile.write(String.format("%d\t%d\t%d\t%f\t%f\n", // \t%d\t%f\t%d\t%d\t%f\t%f\t%d",
                        no, n, k, logLike, mpdBali)); // , components, avgCompSize, maxCompSize,
                        // minCompSize, avgEdgeLen, minEdgeLen, shortEdges));
            } catch (IOException e) {
                // Ignore...
            }
        }
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
        private int alignmentAccepted = 0;
        private double alignmentRatio = Double.NaN;
        private int edgeSamples = 0;
        private int edgeAccepted = 0;
        private double edgeRatio = Double.NaN;
        private int topologySamples = 0;
        private int topologyAccepted = 0;
        private double topologyRatio = Double.NaN;
        private int indelSamples = 0;
        private int indelAccepted = 0;
        private double indelRatio = Double.NaN;
        private int substSamples = 0;
        private int substAccepted = 0;
        private double substRatio = Double.NaN;

        private double R = Double.NaN;
        private double lambda = Double.NaN;
        private double mu = Double.NaN;

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
            drawInfoLine(graphics, "Alignment", 170, alignmentRatio, alignmentAccepted, alignmentSamples);
            drawInfoLine(graphics, "Edge", 200, edgeRatio, edgeAccepted, edgeSamples);
            drawInfoLine(graphics, "Topology", 230, topologyRatio, topologyAccepted, topologySamples);
            drawInfoLine(graphics, "Indel", 260, indelRatio, indelAccepted, indelSamples);
            drawInfoLine(graphics, "Subst", 290, substRatio, substAccepted, substSamples);

            graphics.setFont(titleFont);
            graphics.drawString("Model parameters", 10, 360);

            graphics.setFont(textFont);
            graphics.drawString("R:", 10, 390);
            graphics.drawString("" + R, 100, 390);

            graphics.drawString("λ:", 10, 420);
            graphics.drawString("" + lambda, 100, 420);

            graphics.drawString("μ:", 10, 450);
            graphics.drawString("" + mu, 100, 450);
        }

        private void drawInfoLine(Graphics2D graphics, String text, int y, double ratio, int accepted, int total) {
            graphics.drawString(text + ":", 10, y);
            graphics.drawString(String.format("%.3f",ratio), 170, y);
            graphics.drawString("Accepted/Total: " + accepted + "/" + total, 450, y);
        }

        public void update(Mcmc mcmc, State state) {
            logLike = state.logLike;

            alignmentSamples = mcmc.alignmentSampled;
            alignmentAccepted = mcmc.alignmentAccepted;
            alignmentRatio = (double) mcmc.alignmentAccepted / mcmc.alignmentSampled;

            edgeSamples = mcmc.edgeSampled;
            edgeAccepted = mcmc.edgeAccepted;
            edgeRatio = (double) mcmc.edgeAccepted / mcmc.edgeSampled;

            topologySamples = mcmc.topologySampled;
            topologyAccepted = mcmc.topologyAccepted;
            topologyRatio = (double) mcmc.topologyAccepted / mcmc.topologySampled;

            indelSamples = mcmc.indelSampled;
            indelAccepted = mcmc.indelAccepted;
            indelRatio = (double) mcmc.indelAccepted / mcmc.indelSampled;

            substSamples = mcmc.substSampled;
            substAccepted = mcmc.substAccepted;
            substRatio = (double) mcmc.substAccepted / mcmc.substSampled;

            R = state.indelParams[0];
            lambda = state.indelParams[1];
            mu = state.indelParams[2];
        }
    }
}
