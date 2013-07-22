package statalign.postprocess.plugins;

import statalign.base.State;
import statalign.postprocess.Postprocess;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.Scanner;

public class ParameterHistogram extends Postprocess {
    private HistogramGUI gui = new HistogramGUI();

    private StringBuilder RMeasurements = new StringBuilder();
    private StringBuilder lambdaMeasurements = new StringBuilder();
    private StringBuilder muMeasurements = new StringBuilder();

    private int samples = 0;

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
        samples++;

        RMeasurements.append(state.indelParams[0] + "\n");
        lambdaMeasurements.append(state.indelParams[1] + "\n");
        muMeasurements.append(state.indelParams[2] + "\n");
    }

    @Override
    public void afterLastSample() {
        writeAndGenerateHistogram("R", RMeasurements, 50, 1);
        writeAndGenerateHistogram("lambda", lambdaMeasurements, 50, 0.002);
        writeAndGenerateHistogram("mu", muMeasurements, 50, 0.002);

        gui.setReady();
        gui.repaint();
    }

    private void writeAndGenerateHistogram(String name, StringBuilder measurements, double n, double max) {
        OutputStreamWriter writer = null;
        try {
            File outputFile = new File("plots/" + name + ".png");
            File gnuplotFile = new File("plots/gnuplot.gnuplot");

            //Scanner scanner = new Scanner(gnuplot.getInputStream()).useDelimiter("\\n");
            //writer = new OutputStreamWriter(gnuplot.getOutputStream());
            writer = new FileWriter(gnuplotFile);
            writer.write(
                    "n=" + n + " #number of intervals\n" +
                    "max="+ max +" #max value\n" +
                    "min=0. #min value\n" +
                    "width=(max-min)/n #interval width\n" +
                    "#function used to map a value to the intervals\n" +
                    "hist(x,width)=width*floor(x/width)+width/2.0\n" +
                    "set term png size 470,300 font \"Helvetica, 10\" #output terminal and file\n" +
                    "set output \"" + outputFile.getAbsolutePath() + "\"\n" +
                    "#set xrange [min:max]\n" +
                    "set xrange [0:max]\n" +
                    "set yrange [0:]\n" +
                    "set boxwidth width*0.9\n" +
                    "set style fill solid 0.5 #fillstyle\n" +
                    "set tics out nomirror\n" +
                    "set xlabel \"" + name + "\"\n" +
                    "set ylabel \"Frequency\"\n" +
                    "#count and plot\n" +
                    "plot \"-\" u (hist($1,width)):(" + (double) 1 / samples + ") smooth freq w boxes lc rgb\"green\" notitle\n");

            writer.write(measurements + "e");
            writer.close();


            // System.out.println(scanner.next());

            Process gnuplot = Runtime.getRuntime().exec("/opt/local/bin/gnuplot " + gnuplotFile.getAbsolutePath());
            gnuplot.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                if (writer != null)
                    writer.close();
            } catch (IOException e) {
                // Bad luck...
            }
        }
    }

    public static class HistogramGUI extends JPanel
    {
        private boolean isHistogramReady = false;

        public void setReady() {
            isHistogramReady = true;
        }

        @Override
        public void paintComponent(Graphics graphics) {
            Graphics2D g = (Graphics2D) graphics;

            if (isHistogramReady) {
                graphics.setColor(Color.WHITE);
                graphics.fillRect(0, 0, getWidth(), getHeight());

                try {
                    BufferedImage R = ImageIO.read(new File("plots/R.png"));
                    BufferedImage lambda = ImageIO.read(new File("plots/lambda.png"));
                    BufferedImage mu = ImageIO.read(new File("plots/mu.png"));

                    graphics.drawImage(R, 10, 10, null);
                    graphics.drawImage(lambda, 500, 10, null);
                    graphics.drawImage(mu, 10, 340, null);
                } catch (IOException e) {
                }
            } else {
                Font font = new Font("Andale Mono", Font.BOLD, 36);
                g.setFont(font);

                g.drawString("Waiting for data...", 10, 70);
            }
        }
    }
}
