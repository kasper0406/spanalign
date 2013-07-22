package statalign.postprocess.plugins;

import statalign.base.State;
import statalign.postprocess.Postprocess;
import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.utils.NewickParser;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;

public class SpannoidViewer extends Postprocess {
    private JPanel panel = new SpannoidViewerGUI();

    @Override
    public String getTabName() {
        return "Spannoid Viewer";
    }

    @Override
    public Icon getIcon() {
        return new ImageIcon(ClassLoader.getSystemResource("icons/tree.gif"));
    }

    @Override
    public JPanel getJPanel() {
        displayImage();
        return panel;
    }

    @Override
    public String getTip() {
        return "Displays a spannoid.";
    }

    @Override
    public void setSampling(boolean enabled) {
        // Always sample..
    }

    public void newSample(State state, int no, int total) {
        PrintWriter out = null;
        try {
            out = new PrintWriter("spannoid.dot", "UTF-8");
            writeSpannoidToFile(state, out);
            out.close();

            File curDir = new File(".");
            Process p = Runtime.getRuntime().exec("/opt/local/bin/dot -Tjpg -o"+ curDir.getAbsolutePath() +"/spannoid.jpg < "+ curDir.getAbsolutePath() +"/spannoid.dot");
            p.waitFor();
            displayImage();
        } catch (Exception e) {
            e.printStackTrace();
            // Something went wrong.
            // We don't care!
        }
    }

    private void writeSpannoidToFile(State state, PrintWriter out) {
        String tree = state.getNewickString();
        NewickParser parser = new NewickParser(tree);
        TreeNode root = parser.parse();

        out.write("graph G {\n");

        Queue<TreeNode> queue = new LinkedList<TreeNode>();
        queue.add(root);

        int i = 0;
        Map<TreeNode, Integer> lookup = new HashMap<TreeNode, Integer>();
        while (!queue.isEmpty()) {
            TreeNode current = queue.poll();
            int currentId = i++;
            lookup.put(current, currentId);

            if (current.name == null) {
                // Steiner
                out.write("\"N" + currentId + "\" [label=\"\"];\n");
            } else {
                // Labeled
                out.write("\"N" + currentId + "\" [label=\""+ current.name +"\",fillcolor=black,fontcolor=white,style=\"filled\"];\n");
            }

            if (current.parent != null) {
                int parentId = lookup.get(current.parent);
                out.write("\"N"+ parentId +"\" -- \"N"+ currentId +"\" [label=\"" + current.edgeLength + "\"];\n");
            }

            for (TreeNode child : current.children)
                queue.add(child);
        }

        out.write("}\n");
    }

    private void displayImage() {
        panel.repaint();
    }

    public static class SpannoidViewerGUI extends JPanel {
        @Override
        public void paintComponent(Graphics graphics) {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, getWidth(), getHeight());

            BufferedImage image = null;
            try {
                image = ImageIO.read(new File("spannoid.jpg"));
                graphics.drawImage(image, 0, 0, null);
            } catch (IOException e) {
            }
        }
    }
}
