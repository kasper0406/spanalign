package statalign.postprocess.plugins;

import statalign.base.InputData;
import statalign.base.State;
import statalign.postprocess.Postprocess;
import statalign.postprocess.utils.NewickParser;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;

public class SpannoidViewer extends Postprocess {
    private JPanel panel = new JPanel(new BorderLayout());
    private SpannoidViewerGUI spannoid = new SpannoidViewerGUI();

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

    @Override
    public void beforeFirstSample(InputData input) {
        panel.removeAll();
        JScrollPane scrollPane = new JScrollPane(spannoid);
        panel.add(scrollPane);
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
            BufferedImage image = null;
            try {
                image = ImageIO.read(new File("spannoid.jpg"));
            } catch (IOException e) {
            }
            if (image != null) {
                spannoid.setImage(image);
            }
            panel.repaint();
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

        out.write("graph G { rankdir=LR; \n");

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
                String edgeLen = String.format("%.2f", current.edgeLength);
                out.write("\"N"+ parentId +"\" -- \"N"+ currentId +"\" [label=\"" + edgeLen + "\"];\n");
            }

            for (TreeNode child : current.children)
                queue.add(child);
        }

        out.write("}\n");
    }

    public static class SpannoidViewerGUI extends JPanel {
        private BufferedImage image = null;

        public void setImage(BufferedImage image) {
            this.image = image;
        }

        @Override
        public void paintComponent(Graphics graphics) {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, image.getWidth(), image.getHeight());

            graphics.drawImage(image, 0, 0, null);
        }

        @Override
        public Dimension getMinimumSize(){
            return getPreferredSize();
        }

        @Override
        public Dimension getPreferredSize() {
            if (image == null)
                return new Dimension(500, 300);
            return new Dimension(image.getWidth(), image.getHeight());
        }
    }
}
