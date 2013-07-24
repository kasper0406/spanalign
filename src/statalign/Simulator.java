package statalign;

import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.utils.NewickParser;

import java.io.FileWriter;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Preeyan
 * Date: 22/07/13
 * Time: 12:37
 * To change this template use File | Settings | File Templates.
 */
public class Simulator {
    private static Random random;

    private static double Alpha, Lambda, Mu, R;

    public static void main(String[] args) {
        int seed = 1;
        random = new Random(seed);

        double averageFragmentLength = 1;
        double averageSequenceLength = 100;

        Alpha = 0.02;
        Lambda = 0.01;

        Mu = Lambda / (1 - 1 / (1 + (averageSequenceLength / averageFragmentLength)));
        R = 1 - 1 / (1 + averageFragmentLength);

        NewickParser parser = new NewickParser("(A:1,B:1,(C:1,D:1):1):1;");
        Node root = buildTree(parser.parse());

        root.buildAlignment();

        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter("simulated.fasta");
            String fasta = buildFasta(root);
            fileWriter.append(fasta);
            System.out.print(fasta);
        } catch (Exception e) {
            System.out.println("Something horrible happened.\n" + e);
        } finally {
            if (fileWriter != null) {
                try {
                    fileWriter.close();
                } catch (Exception e) {
                    System.out.println("Something really horrible happened.\n" + e);
                }
            }
        }
    }

    static Node buildTree(TreeNode root) {
        Node result = new Node(root.name, new Model(random, Alpha, Lambda, Mu, R, root.edgeLength));

        for (TreeNode child : root.children) {
            result.getChildren().add(buildTree(child));
        }

        return result;
    }

    static String buildFasta(Node root) {
        StringBuilder fasta = new StringBuilder();

        if (root.getName() != null) {
            fasta.append(">" + root.getName() + "\n");
            fasta.append(root.getSequence() + "\n");
        }

        for (Node child : root.getChildren()) {
            fasta.append(buildFasta(child));
        }

        return fasta.toString();
    }
}

class Model {
    private Random random;

    // Jukes-Cantor
    private double Alpha;

    private double Lambda, Mu, R;
    private double T;

    Model(Random random, double Alpha, double Lambda, double Mu, double R, double T) {
        this.random = random;
        this.Alpha = Alpha;
        this.Lambda = Lambda;
        this.Mu = Mu;
        this.R = R;
        this.T = T;
    }

    private double Beta() {
        double Exp = Math.exp((Lambda - Mu) * T);
        return (1 - Exp) / (Mu - Lambda * Exp);
    }

    char generateNucleotide() {
        // Jukes-Cantor
        return "ACGT".charAt(random.nextInt(4));
    }

    char evolveNucleotide(char parent) {
        // Jukes-Cantor
        if (withProbability((1 + 3 * Math.exp(- 4 * Alpha * T)) / 4)) {
            return parent;
        }
        return "ACGT".replace("" + parent, "").charAt(random.nextInt(3));
    }

    String generateFragment() {
        StringBuilder fragment = new StringBuilder();

        int length = geometric(R);

        for (int i = 1; i <= length; i++) {
            fragment.append(generateNucleotide());
        }

        return fragment.toString();
    }

    String evolveFragment(String parent) {
        StringBuilder child = new StringBuilder();

        for (int i = 0; i < parent.length(); i++) {
            child.append(evolveNucleotide(parent.charAt(i)));
        }

        return child.toString();
    }

    boolean shouldParentFragmentSurvive() {
        return withProbability(Math.exp(- Mu * T));
    }

    int numberOfNewFragmentsIfParentFragmentSurvived() {
        return geometric(Lambda * Beta());
    }

    int numberOfNewFragmentsIfParentFragmentDidNotSurvive() {
        if (withProbability((Mu * Beta()) / (1 - Math.exp(- Mu * T)))) {
            return 0;
        }

        return 1 + numberOfNewFragmentsIfParentFragmentSurvived();
    }

    int numberOfFragmentsAtRootNode() {
        return geometric(Lambda / Mu);
    }

    private int geometric(double p) {
        return (int)Math.floor(Math.log(random.nextDouble()) / Math.log(p));
    }

    private boolean withProbability(double p) {
        return (random.nextDouble() <= p);
    }
}

class Node {
    private List<Node> children;

    private String name;
    private StringBuilder sequence;

    private Model model;

    Node(String name, Model model) {
        this.name = name;
        this.model = model;
        this.children = new ArrayList<Node>();
        sequence = new StringBuilder();
    }

    String getName() {
        return name;
    }

    String getSequence() {
        return sequence.toString();
    }

    List<Node> getChildren() {
        return children;
    }

    void buildAlignment() {
        createInitialFragmentsForChildren();
        createFragments(model.numberOfFragmentsAtRootNode());
    }

    void createInitialFragments() {
        createInitialFragmentsForChildren();
        createFragments(model.numberOfNewFragmentsIfParentFragmentSurvived());
    }

    private void createInitialFragmentsForChildren() {
        for (Node child : children) {
            child.createInitialFragments();
        }

        insertCorrectGaps();
    }

    void parentHasNewFragment(String parentFragment) {
        if (model.shouldParentFragmentSurvive()) {
            inheritFragment(parentFragment);
            createFragments(model.numberOfNewFragmentsIfParentFragmentSurvived());
        } else {
            insertGapRecursively(parentFragment.length());
            createFragments(model.numberOfNewFragmentsIfParentFragmentDidNotSurvive());
        }
    }

    void inheritFragment(String parentFragment) {
        String newFragment = model.evolveFragment(parentFragment);
        insertFragment(newFragment);
    }

    void createFragments(int numberOfFragments) {
        for (int i = 1; i <= numberOfFragments; i++) {
            createFragment();
        }
    }

    void createFragment() {
        String newFragment = model.generateFragment();
        insertFragment(newFragment);
    }

    void insertFragment(String fragment) {
        sequence.append(fragment);

        for (Node child : children) {
            child.parentHasNewFragment(fragment);
        }

        insertCorrectGaps();
    }

    // Before: parent is "A"    and children are "AB"   "AC"   "AD"
    // After:  parent is "A---" and children are "AB--" "A-C-" "A--D"
    private void insertCorrectGaps() {
        int totalSurplusLength = 0;
        Map<Node, Integer> surplusLengthByNode = new HashMap<Node, Integer>();

        for (Node child : children) {
            int surplusLength = child.sequence.length() - sequence.length();
            totalSurplusLength += surplusLength;
            surplusLengthByNode.put(child, surplusLength);
        }

        int lengthToInsertBefore = 0;
        int lengthToInsertAfter = totalSurplusLength;

        for (Node child : children) {
            lengthToInsertAfter -= surplusLengthByNode.get(child);
            child.insertGapRecursively(lengthToInsertBefore, sequence.length());
            child.insertGapRecursively(lengthToInsertAfter);
            lengthToInsertBefore += surplusLengthByNode.get(child);
        }

        insertGap(totalSurplusLength);
    }

    void insertGap(int length) {
        for (int i = 1; i <= length; i++) {
            sequence.append('-');
        }
    }

    void insertGap(int length, int offset) {
        for (int i = 1; i <= length; i++) {
            sequence.insert(offset, '-');
        }
    }

    void insertGapRecursively(int length) {
        insertGap(length);

        for (Node child : children) {
            child.insertGapRecursively(length);
        }
    }

    void insertGapRecursively(int length, int offset) {
        insertGap(length, offset);

        for (Node child : children) {
            child.insertGapRecursively(length, offset);
        }
    }
}