package statalign;

import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.model.subst.plugins.JukesCantor;
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
    private static SubstitutionModel substitutionModel;

    private static double Lambda, Mu, R;
    private static double defaultEdgeLength;

    public static void main(String[] args) {
        int seed = 1;
        String tree = "(((A,B),C),((D,E),F))G;";
        String file = "simulated.fasta";

        random = new Random(seed);
        try {
            substitutionModel = new Dayhoff();
        } catch (Exception e) {
            System.out.println("Something fairly horrible happened.\n" + e);
        }

        double averageFragmentLength = 1;
        double averageSequenceLength = 100;

        defaultEdgeLength = 0.2;
        Lambda = 0.05;
        Mu = Lambda / (1 - 1 / (1 + (averageSequenceLength / averageFragmentLength)));
        R = 1 - 1 / (1 + averageFragmentLength);

        NewickParser parser = new NewickParser(tree);
        Node root = buildTree(parser.parse());

        root.buildAlignment();

        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(file);
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
        double edgeLength = root.edgeLength;
        if (edgeLength == 0) {
            edgeLength = defaultEdgeLength;
        }

        Node result = new Node(root.name, new Model(random, substitutionModel, Lambda, Mu, R, edgeLength));

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

    private SubstitutionModel substitutionModel;
    private double[][] substitutionProb;

    private double Lambda, Mu, R;
    private double T;

    Model(Random random, SubstitutionModel substitutionModel, double Lambda, double Mu, double R, double T) {
        this.random = random;

        this.substitutionModel = substitutionModel;
        substitutionProb = substitutionModel.updateTransitionMatrix(null, T);

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
        return substitutionModel.alphabet[choose(substitutionModel.e)];
    }

    char evolveNucleotide(char parent) {
        int i;
        for (i = 0; i < substitutionModel.alphabet.length; i++) {
            if (substitutionModel.alphabet[i] == parent) {
                break;
            }
        }
        return substitutionModel.alphabet[choose(substitutionProb[i])];
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

    private int choose(double[] p) {
        double u = random.nextDouble();
        int choice;
        for (choice = 0; (u -= p[choice]) > 0; choice++);
        return choice;
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