package statalign.base;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/17/13
 * Time: 11:10 AM
 * To change this template use File | Settings | File Templates.
 */
public class Hypernode {
    static final int ROUNDING = 100; // tells the number of digits in rounding when printing the tree
    static final double SELECTING = 0.5; /*this is the probability for selecting a
    										homologous column not to be changed during topology changing
	 */

    static final double EMPTY_WINDOW = 0.01; /* this is the probability that an empty window will be realigned*/

    Tree owner;

    Hypernode old;

    /**
     * The name of the sequence associated to the vertex.
     * Used only if the vertex is a leaf of the tree.
     */
    public String name;

    /** This reference points to the parent of the vertex */
    public Hypernode parent;
    /**
     * This reference points to the left child of the vertex. If the
     * vertex is a leaf, it is set to null.
     */
    public Hypernode left;
    /**
     * This reference points to the right child of the vertex. If the
     * vertex is a leaf, it is set to null.
     */
    public Hypernode right;

    public List<Hypernode> children;



    int length;					// sequence length
    AlignColumn first;			// first alignment column of Vertex
    AlignColumn last;			// last, virtual alignment column of Vertex (always present)
    String seq;					// original sequence of this Vertex (given for leaves only)

    public boolean isLabeled;           // flag determining if the Hypernode is nor not labeled.

    int winLength;                    // length of window
    AlignColumn winFirst;        // first alignment column of window
    AlignColumn winLast;        // first alignment column past window end
    boolean selected;                // shows if vertex is part of the selected subtree

    /** The length of the edge that connects this vertex with its parent. */
    public double edgeLength;                            // length of edge to parent vertex
    double[][] charTransMatrix;            // precalculated character transition likelihoods (subst. model)
    double[][] charPropTransMatrix;        // precalculated character transition likelihoods for proposals (subst. model)
    double[][] hmm2TransMatrix;            // precalculated state transition likelihoods for 2-seq HMM (indel model)
    double[][] hmm2PropTransMatrix;        // precalculated state transition likelihoods for 2-seq HMM used for proposals (indel model)
    double[][] hmm3TransMatrix;            // precalculated state transition likelihoods for 3-seq HMM (indel model)
    double[][] hmm3RedTransMatrix;    // precalculated st. trans. likelihoods for 3-seq HMM, silent st. removed (indel model)

    /**
     * The log-sum of the Felsenstein's likelihoods of characters that are inserted into the
     * sequence of this vertex.
     */
    public double orphanLogLike;        // log-sum of the likelihood of each orphan column in subtree (incl. this vertex)
    /**
     * The log-sum of the cumulative insertion-deletion loglikelihoods up to this vertex (ie. summed over the
     * subtree below this vertex.).
     */
    public double indelLogLike;

    public int leafCount;


    ////////////           CONSTRUCTORS       ////////////////////////////////////////////

    Hypernode() {
    }



    Hypernode(String descriptor, Hypernode parent) {

        /* descriptor is a string describing the Hypernodes below this Hypernode */

        //System.out.println(descriptor);

        this.parent = parent;
        owner = parent.owner;
        old = new Hypernode();


        if (descriptor.charAt(0) == '(') {
            String leftDescriptor = "";
            int counter = 0;
            int j;
            for (j = 1; counter != 0 || descriptor.charAt(j) != ','; j++) {
                leftDescriptor += descriptor.charAt(j);
                if (descriptor.charAt(j) == '(') {
                    counter++;
                }
                if (descriptor.charAt(j) == ')') {
                    counter--;
                }
            }
            Hypernode leftChild = new Hypernode(leftDescriptor, this);
            left = leftChild;

            String rightDescriptor = "";
            for (j += 1; counter != -1; j++) {
                rightDescriptor += descriptor.charAt(j);
                if (descriptor.charAt(j) == '(') {
                    counter++;
                }
                if (descriptor.charAt(j) == ')') {
                    counter--;
                }
            }
            rightDescriptor = rightDescriptor.substring(0, rightDescriptor.length() - 1);
            Hypernode rightChild = new Hypernode(rightDescriptor, this);
            right = rightChild;
            String tempStringlength = "";
            for (j += 1; j < descriptor.length(); j++) {
                tempStringlength += descriptor.charAt(j);
            }
            //System.out.println(tempStringlength);
            if (tempStringlength.length() > 0) {
                edgeLength = Math.max(Double.parseDouble(tempStringlength), 0.01);
            } else {
                edgeLength = 0.01;
            }
            //  edgeLength = 0.01;
            //else we have the root Hypernode, which does not need an edgeLength
        } else {
            left = null;
            right = null;
            name = "";
            int i = 0;
            while (descriptor.charAt(i) != ':') {
                name += descriptor.charAt(i);
                i++;
            }
            String tempString = "";
            for (i++; i < descriptor.length(); i++) {
                tempString += descriptor.charAt(i);
            }
            //System.out.println(tempString);
            edgeLength = Math.max(Double.parseDouble(tempString), 0.01);

        }

    }

    boolean isLabeled() {
        return this.isLabeled;
    }

    /**
     * This function returns the number of leaves that are below this vertex.
     * @return the number of leaves that are below this vertex.
     */
    int countLeaves() {
        return leafCount = (left == null && right == null) ? 1 : left.countLeaves() + right.countLeaves();
    }


    String print() {

        // print this Hypernode and Hypernodes below in bracket notation
        DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat("0.#####", dfs);

        if (left == null && right == null) {
            String x = df.format(edgeLength, new StringBuffer(), new FieldPosition(1)).toString();

            return name.replaceAll(" ", "") + ":" +
                    x;//.substring(0,Math.min(x.length(),ROUNDING));
        } else {
            if (parent != null) {
                //	String x = Double.toString(edgeLength);
                String x = df.format(edgeLength, new StringBuffer(), new FieldPosition(1)).toString();
                return "(" + left.print() + "," + right.print() + "):" +
                        x;//.substring(0,Math.min(x.length(),ROUNDING));
            } else {
                return "(" + left.print() + "," + right.print() + ");";
            }
        }
    }

    String print(int digits) {

        String pattern = "0.";
        for (int i = 0; i < digits; i++) {
            pattern += "#";
        }
        // print this Hypernode and Hypernodes below in bracket notation
        DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat(pattern, dfs);

        if (left == null && right == null) {
            String x = df.format(edgeLength, new StringBuffer(), new FieldPosition(1)).toString();

            return name.replaceAll(" ", "") + (digits == 0 ? "" : ":" + x);
        } else {
            if (parent != null) {
                //	String x = Double.toString(edgeLength);
                String x = df.format(edgeLength, new StringBuffer(), new FieldPosition(1)).toString();
                return "(" + left.print(digits) + "," + right.print(digits) + ")" + (digits == 0 ? "" : ":" + x);
            } else {
                return "(" + left.print(digits) + "," + right.print(digits) + ");";
            }
        }
    }

    /** this function set selected to false for all AlignColumns */
    void setAllAlignColumnsUnselected() {
        AlignColumn a = first;
        while (a != null) {
            a.selected = false;
            a = a.next;
        }
    }


    void parentNewChild(Hypernode child) {
        child.last.parent = parent.last;
        if (parent.left == this) {
            parent.left = child;
            parent.last.left = child.last;
        } else {
            parent.right = child;
            parent.last.right = child.last;
        }
    }


    void printPointers() {
        AlignColumn c = first;
        AlignColumn p = parent.first;
        while (p != null) {
            System.out.print(p + " ");
            p = p.next;
        }
        System.out.println();
        while (c != null) {
            System.out.print(c.parent + " ");
            c = c.next;
        }
        System.out.println("\n");

        c = first;
        p = parent.first;
        while (p.next != null) {
            System.out.print(p.mostLikely() + "");
            p = p.next;
        }
        System.out.println();
        while (c != null) {
            if (c.parent.next != null) {
                System.out.print(c.parent.mostLikely() + "");
            }
            c = c.next;
        }
        System.out.println();
        c = first;
        while (c != null) {
            if (c.parent.next != null) {
                System.out.print((c.orphan ? "y" : "n"));
            }
            c = c.next;
        }
        System.out.println("\n");


    }

    /** Returns most likely sequence of this vertex or the original sequence for leaves. */
    String sequence() {
        if(seq != null)
            return seq;
        StringBuilder s = new StringBuilder();
        for (AlignColumn a = first; a != last; a = a.next)
            s.append(a.mostLikely());
        return s.toString();
    }


    /**
     * Calculate the sum of branch lengths below on this vertex.
     * Used to generate the prior for the whole state.
     * @return sum of branchlengths below this Hypernode.
     */
    public double calcSumOfEdges() {
        if (left == null) {
            return edgeLength;
        } else return edgeLength + left.calcSumOfEdges() + right.calcSumOfEdges();
    }


    /**
     * Calculate the maximum depth below on this vertex.
     * Used in tree visualisation.
     * @return maximum depth below this vertex.
     */
    public double maxDepth() {
        if (left == null) {
            return edgeLength;
        } else return edgeLength + Math.max(left.maxDepth(), right.maxDepth());
    }


    /** this function checks if the pointers are all right... */
    void checkPointers() {
        //parent
        for (AlignColumn p = first; p != null; p = p.next) {
            if (p.left != null && (p.left.orphan || p.left.parent != p)) {
                throw new Error("Problem is vertex " + this + ":\np is: " + p + " p.left is " + p.left + " p.left.orphan: " + p.left.orphan +
                        " p.left.parent: " + p.left.parent);
            }
            if (p.right != null && (p.right.orphan || p.right.parent != p)) {
                throw new Error("Problem is vertex " + this + ":\np is: " + p + " p.right is " + p.right + " p.right.orphan: " + p.right.orphan +
                        " p.right.parent: " + p.right.parent);
            }
        }
        for (AlignColumn l = left.first; l != null; l = l.next) {
            if (!l.orphan) {
                if (l.parent == null || l.parent.left != l) {
                    throw new Error("Problem in vertex " + this + ":\nl is: " + l + (l.parent == null ? " l does not have a parent" : " l parent is: " + l.parent + " l parent left is: " + l.parent.left));
                }
            }
        }
        for (AlignColumn r = right.first; r != null; r = r.next) {
            if (!r.orphan) {
                if (r.parent == null || r.parent.right != r) {
                    throw new Error("Problem in vertex " + this + ":\nr is: " + r + (r.parent == null ? " r does not have a parent" : " r parent is: " + r.parent + " r parent right is: " + r.parent.right));
                }
            }
        }
    }


//    /**
//     * This function is merely for testing/debugging purposes
//     * @param args Arguments are not used, all input is directly written into the
//     *             function.
//     */
//    public static void main(String[] args) throws IOException {
//        try {
//            Tree tree = new Tree(new String[]{"qqqqqqqqqqqqkkkwwwwwlidwwwwwkkk",
//                    "kkkwwwwwlidwwwwwkkk",
//                    "qqqqqqqqqqqqkkkwwwwwlidwwwwwkkk",
//                    "kkkwwwwwlidwwwwwkkkeqeq",
//                    "kkkwwwwwlidwwwwwkkkddkldkl",
//                    "kkkwwwwwlidwwwwwkkkeqiqii",
//                    "kkkwwwwwlidwwwwwkkkddkidkil",
//                    "kkkwwwwwlidwwwwwkkkeqiq",
//                    "kkkwwwwwlidwwwwwkkkddkldkll",
//                    "kkkwwwwwlidwwwwwkkkddkldkil"},
//                    new String[]{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
//                    new Dayhoff(),
//                    new Blosum62(), "");
//            System.out.println(tree.printedTree());
//            tree.root.calcFelsRecursively();
//            System.out.println(tree.root.orphanLogLike);
//
////			String[] s = tree.printedAlignment("StatAlign");
//            /*String[] s = null;
//            for (int i = 0; i < s.length; i++) {
//                System.out.println(s[i]);
//            }*/
//
//            /////////////////////////////////
//            double[] weights = new double[tree.vertex.length];
//            tree.countLeaves(); // calculates recursively how many leaves we have below this Hypernode
//            for (int j = 0; j < 100; j++) {
//                for (int i = 0; i < weights.length; i++) {
//                    weights[i] = Math.pow(tree.vertex[i].leafCount, Mcmc.LEAFCOUNT_POWER);
//                }
//                int k = Utils.weightedChoose(weights, null);
//                tree.vertex[k].selectSubtree(Mcmc.SELTRLEVPROB, 0);
//                System.out.println("Selected vertices: ");
//                for (int i = 0; i < tree.vertex.length; i++) {
//                    if (tree.vertex[i].selected) {
//                        System.out.print(i + " ");
//                        tree.vertex[i].selected = false;
//                    }
//                }
//                System.out.println("\n");
//
//            }
//        } catch (StoppedException e) {
//            // stopped during tree construction
//        }
//
//    }




}
