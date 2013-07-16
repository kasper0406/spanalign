package statalign.base;

import statalign.base.hmm.HmmNonParam;
import statalign.base.hmm.HmmTkf92;
import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.model.score.SubstitutionScore;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Kimura3;
import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.utils.NewickParser;

import java.util.*;

public class Spannoid extends Stoppable implements ITree {
    private Set<Tree> components = new HashSet<Tree>();

    private Map<Vertex, Integer> labeledVertexIds = new HashMap<Vertex, Integer>();

    private int n;

    /*
     * Description...
     */
    private Map<Integer, Set<Vertex>> componentConnections = new HashMap<Integer, Set<Vertex>>();

    public Spannoid(String[] sequences, String[] names,
                    SubstitutionModel model, SubstitutionScore ss)
            throws StoppedException
    {
        throw new RuntimeException("Not yet implemented.");
    }

    /**
     * Constructs a spannoid from a Newick representation and sequences.
     * It is assumed that the Newick representation is rooted at a leaf node.
     * @param newick
     * @param sequences
     * @param nameMap
     * @param model
     * @param ss
     * @throws StoppedException
     */
    public Spannoid(String newick, String[] sequences, Map<String, Integer> nameMap,
                    SubstitutionModel model, SubstitutionScore ss) throws StoppedException
    {
        n = sequences.length;

        NewickParser parser = new NewickParser(newick);
        TreeNode root = parser.parse();

        int[][][] convertedSequences = convertSequences(sequences, model, ss);
        createComponents(root, model, nameMap, convertedSequences, sequences);
    }

    private void createComponents(TreeNode node, SubstitutionModel model,
                                  Map<String, Integer> nameMap, int[][][] sequences, String[] originalSequences) {
        // TODO: Re-root tree at a leaf automatically
        if (node.name == null)
            throw new RuntimeException("The tree should be rooted at a leaf!");

        Queue<TreeNode> labeledNodesToVisit = new LinkedList<TreeNode>();
        labeledNodesToVisit.add(node);

        while (!labeledNodesToVisit.isEmpty()) {
            TreeNode current = labeledNodesToVisit.poll();

            for (TreeNode child : current.children)
                labeledNodesToVisit.addAll(createComponent(current, child, model, nameMap, sequences, originalSequences));
        }
    }

    /**
     * Create a component.
     * @param startLeaf
     * @param next The nodes adjacent to startLeaf which should be visited.
     * @param model
     * @param nameMap
     * @param sequences
     * @param originalSequences
     * @return The leaf nodes of the tree except startLeaf.
     */
    private Set<TreeNode> createComponent(TreeNode startLeaf, TreeNode next, SubstitutionModel model,
                                          Map<String, Integer> nameMap, int[][][] sequences, String[] originalSequences) {
        Set<TreeNode> labeledNodes = new HashSet<TreeNode>();
        List<Vertex> internalVertices = new LinkedList<Vertex>();
        List<Vertex> leafVertices = new LinkedList<Vertex>();

        Map<TreeNode, Vertex> treeNodeToVertex = new HashMap<TreeNode, Vertex>();

        Tree tree = new Tree();
        tree.substitutionModel = model;
        tree.hmm2 = new HmmTkf92(null);
        tree.hmm3 = new HmmNonParam();

        Vertex rootVertex = new Vertex(tree, 0.0);
        internalVertices.add(rootVertex);
        tree.root = rootVertex;

        // Add the new vertices
        Vertex firstLeaf = new Vertex(tree, next.edgeLength / 2, sequences[nameMap.get(startLeaf.name)],
                startLeaf.name, originalSequences[nameMap.get(startLeaf.name)]);
        labeledVertexIds.put(firstLeaf, nameMap.get(nameMap.get(startLeaf.name)));
        leafVertices.add(firstLeaf);
        firstLeaf.parent = rootVertex;
        addChildToVertex(rootVertex, firstLeaf);
        treeNodeToVertex.put(startLeaf, firstLeaf);
        firstLeaf.selected = true;

        // Collect information for building tree
        Queue<TreeNode> queue = new LinkedList<TreeNode>();
        queue.add(next);

        while (!queue.isEmpty()) {
            TreeNode current = queue.poll();

            double edgeLength = current.edgeLength;
            // Special case for fake root vertex. Edge length should be halved
            if (current == next)
                edgeLength /= 2.;

            Vertex newVertex;
            if (current.name == null) {
                // Steiner node
                newVertex = new Vertex(tree, edgeLength);
                internalVertices.add(newVertex);

                // Add new vertices to visit
                for (TreeNode child : current.children)
                    queue.add(child);
            } else {
                // Leaf node
                newVertex = new Vertex(tree, edgeLength, sequences[nameMap.get(current.name)],
                        current.name, originalSequences[nameMap.get(current.name)]);
                leafVertices.add(newVertex);
                labeledVertexIds.put(newVertex, nameMap.get(current.name));

                labeledNodes.add(current);
            }
            newVertex.name = current.name;
            newVertex.selected = true;

            Vertex parent = (current == next)
                                ? rootVertex // Special case. Add to fake root.
                                : treeNodeToVertex.get(current.parent);
            newVertex.parent = parent;
            addChildToVertex(parent, newVertex);

            treeNodeToVertex.put(current, newVertex);
        }

        tree.vertex = new Vertex[leafVertices.size() + internalVertices.size()];
        tree.names = new String[leafVertices.size()];
        int i = 0;
        for (Vertex vertex : leafVertices) {
            // Do full alignment!
            vertex.fullWin();

            // Ensure transition matrices are updated
            vertex.edgeChangeUpdate();

            tree.names[i] = vertex.name;
            tree.vertex[i] = vertex;
            i++;
        }

        // Iterate backwards through internal vertices.
        // Ensures that they are added bottom up.
        ListIterator<Vertex> iterator = internalVertices.listIterator(internalVertices.size());
        while (iterator.hasPrevious()) {
            Vertex vertex = iterator.previous();

            // Ensure transition matrices are updated
            vertex.edgeChangeUpdate();

            tree.vertex[i++] = vertex;

            // Add empty alignment column to internal nodes, and set 
            makeFakeAlignmentColumn(vertex);

            // Do full alignment!
            vertex.fullWin();
        }

        // Align all steiner/internal nodes
        rootVertex.doRecAlign();
        rootVertex.calcFelsRecursively();
        rootVertex.calcIndelLikeRecursively();

        // Add the component
        components.add(tree);

        return labeledNodes;
    }

    private void addChildToVertex(Vertex parent, Vertex child) {
        child.parent = parent;
        if (parent.left == null)
            parent.left = child;
        else if (parent.right == null)
            parent.right = child;
        else
            throw new RuntimeException("The child is full!");
    }

    private void makeFakeAlignmentColumn(Vertex vertex) {
        AlignColumn fake = new AlignColumn(vertex);
        vertex.first = fake;
        vertex.last = fake;
        vertex.length = 0;
        fake.left = vertex.left.last;
        fake.right = vertex.right.last;
        vertex.left.last.parent = fake;
        vertex.left.last.orphan = false;
        vertex.right.last.parent = fake;
        vertex.right.last.orphan = false;
    }

    // TODO: Is it necessary to pass along a substitution score?
    //       Couldn't it be obtained from the substitution model?
    /**
     * Convert the string input sequences to an array of Felchenstein probabilities.
     *   st. for a given column i, row j is 1 iff seq[i] corresponds to symbol j in
     *   the substitution model.
     * @param sequences The string representation of the sequence.
     * @param model The substitution model used.
     * @param ss The substitution score.
     * @return
     */
    private int[][][] convertSequences(String[] sequences, SubstitutionModel model, SubstitutionScore ss) {
        int[][][] seq = new int[sequences.length][][];
        for (int i = 0; i < sequences.length; i++) {
            int k = 0;
            for (int j = 0; j < sequences[i].length(); j++) {
                int sum = 0;
                char ch = sequences[i].charAt(j);
                for (int l = 0; l < ss.which[ch].length; l++) {
                    sum += ss.which[ch][l];
                }
                if (sum > 0) {
                    k++;
                }
            }
            seq[i] = new int[k][model.e.length];
            k = 0;
            for (int j = 0; j < sequences[i].length(); j++) {
                int sum = 0;
                char ch = sequences[i].charAt(j);
                for (int l = 0; l < ss.which[ch].length; l++) {
                    sum += ss.which[ch][l];
                }
                if (sum > 0) {
                    seq[i][k] = ss.which[ch];
                    k++;
                }
            }
        }
        return seq;
    }

    public double getLogLike() {
        double logLike = 0;
        for (Tree component : components)
            logLike += component.getLogLike();
        return logLike;
    }

    public State getState() {
        Tree firstComponent = components.iterator().next();

        int nodes = countNodes() - components.size() + 1;
        State state = new State(nodes);

        ArrayList<LinkedList<Integer>> children = new ArrayList<LinkedList<Integer>>(nodes);
        for (int i = 0; i < nodes; i++)
            children.add(new LinkedList<Integer>());

        Map<Vertex, Integer> lookup = new HashMap<Vertex, Integer>();

        // int steinerIndex = n;
        int leafCounter = 0;
        int internalCounter = n - components.size() + 1;

        Queue<Vertex> queue = new LinkedList<Vertex>();
        queue.add(firstComponent.vertex[0]); // Add node to start BFS traversal

        while (!queue.isEmpty()) {
            Vertex current = queue.poll();

            boolean visitedBefore = lookup.containsKey(current); // This happens when searching connected components
            if (!visitedBefore) {
                boolean isInternal = !labeledVertexIds.containsKey(current)
                        || componentConnections.get(labeledVertexIds.get(current)).size() > 1;
                if (isInternal)
                    lookup.put(current, internalCounter++);
                else
                    lookup.put(current, leafCounter++);
            }

            if (current.parent != null) {
                if (lookup.containsKey(current.parent)) {
                    state.parent[lookup.get(current)] = lookup.get(current.parent);
                    children.get(lookup.get(current.parent)).add(lookup.get(current));
                } else
                    queue.add(current.parent);
            }
            if (current.left != null) {
                if (lookup.containsKey(current.left)) {
                    state.parent[lookup.get(current)] = lookup.get(current.left);
                    children.get(lookup.get(current.left)).add(lookup.get(current));
                } else
                    queue.add(current.left);
            }
            if (current.right != null) {
                if (lookup.containsKey(current.right)) {
                    state.parent[lookup.get(current)] = lookup.get(current.right);
                    children.get(lookup.get(current.right)).add(lookup.get(current));
                } else
                    queue.add(current.right);
            }

            if (labeledVertexIds.containsKey(current)) {
                // If labeled nodes, look for new connected components
                Set<Vertex> connectedVertices = componentConnections.get(labeledVertexIds.get(current));
                for (Vertex connection : connectedVertices) {
                    if (!lookup.containsKey(connection)) {
                        queue.add(connection);
                        lookup.put(connection, lookup.get(current));
                    }
                }
            }

            state.name[lookup.get(current)] = (current.name == null) ? "" : current.name;
            state.seq[lookup.get(current)] = current.sequence();
        }

        state.parent[lookup.get(firstComponent.vertex[0])] = -1;
        state.root = lookup.get(firstComponent.vertex[0]);

        for (int i = 0; i < nodes; i++) {
            int size = children.get(i).size();
            state.children[i] = new int[size];
            int j = 0;
            for (int child : children.get(i))
                state.children[i][j++] = child;
        }

        // TODO: Alignments should be merged!
        //       For now just make dummy alignment.
        for (int i = 0; i < nodes; i++) {
            int parentLength = (state.parent[i] != -1) ? state.seq[state.parent[i]].length() // Find parent sequence length
                                                       : state.seq[i].length(); // We are the root
            state.align[i] = new int[parentLength];
            for (int j = 0; j < parentLength; j++)
                state.align[i][j] = j;
        }

        // TODO: Fix this!
        state.indelParams = firstComponent.hmm2.params.clone();
        state.substParams = firstComponent.substitutionModel.params.clone();
        state.logLike = getLogLike();

        return state;
    }

    private int countNodes() {
        int nodes = 0;
        for (Tree component : components)
            nodes += component.vertex.length;
        return nodes;
    }

    public static void main(String[] args) throws Exception {
        // String tree = "((B:0.5,C:0.5):0.5)A;";
        // String tree = "((A:0.5,B:0.2):1,(D:1,E:0.2):1)C;";
        String tree = "((D:0.1,(F:0.1,(E:0.2,C:0.05):0.2)B:0.1):0.2)A;";
        String[] seqs = new String[] { "AAGT", "CGATTC", "CCGAAG", "AGACA", "TTGACC", "GTAC" };
        Map<String, Integer> nameMap = new HashMap<String, Integer>();
        for (int i = 0; i < seqs.length; i++)
            nameMap.put("" + (char)((int)'A' + i), i);

        SubstitutionModel model = new Kimura3();
        SubstitutionScore ss = model.attachedScoringScheme;

        Spannoid spannoid = new Spannoid(tree, seqs, nameMap, model, ss);
        System.out.println(String.format("Log-like: %f", spannoid.getLogLike()));

        // Print alignments
        for (Tree component : spannoid.components) {
            System.out.println();
            String[] alignment = component.root.printedMultipleAlignment();
            for (String seq : alignment)
                System.out.println(seq);
        }
    }
}
