package statalign.base;

import statalign.base.hmm.HmmNonParam;
import statalign.base.hmm.HmmTkf92;
import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.model.score.SubstitutionScore;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Kimura3;
import statalign.postprocess.plugins.SpannoidViewer;
import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.utils.NewickParser;
import sun.misc.IOUtils;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.StringWriter;
import java.util.*;

public class Spannoid extends Stoppable implements ITree {
    private int n;
    private List<Tree> components = new ArrayList<Tree>();

    private final String BONPHY_PATH = "/Users/kasper0406/Desktop/bonphy/bonphy.py";

    private SubstitutionModel substitutionModel;

    // TODO: Add more options
    public enum BonphyStrategy {
        TOTAL_LENGTH("Total length", 1),
        CONTRACTED("Contracted length", 2),
        INTERNAL_MOVED("Internal nodes moved", 3);

        private int bonphyOptimizationNumber;
        private String name;

        private BonphyStrategy(String name, int bonphyOptimizationNumber) {
            this.name = name;
            this.bonphyOptimizationNumber = bonphyOptimizationNumber;
        }

        public String getOptimizationNumber() {
            return Integer.toString(bonphyOptimizationNumber);
        }

        public String toString() {
            return name;
        }
    };

    /*
     * Description...
     */
    private Map<Integer, Set<Vertex>> componentConnections = new HashMap<Integer, Set<Vertex>>();
    private Map<Vertex, Integer> labeledVertexIds = new HashMap<Vertex, Integer>();
    private List<Integer> innerBlackNodes = new ArrayList<Integer>();



    private double heat = 1.0d;

    public Spannoid(int componentSize, BonphyStrategy optimizationStrategy,
                    String[] sequences, String[] names,
                    SubstitutionModel model, SubstitutionScore ss)
            throws StoppedException, IOException, InterruptedException {
        this.substitutionModel = model;
        n = sequences.length;
        for (int i = 0; i < n; i++)
            componentConnections.put(i, new HashSet<Vertex>());

        int[][][] convertedSequences = convertSequences(sequences, model, ss);
        String njTree = NJTree.getNJTree(convertedSequences, names, ss);

        Map<String, Integer> nameMap = new HashMap<String, Integer>();
        for (int i = 0; i < names.length; i++)
            nameMap.put(names[i], i);

        Process bonphy = Runtime.getRuntime().exec(BONPHY_PATH + " -k" + componentSize
                + " -s" + optimizationStrategy.getOptimizationNumber());
        OutputStreamWriter output = new OutputStreamWriter(bonphy.getOutputStream());
        output.write(njTree);
        output.close();

        bonphy.waitFor();
        Scanner scanner = new Scanner(bonphy.getInputStream()).useDelimiter("\\n");
        scanner.next(); // Skip first line
        while (scanner.hasNext()) {
            String newickComponent = scanner.next();
            NewickParser parser = new NewickParser(newickComponent);
            TreeNode root = parser.parse();
            root = shrinkTree(root);
            createComponents(root, model, nameMap, convertedSequences, sequences);
        }
        scanner.close();
        setupInnerBlackNodes();
    }

    /**
     * Removes all Steiner nodes of degree 2 from the tree.
     */
    /*
    private TreeNode shrinkTree(TreeNode node) {
        // TODO: This seems overly complicated, since it seems only root node can have the problem.
        //       Should be made simpler!
        //       If it turns out we need to do it for all nodes, this should be generalized instead.
        if (node.children.size() == 2 && node.parent == null) {
            if (node.name != null && !node.name.isEmpty())
                throw new RuntimeException("Invalid format of component!");

            TreeNode left = shrinkTree(node.children.get(0));
            TreeNode right = shrinkTree(node.children.get(1));
            left.parent = null;
            left.addChild(right);
            right.edgeLength = left.edgeLength + right.edgeLength;
            left.edgeLength = node.edgeLength;
            return left;
        }
        for (int i = 0; i < node.children.size(); i++)
            node.children.set(i, shrinkTree(node.children.get(i)));

        return node;
    }
    */

    private TreeNode shrinkTree(TreeNode node) {
        if (node.children.size() == 2 && node.parent == null) {
            if (node.name != null && !node.name.isEmpty())
                throw new RuntimeException("Invalid format of component!");

            TreeNode left = node.children.get(0);
            TreeNode right = node.children.get(1);
            left.parent = null;
            left.addChild(right);
            right.edgeLength = left.edgeLength + right.edgeLength;
            left.edgeLength = node.edgeLength;
            return left;
        }
        return node;
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
        this.substitutionModel = model;
        constructFromNewick(newick, sequences, nameMap, model, ss);
        setupInnerBlackNodes();
    }

    private void constructFromNewick(String newick, String[] sequences, Map<String, Integer> nameMap,
                                     SubstitutionModel model, SubstitutionScore ss) {
        n = sequences.length;
        for (int i = 0; i < n; i++)
            componentConnections.put(i, new HashSet<Vertex>());

        NewickParser parser = new NewickParser(newick);
        TreeNode root = parser.parse();

        int[][][] convertedSequences = convertSequences(sequences, model, ss);
        createComponents(root, model, nameMap, convertedSequences, sequences);

    }

    private void createComponents(TreeNode node, SubstitutionModel model,
                                  Map<String, Integer> nameMap, int[][][] sequences, String[] originalSequences) {
        node = node.rootAtLeaf();

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
        //Vertex firstLeaf = new Vertex(tree, next.edgeLength, sequences[nameMap.get(startLeaf.name)],
        //        startLeaf.name, originalSequences[nameMap.get(startLeaf.name)]);
        labeledVertexIds.put(firstLeaf, nameMap.get(startLeaf.name));
        componentConnections.get(nameMap.get(startLeaf.name)).add(firstLeaf);
        leafVertices.add(firstLeaf);
        addChildToVertex(rootVertex, firstLeaf);
        firstLeaf.parent = rootVertex;
        treeNodeToVertex.put(startLeaf, firstLeaf);
        firstLeaf.selected = true;

        // Collect information for building tree
        Queue<TreeNode> queue = new LinkedList<TreeNode>();
        queue.add(next);

        while (!queue.isEmpty()) {
            TreeNode current = queue.poll();

            double edgeLength = current.edgeLength;

            // Special case for fake root vertex. Edge length should be halved
            if (current == next) {
                edgeLength /= 2;
                // edgeLength = 0;
            }

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
                componentConnections.get(nameMap.get(current.name)).add(newVertex);

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

        tree.vertex = new ArrayList<Vertex>(Collections.<Vertex>nCopies(leafVertices.size() + internalVertices.size(), null));
        //tree.names = new ArrayList<String>(Collections.<String>nCopies(leafVertices.size(), null));
        int i = 0;
        for (Vertex vertex : leafVertices) {
            // Do full alignment!
            vertex.fullWin();

            // Ensure transition matrices are updated
            vertex.edgeChangeUpdate();

            //tree.names.set(i, vertex.name);
            tree.vertex.set(i, vertex);
            i++;
        }

        // Iterate backwards through internal vertices.
        // Ensures that they are added bottom up.
        ListIterator<Vertex> iterator = internalVertices.listIterator(internalVertices.size());
        while (iterator.hasPrevious()) {
            Vertex vertex = iterator.previous();

            // Ensure transition matrices are updated
            vertex.edgeChangeUpdate();

            tree.vertex.set(i++, vertex);

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

    private void setupInnerBlackNodes(){
        innerBlackNodes = new ArrayList<Integer>();
        /// loop through all the nodes and build the inner black nodes list
        for (int i = 0; i < n; i++){
            if(componentConnections.get(i).size() > 1){
                innerBlackNodes.add(i);

            }

        }

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
        if (vertex.right != null) {
            fake.right = vertex.right.last;
            vertex.right.last.parent = fake;
            vertex.right.last.orphan = false;
        }
        if (vertex.left != null) {
            fake.left = vertex.left.last;
            vertex.left.last.parent = fake;
            vertex.left.last.orphan = false;
        }
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

    private double probOfSequence(Vertex vertex) {
        double r = vertex.owner.hmm2.params[0];
        double lambda = vertex.owner.hmm2.params[1];
        double mu = vertex.owner.hmm2.params[2];
        String sequence = vertex.sequence();
        final int n = sequence.length();

        double prob = 0;
        /*
        if (n == 0)
            prob = Math.log(1 - lambda / mu);
        else
            prob = Math.log(1 - lambda / mu) + Math.log(lambda / mu)
                + Math.log(1 - r)
                + (n - 1) * (Math.log((lambda / mu) * (1 - r) + r));
        */

        // TODO: Consider generalizing this!
        prob += n * Math.log((double)1 / 4);

        return prob;
    }

    public double getLogLike() {
        double logLike = 0;
        for (Tree component : components)
            logLike += component.getLogLike();

        for (Set<Vertex> connections : componentConnections.values()) {
            Vertex vertex = connections.iterator().next();
            logLike -= probOfSequence(vertex) * (connections.size() - 1);
        }

        return logLike;
    }

    @Override
    public double getLogPrior() {
        double edgeSum = 0.0;
        for (Tree component : components)
            edgeSum += component.root.calcSumOfEdges();

        return - edgeSum - Math.log(substitutionModel.getPrior())
                - getLambda() - getMu();
    }

    @Override
    public double getOrphanLogLike() {
        double orphanLogLike = 0;
        for (Tree component : components)
            orphanLogLike += component.getOrphanLogLike();

        for (Set<Vertex> connections : componentConnections.values()) {
            Vertex vertex = connections.iterator().next();
            orphanLogLike -= probOfSequence(vertex) * (connections.size() - 1);
        }

        return orphanLogLike;
    }

    private Tree getRepresentant() {
        return components.iterator().next();
    }

    public double getR() {
        // ASSUMPTION: These two params are consistent between components.
        // TODO: Get rid of assumption?
        return getRepresentant().hmm2.params[0];
    }

    public double getLambda() {
        // ASSUMPTION: These two params are consistent between components.
        // TODO: Get rid of assumption?
        return getRepresentant().hmm2.params[1];
    }

    public double getMu() {
        // ASSUMPTION: These two params are consistent between components.
        // TODO: Get rid of assumption?
        return getRepresentant().hmm2.params[2];
    }

    @Override
    public double getHeat() {
        return heat;
    }

    public State getState() {
        Tree firstComponent = components.iterator().next();
        final Vertex rootVertex = firstComponent.vertex.get(0);

        int nodes = countNodes() - components.size() + 1;
        State state = new State(nodes);

        ArrayList<LinkedList<Integer>> children = new ArrayList<LinkedList<Integer>>(nodes);
        for (int i = 0; i < nodes; i++)
            children.add(new LinkedList<Integer>());

        Map<Vertex, Integer> lookup = new HashMap<Vertex, Integer>();

        int labeledCounter = 0;
        int unlabeledCounter = n;

        Queue<Vertex> queue = new LinkedList<Vertex>();
        queue.add(rootVertex); // Add node to start BFS traversal

        while (!queue.isEmpty()) {
            Vertex current = queue.poll();

            boolean visitedBefore = lookup.containsKey(current); // This happens when searching connected components
            if (!visitedBefore) {
                boolean isLabeled = labeledVertexIds.containsKey(current);
                if (isLabeled)
                    lookup.put(current, labeledCounter++);
                else
                    lookup.put(current, unlabeledCounter++);
            }

            if (current.parent != null) {
                if (lookup.containsKey(current.parent)) {
                    state.parent[lookup.get(current)] = lookup.get(current.parent);
                    children.get(lookup.get(current.parent)).add(lookup.get(current));

                    if (state.align[lookup.get(current)] != null)
                        throw new RuntimeException("Alignment set twice!");
                    state.align[lookup.get(current)] = current.getAlign();
                    state.edgeLen[lookup.get(current)] = current.edgeLength;
                } else
                    queue.add(current.parent);
            }
            if (current.left != null) {
                if (lookup.containsKey(current.left)) {
                    state.parent[lookup.get(current)] = lookup.get(current.left);
                    children.get(lookup.get(current.left)).add(lookup.get(current));

                    if (state.align[lookup.get(current)] != null)
                        throw new RuntimeException("Alignment set twice!");
                    state.align[lookup.get(current)] = reverseAlign(current.left.getAlign());
                    state.edgeLen[lookup.get(current)] = current.left.edgeLength;
                } else
                    queue.add(current.left);
            }
            if (current.right != null) {
                if (lookup.containsKey(current.right)) {
                    state.parent[lookup.get(current)] = lookup.get(current.right);
                    children.get(lookup.get(current.right)).add(lookup.get(current));

                    if (state.align[lookup.get(current)] != null)
                        throw new RuntimeException("Alignment set twice!");
                    state.align[lookup.get(current)] = reverseAlign(current.right.getAlign());
                    state.edgeLen[lookup.get(current)] = current.right.edgeLength;
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
            state.felsen[lookup.get(current)] = current.getFelsen();
        }

        state.parent[lookup.get(rootVertex)] = -1;
        state.root = lookup.get(rootVertex);

        // Set alignment for root
        state.align[state.root] = new int[rootVertex.length];
        Arrays.fill(state.align[state.root], 0);

        for (int i = 0; i < nodes; i++) {
            int size = children.get(i).size();
            state.children[i] = new int[size];
            int j = 0;
            for (int child : children.get(i))
                state.children[i][j++] = child;
        }

        state.indelParams = firstComponent.hmm2.params.clone();
        state.substParams = firstComponent.substitutionModel.params.clone();
        state.logLike = getLogLike();

        return state;
    }

    @Override
    public SubstitutionModel getSubstitutionModel() {
        return substitutionModel;
    }

    /**
     * Given two nodes parent and child, and given the alignment of child relative
     * to parent as described in Vertex.getAlign(), this method returns the
     * alignment of parent relative to child.
     * @param toParent Alignment of child relative to parent.
     * @return Alignment of parent relative to child.
     */
    private int[] reverseAlign(int[] toParent) {
        // Find the length of the parent sequence
        int length;
        int lastEntry = toParent[toParent.length - 1];
        if (lastEntry < 0)
            length = -lastEntry;
        else
            length = lastEntry + 1;

        int[] toChild = new int[length];
        int parentPos = 0;
        for (int childPos = 0; childPos < length; childPos++) {
            if (parentPos >= toParent.length) {
                toChild[childPos] = -(parentPos + 1);
                parentPos++;
            } else

            if (toParent[parentPos] == childPos) {
                // Match
                toChild[childPos] = parentPos++;
            } else if (toParent[parentPos] < 0) {
                // Insertion

                // TODO: Consider if this is always okay!
                childPos--;
                parentPos++;
            } else {
                // Deletion
                toChild[childPos] = -(parentPos + 1);
            }
        }

        return toChild;
    }

    private int countNodes() {
        int nodes = 0;
        for (Tree component : components)
            nodes += component.vertex.size();
        return nodes;
    }

    public String printedTree() {
        return getState().getNewickString();
    }

    public static void main(String[] args) throws Exception {
        /*
        String[] seqs = new String[] { "AAGT", "CGATTC", "CCGAAG", "AGACA", "TTGACC", "GTAC" };
        String[] names = new String[] { "A", "B", "C", "D", "E", "F" };

        // String[] seqs = new String[] { "AAGT", "GTAC" };
        // String[] names = new String[] { "A", "B" };

        SubstitutionModel model = new Kimura3();
        SubstitutionScore ss = model.attachedScoringScheme;

        Spannoid spannoid = new Spannoid(3, BonphyStrategy.TOTAL_LENGTH, seqs, names, model, ss);
        System.out.println("Log-like of spannoid: " + spannoid.getLogLike());
        */

        SpannoidUpdater updater = new SpannoidUpdater();

        // String tree = "((B:0.5,C:0.5):0.5)A;";
        // String tree = "((A:0.5,B:0.2):1,(D:1,E:0.2):1)C;";
        //String tree = "((D:0.1,(F:0.1,(E:0.2,C:0.05):0.2)B:0.1):0.2)A;";

        // String tree = "((B:0.1,((C:0.1,D:0.1):0.1,(G:0.1,(E:0.1,F:0.1):0.1):0.1):0.1):0.1)A;";
        // String tree = "(((B:0.1,((C:0.1,D:0.1):0.1,((E:0.1,F:0.1):0.1)G:0.1):0.1):0.1)A;";

        String tree = "((B:1,(C:1,(D:1,(E:1,(F:1,G:1):1):1):1):1):1)A:1;";

        String[] seqs = new String[] { "AAGT", "CGATTC", "CCGAAG", "AG", "TTGACCAAGC", "G", "ACGGT" };
        Map<String, Integer> nameMap = new HashMap<String, Integer>();
        for (int i = 0; i < seqs.length; i++)
            nameMap.put("" + (char)((int)'A' + i), i);

        SubstitutionModel model = new Kimura3();
        SubstitutionScore ss = model.attachedScoringScheme;

        Spannoid spannoid = new Spannoid(tree, seqs, nameMap, model, ss);
        System.out.println(String.format("Log-like: %f", spannoid.getLogLike()));

        // Find vertex named A
        Vertex testVertex = null;
        for (Vertex v : spannoid.components.get(0).vertex) {
            if ("D".equals(v.name)) {
                testVertex = v;
                break;
            }
        }

        SpannoidViewer viewer = new SpannoidViewer();
        viewer.newSample(spannoid.getState(), 0, 0);

        updater.contractEdge(spannoid, testVertex.parent, testVertex);

        viewer.newSample(spannoid.getState(), 0, 0);

        /*
        // Print alignments
        for (Tree component : spannoid.components) {
            System.out.println();
            String[] alignment = component.root.printedMultipleAlignment();
            for (String seq : alignment)
                System.out.println(seq);
        }

        System.out.println();
        System.out.println("Combined alignment:");
        State state = spannoid.getState();
        String[] alignment = state.getLeafAlign();
        // String[] alignment = state.getFullAlign();
        int i = 0;
        for (String seq : alignment) {
            String name = (state.name[i].isEmpty()) ? "-" : state.name[i];
            System.out.println(String.format("%s:\t%s", name, seq));
            i++;
        } */
    }

    public static class SpannoidUpdater extends AbstractUpdater<Spannoid> {
        private static SteinerTreeUpdater updater = new SteinerTreeUpdater();

        @Override
        public void updateR(Spannoid tree, double newR) {
            for (Tree component : tree.components)
                updater.updateR(component, newR);
        }

        @Override
        public void updateLambda(Spannoid tree, double newLambda) {
            for (Tree component : tree.components)
                updater.updateLambda(component, newLambda);
        }

        @Override
        public void updateMu(Spannoid tree, double newMu) {
            for (Tree component : tree.components)
                updater.updateMu(component, newMu);
        }

        @Override
        public void recalcSubstitutionParameters(Spannoid tree) {
            for (Tree component : tree.components)
                updater.recalcSubstitutionParameters(component);
        }

        @Override
        public void revertNNI(Tree tree, AbstractUpdater.NNIResult nni){
            super.revertNNI(tree, nni);
        }

        public Tree getRandomComponent(Spannoid spannoid) {
            int i = Utils.generator.nextInt(spannoid.components.size());
            return spannoid.components.get(i);
        }

        public Vertex getRandomVertex(Spannoid spannoid) {
            Tree component = getRandomComponent(spannoid);
            int j = Utils.generator.nextInt(component.vertex.size());
            return component.vertex.get(j);
        }

        public Vertex getComponentRandomVertex(Tree component){
            int k = Utils.generator.nextInt(component.vertex.size());
            Vertex observed  = component.vertex.get(k);
            while (observed.left != null || observed.right != null){
                k = Utils.generator.nextInt(component.vertex.size());
                observed  = component.vertex.get(k);
            }
            return observed;

        }

        public Tree getRandomNeighbouringComponent(Spannoid spannoid, Vertex source){
            int sourceIndex = spannoid.labeledVertexIds.get(source);
            Set<Vertex> neighborhood = spannoid.componentConnections.get(sourceIndex);
            Vertex[] overlapping_vertices = new Vertex[neighborhood.size() - 1];
            int i = 0;
            for (Vertex v : neighborhood) {
                if (v != source)
                    overlapping_vertices[i++] = v;
            }

            // neighborhood.remove(source);
            int k = Utils.generator.nextInt(overlapping_vertices.length); //// TODO:if neighborhood is empty now ??????
            return overlapping_vertices[k].owner;

        }

        public Vertex getRandomBlack(Spannoid spannoid){
            int j = Utils.generator.nextInt(spannoid.n);
            Set<Vertex> neighborhood  = spannoid.componentConnections.get(j);
            Vertex[] overlapping_vertices = neighborhood.toArray(new Vertex[0]);
            int k = Utils.generator.nextInt(overlapping_vertices.length);
            return overlapping_vertices[k];
        }


        public Vertex getRandomInnerBlack(Spannoid spannoid){
            int j = Utils.generator.nextInt(spannoid.innerBlackNodes.size());
            Integer index = spannoid.innerBlackNodes.get(j);
            Set<Vertex> neighborhood = spannoid.componentConnections.get(index);
            Vertex[] overlapping_vertices = neighborhood.toArray(new Vertex[0]);
            int k = Utils.generator.nextInt(overlapping_vertices.length);
            return overlapping_vertices[k];
        }

        public Vertex getDestinationFromSourceMoveComponent(Spannoid spannoid, Vertex source) {
            List<Vertex> destSet = new ArrayList<Vertex>();
            Set<Vertex> connected = spannoid.componentConnections.get(spannoid.labeledVertexIds.get(source));
            for (Vertex con : connected) {
                if (con != source) {
                    for (Vertex v : con.owner.vertex) {
                        if (v != con && v.name != null)
                            destSet.add(v);
                    }
                }
            }
            int k = Utils.generator.nextInt(destSet.size());
            return destSet.get(k);
        }

        public Vertex getConnection(Spannoid spannoid, Vertex vertex) {
            int id = spannoid.labeledVertexIds.get(vertex);
            Set<Vertex> connections = spannoid.componentConnections.get(id);
            for (Vertex v : connections) {
                if (v != vertex)
                    return v;
            }
            return null;
        }

        public double moveComponent(Spannoid spannoid, Vertex source, Vertex dest){
            source.fullWin();
            source.parent.fullWin();
            double bpp = source.hmm2BackProp();

            // Update internal Spannoid structure
            int vId = spannoid.labeledVertexIds.get(source);
            spannoid.componentConnections.get(vId).remove(source);
            int destId = spannoid.labeledVertexIds.get(dest);
            spannoid.componentConnections.get(destId).add(source);
            spannoid.labeledVertexIds.put(source, destId);
            spannoid.setupInnerBlackNodes();

            // Update Vertex with new information and alignment
            source.seq = dest.seq;
            source.length = dest.length;
            source.name = dest.name;
            source.first = new AlignColumn(source);

            AlignColumn cur = dest.first;
            source.first.seq = cur.seq.clone();
            AlignColumn prev = source.first;
            cur  = cur.next;

            while (cur != dest.last && cur != null){
                AlignColumn actual = new AlignColumn(source);
                actual.seq = cur.seq.clone();

                actual.prev = prev;
                prev.next = actual;
                prev = actual;

                cur = cur.next;
            }

            AlignColumn last = new AlignColumn(source);
            last.prev = prev;
            prev.next = last;

            source.last = last;

            AlignColumn c = source.parent.first;
            AlignColumn n = source.first;
            while (c != source.parent.last) {
                if (source.parent.left == source) {
                    c.left = null;
                } else {
                    c.right = null;
                }

                c = c.next;
                if (n != null)
                    n = n.next;
            }

            source.last.parent = source.parent.last;
            source.last.orphan = false;

            if (source.parent.left == source)
                source.parent.last.left = source.last;
            else
                source.parent.last.right = source.last;

            source.fullWin();
            source.parent.fullWin();
            source.brother().fullWin();

            bpp += source.hmm2AlignWithSave();
            source.calcAllUp();

            return bpp;
        }

        /**
         * Swaps the alignment between child and parent.
         * @param child
         * @param parent Parent needs to be parent of child!
         */
        private void swapAlignment(Vertex child, Vertex parent, boolean addToLeft) {
            // Realign cur and parent
            AlignColumn acNewLeaf = parent.first;
            AlignColumn acNewParent = child.first;
            while (acNewLeaf != parent.last || acNewParent != child.last) {
                if (acNewParent.parent != acNewLeaf) {           // Deletion     (* -)   -> Insertion    (- *)
                    acNewLeaf.parent = acNewParent;
                    acNewLeaf.orphan = true;

                    acNewLeaf = acNewLeaf.next;
                } else if (acNewParent.orphan) {                 // Insertion    (- *)   -> Deletion     (* -)
                    acNewParent.orphan = false;
                    acNewParent.parent = null; // TODO: Should be set somewhere...

                    acNewParent = acNewParent.next;
                } else {                                         // Substitution (* *)   -> Substitution (* *)
                    // AlignColumn oldParent = acNewLeaf.parent;
                    acNewLeaf.parent = acNewParent;
                    acNewParent.parent = null;

                    if (addToLeft)
                        acNewParent.left = acNewLeaf;
                    else
                        acNewParent.right = acNewLeaf;

                    acNewLeaf = acNewLeaf.next;
                    acNewParent = acNewParent.next;
                }
            }

            parent.last.parent = child.last;
            child.last.parent = null;
        }

        private Tree createEmptyComponent(SubstitutionModel model) {
            Tree tree = new Tree();
            tree.substitutionModel = model;
            tree.hmm2 = new HmmTkf92(null);
            tree.hmm3 = new HmmNonParam();

            return tree;
        }

        /**
         * (actual -> guide -> reference) --> (actual -> reference)
         * @param actual The node to be aligned.
         * @param guide The guide node describing alignment between actual and reference.
         * @param reference The reference alignment.
         */
        private void alignAlignment(Vertex actual, Vertex guide, Vertex reference) {
            AlignColumn actualAC = actual.first;
            AlignColumn guideAC = guide.first;
            AlignColumn referenceAC = reference.first;

            // TODO: Update left, right pointers!!!
            while (actualAC != actual.last || guideAC != guide.last || referenceAC != reference.last) {
                if (actualAC.parent == guideAC && guideAC.parent == referenceAC) {
                    // Substitution (* * *) -> Substitution (* *)
                    actualAC.parent = referenceAC;

                    if (referenceAC.left == guideAC)
                        referenceAC.left = actualAC;
                    else
                        referenceAC.right = actualAC;

                    actualAC = actualAC.next;
                    guideAC = guideAC.next;
                    referenceAC = referenceAC.next;
                } else if (actualAC.orphan) {
                    // (* - -) -> (* -)
                    actualAC.parent = referenceAC;

                    actualAC = actualAC.next;
                } else if (actualAC.parent == guideAC && guideAC.orphan) {
                    // (* * -) -> (* -)        """"?
                    actualAC.parent = referenceAC;
                    actualAC.orphan = true;

                    actualAC = actualAC.next;
                    guideAC = guideAC.next;
                } else if (guideAC.parent == referenceAC && actualAC.parent != guideAC) {
                    // (- * *) -> (- *)
                    if (referenceAC.left == guideAC)
                        referenceAC.left = null;
                    else
                        referenceAC.right = null;

                    guideAC = guideAC.next;
                    referenceAC = referenceAC.next;
                } else if (guideAC.parent != referenceAC) {
                    // (- - *) -> (- *)
                    referenceAC = referenceAC.next;
                } else {
                    throw new RuntimeException("We fucked up!");
                }
            }

            actual.last.parent = reference.last;
        }

        private void copyAlignColumn(Vertex copy, Vertex original) {
            copy.name = original.name;
            copy.length = original.length;
            copy.seq = original.seq;

            AlignColumn cur = original.first;

            AlignColumn prev = null;
            while (cur != null) {
                AlignColumn ac = new AlignColumn(copy);
                ac.orphan = cur.orphan;
                ac.parent = cur.parent;
                ac.emptyWindow = cur.emptyWindow;
                ac.left = cur.left;
                ac.right = cur.right;
                if (cur.seq != null)
                    ac.seq = cur.seq.clone();

                ac.prev = prev;
                if (prev == null)
                    copy.first = ac;
                else
                    prev.next = ac;

                prev = ac;
                cur = cur.next;
            }
            copy.last = prev;
        }

        private void dfsCopyVertex(Vertex cur, Vertex prev, List<Vertex> vertices) {
            boolean labeled = cur.seq != null && !cur.seq.isEmpty();
            if (labeled)
                vertices.add(cur);

            if (cur.parent != null && cur.parent != prev)
                dfsCopyVertex(cur.parent, cur, vertices);
            if (cur.left != null && cur.left != prev)
                dfsCopyVertex(cur.left, cur, vertices);
            if (cur.right != null && cur.right != prev)
                dfsCopyVertex(cur.right, cur, vertices);

            if (!labeled)
                vertices.add(cur);
        }

        /**
         * Swaps the path starting from the labeled root and down guided by the direction array.
         * @param labelRoot
         * @param direction 0 -> left, 1 -> right
         */
        private void swapPath(Vertex labelRoot, boolean[] direction) {
            if (direction.length == 0) return;

            Vertex parent = labelRoot;
            Vertex child = direction[0] ? parent.right : parent.left;

            for (int i = 0; i < direction.length; i++) {
                Vertex nextChild = direction[i+1] ? child.right : child.left;
                swapAlignment(child, parent, !direction[i]);

                parent.parent = child;
                if (direction[i])
                    child.right = parent;
                else
                    child.left = parent;

                parent = child;
                child = nextChild;
            }
        }

        /**
         *
         * @param component
         * @param where Where to place the new fake root.
         */
        private void rerootComponent(Tree component, Vertex where) {
            List<Boolean> directions = new LinkedList<Boolean>();
            Vertex oldParent = where;
            while (oldParent.parent != null) {
                if (oldParent.parent.left == oldParent)
                    directions.add(false);
                else if (oldParent.parent.right == oldParent)
                    directions.add(true);
                else
                    throw new RuntimeException("Something horrible happened!");

                oldParent = oldParent.parent;
            }
            Collections.reverse(directions); // Get directions ordered from (current) parent to child.

            // Convert from list to array representation.
            // This is required since Java auto-unboxing is crap.
            int i = 0;
            boolean[] directionsArray = new boolean[directions.size()];
            for (boolean b : directions)
                directionsArray[i++] = b;

            Vertex leftSubtree = where.parent;

            // oldParent is labeled root node
            makeFakeAlignment(oldParent.left);
            oldParent.fullWin();
            oldParent.left.fullWin();
            oldParent.updateHmmMatrices(); // TODO: Consider if updating HMM matrices is necessary
            oldParent.left.updateHmmMatrices();
            oldParent.left.hmm2AlignWithSave();

            // Swap the alignment and pointers from the root to the 'where' vertex.
            swapPath(oldParent, directionsArray);

            Vertex newRoot = new Vertex(component, 0.0);
            newRoot.left = leftSubtree;
            newRoot.right = where;

            leftSubtree.parent = newRoot;
            where.parent = newRoot;

            component.root = newRoot;
            component.vertex.add(newRoot);

            drawNewAlignment(newRoot);
        }

        private void makeFakeAlignment(Vertex vertex) {
            AlignColumn fake = new AlignColumn(vertex);
            vertex.first = fake;
            vertex.last = fake;

            // Null check are required because of labeled root node at intermediate representation.
            if (vertex.left != null) {
                AlignColumn ac = vertex.left.first;
                while (ac != vertex.left.last) {
                    ac.orphan = true;
                    ac.parent = fake;
                    ac = ac.next;
                }
                ac.parent = fake;
                fake.left = ac;
            }

            if (vertex.right != null) {
                AlignColumn ac = vertex.right.first;
                while (ac != vertex.right.last) {
                    ac.orphan = true;
                    ac.parent = fake;
                    ac = ac.next;
                }
                ac.parent = fake;
                fake.right = ac;
            }
        }

        private void drawNewAlignment(Vertex vertex) {
            makeFakeAlignment(vertex);

            vertex.left.fullWin();
            vertex.right.fullWin();
            vertex.fullWin();

            // TODO: Consider if this is necessary
            vertex.updateHmmMatrices();
            vertex.updateHmmMatrices();
            vertex.updateHmmMatrices();

            vertex.hmm3AlignWithSave();
        }

        /**
         * @param spannoid The Spannoid
         * @param steiner A Steiner node in a component of the Spannoid.
         * @param labeled A labeled node adjacent to steiner which should be contracted onto steiner.
         * @return The back-proposal of the move.
         */
        public double contractEdge(Spannoid spannoid, Vertex steiner, Vertex labeled) {
            /*
             * Parent Component
             */
            Tree parentComponent = createEmptyComponent(spannoid.getSubstitutionModel());
            Vertex parentLeaf = new Vertex(parentComponent, steiner.edgeLength + labeled.edgeLength / 2);
            copyAlignColumn(parentLeaf, labeled);
            alignAlignment(parentLeaf, steiner, steiner.parent);

            // Copy vertices to new component
            parentComponent.vertex.add(parentLeaf);
            dfsCopyVertex(steiner.parent, steiner, parentComponent.vertex);

            parentLeaf.parent = steiner.parent;
            if (steiner.parent.left == steiner) {
                steiner.parent.left = parentLeaf;
            } else {
                steiner.parent.right = parentLeaf;
            }

            // Find the root
            Vertex root = parentLeaf;
            while (root.parent != null)
                root = root.parent;
            parentComponent.root = root;

            /*
             * Child Component
             */
            Vertex subTree = labeled.brother();
            Tree childComponent = createEmptyComponent(spannoid.getSubstitutionModel());
            Vertex childLeaf = new Vertex(childComponent, subTree.edgeLength + labeled.edgeLength / 2);
            childLeaf.left = subTree;
            subTree.parent = childLeaf;
            copyAlignColumn(childLeaf, labeled);
            dfsCopyVertex(subTree, labeled.parent, childComponent.vertex);

            rerootComponent(childComponent, subTree.right.left); // TODO: Root at random vertex

            childLeaf.left = null;

            // Update Spannoid information
            spannoid.components.remove(steiner.owner);
            spannoid.components.add(parentComponent);
            spannoid.components.add(childComponent);

            int id = spannoid.labeledVertexIds.get(labeled);
            Set<Vertex> connections = spannoid.componentConnections.get(id);
            connections.remove(labeled);
            connections.add(parentLeaf);
            connections.add(childLeaf);

            spannoid.labeledVertexIds.remove(labeled);
            spannoid.labeledVertexIds.put(parentLeaf, id);
            spannoid.labeledVertexIds.put(childLeaf, id);

            spannoid.setupInnerBlackNodes();

            return 0.0;
        }

        private void foobar(Vertex v) {
            AlignColumn ac = v.first;
            AlignColumn p = v.parent.first;

            while (ac != v.last || p != v.parent.last) {
                if (ac.parent == p) { // Substitution
                    ac.parent = p;
                    p.left = ac;
                } else if (ac.orphan) { // Insertion (- *) -> (* -)

                } else { // Deletion

                }
            }
        }



        /**
         * @param spannoid The Spannoid
         * @param steiner A Steiner node in a component of the Spannoid.
         * @param labeled A labeled node adjacent to steiner which should be contracted onto steiner.
         * @return The back-proposal of the move.
         */
        /*
        public double contractEdge(Spannoid spannoid, Vertex steiner, Vertex labeled) {
            Tree component = steiner.owner;

            // Remove the old component
            spannoid.components.remove(component);

            // Create <= 3 new components


            return 0.0; // TODO: Find correct bpp
        }       */

        private Set<Tree> splitComponent(Spannoid spannoid, Vertex steiner, Vertex labeled) {
            // TODO: Add some code...

            return null;
        }

        private Vertex traverse(Spannoid spannoid, Tree newComponent, Vertex prev, Vertex cur) {
            if (cur == null) return null;

            if (cur.left == null && cur.right == null) {
                // Leaf node
                return cur;
            }

            Vertex p = null, l = null, r = null;
            if (cur.parent != prev)
                p = traverse(spannoid, newComponent, cur, cur.parent);
            if (cur.left != prev)
                l = traverse(spannoid, newComponent, cur, cur.left);
            if (cur.right != prev)
                r = traverse(spannoid, newComponent, cur, cur.right);

            cur.owner = newComponent;

            if (cur.parent == null) { // Fake root vertex
                l.parent = prev;

                // Align to new parent


                return l;
            }

            // Reassign labeles s.t. only l and r are specified.
            if (p == null) {
                // The entire subtree is as it should be. Just update some pointers.
                return cur;
            } else {
                if (l == null)
                    cur.left = p;
                else if (r == null)
                    cur.right = p;
                else
                    throw new RuntimeException("Should not happen!");

                // Parent goes to left side
                cur.parent.parent = cur;
                cur.parent = prev;

                // Realign cur and parent
                AlignColumn acNewLeaf = p.first;
                AlignColumn acNewParent = cur.first;
                while (acNewLeaf != p.last || acNewParent != cur.last) {
                    if (acNewParent.parent != acNewLeaf) {           // Deletion     (* -)   -> Insertion    (- *)
                        AlignColumn ins = new AlignColumn(p);
                        acNewLeaf.next.prev = ins;
                        ins.next = acNewLeaf.next;
                        acNewLeaf.next = ins;
                        ins.prev = acNewLeaf;

                        acNewLeaf = ins.next;
                    } else if (acNewParent.orphan) {                 // Insertion    (- *)   -> Deletion     (* -)
                        AlignColumn tmp = acNewParent.prev;
                        tmp.next = acNewParent.next;
                        acNewParent.next.prev = tmp;

                        acNewParent = acNewParent.next;
                    } else {                                         // Substitution (* *)   -> Substitution (* *)
                        AlignColumn tmp = acNewLeaf.parent;
                        acNewLeaf.parent = acNewParent;
                        acNewParent.parent = tmp;

                        acNewLeaf = acNewLeaf.next;
                        acNewParent = acNewParent.next;
                    }
                }

                // Recalculate felsen
                p.calcFelsen();
                p.calcOrphan();
                cur.calcFelsen();
                cur.calcOrphan();

                return cur;
            }
        }

        /**
         * Return a dot representation of the AlignColumns of Vertex v relative to its parent.
         */
        public String printAlignColumns(Vertex v) {
            StringBuilder builder = new StringBuilder();
            builder.append("digraph ac {");
            builder.append("rankdir=LR; node [shape=Mrecord]; edge [tailclip=false];");

            int nameIndex = 0;
            Map<AlignColumn, String> nameMap = new HashMap<AlignColumn, String>();

            builder.append("{ rank=same; ");

            AlignColumn ac = v.first;
            while (ac != null) {
                String name;
                if (!nameMap.containsKey(ac)) {
                    name = "AC" + nameIndex++;
                    nameMap.put(ac, name);
                } else {
                    name = nameMap.get(ac);
                }

                char symbol = '-';
                if (ac != v.last) {
                    if (ac.seq == null)
                        symbol = '*';
                    else
                        symbol = ac.mostLikely();
                }

                String orphan = "";
                if (ac.orphan)
                    orphan = ",fillcolor=snow2,style=filled";

                builder.append("\"" + name + "\" [label=\"{ <prev> | <data> " + symbol + " | <next> }\"" + orphan + "];");

                ac = ac.next;
            }

            builder.append("}");
            builder.append("{ rank=same; ");

            AlignColumn parentAC = v.parent.first;
            while (parentAC != null) {
                String name;
                if (!nameMap.containsKey(parentAC)) {
                    name = "parentAC" + nameIndex++;
                    nameMap.put(parentAC, name);
                } else {
                    name = nameMap.get(parentAC);
                }

                char symbol = '-';
                if (parentAC != v.parent.last) {
                    if (parentAC.seq == null)
                        symbol = '*';
                    else
                        symbol = parentAC.mostLikely();
                }

                String orphan = "";
                if (parentAC.orphan)
                    orphan = ",fillcolor=snow2,style=filled";

                builder.append("\"" + name + "\" [label=\"{ <prev> | <data> " + symbol + " | <next> }\"" + orphan + "];");

                parentAC = parentAC.next;
            }
            builder.append("}");

            ac = v.first;
            while (ac != null) {
                if (ac.prev!= null) {
                    builder.append("\"" + nameMap.get(ac.prev) + "\":next:c -> \"" + nameMap.get(ac) + "\":data [arrowhead=vee, arrowtail=dot, dir=both];");
                }

                if (ac.next != null) {
                    builder.append("\"" + nameMap.get(ac.next) + "\":prev:c -> \"" + nameMap.get(ac) + "\":data [arrowhead=vee, arrowtail=dot, dir=both];");
                }

                if (ac.parent != null) {
                    builder.append("\"" + nameMap.get(ac) + "\":data -> \"" + nameMap.get(ac.parent) + "\":data [arrowhead=vee, arrowtail=none, dir=both];");
                }

                ac = ac.next;
            }

            parentAC = v.first.parent;
            while (parentAC != null) {
                if (parentAC.prev != null) {
                    builder.append("\"" + nameMap.get(parentAC.prev) + "\":next:c -> \"" + nameMap.get(parentAC) + "\":data [arrowhead=vee, arrowtail=dot, dir=both];");
                }

                if (parentAC.next != null) {
                    builder.append("\"" + nameMap.get(parentAC.next) + "\":prev:c -> \"" + nameMap.get(parentAC) + "\":data [arrowhead=vee, arrowtail=dot, dir=both];");
                }

                parentAC = parentAC.next;
            }

            builder.append("}");

            return builder.toString();
        }
    }
}
