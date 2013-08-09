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

import java.io.IOException;
import java.io.OutputStreamWriter;
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
     * Map from a given sequence ID to associated vertices.
     */
    private Map<Integer, Set<Vertex>> componentConnections = new HashMap<Integer, Set<Vertex>>();

    /*
     * Map from a given vertex to its associated ID.
     */
    private Map<Vertex, Integer> labeledVertexIds = new HashMap<Vertex, Integer>();

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
    }

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
        int i = 0;
        for (Vertex vertex : leafVertices) {
            // Do full alignment!
            vertex.fullWin();

            // Ensure transition matrices are updated
            vertex.edgeChangeUpdate();

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

            // TODO: Hack for making us able to check alignments!
            Vertex tmp = vertex.parent;
            vertex.parent = null;
            vertex.checkPointers();
            vertex.parent = tmp;

            // Do full alignment!
            vertex.fullWin();
            vertex.left.fullWin();
            vertex.right.fullWin();

            vertex.hmm3AlignWithRecalc();
        }

        // Align all steiner/internal nodes
        rootVertex.doRecAlign();
        rootVertex.calcFelsRecursively();
        rootVertex.calcIndelLikeRecursively();

        for (Vertex v : tree.vertex)
            v.checkPointers();

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
        fake.orphan = false;

        vertex.first = fake;
        vertex.last = fake;
        vertex.length = 0;
        if (vertex.right != null) {
            for (AlignColumn ac = vertex.right.first; ac != vertex.right.last; ac = ac.next) {
                ac.parent = fake;
                ac.orphan = true;
            }

            fake.right = vertex.right.last;
            vertex.right.last.parent = fake;
            vertex.right.last.orphan = false;
        }
        if (vertex.left != null) {
            for (AlignColumn ac = vertex.left.first; ac != vertex.left.last; ac = ac.next) {
                ac.parent = fake;
                ac.orphan = true;
            }

            fake.left = vertex.left.last;
            vertex.left.last.parent = fake;
            vertex.left.last.orphan = false;
        }
    }

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
        double logProb = 0;

        for (AlignColumn column = vertex.first; column != vertex.last; column = column.next) {
            logProb += Math.log(Utils.calcEmProb(column.seq, substitutionModel.e));
        }

        return logProb;
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
        return components.get(0);
    }

    public double getR() {
        // ASSUMPTION: These two params are consistent between components.
        // TODO: Get rid of assumption?
        return getRepresentant().getR();
    }

    public double getLambda() {
        // ASSUMPTION: These two params are consistent between components.
        // TODO: Get rid of assumption?
        final double lambda = getRepresentant().getLambda();

        // Consistenct check!
        for (Tree component : components) {
            if (lambda != component.getLambda())
                throw new RuntimeException();
        }

        return lambda;
    }

    public double getMu() {
        // ASSUMPTION: These two params are consistent between components.
        // TODO: Get rid of assumption?
        final double mu = getRepresentant().getMu();

        // Consistenct check!
        for (Tree component : components) {
            if (mu!= component.getMu())
                throw new RuntimeException();
        }

        return mu;
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
                    state.align[lookup.get(current)] = current.left.getReverseAlign();
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
                    state.align[lookup.get(current)] = current.right.getReverseAlign();
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

    private int countNodes() {
        int nodes = 0;
        for (Tree component : components)
            nodes += component.vertex.size();
        return nodes;
    }

    public String printedTree() {
        return getState().getNewickString();
    }

    private static Vertex getVertexByName(String name, int skip, Spannoid spannoid) {
        for (Tree component : spannoid.components) {
            for (Vertex v : component.vertex) {
                if (v.name != null && name.equals(v.name)) {
                    if (skip-- == 0)
                        return v;
                }
            }
        }
        return null;
    }

    public Set<Vertex> getLabelledVertices() {
        return labeledVertexIds.keySet();
    }

    public List<Vertex> getInnerLabelledVertices() {
        List<Vertex> result = new ArrayList<Vertex>();

        for (int i = 0; i < n; i++) {
            Set<Vertex> neighborhood = componentConnections.get(i);
            if (neighborhood.size() > 1) {
                result.addAll(neighborhood);
            }
        }

        return result;
    }

    public List<Vertex> getLabelledVerticesForContraction() {
        List<Vertex> nodes = new ArrayList<Vertex>();
        for (Vertex v : getLabelledVertices()) {
            if (v.owner.vertex.size() > 3) {
                nodes.add(v);
            }
        }
        return nodes;
    }

    public Set<Vertex> getNeighbourhood(Vertex vertex) {
        int vertexId = labeledVertexIds.get(vertex);
        return componentConnections.get(vertexId);
    }

    public static void main(String[] args) throws Exception {
        String tree = "((B:1,(C:1,(D:1,(E:1,(F:1,G:1):1):1):1):1):1)A:1;";

        String[] seqs = new String[] { "AAGT", "CGATTC", "CCGAAG", "AG", "TTGACCAAGC", "G", "ACGGT" };
        Map<String, Integer> nameMap = new HashMap<String, Integer>();
        for (int i = 0; i < seqs.length; i++)
            nameMap.put("" + (char)((int)'A' + i), i);

        SubstitutionModel model = new Kimura3();
        SubstitutionScore ss = model.attachedScoringScheme;

        Spannoid spannoid = new Spannoid(tree, seqs, nameMap, model, ss);
        System.out.println(String.format("Log-like: %f", spannoid.getLogLike()));

        SpannoidViewer viewer = new SpannoidViewer();
        SpannoidUpdater updater = new SpannoidUpdater();

        viewer.newSample(spannoid.getState(), 0, 0);
        SpannoidUpdater.ContractEdgeResult contraction = updater.contractEdge(spannoid, getVertexByName("C", 0, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);
        contraction = updater.contractEdge(spannoid, getVertexByName("B", 0, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);

        contraction = updater.contractEdge(spannoid, getVertexByName("E", 0, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);
        contraction = updater.contractEdge(spannoid, getVertexByName("G", 0, spannoid));
        // contraction = updater.contractEdge(spannoid, getVertexByName("D", 0, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);

        SpannoidUpdater.ExpandEdgeResult expansion = updater.expandEdge(spannoid, getVertexByName("C", 0, spannoid), getVertexByName("C", 1, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);
        updater.revertEdgeExpansion(expansion);
        updater.expandEdge(spannoid, getVertexByName("C", 0, spannoid), getVertexByName("C", 1, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);
        updater.expandEdge(spannoid, getVertexByName("E", 0, spannoid), getVertexByName("E", 1, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);
        updater.expandEdge(spannoid, getVertexByName("G", 0, spannoid), getVertexByName("G", 1, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);
        //updater.expandEdge(spannoid, getVertexByName("D", 0, spannoid), getVertexByName("D", 1, spannoid));
        //viewer.newSample(spannoid.getState(), 0, 0);
        updater.expandEdge(spannoid, getVertexByName("B", 0, spannoid), getVertexByName("B", 1, spannoid));
        viewer.newSample(spannoid.getState(), 0, 0);
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

        private void swapAlignment(Vertex child, Vertex parent, Direction direction) {
            swapAlignmentWithoutRecalc(child, parent, direction);
            parent.calcAllUp();
        }

        /**
         * Swaps the alignment between child and parent.
         * @param child
         * @param parent Parent needs to be parent of child!
         */
        private void swapAlignmentWithoutRecalc(Vertex child, Vertex parent, Direction direction) {
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
                    acNewParent.parent = null;

                    if (direction == Direction.LEFT)
                        acNewParent.left = null;
                    else
                        acNewParent.right = null;

                    acNewParent = acNewParent.next;
                } else {                                         // Substitution (* *)   -> Substitution (* *)
                    // AlignColumn oldParent = acNewLeaf.parent;
                    acNewLeaf.parent = acNewParent;
                    acNewParent.parent = null;

                    acNewParent.orphan = acNewLeaf.orphan;
                    acNewLeaf.orphan = false;

                    if (direction == Direction.LEFT)
                        acNewParent.left = acNewLeaf;
                    else
                        acNewParent.right = acNewLeaf;

                    acNewLeaf = acNewLeaf.next;
                    acNewParent = acNewParent.next;
                }
            }

            parent.last.orphan = child.last.orphan;
            child.last.orphan = false;

            // Reset pointers from previous parent to child.
            for (AlignColumn ac = parent.first; ac != null; ac = ac.next) {
                if (parent.left == child)
                    ac.left = null;
                else if (parent.right == child)
                    ac.right = null;
            }

            if (parent.left == child)
                parent.left = null;
            else if (parent.right == child)
                parent.right = null;

            // Change tree topology
            parent.last.parent = child.last;
            child.last.parent = null;
            parent.parent = child;
            child.parent = null;
            if (direction == Direction.LEFT) {
                child.last.left = parent.last;
                child.left = parent;
            } else {
                child.last.right = parent.last;
                child.right = parent;
            }

            parent.edgeLength = child.edgeLength;

            parent.checkPointers();

            parent.edgeChangeUpdate();
            child.edgeChangeUpdate();

            child.checkPointers();
        }

        private Tree createEmptyComponent(SubstitutionModel model,
                                          double R, double lambda, double mu) {
            Tree tree = new Tree();
            tree.substitutionModel = model;
            tree.hmm2 = new HmmTkf92(null);
            tree.hmm3 = new HmmNonParam();

            tree.hmm2.params = new double[] { R, lambda, mu };

            return tree;
        }

        private void alignAlignmentWithTopologyUpdate(Vertex actual, Vertex guide, Vertex reference) {
            alignAlignment(actual, guide, reference);

            actual.parent = reference;
            if (reference.left == guide)
                reference.left = actual;
            else if (reference.right == guide)
                reference.right = actual;
            else
                throw new RuntimeException();

            actual.edgeChangeUpdate();
            reference.edgeChangeUpdate();

            // TODO: Update edge lengths?
        }

        /**
         * (actual <-- guide <-- reference) --> (actual <-- reference)
         * @param actual The node to be aligned.
         * @param guide The guide node describing alignment between actual and reference.
         * @param reference The reference alignment.
         */
        private void alignAlignment(Vertex actual, Vertex guide, Vertex reference) {
            AlignColumn actualAC = actual.first;
            AlignColumn guideAC = guide.first;
            AlignColumn referenceAC = reference.first;

            while (actualAC != actual.last || guideAC != guide.last || referenceAC != reference.last) {
                if (actualAC.parent == guideAC && guideAC.parent == referenceAC && !actualAC.orphan && !guideAC.orphan) {
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
                    // (* * -) -> (* -)
                    actualAC.parent = referenceAC;
                    actualAC.orphan = true;

                    actualAC = actualAC.next;
                    guideAC = guideAC.next;
                } else if (guideAC.parent == referenceAC && actualAC.parent != guideAC && !guideAC.orphan) {
                    // (- * *) -> (- *)
                    if (referenceAC.left == guideAC)
                        referenceAC.left = null;
                    else
                        referenceAC.right = null;

                    guideAC = guideAC.next;
                    referenceAC = referenceAC.next;
                } else if (actualAC.parent != guideAC && guideAC.orphan) {
                    // (- * -) -> ignore
                    guideAC = guideAC.next;
                } else if (guideAC.parent != referenceAC) {
                    // (- - *) -> (- *)
                    referenceAC = referenceAC.next;
                } else {
                    throw new RuntimeException("We fucked up!");
                }
            }

            actual.last.parent = reference.last;
            if (reference.last.left == guide.last)
                reference.last.left = actual.last;
            else
                reference.last.right = actual.last;
        }

        private void copyAlignmentColumn(Vertex copy, Vertex original) {
            copy.length = original.length;

            AlignColumn cur = original.first;

            AlignColumn prev = null;
            while (cur != null) {
                AlignColumn ac = new AlignColumn(copy);
                ac.orphan = cur.orphan;
                ac.parent = cur.parent;
                ac.emptyWindow = cur.emptyWindow;
                ac.selected = cur.selected;
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

        private void copyVertex(Vertex copy, Vertex original) {
            copy.name = original.name;
            copy.seq = original.seq;

            copyAlignmentColumn(copy, original);
        }

        private void dfsMoveVertex(Vertex cur, Vertex prev, List<Vertex> labeled, List<Vertex> unlabeled, Tree newComponent) {
            boolean labeledNode = cur.seq != null && !cur.seq.isEmpty();
            cur.owner = newComponent;

            if (labeledNode)
                labeled.add(cur);
            else
                unlabeled.add(cur);

            if (cur.parent != null && cur.parent != prev)
                dfsMoveVertex(cur.parent, cur, labeled, unlabeled, newComponent);
            if (cur.left != null && cur.left != prev)
                dfsMoveVertex(cur.left, cur, labeled, unlabeled, newComponent);
            if (cur.right != null && cur.right != prev)
                dfsMoveVertex(cur.right, cur, labeled, unlabeled, newComponent);
        }

        private void dfsMoveVertex(Vertex cur, Vertex prev, Tree newComponent) {
            List<Vertex> labeled = new LinkedList<Vertex>();
            List<Vertex> unlabeled = new LinkedList<Vertex>();

            dfsMoveVertex(cur, prev, labeled, unlabeled, newComponent);
            newComponent.vertex.addAll(labeled);
            newComponent.vertex.addAll(unlabeled);
        }

        private void changeOwnerVertices(Tree newOwner, Vertex root) {
            if (root.left != null && root.right != null) {
                changeOwnerVertices(newOwner, root.left);
                changeOwnerVertices(newOwner, root.right);
            }

            root.owner = newOwner;
        }

        public enum Direction {
            LEFT, RIGHT;

            public Direction flip() {
                if (this == LEFT)
                    return RIGHT;
                else
                    return LEFT;
            }

            public Vertex getChild(Vertex parent) {
                if (this == LEFT)
                    return parent.left;
                else
                    return parent.right;
            }
        }

        /**
         * Swaps the path starting from the labeled root and down guided by the direction array.
         */
        private void swapPath(Vertex labelRoot, Direction[] direction) {
            if (direction.length == 0) return;

            Vertex parent = labelRoot;
            Vertex child = direction[0].getChild(parent);

            for (int i = 0; i < direction.length - 1; i++) {
                Vertex nextChild = direction[i+1].getChild(child);

                swapAlignmentWithoutRecalc(child, parent, direction[i+1]);

                parent = child;
                child = nextChild;
            }

            child.edgeChangeUpdate();
            // child.checkPointers();

            labelRoot.calcAllUp();
        }

        /**
         *
         * @param where Where to place the new fake root.
         */
        private Vertex rerootComponent(Vertex where, boolean doAlign) {
            return rerootComponent(where, null, doAlign);
        }

        private Vertex rerootComponent(Vertex where, MuDouble p, boolean doAlign) {
            final double edgeSum = where.edgeLength;

            List<Direction> directions = new LinkedList<Direction>();
            Vertex oldParent = where;
            while (oldParent.parent != null) {
                if (oldParent.parent.left == oldParent)
                    directions.add(Direction.LEFT);
                else if (oldParent.parent.right == oldParent)
                    directions.add(Direction.RIGHT);
                else
                    throw new RuntimeException("Something horrible happened!");

                oldParent = oldParent.parent;
            }
            Collections.reverse(directions); // Get directions ordered from (current) parent to child.
            Direction[] directionsArray = directions.toArray(new Direction[0]);

            Vertex newRoot = new Vertex(where.owner, 1);
            if (where == oldParent) {
                newRoot.left = oldParent.left;
                newRoot.right = oldParent;
                newRoot.right.edgeLength = 0;

                oldParent.left.parent = newRoot;
                oldParent.parent = newRoot;

                oldParent.left = null;
                oldParent.right = null;
            } else {
                final Vertex brotherSubtree = where.parent;

                // Swap the alignment and pointers from the root to the 'where' vertex.
                swapPath(oldParent, directionsArray);

                newRoot.left = brotherSubtree;
                newRoot.right = where;

                brotherSubtree.parent = newRoot;
                where.parent = newRoot;
            }

            oldParent.left = null;
            for (AlignColumn ac = oldParent.first; ac != null; ac = ac.next)
                ac.left = null;

            double bpp = 0;

            double[] edgeSplit = splitEdgeLength(edgeSum);
            newRoot.left.edgeLength = edgeSplit[0];
            newRoot.right.edgeLength = edgeSplit[1];

            // choose how to split the edge
            bpp -= -Math.log(edgeSum - 0.01);

            // choose new alignment at root
            if (doAlign) {
                bpp -= -drawNewAlignment(newRoot);
            } else {
                makeFakeAlignment(newRoot);
            }

            newRoot.left.calcOrphan();
            newRoot.right.calcOrphan();
            newRoot.left.calcAllUp();

            newRoot.left.checkPointers();
            newRoot.right.checkPointers();
            newRoot.checkPointers();

            if (p != null) {
                p.value = bpp;
            }

            oldParent.edgeChangeUpdate();

            return newRoot;
        }

        private double[] splitEdgeLength(double edgeSum) {
            final double e1 = Utils.generator.nextDouble() * (edgeSum - 0.01);
            final double e2 = edgeSum - 0.01 - e1;
            return new double[] { 0.01 + e1, 0.01 + e2 };
        }

        private void makeFakeAlignment(Vertex vertex) {
            AlignColumn fake = new AlignColumn(vertex);
            fake.orphan = false;

            vertex.first = fake;
            vertex.last = fake;
            vertex.length = 0;

            // Null check are required because of labeled root node at intermediate representation.
            if (vertex.left != null) {
                for (AlignColumn ac = vertex.left.first; ac != vertex.left.last; ac = ac.next) {
                    ac.orphan = true;
                    ac.parent = fake;
                }
                vertex.left.last.parent = fake;
                vertex.left.last.orphan = false;
                fake.left = vertex.left.last;
            }

            if (vertex.right != null) {
                for (AlignColumn ac = vertex.right.first; ac != vertex.right.last; ac = ac.next) {
                    ac.orphan = true;
                    ac.parent = fake;
                }
                vertex.right.last.parent = fake;
                vertex.right.last.orphan = false;
                fake.right = vertex.right.last;
            }
        }

        private double drawNewAlignment(Vertex vertex) {
            double bpp;

            makeFakeAlignment(vertex);
            vertex.left.fullWin();
            vertex.right.fullWin();
            vertex.fullWin();

            vertex.checkPointers();

            vertex.left.edgeChangeUpdate();
            vertex.right.edgeChangeUpdate();
            vertex.edgeChangeUpdate();

            bpp = vertex.hmm3AlignWithRecalc();
            vertex.checkPointers();

            return bpp;
        }

        private void backupTree(Vertex v) {
            if (v.parent == null) {
                // We are the parent. We store alignment.
                v.topologyBackup = new Vertex();

                v.topologyBackup.name = v.name;
                v.topologyBackup.length = v.length;
                v.topologyBackup.seq = v.seq;

                v.topologyBackup.first = v.first;
                v.topologyBackup.last = v.last;
                copyVertex(v, v.topologyBackup);
            }

            Map<AlignColumn, AlignColumn> oldToNewMap = new HashMap<AlignColumn, AlignColumn>();
            AlignColumn oldAC = v.topologyBackup.first;
            for (AlignColumn newAC = v.first; newAC != null; newAC = newAC.next) {
                oldToNewMap.put(oldAC, newAC);
                oldAC = oldAC.next;
            }

            if (v.left != null && v.right != null) {
                v.left.topologyBackup = new Vertex();
                v.left.topologyBackup.name = v.left.name;
                v.left.topologyBackup.length = v.left.length;
                v.left.topologyBackup.seq = v.left.seq;

                v.left.topologyBackup.first = v.left.first;
                v.left.topologyBackup.last = v.left.last;
                copyVertex(v.left, v.left.topologyBackup);

                v.right.topologyBackup = new Vertex();
                v.right.topologyBackup.name = v.right.name;
                v.right.topologyBackup.length = v.right.length;
                v.right.topologyBackup.seq = v.right.seq;

                v.right.topologyBackup.first = v.right.first;
                v.right.topologyBackup.last = v.right.last;
                copyVertex(v.right, v.right.topologyBackup);

                backupTree(v.left);
                backupTree(v.right);

                for (AlignColumn ac = v.left.first; ac != null; ac = ac.next) {
                    if (ac.parent != null) {
                        if (ac.parent.left != null)
                            oldToNewMap.get(ac.parent).left = ac;
                        ac.parent = oldToNewMap.get(ac.parent);
                    }
                }

                for (AlignColumn ac = v.right.first; ac != null; ac = ac.next) {
                    if (ac.parent != null) {
                        if (ac.parent.right != null)
                            oldToNewMap.get(ac.parent).right = ac;
                        ac.parent = oldToNewMap.get(ac.parent);
                    }
                }
            }

            v.topologyBackup.name = v.name;
            v.topologyBackup.length = v.length;
            v.topologyBackup.edgeLength = v.edgeLength;
            v.topologyBackup.leafCount = v.leafCount;

            v.topologyBackup.parent = v.parent;
            v.topologyBackup.right = v.right;
            v.topologyBackup.left = v.left;

            v.topologyBackup.winFirst = v.winFirst;
            v.topologyBackup.winLast = v.winLast;
            v.topologyBackup.winLength = v.winLength;

            v.winFirst = oldToNewMap.get(v.topologyBackup.winFirst);
            v.winLast = oldToNewMap.get(v.topologyBackup.winLast);

            v.topologyBackup.hmm2TransMatrix = v.hmm2TransMatrix;
            v.topologyBackup.hmm2PropTransMatrix = v.hmm2PropTransMatrix;
            v.topologyBackup.hmm3RedTransMatrix = v.hmm3RedTransMatrix;
            v.topologyBackup.hmm3TransMatrix = v.hmm3TransMatrix;
        }

        private void restoreVertex(Vertex v) {
            if (v.topologyBackup == null) {
                throw new RuntimeException("Trying to restore topology from non-backuped node.");
            }

            v.name = v.topologyBackup.name;
            v.length = v.topologyBackup.length;
            v.edgeLength = v.topologyBackup.edgeLength;
            v.leafCount = v.topologyBackup.leafCount;

            v.parent = v.topologyBackup.parent;
            v.right = v.topologyBackup.right;
            v.left = v.topologyBackup.left;

            v.first = v.topologyBackup.first;
            v.last = v.topologyBackup.last;
            v.winFirst = v.topologyBackup.winFirst;
            v.winLast = v.topologyBackup.winLast;
            v.winLength = v.topologyBackup.winLength;

            v.hmm2TransMatrix = v.topologyBackup.hmm2TransMatrix;
            v.hmm2PropTransMatrix = v.topologyBackup.hmm2PropTransMatrix;
            v.hmm3RedTransMatrix = v.topologyBackup.hmm3RedTransMatrix;
            v.hmm3TransMatrix = v.topologyBackup.hmm3TransMatrix;
        }

        public static class ContractEdgeResult extends MCMCResult {
            public Spannoid spannoid;
            public Vertex up;
            public Vertex down;
            public Vertex contractedVertex;
            public Vertex steiner;
            public boolean specialCase;

            public ContractEdgeResult(Spannoid spannoid, Vertex contractedVertex, Vertex steiner,
                                      Vertex parentTree, Vertex childTree, boolean specialCase) {
                this.spannoid = spannoid;
                this.up = parentTree;
                this.down = childTree;
                this.contractedVertex = contractedVertex;
                this.steiner = steiner;
                this.specialCase = specialCase;
            }
        }

        public int countNodes(Vertex v) {
            int res = 1;
            if (v.left != null)
                res += countNodes(v.left);
            if (v.right != null)
                res += countNodes(v.right);
            return res;
        }

        private void updateAlignColumnReferences(Vertex cur, Vertex old) {
            AlignColumn oldAC = old.first;
            for (AlignColumn ac = cur.first; ac != null; ac = ac.next) {
                if (ac.parent.left == oldAC)
                    ac.parent.left = ac;
                if (ac.parent.right == oldAC)
                    ac.parent.right = ac;

                oldAC = oldAC.next;
            }
        }

        public ContractEdgeResult contractEdge(Spannoid spannoid, Vertex labeled) {
            final Vertex originalSteiner = labeled.parent;
            Vertex steiner = originalSteiner;
            final Tree originalComponent = labeled.owner;

            double bpp = 0;

            // Backup the entire tree
            backupTree(originalComponent.root);

            // SPECIAL CASE: Steiner is the fake root vertex!
            final boolean specialCase = (steiner.parent == null);

            if (specialCase) {
                Vertex brother = labeled.brother();

                brother.fullWin();
                brother.left.fullWin();
                brother.right.fullWin();
                bpp += brother.hmm3BackProp();
                steiner.fullWin();
                steiner.left.fullWin();
                steiner.right.fullWin();
                bpp += steiner.hmm3BackProp();

                // Swap alignment of labeled and fake root.
                swapAlignment(labeled, steiner, Direction.LEFT);
                // Make alignment of labeled --> brother
                alignAlignment(brother, steiner, labeled);
                labeled.edgeLength += brother.edgeLength - 0.01;
                bpp += -Math.log(labeled.edgeLength - 0.01); // checked

                // Forget steiner node
                labeled.parent = null;
                brother.parent = labeled;
                labeled.left = brother;

                labeled.checkPointers();
                brother.checkPointers();

                // choose root for up component
                List<Vertex> choices = getSubtreeNodes(brother.left);
                Vertex newRootPosition = choices.get(Utils.generator.nextInt(choices.size()));
                bpp -= -Math.log(choices.size());

                MuDouble p = new MuDouble();
                Vertex newRoot = rerootComponent(newRootPosition, p, true); // brother.left is okay, since otherwise component would be non-contractible!
                bpp -= p.value; // checked
                labeled.owner.root = newRoot;
                labeled.owner.vertex.add(newRoot);

                brother.edgeChangeUpdate();

                steiner = labeled.parent;
            } else {
                steiner.parent.fullWin();
                steiner.fullWin();
                steiner.left.fullWin();
                steiner.right.fullWin();
                bpp += steiner.hmm2BackProp();
                bpp += steiner.hmm3BackProp();
            }

            final double R = originalComponent.getR(),
                         lambda = originalComponent.getLambda(),
                         mu = originalComponent.getMu();
            Tree upComponent = createEmptyComponent(originalComponent.getSubstitutionModel(), R, lambda, mu);
            Tree downComponent = createEmptyComponent(originalComponent.getSubstitutionModel(), R, lambda, mu);

            // Create copy of nodes, with edge length corresponding to original case.
            // (We ignore edge length of labeled node)
            Vertex up = new Vertex(upComponent, steiner.edgeLength);
            Vertex down = new Vertex(downComponent, labeled.brother().edgeLength);

            // Handle up component
            upComponent.root = originalComponent.root;

            copyVertex(up, labeled);

            // Update steiner ACs to point to up instead of labeled
            updateAlignColumnReferences(up, labeled);
            alignAlignmentWithTopologyUpdate(up, steiner, steiner.parent);

            up.checkPointers();
            steiner.parent.checkPointers();

            dfsMoveVertex(up, null, upComponent);

            up.calcIndelLogLikeUp();

            // Handle down component
            Vertex brother = labeled.brother();

            copyVertex(down, labeled);

            // Make (non-trivial and fast) alignment of down and brother
            updateAlignColumnReferences(down, up);

            down.parent = steiner;
            steiner.parent = null;
            for (AlignColumn ac = steiner.first; ac != null; ac = ac.next)
                ac.parent = null;

            if (steiner.left == labeled)
                steiner.left = down;
            else if (steiner.right == labeled)
                steiner.right = down;
            else
                throw new RuntimeException();

            down.checkPointers();
            steiner.checkPointers();

            swapAlignmentWithoutRecalc(down, steiner, Direction.LEFT);
            alignAlignmentWithTopologyUpdate(brother, steiner, down);

            down.checkPointers();
            brother.checkPointers();

            dfsMoveVertex(down, null, downComponent);

            // choose root for down component
            List<Vertex> choices = getSubtreeNodes(brother);
            choices.add(down);
            final Vertex rootAt = choices.get(Utils.generator.nextInt(choices.size()));
            bpp -= -Math.log(choices.size());

            MuDouble p = new MuDouble();
            Vertex newRoot = rerootComponent(rootAt, p, true);
            bpp -= p.value; // checked
            downComponent.vertex.add(newRoot);
            downComponent.root = newRoot;

            down.checkPointers();
            down.parent.checkPointers();
            brother.checkPointers();
            newRoot.checkPointers();

            // Update Spannoid information
            final int before = spannoid.components.size();
            spannoid.components.remove(steiner.owner);
            spannoid.components.add(upComponent);
            spannoid.components.add(downComponent);
            if (spannoid.components.size() != before + 1)
                throw new RuntimeException("Components not added/removed correctly!");

            int id = spannoid.labeledVertexIds.get(labeled);
            Set<Vertex> connections = spannoid.componentConnections.get(id);
            connections.remove(labeled);
            connections.add(up);
            connections.add(down);

            spannoid.labeledVertexIds.remove(labeled);
            spannoid.labeledVertexIds.put(up, id);
            spannoid.labeledVertexIds.put(down, id);

            upComponent.root.calcFelsRecursively();
            upComponent.root.calcIndelLikeRecursively();
            downComponent.root.calcFelsRecursively();
            downComponent.root.calcIndelLikeRecursively();

            checkSpannoid(spannoid);

            ContractEdgeResult result = new ContractEdgeResult(spannoid, labeled, originalSteiner, up, down, specialCase);
            result.bpp = bpp;

            return result;
        }

        public void revertEdgeContraction(ContractEdgeResult contraction)
        {
            Tree originalComponent = contraction.contractedVertex.owner;

            if (contraction.steiner.topologyBackup.parent == null) {
                // Handle special case
                originalComponent.vertex.remove(originalComponent.root);
                originalComponent.root = contraction.steiner;
            }

            for (Vertex v : originalComponent.vertex)
                restoreVertex(v);

            // Revert the owner of the vertices
            changeOwnerVertices(originalComponent, originalComponent.root);

            // TODO: Remove this check!
            for (Vertex v : originalComponent.vertex)
                v.checkPointers();

            if (originalComponent.vertex.size() != countNodes(originalComponent.root))
                throw new RuntimeException("Wrong number of nodes!");

            originalComponent.root.calcFelsRecursively();
            originalComponent.root.calcIndelLikeRecursively();

            // Update spannoid information
            Spannoid spannoid = contraction.spannoid;
            final int before = spannoid.components.size();
            spannoid.components.remove(contraction.up.owner);
            spannoid.components.remove(contraction.down.owner);
            spannoid.components.add(originalComponent);
            if (spannoid.components.size() != before - 1)
                throw new RuntimeException("Components not added/removed correctly!");

            int id = spannoid.labeledVertexIds.get(contraction.up);
            Set<Vertex> connections = spannoid.componentConnections.get(id);
            connections.remove(contraction.up);
            connections.remove(contraction.down);
            connections.add(contraction.contractedVertex);

            spannoid.labeledVertexIds.remove(contraction.up);
            spannoid.labeledVertexIds.remove(contraction.down);
            spannoid.labeledVertexIds.put(contraction.contractedVertex, id);
        }

        private List<Vertex> getSubtreeNodes(Vertex root) {
            List<Vertex> nodes = new ArrayList<Vertex>();
            collectNodes(root, nodes);
            return nodes;
        }

        private void collectNodes(Vertex v, List<Vertex> nodes) {
            nodes.add(v);
            if (v.left != null && v.right != null) {
                collectNodes(v.left, nodes);
                collectNodes(v.right, nodes);
            }
        }

        private void dummyAlignToParent(Vertex v) {
            AlignColumn fake = new AlignColumn(v);
            fake.orphan = false;
            fake.parent = v.parent.last;
            v.first = fake;
            v.last = fake;

            for (AlignColumn ac = v.parent.first; ac != v.parent.last; ac = ac.next) {
                if (v.parent.left == v) {
                    ac.left = null;
                } else if (v.parent.right == v) {
                    ac.right = null;
                } else {
                    throw new RuntimeException();
                }
            }
            if (v.parent.left == v)
                v.parent.last.left = fake;
            else if (v.parent.right == v)
                v.parent.last.right = fake;
            else
                throw new RuntimeException();

            v.fullWin();
            v.parent.fullWin();
            v.edgeChangeUpdate();
            v.parent.edgeChangeUpdate();

            v.checkPointers();
        }

        public class ExpandEdgeResult extends MCMCResult {
            Spannoid spannoid;
            Vertex up;
            Vertex down;
            Vertex labeled;

            public ExpandEdgeResult(Spannoid spannoid, Vertex up, Vertex down, Vertex labeled) {
                this.spannoid = spannoid;
                this.up = up;
                this.down = down;
                this.labeled = labeled;
            }
        }

        public ExpandEdgeResult expandEdge(Spannoid spannoid, Vertex up, Vertex down)
        {
            double bpp = 0;

            int sizeOfUp = up.owner.vertex.size() - 2;
            int sizeOfDown = down.owner.vertex.size() - 2;

            Vertex upRoot = up.owner.root;
            Vertex downRoot = down.owner.root;

            int choice = Utils.weightedChoose(new int[] {sizeOfUp, sizeOfDown, 1});

            boolean removeUpRoot = false, removeDownRoot = false;

            boolean rootAtExpandedEdge = false;
            switch (choice) {
                case 0:
                    removeUpRoot = false;
                    removeDownRoot = true;

                    bpp -= Math.log(sizeOfUp);
                    break;
                case 1:
                    removeUpRoot = true;
                    removeDownRoot = false;

                    Vertex temp = up;
                    up = down;
                    down = temp;

                    bpp -= Math.log(sizeOfDown);
                    break;
                case 2:
                    removeUpRoot = true;
                    removeDownRoot = true;

                    rootAtExpandedEdge = true;

                    // bpp -= Math.log(1);
                    break;
            }
            bpp -= -Math.log(sizeOfUp + sizeOfDown + 1);

            if (removeUpRoot) {
                // backproposal for placing the up root
                bpp += -Math.log(sizeOfUp);
                bpp += -Math.log((upRoot.left.edgeLength - 0.01)
                        +(upRoot.right.edgeLength - 0.01));

                upRoot.fullWin();
                upRoot.left.fullWin();
                upRoot.right.fullWin();
                bpp += upRoot.hmm3BackProp();
            }
            if (removeDownRoot) {
                // backproposal for placing the down root
                bpp += -Math.log(sizeOfDown);
                bpp += -Math.log((downRoot.left.edgeLength - 0.01)
                        +(downRoot.right.edgeLength - 0.01));

                downRoot.fullWin();
                downRoot.left.fullWin();
                downRoot.right.fullWin();
                bpp += downRoot.hmm3BackProp();
            }

            backupTree(up.owner.root);
            backupTree(down.owner.root);

            final double R = spannoid.getR(), lambda = spannoid.getLambda(), mu = spannoid.getMu();
            Tree newComponent = createEmptyComponent(spannoid.getSubstitutionModel(), R, lambda, mu);

            // choose edge length
            double labeledEdgeLength = 0.01 - Math.log(Utils.generator.nextDouble());
            bpp -= -(labeledEdgeLength - 0.01);

            Vertex labeled = new Vertex(newComponent, labeledEdgeLength);
            copyVertex(labeled, up);

            Vertex steiner = new Vertex(newComponent, up.edgeLength);

            /*
             * Create connection to labeled node (left)
             */
            if (!rootAtExpandedEdge) {
                // Keep the root in up-component.
                steiner.parent = up.parent;
                if (up.parent.left == up) {
                    up.parent.left = steiner;
                }
                else if (up.parent.right == up) {
                    up.parent.right = steiner;
                }
                else {
                    throw new RuntimeException();
                }
                dummyAlignToParent(steiner);

                newComponent.root = up.owner.root;

                steiner.left = labeled;
                labeled.parent = steiner;
                for (AlignColumn ac = labeled.first; ac != labeled.last; ac = ac.next) {
                    ac.parent = steiner.last;
                    ac.orphan = true;

                    if (ac.left != null || ac.right != null) {
                        throw new RuntimeException();
                    }
                }
                labeled.last.parent = steiner.last;
                labeled.last.orphan = false;
                if (labeled.last.left != null || labeled.last.right != null) {
                    throw new RuntimeException();
                }
                steiner.last.left = labeled.last;
                labeled.fullWin();
                labeled.parent.fullWin();
                labeled.edgeChangeUpdate();
                labeled.parent.edgeChangeUpdate();

                labeled.checkPointers();
                steiner.checkPointers();
            } else {
                // Put root between steiner and labeled.
                Vertex newRoot = new Vertex(newComponent, 0);
                newComponent.root = newRoot;

                newRoot.left = labeled;
                labeled.parent = newRoot;
                newRoot.right = steiner;
                steiner.parent = newRoot;

                bpp -= -Math.log(labeledEdgeLength - 0.01); // checked
                double[] newEdgeLengths = splitEdgeLength(labeledEdgeLength);
                labeled.edgeLength = newEdgeLengths[0];
                steiner.edgeLength = newEdgeLengths[1];
                labeled.edgeChangeUpdate();
                steiner.edgeChangeUpdate();

                makeFakeAlignment(steiner);
                getRidOfRooting(up, steiner, Direction.LEFT);
            }

            getRidOfRooting(down, steiner, Direction.RIGHT);

            steiner.parent.fullWin();
            steiner.fullWin();
            steiner.left.fullWin();
            steiner.right.fullWin();

            // proposals for new alignment
            bpp -= -steiner.hmm3AlignWithRecalc(); //checked
            if (!rootAtExpandedEdge) {
                bpp -= -steiner.hmm2AlignWithRecalc(); //checked
            } else {
                bpp -= -drawNewAlignment(newComponent.root); //checked
            }


            dfsMoveVertex(newComponent.root, null, newComponent);
            changeOwnerVertices(newComponent, newComponent.root);

            newComponent.root.calcFelsRecursively();
            newComponent.root.calcIndelLikeRecursively();

            /*
             * Update spannoid information
             */
            final int before = spannoid.components.size();
            spannoid.components.remove(up.owner);
            spannoid.components.remove(down.owner);
            spannoid.components.add(newComponent);
            if (spannoid.components.size() != before - 1)
                throw new RuntimeException("Components not added/removed correctly!");

            int id = spannoid.labeledVertexIds.get(down);
            Set<Vertex> connections = spannoid.componentConnections.get(id);
            connections.remove(down);
            connections.remove(up);
            connections.add(labeled);

            spannoid.labeledVertexIds.remove(down);
            spannoid.labeledVertexIds.remove(up);
            spannoid.labeledVertexIds.put(labeled, id);

            checkSpannoid(spannoid);

            ExpandEdgeResult result = new ExpandEdgeResult(spannoid, up, down, labeled);
            result.bpp = bpp;

            return result;
        }

        private void rerootAtLeafNoRootAlign(Vertex vertex) {
            Vertex otherSubtree = null;
            Vertex fakeRoot = vertex.owner.root;
            final Direction direction = getDirectionFromRoot(vertex);

            if (direction.getChild(fakeRoot) != vertex) { // If root is not in correct position, do re-rooting
                // Save the old subtree s.t. references can be used for swapping.
                switch (direction) {
                    case LEFT:
                        otherSubtree = fakeRoot.right;
                        fakeRoot.right = null;
                        for (AlignColumn ac = fakeRoot.first; ac != null; ac = ac.next)
                            ac.right = null;
                        break;

                    case RIGHT:
                        otherSubtree = fakeRoot.left;
                        fakeRoot.left = null;
                        for (AlignColumn ac = fakeRoot.first; ac != null; ac = ac.next)
                            ac.left = null;
                        break;

                    default: throw new RuntimeException();
                }

                // Reroot at hte leaf
                rerootComponent(vertex, false);

                // Reinsert otherSubtree into the place where fakeRoot was located before.
                otherSubtree.parent = fakeRoot.parent;
                if (fakeRoot.parent.left == fakeRoot)
                    fakeRoot.parent.left = otherSubtree;
                else if (fakeRoot.parent.right == fakeRoot)
                    fakeRoot.parent.right = otherSubtree;
                else
                    throw new RuntimeException();
                alignAlignment(otherSubtree, fakeRoot, fakeRoot.parent);
                otherSubtree.calcAllUp();
            }
        }

        private void getRidOfRooting(Vertex leafOfTree, Vertex whereToPlace, Direction direction) {
            rerootAtLeafNoRootAlign(leafOfTree);

            Vertex brother = leafOfTree.brother();
            if (direction == Direction.LEFT) {
                whereToPlace.left = brother;
            } else {
                whereToPlace.right = brother;
            }
            brother.parent = whereToPlace;
            for (AlignColumn ac = brother.first; ac != brother.last; ac = ac.next) {
                ac.parent = whereToPlace.last;
                ac.orphan = true;
            }
            brother.last.parent = whereToPlace.last;
            brother.last.orphan = false;
            if (direction == Direction.LEFT) {
                whereToPlace.last.left = brother.last;
            } else {
                whereToPlace.last.right = brother.last;
            }
            brother.edgeLength = leafOfTree.edgeLength;
            brother.fullWin();
            brother.parent.fullWin();
            brother.edgeChangeUpdate();
            brother.parent.edgeChangeUpdate();
        }

        public void revertEdgeExpansion(ExpandEdgeResult expansion) {
            final Tree upComponent = expansion.up.owner;
            final Tree downComponent = expansion.down.owner;

            for (Vertex v : upComponent.vertex)
                restoreVertex(v);
            for (Vertex v : downComponent.vertex)
                restoreVertex(v);

            // Revert the owner of the vertices
            changeOwnerVertices(upComponent, upComponent.root);
            changeOwnerVertices(downComponent, downComponent.root);

            // TODO: Remove this check!
            for (Vertex v : upComponent.vertex)
                v.checkPointers();
            for (Vertex v : downComponent.vertex)
                v.checkPointers();

            // TODO: Consider if these are needed!
            upComponent.root.calcFelsRecursively();
            upComponent.root.calcIndelLikeRecursively();
            downComponent.root.calcFelsRecursively();
            downComponent.root.calcIndelLikeRecursively();

            /*
             * Update spannoid information
             */
            final Spannoid spannoid = expansion.spannoid;
            final int before = spannoid.components.size();
            spannoid.components.remove(expansion.labeled.owner);
            spannoid.components.add(upComponent);
            spannoid.components.add(downComponent);
            if (spannoid.components.size() != before + 1)
                throw new RuntimeException("Components not added/removed correctly!");

            int id = spannoid.labeledVertexIds.get(expansion.labeled);
            Set<Vertex> connections = spannoid.componentConnections.get(id);
            connections.remove(expansion.labeled);
            connections.add(expansion.up);
            connections.add(expansion.down);

            spannoid.labeledVertexIds.remove(expansion.labeled);
            spannoid.labeledVertexIds.put(expansion.up, id);
            spannoid.labeledVertexIds.put(expansion.down, id);

            checkSpannoid(spannoid);
        }

        private Direction getDirectionFromRoot(Vertex child) {
            Direction lastDirection = null;
            Vertex root = child.owner.root;
            while (child != root) {
                if (child.parent.left == child)
                    lastDirection = Direction.LEFT;
                else if (child.parent.right == child)
                    lastDirection = Direction.RIGHT;
                else
                    throw new RuntimeException();

                child = child.parent;
            }
            return lastDirection;
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

        public Vertex getRandomInnerBlack(Spannoid spannoid){
            List<Vertex> innerBlackNodes = spannoid.getInnerLabelledVertices();

            if (innerBlackNodes.size() == 0)
                return null;

            int j = Utils.generator.nextInt(innerBlackNodes.size());
            return innerBlackNodes.get(j);
        }

        public Vertex[] getRandomNodesForExpansion(Spannoid spannoid) {
            Vertex innerBlack = getRandomInnerBlack(spannoid);
            if (innerBlack == null)
                return null;

            int id = spannoid.labeledVertexIds.get(innerBlack);
            Vertex[] connections = spannoid.componentConnections.get(id).toArray(new Vertex[0]);

            int k = Utils.generator.nextInt(connections.length);
            int j;
            do {
                j = Utils.generator.nextInt(connections.length);
            } while (j == k);

            return new Vertex[] { connections[k], connections[j] };
        }

        private int countLabeledNodes(Tree tree) {
            final int nodes = tree.vertex.size();
            return (nodes + 1) / 2;
        }

        public List<Vertex[]> getNodesEligibleForExpansion(Spannoid spannoid, final int k) {
            List<Vertex[]> possiblePairs = new ArrayList<Vertex[]>();

            for (Set<Vertex> connections : spannoid.componentConnections.values()) {
                if (connections.size() <= 1) continue;

                for (Vertex v1 : connections) {
                    for (Vertex v2 : connections) {
                        if (v1 == v2) continue;

                        if (countLabeledNodes(v1.owner) + countLabeledNodes(v2.owner) - 1 <= k)
                            possiblePairs.add(new Vertex[] { v1, v2 });
                    }
                }
            }

            return possiblePairs;
        }

        public void checkSpannoid(Spannoid spannoid) {
            spannoid.getR();
            spannoid.getLambda();
            spannoid.getMu();

            for (Tree component : spannoid.components) {
                for (Vertex v : component.vertex)
                    v.checkPointers();
            }
        }
    }

    static class Transplanter extends MCMCMove<Spannoid, TransplantResult> {
        public Transplanter(Spannoid tree) {
            super(tree);
        }

        private boolean canSample() {
            return tree.getInnerLabelledVertices().size() > 0;
        }

        public boolean sample() {
            if (!canSample())
                return false;
            return super.sample();
        }

        @Override
        protected TransplantResult jump() {
            TransplantResult result = new TransplantResult();

            // choose source vertex
            List<Vertex> sources = tree.getInnerLabelledVertices();
            result.source = sources.get(Utils.generator.nextInt(sources.size()));
            result.bpp -= Math.log(sources.size());

            // choose destination vertex
            List<Vertex> destinations = getDestinations(result.source);
            result.destination = destinations.get(Utils.generator.nextInt(destinations.size()));
            result.bpp -= Math.log(destinations.size());

            // where to move the source vertex back to in case of rejection
            result.prev = getArbitraryNeighbour(result.source);

            result.bpp += moveSubtree(result.source, result.destination);

            result.bpp += Math.log(tree.getInnerLabelledVertices().size());
            result.bpp += Math.log(getDestinations(result.destination).size());

            return result;
        }

        @Override
        protected void restore(TransplantResult result) {
            restoreSubtree(result.source, result.prev);
        }

        private Vertex getArbitraryNeighbour(Vertex vertex) {
            int vertexId = tree.labeledVertexIds.get(vertex);
            Set<Vertex> neighbourhood = tree.componentConnections.get(vertexId);
            for (Vertex neighbour : neighbourhood) {
                if (neighbour != vertex) {
                    return neighbour;
                }
            }

            // should never reach this
            return null;
        }

        private List<Vertex> getDestinations(Vertex source) {
            List<Vertex> destinations = new ArrayList<Vertex>();

            int sourceId = tree.labeledVertexIds.get(source);
            Set<Vertex> neighbourhood = tree.componentConnections.get(sourceId);

            for (Vertex neighbour : neighbourhood) {
                if (neighbour == source) {
                    continue;
                }

                for (Vertex v : neighbour.owner.vertex) {
                    if (v != neighbour && v.name != null) {
                        destinations.add(v);
                    }
                }
            }

            return destinations;
        }

        private double moveSubtree(Vertex source, Vertex dest){
            // update internal Spannoid structure
            moveComponent(source, dest);

            double bpp = 0;

            source.fullWin();
            source.parent.fullWin();
            bpp += source.hmm2BackProp();

            // update vertex information
            source.seq = dest.seq;
            source.length = dest.length;
            source.name = dest.name;

            // save old alignment
            source.old.first = source.first;
            source.old.last = source.last;

            // create new alignment columns for new sequence
            source.first = new AlignColumn(source);
            source.first.seq = dest.first.seq.clone();
            source.first.parent = source.parent.last;
            AlignColumn prev = source.first;
            for (AlignColumn cur = dest.first.next; cur != dest.last; cur = cur.next) {
                AlignColumn actual = new AlignColumn(source);
                actual.seq = cur.seq.clone();
                actual.parent = source.parent.last;

                actual.prev = prev;
                prev.next = actual;

                prev = actual;
            }
            AlignColumn last = new AlignColumn(source);
            last.parent = source.parent.last;
            last.orphan = false;
            last.prev = prev;
            prev.next = last;
            source.last = last;

            // destroy parent alignment pointers
            for (AlignColumn c = source.parent.first; c != source.parent.last; c = c.next) {
                if (source.parent.left == source) {
                    c.left = null;
                } else {
                    c.right = null;
                }
            }

            if (source.parent.left == source) {
                source.parent.last.left = source.last;
            } else {
                source.parent.last.right = source.last;
            }

            source.fullWin();
            source.parent.fullWin();
            bpp += source.hmm2AlignWithRecalc();

            source.calcAllUp();

            return bpp;
        }

        private void restoreSubtree(Vertex source, Vertex prev){
            // update internal Spannoid structure
            moveComponent(source, prev);

            // restore old vertex information
            source.seq = prev.seq;
            source.length = prev.length;
            source.name = prev.name;

            // restore old alignment
            source.first = source.old.first;
            source.last = source.old.last;

            boolean isLeft = (source.parent.left == source);

            AlignColumn c = source.first;
            AlignColumn p = source.parent.first;
            while (p != null) {
                if (c.parent != p) {
                    if (isLeft) {
                        p.left = null;
                    } else {
                        p.right = null;
                    }
                    p = p.next;
                } else if (c.orphan) {
                    c = c.next;
                } else {
                    if (isLeft) {
                        p.left = c;
                    } else {
                        p.right = c;
                    }
                    p = p.next;
                    c = c.next;
                }
            }

            source.calcAllUp();
        }

        private void moveComponent(Vertex source, Vertex destination) {
            int sourceId = tree.labeledVertexIds.get(source);
            int destinationId = tree.labeledVertexIds.get(destination);

            tree.componentConnections.get(sourceId).remove(source);
            tree.componentConnections.get(destinationId).add(source);

            tree.labeledVertexIds.put(source, destinationId);
        }
    }
}

class TransplantResult extends MCMCResult {
    Vertex source, destination, prev;
}
