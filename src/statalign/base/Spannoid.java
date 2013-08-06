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
     * Description...
     */
    private Map<Integer, Set<Vertex>> componentConnections = new HashMap<Integer, Set<Vertex>>();
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

    public static void main(String[] args) throws Exception {
        String[] seqs = new String[] { "AAGT", "CGATTC", "CCGAAG", "AGACA", "TTGACC", "GTAC" };
        String[] names = new String[] { "A", "B", "C", "D", "E", "F" };

        // String[] seqs = new String[] { "AAGT", "GTAC" };
        // String[] names = new String[] { "A", "B" };

        SubstitutionModel model = new Kimura3();
        SubstitutionScore ss = model.attachedScoringScheme;

        Spannoid spannoid = new Spannoid(3, BonphyStrategy.TOTAL_LENGTH, seqs, names, model, ss);
        System.out.println("Log-like of spannoid: " + spannoid.getLogLike());

        /*
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
        }
        */
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
            revertNNI(tree, nni);

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



    }

    static class Transplanter extends MCMCMove<Spannoid, TransplantResult> {
        public Transplanter(Spannoid tree) {
            super(tree);
        }

        @Override
        protected TransplantResult jump() {
            TransplantResult result = new TransplantResult();

            // choose source vertex
            List<Vertex> sources = getSources();
            result.source = sources.get(Utils.generator.nextInt(sources.size()));
            result.bpp -= Math.log(sources.size());

            // choose destination vertex
            List<Vertex> destinations = getDestinations(result.source);
            result.destination = destinations.get(Utils.generator.nextInt(destinations.size()));
            result.bpp -= Math.log(destinations.size());

            // where to move the source vertex back to in case of rejection
            result.prev = getArbitraryNeighbour(result.source);

            result.bpp += moveSubtree(result.source, result.destination);

            result.bpp += Math.log(getSources().size());
            result.bpp += Math.log(getDestinations(result.destination).size());

            return result;
        }

        @Override
        protected void restore(TransplantResult result) {
            restoreSubtree(result.source, result.prev);
        }

        private List<Vertex> getSources(){
            List<Vertex> sources = new ArrayList<Vertex>();

            for (int i = 0; i < tree.n; i++) {
                Set<Vertex> neighborhood = tree.componentConnections.get(i);
                if (neighborhood.size() > 1) {
                    sources.addAll(neighborhood);
                }
            }

            return sources;
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
            bpp += source.hmm2Align();

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