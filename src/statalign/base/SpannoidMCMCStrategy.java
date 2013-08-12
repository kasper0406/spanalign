package statalign.base;

import java.util.List;
import java.util.Set;

public class SpannoidMCMCStrategy extends AbstractTreeMCMCStrategy<Spannoid, Spannoid.SpannoidUpdater> {
    private Transplanter transplanter;
    private final int componentSize;

    public SpannoidMCMCStrategy(Spannoid spannoid, int componentSize) {
        super(spannoid, new Spannoid.SpannoidUpdater());
        transplanter = new Transplanter(spannoid);

        this.componentSize = componentSize;
    }

    @Override
    public ITree getTree() {
        return tree;
    }

    @Override
    public boolean sampleEdge() {
        updater.checkSpannoid(tree);

        // Update random edge in random component
        // TODO: Maybe take fake root element into account?!
        Vertex vertex = updater.getRandomVertex(tree);
        boolean res = sampleEdge(vertex);

        updater.checkSpannoid(tree);
        return res;
    }

    private boolean sampleInnerTopology() {
        updater.checkSpannoid(tree);

        Tree component = updater.getRandomComponent(tree);
        boolean res = sampleTopology(component);
        updater.checkSpannoid(tree);
        return res;
    }

    private boolean sampleContract() {
        updater.checkSpannoid(tree);

        double oldLogLike = tree.getLogLike();

        List<Vertex> choices = tree.getLabelledVerticesForContraction();
        if (choices.size() == 0) {
            return false;
        }

        // choose a node and edge to contract
        Vertex labeled = choices.get(Utils.generator.nextInt(choices.size()));
        double bpp = -Math.log(choices.size());

        // backproposal of choosing the contracted edge length
        bpp += -(labeled.edgeLength - 0.01);

        Spannoid.SpannoidUpdater.ContractEdgeResult contraction = updater.contractEdge(tree, labeled);
        bpp += contraction.bpp;

        // backproposal of choosing the new node for expansion
        List<Vertex> innerBlackNodes = tree.getInnerLabelledVertices();
        Set<Vertex> neighbourhood = tree.getNeighbourhood(contraction.up);
        bpp += Math.log(2);
        // bpp += -Math.log((neighbourhood.size() - 1) * innerBlackNodes.size());
        bpp += -Math.log(updater.getNodesEligibleForExpansion(tree, componentSize).size());

        // backproposal of placing the root in this component
        int sizeOfUp = contraction.up.owner.vertex.size() - 2;
        int sizeOfDown = contraction.down.owner.vertex.size() - 2;
        bpp += -Math.log(1 + sizeOfDown + sizeOfUp);

        double newLogLike = tree.getLogLike();

        if (Math.log(Utils.generator.nextDouble()) <= bpp + (newLogLike - oldLogLike)) {
            return true;
        } else {
            updater.revertEdgeContraction(contraction);
            return false;
        }
    }

    private boolean sampleExpand() {
        updater.checkSpannoid(tree);

        double oldLogLike = tree.getLogLike();

        // Vertex[] nodes = updater.getRandomNodesForExpansion(tree);
        List<Vertex[]> possibleExpansions = updater.getNodesEligibleForExpansion(tree, componentSize);
        if (possibleExpansions.isEmpty())
            return false;

        final Vertex[] nodes = possibleExpansions.get(Utils.generator.nextInt(possibleExpansions.size()));

        double bpp = 0;

        // proposal of choosing this node and edge
        List<Vertex> innerBlackNodes = tree.getInnerLabelledVertices();
        Set<Vertex> neighbourhood = tree.getNeighbourhood(nodes[0]);
        bpp -= Math.log(2);
        // bpp -= -Math.log((neighbourhood.size() - 1) * innerBlackNodes.size());
        bpp -= -Math.log(possibleExpansions.size());

        Spannoid.SpannoidUpdater.ExpandEdgeResult expansion = updater.expandEdge(tree, nodes[0], nodes[1]);
        bpp += expansion.bpp;

        // backproposal of choosing the new edge for contraction
        bpp += Math.log(tree.getLabelledVerticesForContraction().size());

        double newLogLike = tree.getLogLike();

        if (Math.log(Utils.generator.nextDouble()) <= bpp + (newLogLike - oldLogLike)) {
            return true;
        } else {
            updater.revertEdgeExpansion(expansion);
            return false;
        }
    }

    private boolean sampleSplitMergeComponents() {
        switch (Utils.generator.nextInt(2)) {
            case 0:
                return sampleContract();
            case 1:
                return sampleExpand();
        }
        return false;
    }

    @Override
    public boolean sampleTopology() {
        updater.checkSpannoid(tree);

        boolean res = false;
        int j = Utils.generator.nextInt(3);
        switch (j) {
            case 0:
                res = sampleInnerTopology();
                break;
            case 1:
                res = transplanter.sample();
                break;
            case 2:
                res = sampleSplitMergeComponents();
                break;
        }
        updater.checkSpannoid(tree);
        return res;
    }

    @Override
    public boolean sampleAlignment() {
        updater.checkSpannoid(tree);

        // TODO: Consider picking component in which to resample alignment flipping some coin distributed by component sizes.
        Tree component = updater.getRandomComponent(tree);
        boolean res = sampleAlignment(component);

        updater.checkSpannoid(tree);
        return res;
    }
}
