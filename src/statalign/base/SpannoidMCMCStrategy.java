package statalign.base;

/**
 * !!!ATTENTION!!!
 *
 * TOPOLOGY CHANGE ACCEPTANCE RATE PROBABLY NOT CORRECT!
 *
 * !!!ATTENTION!!!
 */

public class SpannoidMCMCStrategy extends AbstractTreeMCMCStrategy<Spannoid, Spannoid.SpannoidUpdater> {
    public SpannoidMCMCStrategy(Spannoid spannoid) {
        super(spannoid, new Spannoid.SpannoidUpdater());
    }

    @Override
    public ITree getTree() {
        return tree;
    }

    @Override
    public boolean sampleEdge() {
        // Update random edge in random component
        // TODO: Maybe take fake root element into account?!
        Vertex vertex = updater.getRandomVertex(tree);
        return sampleEdge(vertex);
    }


    private boolean sampleInnerTopology() {
        Tree component = updater.getRandomComponent(tree);
        return sampleTopology(component);
    }

    private boolean sampleMoveComponents(){
        double oldLogLi = tree.getLogLike();

        Vertex source = updater.getRandomInnerBlack(tree);
        Vertex prev = updater.getConnection(tree, source);
        Vertex dest = updater.getDestinationFromSourceMoveComponent(tree, source);

        double bpp = updater.moveComponent(tree, source, dest);

        double newLogLi = tree.getLogLike();

        // TODO: INCORRECT ACCEPTANCE RATE
        if (Math.log(Utils.generator.nextDouble()) < bpp
                + (newLogLi - oldLogLi) * tree.getHeat()) {
            return true;
        } else {
            updater.moveComponent(tree, source, prev);

            return false;
        }
    }

    private boolean sampleContract() {
        Vertex labeled = updater.getLabeledNodeForContractions(tree);

        Spannoid.SpannoidUpdater.ContractEdgeResult contraction = updater.contractEdge(tree, labeled);

        if (Utils.generator.nextDouble() <= 0.5) {
            return true;
        } else {
            updater.revertEdgeContraction(contraction);
            return false;
        }
    }

    private boolean sampleExpand() {
        return false;
    }

    private boolean sampleSplitMergeComponents() {
        return sampleContract();

        /*
        switch (Utils.generator.nextInt(2)) {
            case 0:
                return sampleContract();
            case 1:
                return sampleExpand();
        }
        return false;
        */
    }

    @Override
    public boolean sampleTopology() {
        return sampleContract();

        /*
        int j = Utils.generator.nextInt(3);
        switch (j) {
            case 0:
                return sampleInnerTopology();
            case 1:
                return sampleMoveComponents();
            case 2:
                return sampleSplitMergeComponents();
        }
        return false;                         */
    }

    @Override
    public boolean sampleAlignment() {
        // TODO: Consider picking component in which to resample alignment flipping some coin distributed by component sizes.
        Tree component = updater.getRandomComponent(tree);
        return sampleAlignment(component);
    }
}
