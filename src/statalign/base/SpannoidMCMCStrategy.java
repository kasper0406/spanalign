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

    private boolean sampleMoveComponents(){
        updater.checkSpannoid(tree);

        double oldLogLi = tree.getLogLike();

        Vertex source = updater.getRandomInnerBlack(tree);
        Vertex prev = updater.getConnection(tree, source);
        Vertex dest = updater.getDestinationFromSourceMoveComponent(tree, source);

        double bpp = updater.moveSubtree(tree, source, dest);

        double newLogLi = tree.getLogLike();

        if (Math.log(Utils.generator.nextDouble()) < bpp
                + (newLogLi - oldLogLi) * tree.getHeat()) {
            updater.checkSpannoid(tree);
            return true;
        } else {
            updater.restoreSubtree(tree, source, prev);

            updater.checkSpannoid(tree);

            return false;
        }
    }

    private boolean sampleContract() {
        updater.checkSpannoid(tree);

        Vertex labeled = updater.getLabeledNodeForContractions(tree);
        if (labeled == null)
            return false;

        Spannoid.SpannoidUpdater.ContractEdgeResult contraction = updater.contractEdge(tree, labeled);

        if (Utils.generator.nextDouble() <= 0.005) {
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

        boolean res = sampleContract();
        updater.checkSpannoid(tree);
        return res;

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
        updater.checkSpannoid(tree);

        // TODO: Consider picking component in which to resample alignment flipping some coin distributed by component sizes.
        Tree component = updater.getRandomComponent(tree);
        boolean res = sampleAlignment(component);

        updater.checkSpannoid(tree);
        return res;
    }
}
