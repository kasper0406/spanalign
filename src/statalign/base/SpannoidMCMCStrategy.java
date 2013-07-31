package statalign.base;

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

        double bpp = updater.moveSubtree(tree, source, dest);

        double newLogLi = tree.getLogLike();

        if (Math.log(Utils.generator.nextDouble()) < bpp
                + (newLogLi - oldLogLi) * tree.getHeat()) {
            return true;
        } else {
            updater.restoreSubtree(tree, source, prev);

            return false;
        }
    }

    /*
    @Override
    public boolean sampleIndelParameter() {
        return super.sampleIndelParameter();    //To change body of overridden methods use File | Settings | File Templates.
    }
    */

    @Override
    public boolean sampleTopology() {
        int j = Utils.generator.nextInt(3);
        switch (j) {
            case 0:
                // return sampleInnerTopology();
                return false;
            case 1:
                return sampleMoveComponents();
            default:
                return false;
        }
    }

    @Override
    public boolean sampleAlignment() {
        // TODO: Consider picking component in which to resample alignment flipping some coin distributed by component sizes.
        Tree component = updater.getRandomComponent(tree);
        return sampleAlignment(component);
    }
}
