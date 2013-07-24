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
        Vertex source = updater.getRandomInnerBlack(tree);
        Tree adj = updater.getRandomNeighbouringComponent(tree, source);
        Vertex dest = updater.getComponentRandomVertex(adj);
        updater.moveComponent(tree, source, dest);
        return true;

    }


    @Override
    public boolean sampleIndelParameter() {
        return super.sampleIndelParameter();    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public boolean sampleTopology() {
        int j = Utils.generator.nextInt(3);
        switch (j) {
            case 0:
                return sampleInnerTopology();
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
