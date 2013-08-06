package statalign.base;

public class SpannoidMCMCStrategy extends AbstractTreeMCMCStrategy<Spannoid, Spannoid.SpannoidUpdater> {
    private Spannoid.Transplanter transplanter;

    public SpannoidMCMCStrategy(Spannoid spannoid) {
        super(spannoid, new Spannoid.SpannoidUpdater());
        transplanter = new Spannoid.Transplanter(spannoid);
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
                return transplanter.sample();
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
