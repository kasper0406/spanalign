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

    @Override
    public boolean sampleTopology() {
        int cases = Utils.generator.nextInt(3);
        switch(cases){
            case  0:
                System.out.println("Changing inner component topology");
                Tree component = updater.getRandomComponent(tree);
                return sampleTopology(component);
            case 1:
                 System.out.println("Changing outer component topology");


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
