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
        Vertex vertex = tree.getRandomVertex();
        return sampleEdge(vertex);
    }

    @Override
    public boolean sampleTopology() {
        return false;
    }

    @Override
    public boolean sampleAlignment() {
        return false;
    }

    @Override
    public boolean sampleSubstParameter() {
        return false;
    }
}
