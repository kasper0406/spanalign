package statalign.base;

public class SteinerTreeMCMCStrategy extends AbstractTreeMCMCStrategy<Tree, SteinerTreeUpdater> {
    public SteinerTreeMCMCStrategy(Tree tree) {
        super(tree, new SteinerTreeUpdater());

        this.tree = tree;
    }

    @Override
    public boolean sampleEdge() {
        int i = Utils.generator.nextInt(tree.vertex.size());
        return sampleEdge(tree.vertex.get(i));
    }

    @Override
    public boolean sampleTopology() {
        return sampleTopology(tree);
    }

    @Override
    public boolean sampleAlignment() {
        return sampleAlignment(tree);
    }
}
