package statalign.base;

public class SteinerTreeUpdater implements ITreeUpdater<Tree> {
    private void updateParams(Tree tree, double newR, double newLambda, double newMu)
    {
        tree.hmm2.params[0] = newR;
        tree.hmm2.params[1] = newLambda;
        tree.hmm2.params[2] = newMu;

        for (int i = 0; i < tree.vertex.length; i++) {
            tree.vertex[i].updateHmmMatrices();
        }
        tree.root.calcIndelLikeRecursively();
    }

    public void updateR(Tree tree, double newR)
    {
        updateParams(tree, newR, tree.getLambda(), tree.getMu());
    }

    public void updateLambda(Tree tree, double newLambda)
    {
        updateParams(tree, tree.getR(), newLambda, tree.getMu());
    }

    public void updateMu(Tree tree, double newMu)
    {
        updateParams(tree, tree.getR(), tree.getLambda(), newMu);
    }

    @Override
    public void recalcSubstitutionParameters(Tree tree) {
        for (int i = 0; i < tree.vertex.length; i++) {
            tree.vertex[i].updateTransitionMatrix();
        }
        tree.root.calcFelsRecursively();
    }
}
