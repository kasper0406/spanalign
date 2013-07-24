package statalign.base;

public class SteinerTreeUpdater extends AbstractUpdater<Tree>  {
    private void updateParams(Tree tree, double newR, double newLambda, double newMu)
    {
        tree.hmm2.params[0] = newR;
        tree.hmm2.params[1] = newLambda;
        tree.hmm2.params[2] = newMu;

        for (Vertex v : tree.vertex)
            v.updateHmmMatrices();
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
        for (Vertex v : tree.vertex)
            v.updateTransitionMatrix();
        tree.root.calcFelsRecursively();
    }

    @Override
    public void revertNNI(Tree tree, AbstractUpdater.NNIResult nni){
        revertNNI(tree, nni);

    }



}
