package statalign.base;

public class SteinerTreeUpdater implements ITreeUpdater<Tree> {
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

    public class NNIResult {
        double bpp;
        Vertex nephew;
        Vertex uncle;

        public NNIResult(double bpp, Vertex nephew, Vertex uncle) {
            this.bpp = bpp;
            this.nephew = nephew;
            this.uncle = uncle;
        }
    }

    public NNIResult performNNI(Tree tree) {
        int vnum = tree.vertex.size();

        if (vnum <= 3)
            return new NNIResult(Double.NEGATIVE_INFINITY, null, null);

        int vertId, rnd = Utils.generator.nextInt(vnum - 3);
        vertId = tree.getTopVertexId(rnd);
        if (vertId != -1) {
            int lastId[] = new int[3], num = 0, newId = vertId;

            for (int i = vnum - 3; i < vnum; i++) {
                int id = tree.getTopVertexId(i);
                if (id == -1)
                    lastId[num++] = i;
                else if (id < vertId)
                    newId--;
            }
            rnd = lastId[newId];
        }
        Vertex nephew = tree.vertex.get(rnd);
        Vertex uncle = nephew.parent.brother();

        double bpp = nephew.fastSwapWithUncle();

        return new NNIResult(bpp, nephew, uncle);
    }

    public void revertNNI(Tree tree, NNIResult nni) {
        if (nni.bpp == Double.NEGATIVE_INFINITY)
            return;

        nni.uncle.fastSwapBackUncle();
    }
}
