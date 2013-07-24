package statalign.base;

/**
 * Created with IntelliJ IDEA.
 * User: aldo
 * Date: 24/07/13
 * Time: 09:16
 * To change this template use File | Settings | File Templates.
 */
public abstract class AbstractUpdater<T extends ITree> implements ITreeUpdater<T>{



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
