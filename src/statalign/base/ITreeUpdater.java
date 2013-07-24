package statalign.base;

public interface ITreeUpdater<T extends ITree> {
    void updateR(T tree, double newR);
    void updateLambda(T tree, double newLambda);
    void updateMu(T tree, double newMu);
    AbstractUpdater.NNIResult performNNI(Tree tree);
    void revertNNI(Tree tree, AbstractUpdater.NNIResult nni);


    void recalcSubstitutionParameters(T tree);


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
}