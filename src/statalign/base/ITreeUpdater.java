package statalign.base;

public interface ITreeUpdater<T extends ITree> {
    void updateR(T tree, double newR);
    void updateLambda(T tree, double newLambda);
    void updateMu(T tree, double newMu);

    void recalcSubstitutionParameters(T tree);
}
