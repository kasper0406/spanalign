package statalign.base;

public abstract class MCMCMove<T extends ITree, R extends MCMCResult> {
    protected T tree;

    public MCMCMove(T tree) {
        this.tree = tree;
    }

    public boolean sample() {
        double oldLogLikelihood = tree.getLogLike();

        R result = jump();

        double newLogLikelihood = tree.getLogLike();

        boolean accepted = Math.log(Utils.generator.nextDouble())
                < result.bpp + (newLogLikelihood - oldLogLikelihood) * tree.getHeat();

        if (!accepted) {
            restore(result);
        }

        return accepted;
    }

    protected abstract R jump();
    protected abstract void restore(R result);
}

class MCMCResult {
    double bpp = 0;
}