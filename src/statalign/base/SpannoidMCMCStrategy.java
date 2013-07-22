package statalign.base;

public class SpannoidMCMCStrategy implements MCMCStrategy {
    private Spannoid spannoid;

    public SpannoidMCMCStrategy(Spannoid spannoid) {
        this.spannoid = spannoid;
    }

    @Override
    public ITree getTree() {
        return spannoid;
    }

    @Override
    public boolean sampleEdge() {
        return false;
    }

    @Override
    public boolean sampleTopology() {
        return false;
    }

    @Override
    public boolean sampleIndelParameter() {
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
