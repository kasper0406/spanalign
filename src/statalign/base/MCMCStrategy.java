package statalign.base;

public interface MCMCStrategy {
    public ITree getTree();

    boolean sampleEdge();
    boolean sampleTopology();
    boolean sampleIndelParameter();
    boolean sampleAlignment();
    boolean sampleSubstParameter();
}
