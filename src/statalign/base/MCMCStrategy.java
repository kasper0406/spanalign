package statalign.base;

public interface MCMCStrategy {
    boolean sampleEdge();
    boolean sampleTopology();
    boolean sampleIndelParameter();
    boolean sampleAlignment();
    boolean sampleSubstParameter();
}
