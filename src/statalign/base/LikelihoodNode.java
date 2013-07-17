package statalign.base;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/17/13
 * Time: 3:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class LikelihoodNode extends IndelNode {

    public LikelihoodNode parent;
    public LikelihoodNode left;
    public LikelihoodNode right;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /** Calculates Felsenstein and indel likelihoods up to root, starting from `parent' */
    void calcAllUp() {
        for (LikelihoodNode v = parent; v != null; v = v.parent) {
            v.calcFelsen();
            v.calcOrphan();
            v.calcIndelLogLike();
        }
    }

    void compareLike(LikelihoodNode vert) {
        if (left != null && right != null && vert.left != null && vert.right != null) {
            left.compareLike(vert.left);
            right.compareLike(vert.right);
        }
        System.out.println("this.indel=" + indelLogLike + " tree.indel=" + vert.indelLogLike);
        System.out.println("indel " + (Math.abs(indelLogLike - vert.indelLogLike) < 1e-8 ? "ok" : "error"));
        System.out.println("this.orphan=" + orphanLogLike + " tree.orphan=" + vert.orphanLogLike);
        System.out.println("orphan " + (Math.abs(orphanLogLike - vert.orphanLogLike) < 1e-8 ? "ok" : "error"));
    }



}
