package statalign.base;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/17/13
 * Time: 3:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class FelsensteinNode  extends Hypernode{



    public FelsensteinNode parent;
    public FelsensteinNode left;
    public FelsensteinNode right;



    void updateTransitionMatrix() {
        //	System.out.println("owner: "+owner);
        //System.out.println("");
        charTransMatrix = owner.substitutionModel.updateTransitionMatrix(charTransMatrix, edgeLength);
        charPropTransMatrix = new double[charTransMatrix.length][];


    }





    /**
     * This function calculates the Felsenstein likelihood. The result is stored
     * in the orphanLogLike at the root.
     */
    void calcFelsRecursively() {
        if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.calcFelsRecursively();
            //System.out.println("calling the right child");
            right.calcFelsRecursively();
        }
        calcFelsen();
        calcOrphan();
    }




    /** Calculates Felsenstein likelihoods of `this' */
    void calcFelsen() {
        if (left != null && right != null) {
            AlignColumn p;
            AlignColumn l = left.first;
            AlignColumn r = right.first;
            double[] fel1, fel2;
            int columnIndex = 0;
            for (p = first; p != last; p = p.next) {
                fel1 = null;
                fel2 = null;
                while (l != left.last && l.orphan)
                    l = l.next;
                while (r != right.last && r.orphan)
                    r = r.next;
                if (l.parent == p) {
                    fel1 = l.seq;
                    l = l.next;
                }
                if (r.parent == p) {
                    fel2 = r.seq;
                    r = r.next;
                }
                Utils.calcFelsen(p.seq, fel1, left.charTransMatrix, fel2, right.charTransMatrix);
                // if sequence at this Hypernode is known, zero out incompatible likelihoods
                if (isLabeled()) {
                    for (int i = 0; i < p.seq.length; i++) {
                        if (owner.substitutionModel.alphabet[i]!=seq.charAt(columnIndex)) {
                            p.seq[i] = 0;
                        }
                    }
                }
                columnIndex++;
            }
        }
    }

    /** Calculates Felsenstein likelihoods of `this' */
    void calcFelsenWithCheck() {
        if (left != null && right != null) {
            AlignColumn p;
            AlignColumn l = left.first;
            AlignColumn r = right.first;
            double[] fel1, fel2;
            for (p = first; p != last; p = p.next) {
                fel1 = null;
                fel2 = null;
                while (l != left.last && l.orphan)
                    l = l.next;
                while (r != right.last && r.orphan)
                    r = r.next;
                if (l.parent == p) {
                    fel1 = l.seq;
                    l = l.next;
                }
                if (r.parent == p) {
                    fel2 = r.seq;
                    r = r.next;
                }
                boolean match = Utils.calcFelsenWithCheck(p.seq, fel1, left.charTransMatrix, fel2, right.charTransMatrix);
                if (!match) {
                    throw new Error("Felsensteins do not match!");
                }
            }
        }
    }




    /**
     * Calculates the sum of orphan likelihoods in the (inclusive) subtree of `this'
     * which is then stored in `this.orphanLogLikge' (this implies that alignment to
     * parent must exists at the time of calling)
     * Saves previous orphan likelihood into `old' so it shouldn't be called twice in a row.
     * (But does not save previous Felsensteins of `this' in AlignColumn's `seq'.)
     */
    void calcOrphan() {
        old.orphanLogLike = orphanLogLike;

        //orphan likelihoods
        orphanLogLike = 0.0;

        for (AlignColumn actual = first; actual != last; actual = actual.next) {
            if (actual.parent == null || actual.orphan) {
                orphanLogLike += Math.log(Utils.calcEmProb(actual.seq, owner.substitutionModel.e));
            }
        }

        if (left != null && right != null)
            orphanLogLike += left.orphanLogLike + right.orphanLogLike;
    }



    /** Returns the Felsenstein loglikelihood vector for the sequence in this vertex. */
    double[][] getFelsen() {
        double[][] felsen = new double[length][];

        int i = 0;
        for (AlignColumn c = first; c != last; c = c.next) {
            felsen[i++] = c.seq.clone();
        }
        return felsen;
    }




}
