package statalign.base;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/17/13
 * Time: 3:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class IndelNode extends FelsensteinNode {


  
   double[][] hmm2TransMatrix;            // precalculated state transition likelihoods for 2-seq HMM (indel model)
   double[][] hmm2PropTransMatrix;        // precalculated state transition likelihoods for 2-seq HMM used for proposals (indel model)

    /**
     * The log-sum of the cumulative insertion-deletion loglikelihoods up to this vertex (ie. summed over the
     * subtree below this vertex.).
     */
    public double indelLogLike;





    void calcIndelLikeRecursively() {
        if (left != null && right != null) {
            left.calcIndelLikeRecursively();
            right.calcIndelLikeRecursively();
        }
        calcIndelLogLike();
    }

    void calcIndelRecursivelyWithCheck() {
        if (left != null && right != null) {
            left.calcIndelRecursivelyWithCheck();
            right.calcIndelRecursivelyWithCheck();
        }
        calcIndelLikeWithCheck();
    }

    void calcIndelLikeWithCheck() {
        double newindelLogLike = 0.0;
        if (left != null && right != null) {
            newindelLogLike += left.calcIndelLogUpWithCheck();
            newindelLogLike += right.calcIndelLogUpWithCheck();
        }
        if (Math.abs(indelLogLike - newindelLogLike) > 1e-5) {
            throw new Error("problem: fast indel: " + indelLogLike + " slow indel: " + newindelLogLike);
        }
    }


    double calcIndelLogUpWithCheck() {
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        final int emitPatt2State[] = owner.hmm2.getEmitPatt2State();

        double indelLogLikeUp = indelLogLike;

        if (parent != null) {
            AlignColumn c = first, p = parent.first;
            int prevk = START, k;

            while (c != last || p != parent.last) {
                if (c.parent != p) {                            // deletion (* -), pattern code 2
                    k = emitPatt2State[2];
                    p = p.next;
                } else if (c.orphan) {                        // insertion (- *), pattern code 1
                    k = emitPatt2State[1];
                    c = c.next;
                } else {                                                // substitution (* *), pattern code 3
                    k = emitPatt2State[3];
                    p = p.next;
                    c = c.next;
                }
                indelLogLikeUp += hmm2TransMatrix[prevk][k];
                prevk = k;
            }

            indelLogLikeUp += hmm2TransMatrix[prevk][END];
        }

        return indelLogLikeUp;
    }


    /**
     * Function to calculate the indel logLikelihood of the subtree of `this' including
     * the indel logLikelihood of the alignment between `this.left'&`this' and `this.right'&`this'.
     * Saves previous logLikelihood into `old' Vertex, so it shouldn't be called twice in a row.
     * Assumes it has been called previously on both `left' and `right'.
     * Result is stored in `this'.
     */
    void calcIndelLogLike() {
        old.indelLogLike = indelLogLike;
        indelLogLike = 0.0;
        if (left != null && right != null) {
            indelLogLike += left.calcIndelLogLikeUp();
            indelLogLike += right.calcIndelLogLikeUp();
        }
    }

    /**
     * Function to calculate the indel log-likelihood of the subtree of `this' plus
     * the indel log-likelihood of the alignment between `this' & `this.parent'
     * Assumes log-likelihood of subtree is already pre-calculated and stored in `this'
     * @return counted log-likelihood
     */
    double calcIndelLogLikeUp() {
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        final int emitPatt2State[] = owner.hmm2.getEmitPatt2State();

        double indelLogLikeUp = indelLogLike;

        //System.out.println("--------------------------------------------------");
        //printPointers();

        if (parent != null) {
            AlignColumn c = first, p = parent.first;
            int prevk = START, k;

            while (c != last || p != parent.last) {
                if (c.parent != p) {                            // deletion (* -), pattern code 2
                    k = emitPatt2State[2];
                    p = p.next;
                } else if (c.orphan) {                        // insertion (- *), pattern code 1
                    k = emitPatt2State[1];
                    c = c.next;
                } else {                                                // substitution (* *), pattern code 3
                    k = emitPatt2State[3];
                    p = p.next;
                    c = c.next;
                }
                indelLogLikeUp += hmm2TransMatrix[prevk][k];
                prevk = k;
            }

            indelLogLikeUp += hmm2TransMatrix[prevk][END];
        }

        return indelLogLikeUp;
    }

    ////////////////////// HMMS ////////////////////////

    void updateHmm2Matrix() {
        hmm2TransMatrix = owner.hmm2.preCalcTransMatrix(hmm2TransMatrix, edgeLength);

        // This is commented out for now. Uncomment to use the rescaling of the transition matrices.
        // Do the same for {@link #updateHmm3Matrix()}.

//        double heat = owner.heat;
        hmm2PropTransMatrix = new double[hmm2TransMatrix.length][];
        for (int i = 0; i < hmm2TransMatrix.length; i++) {
            hmm2PropTransMatrix[i] = hmm2TransMatrix[i].clone();

        }


        }




///////////////////////////////////// HMM management //////////////////////////////

    void updateHmmMatrices() {
        if (parent != null)
            updateHmm2Matrix();

    }

    void edgeChangeUpdate() {
        if (parent != null) {
            updateTransitionMatrix();
            updateHmm2Matrix();

        }
    }


}
