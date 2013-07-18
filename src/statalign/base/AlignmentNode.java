package statalign.base;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/17/13
 * Time: 3:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentNode extends LikelihoodNode {

    public AlignmentNode parent;
    public AlignmentNode left;
    public AlignmentNode right;
    public AlignmentNode old;

    double[][] hmm3TransMatrix;            // precalculated state transition likelihoods for 3-seq HMM (indel model)
    double[][] hmm3RedTransMatrix;    // precalculated st. trans. likelihoods for 3-seq HMM, silent st. removed (indel model)


    int winLength;                    // length of window
    AlignColumn winFirst;        // first alignment column of window
    AlignColumn winLast;        // first alignment column past window end
    boolean selected;                // shows if vertex is part of the selected subtree




    AlignmentNode() {
    }


    AlignmentNode brother() {
        return parent.left == this ? parent.right : parent.left;
    }

    void updateTransitionMatrix() {
        //	System.out.println("owner: "+owner);
        //System.out.println("");
        charTransMatrix = owner.substitutionModel.updateTransitionMatrix(charTransMatrix, edgeLength);
        charPropTransMatrix = new double[charTransMatrix.length][];


    }

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

    void updateHmm3Matrix() {
        if (left != null && right != null) {
            hmm3TransMatrix = owner.hmm3.preCalcTransMatrix(hmm3TransMatrix, left.edgeLength, right.edgeLength);
            hmm3RedTransMatrix = owner.hmm3.preCalcRedTransMatrix(hmm3RedTransMatrix, hmm3TransMatrix);


        }
    }

    void updateHmmMatrices() {
        if (parent != null)
            updateHmm2Matrix();
        updateHmm3Matrix();
    }

    void edgeChangeUpdate() {
        if (parent != null) {
            updateTransitionMatrix();
            updateHmm2Matrix();
            parent.updateHmm3Matrix();
        }
    }



    void fullWin() {
        winFirst = first;
        winLast = last;
        winLength = length;
    }




    private double[][][] hmm3ProbMatrix() {
        double[] equDist = owner.substitutionModel.e;
        int hmm3Parent[] = owner.hmm3.getStateEmit()[0];
        int hmm3Left[] = owner.hmm3.getStateEmit()[1];
        int hmm3Right[] = owner.hmm3.getStateEmit()[2];
        int leftLen = left.winLength, rightLen = right.winLength;
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();

        double probMatrix[][][];                                                                // DP matrix used for 3-seq HMM alignment
        probMatrix = new double[leftLen + 1][rightLen + 1][END];        // don't reserve space for end state

        double emissionProb, felsen[] = new double[equDist.length], tr;
        AlignColumn l = null, r = null;
        int i, j, k, previ, prevj, prevk;

        /* i: left prefix length, j: right prefix length, k: state
	   l: left child's actual alignment column, r: that of right child */
        for (i = 0; i <= leftLen; i++) {
            for (j = 0; j <= rightLen; j++) {
                probMatrix[i][j][START] = (i == 0 && j == 0) ? 0.0 : Utils.log0;        // update matrix for start state
                for (k = START + 1; k < END; k++) {
                    previ = i - hmm3Left[k];
                    prevj = j - hmm3Right[k];
                    if (previ >= 0 && prevj >= 0) {
                        if (hmm3Parent[k] != 0) {
                            // there's a parent character in alignment column, but still can be an (even double-) deletion
                            Utils.calcFelsen(felsen, hmm3Left[k] != 0 ? l.seq : null, left.charTransMatrix, hmm3Right[k] != 0 ? r.seq : null, right.charTransMatrix);
                            emissionProb = Utils.calcEmProb(felsen, equDist);
                        } else {
                            // no parent, it's an insertion on either of (but never both) branches
                            emissionProb = Utils.calcEmProb(hmm3Left[k] != 0 ? l.seq : r.seq, equDist);
                        }
                        if (previ == 0 && prevj == 0)
                            tr = hmm3RedTransMatrix[START][k];
                        else {
                            for (tr = Utils.log0, prevk = START + 1; prevk < END; prevk++)
                                tr = Utils.logAdd(tr, probMatrix[previ][prevj][prevk] + hmm3RedTransMatrix[prevk][k]);
                        }
                        probMatrix[i][j][k] = Math.log(emissionProb) + tr;
                    } else
                        probMatrix[i][j][k] = Utils.log0;
                }
                r = (j > 0) ? r.next : right.winFirst;
            }
            l = (i > 0) ? l.next : left.winFirst;
        }

        return probMatrix;
    }

    /**
     * Computes backproposal of alignment between this & this.left & this.right, window sizes matter.
     * @return log of backproposal probability
     */
    double hmm3BackProp() {
        int emitPatt2State[] = owner.hmm3.getEmitPatt2State();
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();
        final int SILENT = owner.hmm3.getSilent();

        AlignColumn l, r, p;
        double retVal = 0.0;
        int k, previ, prevj, prevk;
        int pattCode;                                    // alignment column pattern's binary code

        double probMatrix[][][] = hmm3ProbMatrix();

        /* backproposal calculation */

        previ = prevj = 0;                        // previous prefix lengths (0 at the beginning)
        prevk = START;                                // previous state (start state at the beginning)
        l = left.winFirst;                        // starting alignment columns (window starts)
        r = right.winFirst;
        p = winFirst;
        int silentNum = 0;                        // number of transitions through silent state since last non-silent state

        while (prevk != END) {
            if (l == left.winLast && r == right.winLast && p == winLast)        // go to end state (virtual pattern code: 8)
                pattCode = 8;
            else if (l.parent == p && l.orphan)                    // left insertion (- * -), emission pattern binary code: 2
                pattCode = 2;
            else if (r.parent == p && r.orphan)                    // right insertion (- - *), emission pattern binary code: 1
                pattCode = 1;
            else if (l.parent != p && r.parent != p) {        // double deletion (* - -), silent state, emission pattern binary code: 4
                silentNum++;
                p = p.next;
                continue;
            } else {
                pattCode = 4;                            // there is a parent in all the remaining cases (* ? ?)
                if (l.parent == p)
                    pattCode += 2;                    // could be left/right deletion (* - *)/(* * -) with codes 5/6 or
                if (r.parent == p)
                    pattCode += 1;                    // "full" alignment column (* * *), binary code 7
            }

            k = emitPatt2State[pattCode];            // real state number for column

            double bp = Utils.log0;
            for (int bpPrevk = START; bpPrevk < END; bpPrevk++)
                bp = Utils.logAdd(bp, probMatrix[previ][prevj][bpPrevk] + hmm3RedTransMatrix[bpPrevk][k]);
            bp = probMatrix[previ][prevj][prevk] - bp;
            if (silentNum == 0)            // non-reducated transition, skipping silent state
                bp += hmm3TransMatrix[prevk][k];
            else                                        // non-reducated transitions, passing through silent state at least once
                bp += hmm3TransMatrix[prevk][SILENT] + hmm3TransMatrix[SILENT][k] + (silentNum - 1) * hmm3TransMatrix[SILENT][SILENT];
            if (bp > 1e-5) {
                //System.out.println("Pozitiv!");
            }

            retVal += bp;

            silentNum = 0;
            prevk = k;
            if ((pattCode & 4) != 0)                // make a parent's step
                p = p.next;
            if ((pattCode & 2) != 0) {            // make a left child's step
                l = l.next;
                previ++;
            }
            if ((pattCode & 1) != 0) {            // make a right child's step
                r = r.next;
                prevj++;
            }
        }

        return retVal;
    }

    protected void saveWin() {
        old.winFirst = winFirst;
        old.winLast = winLast.prev;                // temporary reference to last (real) window element (but could be null)
        old.winLength = winLength;
        old.first = first;
        old.length = length;
        old.orphanLogLike = orphanLogLike;
        old.indelLogLike = indelLogLike;
    }

    /**
     * Toggles parent alignment pointers of `this' between parent's old and new sequence (use `toNew' to select).
     * When called, parent's AlignColumn pointers must always point to the new sequence, regardless of `toNew'.
     */
    void toggleUp(boolean toNew) {
        AlignColumn dp = toNew ? parent.winFirst : parent.old.winFirst;        // "destination" parent
        AlignColumn c = first;
        boolean inWin = false;

        AlignColumn p = parent.first;
        while (p != parent.winLast) {
            if (p == parent.winFirst) {
                if (toNew)
                    p = parent.old.winFirst;
                inWin = true;
            }
            if (c.parent != p) {
                if (inWin)
                    dp = dp.next;
                p = p.next;
            } else {
                if (inWin)
                    c.parent = dp;
                c = c.next;
            }
        }
    }

    /**
     * Samples a new alignment between `this' &amp; `this.left' &amp; `this.right', taking window sizes
     * into account. Creates new ancestor sequence, updates all Vertex-stored suppl. data, <b>does not save</b> old
     * data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    double hmm3Align() {
        int hmm3Parent[] = owner.hmm3.getStateEmit()[0];
        int hmm3Left[] = owner.hmm3.getStateEmit()[1];
        int hmm3Right[] = owner.hmm3.getStateEmit()[2];
        int leftLen = left.winLength, rightLen = right.winLength;
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();
        final int SILENT = owner.hmm3.getSilent();

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm3ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];            // no need to have an element for end state

        /* stochastic traceback */

        old.winLength = winLength;
        winLength = 0;
        //	AlignColumn wfp = winFirst.prev;				// must save this in case winFirst == winLast


        previ = leftLen;
        prevj = rightLen;
        AlignColumn l = left.winLast, r = right.winLast;
        AlignColumn p = winLast;

        if (p == null) {
            throw new Error("Starting p is null!!!");
        }

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm3RedTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);

            if (Utils.chooseOne(Math.exp(hmm3TransMatrix[prevk][k] - hmm3RedTransMatrix[prevk][k]), retVal) == 0) {
                do {
                    p = new AlignColumn(p, true);
                    winLength++;
                } while (Utils.chooseOne(Math.exp(hmm3TransMatrix[SILENT][SILENT]), retVal) == 1);
            }
            if (hmm3Parent[prevk] != 0) {
                p = new AlignColumn(p, true);
                winLength++;
            }
            if (hmm3Left[prevk] != 0) {
                l = l.prev;
                l.parent = p;
                if (l.parent == null) {
                    throw new Error("l.parent is null!!!");
                }
                l.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.left = l;
                previ--;
            }
            if (hmm3Right[prevk] != 0) {
                r = r.prev;
                r.parent = p;
                if (r.parent == null) {
                    throw new Error("r.parent is null!!!");
                }
                r.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.right = r;
                prevj--;
            }
        }
        p.setWinFirst(null);
        length += winLength - old.winLength;



        return retVal.value;
    }


    /**
     * Samples a new alignment between `this' &amp; `this.left' &amp; `this.right', taking window sizes
     * into account. Creates new ancestor sequence, updates all Vertex-stored suppl. data, saving old
     * data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    double hmm3AlignWithSave() {
        int hmm3Parent[] = owner.hmm3.getStateEmit()[0];
        int hmm3Left[] = owner.hmm3.getStateEmit()[1];
        int hmm3Right[] = owner.hmm3.getStateEmit()[2];
        int leftLen = left.winLength, rightLen = right.winLength;
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();
        final int SILENT = owner.hmm3.getSilent();

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm3ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];            // no need to have an element for end state

        /* stochastic traceback */

        saveWin();
        winLength = 0;
        AlignColumn wfp = winFirst.prev;                // must save this in case winFirst == winLast

        AlignColumn ol = left.winLast, or = right.winLast;
        boolean saveLeft = false, saveRight = false;
        if (left.left == null || left.right == null || !left.left.selected || !left.right.selected) {
            left.saveWin();
            saveLeft = true;
        }
        if (right.left == null || right.right == null || !right.left.selected || !right.right.selected) {
            right.saveWin();
            saveRight = true;
        }

        previ = leftLen;
        prevj = rightLen;
        AlignColumn l = left.winLast, r = right.winLast;
        AlignColumn p = winLast;

        if (p == null) {
            throw new Error("Starting p is null!!!");
        }

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm3RedTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);

            if (Utils.chooseOne(Math.exp(hmm3TransMatrix[prevk][k] - hmm3RedTransMatrix[prevk][k]), retVal) == 0) {
                do {
                    p = new AlignColumn(p, true);
                    winLength++;
                } while (Utils.chooseOne(Math.exp(hmm3TransMatrix[SILENT][SILENT]), retVal) == 1);
            }
            if (hmm3Parent[prevk] != 0) {
                p = new AlignColumn(p, true);
                winLength++;
            }
            if (hmm3Left[prevk] != 0) {
                if (saveLeft) {
                    ol = ol.prev;
                    l = new AlignColumn(l, false);
                    l.saveBoth(ol);
                } else {
                    l = l.prev;
                }
                l.parent = p;
                if (l.parent == null) {
                    throw new Error("l.parent is null!!!");
                }
                l.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.left = l;
                previ--;
            }
            if (hmm3Right[prevk] != 0) {
                if (saveRight) {
                    or = or.prev;
                    r = new AlignColumn(r, false);
                    r.saveBoth(or);
                } else {
                    r = r.prev;
                }
                r.parent = p;
                if (r.parent == null) {
                    throw new Error("r.parent is null!!!");
                }
                r.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.right = r;
                prevj--;
            }
        }
        p.setWinFirst(wfp);
        length += winLength - old.winLength;

        //		System.out.println("saveLeft: "+saveLeft+" saveRight: "+saveRight+" length: "+length);

        if (saveLeft) {
            l.setWinFirst(left.winFirst.prev);
            if (left.left != null && left.right != null) {
                left.left.toggleUp(true);
                left.right.toggleUp(true);
            }
        }
        if (saveRight) {
            r.setWinFirst(right.winFirst.prev);
            if (right.left != null && right.right != null) {
                right.left.toggleUp(true);
                right.right.toggleUp(true);
            }
        }

        //	left.printPointers();
        //right.printPointers();

        left.calcOrphan();
        right.calcOrphan();
        calcFelsen();
        calcIndelLogLike();

        return retVal.value;
    }

    private double[][][] hmm2ProbMatrix() {
        double[] equDist = owner.substitutionModel.e;
        int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        int hmm2Child[] = owner.hmm2.getStateEmit()[1];
        int parentLen = parent.winLength, childLen = winLength;
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        AlignmentNode brother = brother();

        double probMatrix[][][];                                                                    // DP matrix used for 2-seq HMM alignment
        probMatrix = new double[parentLen + 1][childLen + 1][END];        // don't reserve space for end state

        double emissionProb, felsen[] = new double[equDist.length], tr;
        AlignColumn c = null;                // child
        AlignColumn p = null;                // parent
        AlignColumn b;                            // brother
        int i, j, k, previ, prevj, prevk;

        /* i: parent prefix length, j: child prefix length, k: state  */
        for (i = 0; i <= parentLen; i++) {
            for (j = 0; j <= childLen; j++) {
                probMatrix[i][j][START] = (i == 0 && j == 0) ? 0.0 : Utils.log0;        // update matrix for start state
                for (k = START + 1; k < END; k++) {
                    previ = i - hmm2Parent[k];
                    prevj = j - hmm2Child[k];
                    if (previ >= 0 && prevj >= 0) {
                        if (hmm2Parent[k] != 0) {                // parent present: substitution (* *) or deletion (* -)
                            b = (this == parent.left) ? p.right : p.left;
                            Utils.calcFelsen(felsen, hmm2Child[k] != 0 ? c.seq : null, charTransMatrix, b != null ? b.seq : null, brother.charTransMatrix);
                            emissionProb = Utils.calcEmProb(felsen, equDist);
                        } else {                    // insertion (- *)
                            emissionProb = Utils.calcEmProb(c.seq, equDist);
                        }
                        if (previ == 0 && prevj == 0)
                            tr = hmm2PropTransMatrix[START][k];
                        else {
                            for (tr = Utils.log0, prevk = START + 1; prevk < END; prevk++)
                                tr = Utils.logAdd(tr, probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k]);
                        }
                        probMatrix[i][j][k] = Math.log(emissionProb) + tr;
                    } else
                        probMatrix[i][j][k] = Utils.log0;
                }
                c = (j > 0) ? c.next : winFirst;
            }
            p = (i > 0) ? p.next : parent.winFirst;
        }

        return probMatrix;
    }

    /**
     * Computes backproposal of alignment between `this' & `this.parent', window sizes matter.
     * @return log of backproposal probability
     */
    double hmm2BackProp() {
        int emitPatt2State[] = owner.hmm2.getEmitPatt2State();
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();

        AlignColumn c, p;
        double retVal = 0.0;
        int k, previ, prevj, prevk;
        int pattCode;                                    // alignment column pattern's binary code

        double probMatrix[][][] = hmm2ProbMatrix();

        /* backproposal calculation */

        previ = prevj = 0;                        // previous prefix lengths (0 at the beginning)
        prevk = START;                                // previous state (start state at the beginning)
        c = winFirst;                                    // starting alignment columns (window starts)
        p = parent.winFirst;

        while (prevk != END) {
            if (c == winLast && p == parent.winLast)            // go to end state (virtual pattern code 4)
                pattCode = 4;
            else if (c.parent != p)                                             // deletion (* -), pattern code 2
                pattCode = 2;
            else if (c.orphan)                                                        // insertion (- *), pattern code 1
                pattCode = 1;
            else                                                                                // substitution (* *), pattern code 3
                pattCode = 3;

            k = emitPatt2State[pattCode];                                // real state number for column

            double bp = Utils.log0;
            for (int bpPrevk = START; bpPrevk < END; bpPrevk++)
                bp = Utils.logAdd(bp, probMatrix[previ][prevj][bpPrevk] + hmm2PropTransMatrix[bpPrevk][k]);
            bp = probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k] - bp;

            retVal += bp;

            if ((pattCode & 2) != 0) {            // make a parent's step
                p = p.next;
                previ++;
            }
            if ((pattCode & 1) != 0) {            // make a child's step
                c = c.next;
                prevj++;
            }
            prevk = k;
        }

        return retVal;
    }

    /**
     * Samples a new alignment between `this' &amp; `this.parent', taking window sizes into account.
     * <b>Does not update</b> any Vertex-stored suppl. data, <b>does not save</b> old data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    double hmm2Align() {
        int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        int hmm2Child[] = owner.hmm2.getStateEmit()[1];
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        int parentLen = parent.winLength, childLen = winLength;
        boolean isLeft = parent.left == this;
        //System.out.println("Aligning "+(isLeft ? "left" : "right")+" in hmm2Align");

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm2ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];                        // no need to have an element for end state

        /* stochastic traceback */

        //	parent.saveWin();
        parent.winLength = 0;

        //	AlignColumn oc = winLast;				       			// original child
        //	AlignColumn op = parent.winLast;		     				// original parent

        previ = parentLen;
        prevj = childLen;
        AlignColumn c = winLast;                                   // child
        AlignColumn p = parent.winLast;                                  // parent

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);

            if (hmm2Parent[prevk] != 0) {
                p = p.prev;
                /*
		  if(isLeft){
		  p.left = null;
		  }
		  else{
		  p.right = null;
		  }
				 */
                previ--;
                parent.winLength++;
            }
            if (hmm2Child[prevk] != 0) {
                c = c.prev;
                c.parent = p;
                c.orphan = hmm2Parent[prevk] == 0;
                if (hmm2Parent[prevk] != 0) {
                    if (isLeft)
                        p.left = c;
                    else
                        p.right = c;
                }
                prevj--;
            }
        }

        //	p.setWinFirst(parent.winFirst.prev);
        //	assert (parent.winLength == parent.old.winLength) : "HMM2 alignment error";

        //for(c = winFirst.prev; c != null && c.parent == parent.old.winFirst; c = c.prev)
        //  c.parent = parent.winFirst;

        //calcOrphan();
        //parent.calcFelsen();
        //parent.calcOrphan();
        //parent.calcIndelLogLike();

        return retVal.value;
    }

    /**
     * Samples a new alignment between `this' &amp; `this.parent', taking window sizes into account.
     * Updates all Vertex-stored suppl. data, saving old data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    public double hmm2AlignWithSave() {
        int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        int hmm2Child[] = owner.hmm2.getStateEmit()[1];
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        int parentLen = parent.winLength, childLen = winLength;
        boolean isLeft = parent.left == this;

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm2ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];                        // no need to have an element for end state

        /* stochastic traceback */

        parent.saveWin();
        parent.winLength = 0;

        AlignColumn oc = winLast;                                            // original child
        AlignColumn op = parent.winLast;                            // original parent
        boolean saveChild = false;
        if (left == null || right == null || !left.selected || !right.selected) {
            saveWin();
            saveChild = true;
        }

        previ = parentLen;
        prevj = childLen;
        AlignColumn c = winLast;                                            // child
        AlignColumn p = parent.winLast;                                // parent

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);

            if (hmm2Parent[prevk] != 0) {
                op = op.prev;
                p = new AlignColumn(p, true);
                if (isLeft)
                    p.saveRight(op);
                else
                    p.saveLeft(op);
                p.saveParent(op);
                previ--;
                parent.winLength++;
            }
            if (hmm2Child[prevk] != 0) {
                if (saveChild) {
                    oc = oc.prev;
                    c = new AlignColumn(c, false);
                    c.saveBoth(oc);
                } else {
                    c = c.prev;
                }
                c.parent = p;
                c.orphan = hmm2Parent[prevk] == 0;
                if (hmm2Parent[prevk] != 0) {
                    if (isLeft)
                        p.left = c;
                    else
                        p.right = c;
                }
                prevj--;
            }
        }

        p.setWinFirst(parent.winFirst.prev);
        if (isLeft)
            parent.right.toggleUp(true);
        else
            parent.left.toggleUp(true);
        assert (parent.winLength == parent.old.winLength) : "HMM2 alignment error";

        if (saveChild) {
            c.setWinFirst(winFirst.prev);
            if (left != null && right != null) {
                left.toggleUp(true);
                right.toggleUp(true);
            }
        }

        for (c = winFirst.prev; c != null && c.parent == parent.old.winFirst; c = c.prev)
            c.parent = parent.winFirst;

        calcOrphan();
        parent.calcFelsen();
        parent.calcOrphan();
        parent.calcIndelLogLike();

        return retVal.value;
    }



    /** Selects a subtree */
    void selectSubtree(double[] weights, int level) {
        selected = true;
        // continue below with prescribed probability
        if (left != null && right != null) {
            if (Utils.generator.nextDouble() < weights[level]) {
                left.selectSubtree(weights, level + 1);
                right.selectSubtree(weights, level + 1);
            } else {
                left.selected = false;
                right.selected = false;
            }
        }
    }

    /** Marks `this' as the last selected Vertex of its subtree. */
    void lastSelected() {
        selected = true;
        if (left != null && right != null) {
            left.selected = false;
            right.selected = false;
        }
    }

    /** This function cuts out a window and realigns in the selected subtree. */
    double selectAndResampleAlignment() {

        // this code below checks pointer integrity...
        //    for(int i = 0; i < owner.vertex.length - 1; i++){
        //	owner.vertex[i].calcIndelLogLikeUp();//
        //}
        // select the beginning and end of the alignment

        // select the beginning and end of the window
        MuDouble p = new MuDouble(1.0);
        winLength = Utils.linearizerWeight(length, p);
        int b = (length - winLength == 0 ? 0 : Utils.generator.nextInt(length - winLength));
        AlignColumn actualAC = first;
        for (int i = 0; i < b; i++) {
            actualAC = actualAC.next;
        }
        winFirst = actualAC;
        for (int i = 0; i < winLength; i++) {
            actualAC = actualAC.next;
        }
        winLast = actualAC;

        double bpp = -Math.log(p.value * (length - winLength == 0 ? 1 : length - winLength));

        //	System.out.println("\nbpp after window select: "+bpp);

        // select window down
        if (left != null && right != null) {
            left.selectWindow();
            right.selectWindow();
        }

        selectWindowUp();

        // compute alignment backproposal
        bpp += doRecBackprop();
        if (parent != null)
            bpp += hmm2BackProp();

        //System.out.println("bpp after backproposal: "+bpp);

        // align the sequences
        double bppProp = doRecAlign();
        if (parent != null) {
            bppProp += hmm2AlignWithSave();
            parent.calcAllUp();
        } else {
            calcOrphan();
        }

        bpp += bppProp;
        bpp += Math.log(Utils.linearizerWeightProb(length, winLength) * (length - winLength));

        // 	System.out.print(" Prop: "+bppProp+" it's doublecheck: "+bppBack+" bpp: "+bpp+" ");

        return bpp;
    }


    void selectWindowUp() {
        if (parent != null) {
            if (winFirst.prev == null) {
                parent.winFirst = parent.first;
            } else {
                parent.winFirst = (winFirst.prev.orphan ? winFirst.prev.parent : winFirst.prev.parent.next);
            }
            parent.winLast = winLast.parent;

            parent.winLength = 0;
            for (AlignColumn actualAC = parent.winFirst; actualAC != parent.winLast; actualAC = actualAC.next)
                parent.winLength++;
        }
    }

    double doRecAlign() {
        if (left != null && right != null && left.selected && right.selected) {
            double ret = left.doRecAlign() + right.doRecAlign();
            ret += hmm3AlignWithSave();
            return ret;
        }
        return 0.0;
    }

    /** Computes recursively the log of the backproposal probability of the selected subtree. */
    double doRecBackprop() {
        if (left != null && right != null && left.selected && right.selected) {
            double ret = left.doRecBackprop() + right.doRecBackprop();
            ret += hmm3BackProp();
            return ret;
        }
        return 0.0;
    }

    /**
     * Restores all the changes an alignment resampling on the currently selected subtree has produced.
     * Must be called on selected subtree root.
     */
    void alignRestore() {
        doRecRestore();
        if (parent != null) {
            if (parent.left == this)
                parent.right.toggleUp(false);
            else
                parent.left.toggleUp(false);
            for (AlignColumn c = winFirst.prev; c != null && c.parent == parent.winFirst; c = c.prev)
                c.parent = parent.old.winFirst;
            parent.doRestore();
            if (parent.parent != null) {
                boolean isLeft = parent.parent.left == parent;
                for (AlignColumn p = parent.winFirst; p != parent.winLast; p = p.next) {
                    if (!p.orphan) {
                        if (isLeft)
                            p.parent.left = p;
                        else
                            p.parent.right = p;
                    }
                }
            }
            parent.calcAllUp();
        }
    }

    /**
     * Restores old sequences in subtree of `this' recursively and up-alignments in non-selected AlignmentNodes.
     * Assumes `this' is selected.
     */
    void doRecRestore() {
        if (left != null && right != null) {
            if (left.selected && right.selected) {
                left.doRecRestore();
                right.doRecRestore();
            } else {
                left.toggleUp(false);
                right.toggleUp(false);
            }
        }
        doRestore();
    }

    /**
     * Restores old sequence into `this', pointed to by `old.winFirst' & `old.winLast'.
     * Also restores orphan and indel likelihoods.
     */
    void doRestore() {
        winFirst = old.winFirst;
        winLength = old.winLength;
        length = old.length;
        winLast.prev = old.winLast;
        if (winFirst.prev != null)
            winFirst.prev.next = winFirst;
        else
            first = winFirst;
        orphanLogLike = old.orphanLogLike;
        indelLogLike = old.indelLogLike;
    }



    /** This function recursively selects the boundaries of the sequence that is to be resampled. */
    void selectWindow() {
        if (!selected) {
            return;
        }
        AlignColumn p = parent.first;
        AlignColumn c = first;
        while (p != parent.winFirst) {
            if (c.parent != p) {
                p = p.next;
            } else {
                c = c.next;
            }
        }
        winFirst = c;
        winLength = 0;
        //end of the window
        while (p != parent.winLast || (c.parent == p && c.orphan)) {
            //if(p == parent.winLast){
            //		found = true;
            //	}
            if (c.parent != p) {
                p = p.next;
            } else {
                c = c.next;
                winLength++;
            }
        }
        winLast = c;

        if (left != null && right != null) {
            left.selectWindow();
            right.selectWindow();
        }
    }


    /** This function iteratively aligns the windows for topology change */
    double alignAllWindows() {
        double ret = 0.0;

        length = 0;
        AlignColumn p = last;// pf = last;
        AlignColumn l = left.last, lf = left.last;
        AlignColumn r = right.last, rf = right.last;
        if (left != null) {

            checkPointers();
        }

        p = last; //pf = last;
        l = left.last;
        lf = left.last;
        r = right.last;
        rf = right.last;


        while (l != null) {
            //    System.out.println("****Start finding the left window: ");
            lf = lf.findWindowStart();
            //	    System.out.println("****Finished finding the left window: ");
            rf = rf.findWindowStart();
            winLength = 0;
            left.winFirst = lf;
            left.winLast = l;
            right.winFirst = rf;
            right.winLast = r;
            winLast = p;
            winFirst = p;
            if (p.emptyWindow || true) {

                ret += hmm3Align();
                p = winFirst;
            }

            //	    }
            l = lf.prev;
            lf = lf.prev;
            r = rf.prev;
            rf = rf.prev;
            if (l != null) {
                p = new AlignColumn(p, true);
                length++;
                p.left = l;
                p.right = r;
                l.parent = p;
                l.orphan = false;
                r.parent = p;
                r.orphan = false;
                p.selected = l.selected;
                p.emptyWindow = l.emptyWindow;
            }
            first = p;
        }
        calcFelsen();

        if (left != null) {

            checkPointers();
            //checking pointer integrity
            for (AlignColumn c = left.first; c != null; c = c.next) {
                p = first;
                while (c.parent != p && p != null) {
                    p = p.next;
                }
                if (p == null) {
                    throw new Error("children does not have a parent!!!" + this + " " + this.print());
                }
            }
            for (AlignColumn c = right.first; c != null; c = c.next) {
                p = first;
                while (c.parent != p && p != null) {
                    p = p.next;
                }
                if (p == null) {
                    throw new Error("children does not have a parent!!!" + this + " " + this.print());
                }
            }


        }


        //System.out.println("\n\n************ End of alignAllWindows ***************\n\n");

        return ret;
    }





    /**
     * Returns the alignment of this {@link Vertex} to its parent as an integer array:
     * for each character position <i>i</i> in this sequence it tells the character index
     * in the parent that <i>i</i> is aligned to or <i>-j-1</i> if it is not aligned
     * (gap: -) and the next character in the parent has index <i>j</i> (indices start from 0).
     * If this vertex does not have a parent an all-zero array is returned that has the same
     * length as the sequence in this vertex.
     */
    int[] getAlign() {
        int[] align = new int[length];
        AlignmentNode vp = parent;
        if (vp != null) {
            AlignColumn c = first, p = vp.first;
            int cn = 0, pn = 0;
            while (c != last || p != vp.last) {
                if (c.parent != p) {            // deletion (* -)
                    pn++;
                    p = p.next;
                } else if (c.orphan) {        // insertion (- *)
                    align[cn] = -pn - 1;
                    cn++;
                    c = c.next;
                } else {                    // substitution (* *)
                    align[cn] = pn;
                    pn++;
                    p = p.next;
                    cn++;
                    c = c.next;
                }
            }
        }
        return align;
    }


    String[] printedAlignment() {
        //	try{
        //	printPointers();
        if (parent == null) {
            String[] s = new String[2];
            AlignColumn a = first;
            s[0] = "";
            s[1] = "";
            while (a != last) {
                s[0] += "-";
                s[1] += a.mostLikely();
                a = a.next;
            }
            return s;
        }
        //	System.out.println("The parent !=null!!!");
        String[] s = new String[2];
        s[0] = "";
        s[1] = "";
        AlignColumn p = parent.first;
        AlignColumn c = first;
        while (p != parent.last || c != last) {
            // System.out.println(".");
            if (c.parent != p) {
                if (c.prev == null || (c.prev.parent != p || c.prev.orphan)) {
                    //	s[0] += "*";
                    s[0] += p.mostLikely();
                    s[1] += "-";
                }
                p = p.next;
                //	System.out.println("We increased p");
            } else {
                if (c.orphan) {
                    s[0] += "-";
                } else {
                    s[0] += p.mostLikely();
                }
                s[1] += c.mostLikely();
                c = c.next;
                //	System.out.println("We increased c");
            }
        }

        //System.out.println("Alignment here: "+print()+"\n"+s[0]+"\n"+s[1]);

        return s;
    }

    //////////////////////////////////////////////////
    String[] printedMultipleAlignment() {
        if (left == null && right == null) {
            String n = new String(name);
            String s = "";
            String s1 = printedAlignment()[1];
            for (int i = 0; i < s1.length(); i++) {
                if (s1.charAt(i) != '-') {
                    s += s1.charAt(i);
                }
            }
            return new String[]{n + "\t" + s};
        }
        String[] d1 = left.printedMultipleAlignment();
        String[] d2 = right.printedMultipleAlignment();

        String[] e1 = left.printedAlignment();
        String[] e2 = right.printedAlignment();


        //		System.out.println("Starting the calculation...");
        int tabPosition = d1[0].indexOf('\t');
        String[] s = new String[d1.length + d2.length + 1];
        s[0] = "";
        for (int i = 0; i < tabPosition; i++) {
            s[0] += " ";
        }
        s[0] += "\t";
        for (int i = 0; i < d1.length; i++) {
            s[i + 1] = (String) d1[i].subSequence(0, tabPosition + 1);
        }
        for (int i = d1.length + 1; i < s.length; i++) {
            s[i] = (String) d2[i - d1.length - 1].subSequence(0, tabPosition + 1);
        }
        int x = 0;
        int x1 = tabPosition + 1;
        int y = 0;
        int y1 = tabPosition + 1;
        //	System.out.println("before the while cycle...");
        while (x < e1[0].length() || y < e2[0].length()) {
            while (x < e1[0].length() && e1[0].charAt(x) == '-') {
                while (d1[0].charAt(x1) == '-') {
                    s[0] += "-";
                    for (int i = 0; i < d1.length; i++) {
                        s[i + 1] += d1[i].charAt(x1);
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += "-";
                    }
                    x1++;
                }
                s[0] += "-";
                for (int i = 0; i < d1.length; i++) {
                    s[i + 1] += d1[i].charAt(x1);
                }
                for (int i = d1.length + 1; i < s.length; i++) {
                    s[i] += "-";
                }
                x++;
                x1++;
            }
            while (y < e2[0].length() && e2[0].charAt(y) == '-') {
                while (d2[0].charAt(y1) == '-') {
                    s[0] += "-";
                    for (int i = 1; i < d1.length + 1; i++) {
                        s[i] += "-";
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += d2[i - d1.length - 1].charAt(y1);
                    }
                    y1++;
                }
                s[0] += "-";
                for (int i = 1; i < d1.length + 1; i++) {
                    s[i] += "-";
                }
                for (int i = d1.length + 1; i < s.length; i++) {
                    s[i] += d2[i - d1.length - 1].charAt(y1);
                }
                y++;
                y1++;
            }
            if (x < e1[0].length() && y < e2[0].length()) {
                while (x1 < d1[0].length() && d1[0].charAt(x1) == '-') {
                    s[0] += "-";
                    for (int i = 0; i < d1.length; i++) {
                        s[i + 1] += d1[i].charAt(x1);
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += "-";
                    }
                    x1++;
                }
                while (y1 < d2[0].length() && d2[0].charAt(y1) == '-') {
                    s[0] += "-";
                    for (int i = 1; i < d1.length + 1; i++) {
                        s[i] += "-";
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += d2[i - d1.length - 1].charAt(y1);
                    }
                    y1++;
                }
                //	s[0] += "*";
                s[0] += e1[0].charAt(x);
                for (int i = 0; i < d1.length; i++) {
                    s[i + 1] += (e1[1].charAt(x) != '-' ? d1[i].charAt(x1) : "-");
                }
                x1 += (e1[1].charAt(x) != '-' ? 1 : 0);
                x++;
                for (int i = d1.length + 1; i < s.length; i++) {
                    s[i] += (e2[1].charAt(y) != '-' ? d2[i - d1.length - 1].charAt(y1) : "-");
                }
                y1 += (e2[1].charAt(y) != '-' ? 1 : 0);
                y++;
            }
        }
        for (x1 = x1 + 0; x1 < d1[0].length(); x1++) {
            s[0] += "-";
            for (int i = 0; i < d1.length; i++) {
                s[i + 1] += d1[i].charAt(x1);
            }
            for (int i = d1.length + 1; i < s.length; i++) {
                s[i] += "-";
            }
        }
        for (y1 = y1 + 0; y1 < d2[0].length(); y1++) {
            s[0] += "-";
            for (int i = 0; i < d1.length; i++) {
                s[i + 1] += "-";
            }
            for (int i = d1.length + 1; i < s.length; i++) {
                s[i] += d2[i - d1.length - 1].charAt(y1);
            }
        }


        return s;
    }


    void printToScreenAlignment() {
        String[] s = printedAlignment();
        System.out.println(s[0] + "\n" + s[1]);
    }



}
