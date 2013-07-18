package statalign.base;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/17/13
 * Time: 3:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class TopologyNode extends AlignmentNode{





    /**
     * Swaps `this' with its uncle, hmm3Align'ing parent and grandparent, the latter is hmm2Align'ed as well.
     * Assumes `this' has a non-null grandparent.
     * @return log-quotient of backproposal and proposal
     */
    double swapWithUncle1() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;
        double ret = 0.0;

        fullWin();
        brother().fullWin();
        uncle.fullWin();
        parent.fullWin();
        grandpa.fullWin();
        if (grandpa.parent != null)
            grandpa.parent.fullWin();

        lastSelected();
        brother().lastSelected();
        uncle.lastSelected();
        parent.selected = true;
        grandpa.selected = true;

        ret += parent.hmm3BackProp();
        ret += grandpa.hmm3BackProp();
        if (grandpa.parent != null)
            ret += grandpa.hmm2BackProp();

        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        ret += uncle.parent.hmm3AlignWithSave();
        ret += grandpa.hmm3AlignWithSave();
        if (grandpa.parent != null) {
            ret += grandpa.hmm2AlignWithSave();
            grandpa.parent.calcAllUp();
        } else {
            grandpa.calcOrphan();
        }

        return ret;
    }

    /**
     * Restores the exact state just before the call of `swapWithUncle'. Must be called on ex-uncle.
     * Assumes `this' has a non-null grandparent.
     */
    void swapBackUncle1() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;

        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        grandpa.alignRestore();
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Swaps `this' with its uncle, hmm3Align'ing parent and grandparent, the latter is hmm2Align'ed as well.
     * Assumes `this' has a non-null grandparent.
     * @return log-quotient of backproposal and proposal
     */
    double fastSwapWithUncle() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;
        double ret = 0.0;
        Vertex starter;

        //System.out.println("fast swap here"+grandpa.print());

        lastSelected();
        brother().lastSelected();
        uncle.lastSelected();
        parent.selected = true;
        grandpa.selected = true;

        setAllAlignColumnsUnselected();
        brother().setAllAlignColumnsUnselected();
        uncle.setAllAlignColumnsUnselected();
        parent.setAllAlignColumnsUnselected();
        grandpa.setAllAlignColumnsUnselected();

        if (grandpa.parent != null) {
            grandpa.parent.selected = true;
            grandpa.brother().selected = false;
            grandpa.parent.setAllAlignColumnsUnselected();
            starter = grandpa.parent;
            //  System.out.println("has grandgrandpa");
        } else {
            starter = grandpa;
            //System.out.println("does not have grandgrandpa");
        }

        ret += starter.selectAnchors();
        //System.out.println("RET after selecting the anchors: "+ret);
        //System.out.println("Backproposing this would be:    "+starter.backproposeAnchors());

        //backproposals
        //for parent
        AlignColumn p = parent.last, pf = parent.last;
        AlignColumn t = last, tf = last;
        AlignColumn b = brother().last, bf = brother().last;
        while (p != null) {
            // System.out.println(".");
            pf = pf.findWindowStart();
            tf = tf.findWindowStart();
            bf = bf.findWindowStart();
            if (p.prev != pf || p.emptyWindow || true) {
                parent.winFirst = pf;
                parent.winLast = p;
                winFirst = tf;
                winLast = t;
                brother().winFirst = bf;
                brother().winLast = b;
                ret += parent.hmm3BackProp();
            }
            p = pf.prev;
            pf = pf.prev;
            t = tf.prev;
            tf = tf.prev;
            b = bf.prev;
            bf = bf.prev;
        }
        //for grandpa
        AlignColumn g = grandpa.last, gf = grandpa.last;
        p = parent.last;
        pf = parent.last;
        AlignColumn u = uncle.last, uf = uncle.last;
        while (g != null) {
//			  System.out.println(".");
            gf = gf.findWindowStart();
            pf = pf.findWindowStart();
            uf = uf.findWindowStart();
            if (g.prev != gf || g.emptyWindow || true) {
                grandpa.winFirst = gf;
                grandpa.winLast = g;
                parent.winFirst = pf;
                parent.winLast = p;
                uncle.winFirst = uf;
                uncle.winLast = u;
                ret += grandpa.hmm3BackProp();
            }
            g = gf.prev;
            gf = gf.prev;
            p = pf.prev;
            pf = pf.prev;
            u = uf.prev;
            uf = uf.prev;

        }
        // if there is a grand-grandpa...
        if (grandpa.parent != null) {
            g = grandpa.last;
            gf = grandpa.last;
            AlignColumn gg = grandpa.parent.last, ggf = grandpa.parent.last;
            while (gg != null) {
                ggf = ggf.findWindowStart();
                gf = gf.findWindowStart();
                if (gg.prev != ggf || gg.emptyWindow || true) {
                    grandpa.parent.winFirst = ggf;
                    grandpa.parent.winLast = gg;
                    grandpa.winFirst = gf;
                    grandpa.winLast = g;
                    ret += grandpa.hmm2BackProp();
                }
                g = gf.prev;
                gf = gf.prev;
                gg = ggf.prev;
                ggf = ggf.prev;
            }
        }


        //saving old alignments
        copySequence();
        brother().copySequence();
        uncle.copySequence();
        parent.fullWin();
        parent.saveWin();
        grandpa.fullWin();
        grandpa.saveWin();
        if (grandpa.parent != null) {
            grandpa.copyGreatgrandpaSequence();
        }

        //NNI
        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        //new parent sequence
        ret += uncle.parent.alignAllWindows();
        uncle.parent.checkPointers();


        //new grandpa
        ret += grandpa.alignAllWindows();

        grandpa.checkPointers();


        //if there is a greatgrandpa, then we align it...
        if (grandpa.parent != null) {
            //System.out.println("Greatgrandpa in fastSwap: "+grandpa.parent);
            AlignColumn gg = grandpa.parent.last, ggf = gg;
            g = grandpa.last;
            gf = g;
            while (g != null) {
                gf = gf.findWindowStart();
                ggf = ggf.findWindowStart();
                grandpa.winLast = g;
                grandpa.winFirst = gf;
                grandpa.parent.winLast = gg;
                grandpa.parent.winFirst = ggf;
                if (gg.emptyWindow || true) {
                    ret += grandpa.hmm2Align();
                }
                g = gf.prev;
                gf = g;
                grandpa.parent.first = ggf;
                gg = ggf.prev;
                ggf = gg;
                if (g != null) {
                    g.parent = gg;
                    g.orphan = false;
                    if (grandpa.parent.left == grandpa) {
                        gg.left = g;
                    } else {
                        gg.right = g;
                    }
                }
            }

            //	    System.out.println("RET after aligning everything: "+ret);
            //System.out.println("Calculating Felsenstein, called from fastSwap");
            grandpa.parent.calcFelsen();
            grandpa.parent.checkPointers();
            //System.out.println("Aligned greatgrandpa, pointers:");


        }


        calcOrphan();
        uncle.brother().calcOrphan();
        uncle.calcOrphan();
        uncle.calcAllUp();


        //And finally: the window selecting backproposals...
        //	System.out.println((grandpa.parent == starter ? "grandpa.parent == starter" : "grandpa.parent != starter"));
        ret += starter.backproposeAnchors();

        //	System.out.println("RET after backproposing anchors (final value) "+ret);

        return ret;
    }
    /**
     * Restores the exact state just before the call of `swapWithUncle'. Must be called on ex-uncle.
     * Assumes `this' has a non-null grandparent.
     */
    void fastSwapBackUncle() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;

        fullWin();
        brother().fullWin();
        parent.fullWin();
        uncle.fullWin();
        grandpa.fullWin();
        if (grandpa.parent != null) {
            grandpa.parent.fullWin();
        }


        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;


        //	System.out.println("Starting alignment restauration...");
        grandpa.alignRestore();
        //System.out.println("End alignment restauration...");
    }

    /**
     * This function selects anchor homologous columns that are not changed during topology change
     * returns with the -logarithm of proposal probability!!!
     */
    double selectAnchors() {
        double ret = 0.0;

        AlignColumn actual = first;
        while (actual != null) {
            if (actual.isHomologous()) {
                if (Utils.generator.nextDouble() < SELECTING) {
                    actual.selectDown();
                    ret -= Math.log(SELECTING);
                    // System.out.println("Homologous, selected!");
                    if (actual.isEmptyWindow()) {
                        if (Utils.generator.nextDouble() < EMPTY_WINDOW) {
                            actual.setEmptyWindowDown(true);
                            ret -= Math.log(EMPTY_WINDOW);
                        } else {
                            actual.setEmptyWindowDown(false);
                            ret -= Math.log(1.0 - EMPTY_WINDOW);
                        }
                    } else {// we need this since later on we will not know if it is not an empty window...
                        actual.setEmptyWindowDown(true);
                    }
                } else {
                    actual.selected = false;
                    ret -= Math.log(1.0 - SELECTING);
                    // System.out.println("Homologous, not selected!");
                }
            } else {
                actual.selected = false;
                //	System.out.println("Not homologous!");
            }
            actual = actual.next;
            //System.out.println("...");
        }

        return ret;
    }

    /**
     * This function selects anchor homologous columns that are not changed during topology change
     * returns with the logarithm of proposal probability!!!
     */
    double backproposeAnchors() {
        double ret = 0.0;

        AlignColumn actual = first;
        //	while(actual.next != null){
        while (actual != null) {
            if (actual.isHomologous()) {
                if (actual.selected) {
                    ret += Math.log(SELECTING);
                    //    System.out.println("Homologous, selected!");
                    if (actual.isEmptyWindow()) {
                        if (actual.emptyWindow) {
                            ret += Math.log(EMPTY_WINDOW);
                        } else {
                            ret += Math.log(1.0 - EMPTY_WINDOW);
                        }
                    }
                } else {
                    ret += Math.log(1.0 - SELECTING);
                    //    System.out.println("Homologous, not selected!");
                }
            }
            //    else{
            //	System.out.println("It is not homologous!!!");
            //    }
            actual = actual.next;
            //System.out.println("...");
        }

        return ret;
    }

    /** this function set selected to false for all AlignColumns */
    void setAllAlignColumnsUnselected() {
        AlignColumn a = first;
        while (a != null) {
            a.selected = false;
            a = a.next;
        }
    }

    /** This function copies the sequence to the old vertex */
    void copySequence() {

        fullWin();
        saveWin();
        AlignColumn a = last;
        AlignColumn prev = last;
        while (a.prev != null) {
            a = a.prev;
            prev = new AlignColumn(prev, false);
            prev.saveBoth(a);
        }
        winFirst = prev;
        first = prev;

        if (left != null) {
            //winFirst = first;
            //winLast = last;
            //old.winFirst = old.first;
            //old.winLast = old.last;
            left.toggleUp(true);
            //System.out.println("Left is ready");
            right.toggleUp(true);
        }
    }

    /**
     * This function is used for special copying the geatgrandpa's sequence in topology change.
     * The function is called to the grandpa, whose parent is the greatgrandpa
     */
    void copyGreatgrandpaSequence() {
        parent.fullWin();
        parent.saveWin();
        AlignColumn a = parent.last;
        AlignColumn prev = parent.last;
        boolean isLeft = (parent.left == this);
        while (a.prev != null) {
            a = a.prev;
            prev = new AlignColumn(prev, true);
            if (isLeft) {
                prev.saveRight(a);
            } else {
                prev.saveLeft(a);
            }
            prev.saveParent(a);
        }
        parent.winFirst = prev;
        parent.first = prev;
        if (parent.left == this) {
            parent.right.toggleUp(true);
            //	    System.out.println("Toggling up right!!!");
        } else {
            parent.left.toggleUp(true);
            //System.out.println("Toggling up left!!!");
        }


    }
    // This function was modified so it casts the TopologyNodes into Vertex instances
    void parentNewChild(TopologyNode child) {
        child.last.parent = parent.last;
        if (parent.left == this) {
            parent.left = (Vertex) child;
            parent.last.left = child.last;
        } else {
            parent.right = (Vertex) child;
            parent.last.right = child.last;
        }
    }

}


