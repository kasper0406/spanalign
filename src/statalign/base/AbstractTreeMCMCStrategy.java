package statalign.base;

public abstract class AbstractTreeMCMCStrategy<T extends ITree, Updater extends ITreeUpdater<T>> implements MCMCStrategy {
    protected T tree;
    protected Updater updater;

    final static double LEAFCOUNT_POWER = 1.0;
    final static double SELTRLEVPROB[] = { 0.9, 0.6, 0.4, 0.2, 0 };

    public AbstractTreeMCMCStrategy(T tree, Updater updater) {
        this.tree = tree;
        this.updater = updater;
    }

    public ITree getTree() {
        return tree;
    }

    public boolean sampleIndelParameter() {
        boolean accepted = false;
        switch (Utils.generator.nextInt(3)) {
            case 0: {
                // System.out.print("Indel param R: ");
                double oldR = tree.getR();
                double oldLogLikelihood = tree.getLogLike();
                double newR;
                while ((newR = oldR + Utils.generator.nextDouble()
                        * Utils.R_SPAN - Utils.R_SPAN / 2.0) <= 0.0
                        || newR >= 1.0)
                    ;
                updater.updateR(tree, newR);
                double newLogLikelihood = tree.getLogLike();
                if (Utils.generator.nextDouble() < Math
                        .exp((newLogLikelihood - oldLogLikelihood) * tree.getHeat())
                        * (Math.min(1.0 - oldR, Utils.R_SPAN / 2.0) + Math.min(
                        oldR, Utils.R_SPAN / 2.0))
                        / (Math.min(1.0 - tree.getR(), Utils.R_SPAN / 2.0) + Math
                        .min(tree.getR(), Utils.R_SPAN / 2.0))) {
                    // accept, do nothing
                    // System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
                    accepted = true;
                } else {
                    // restore
                    updater.updateR(tree, oldR);
                    // System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
                }

                break;
            }
            case 1: {
                // ///////////////////////////////////////////////
                // System.out.print("Indel param Lambda: ");
                double oldLambda = tree.getLambda();
                double oldLogLikelihood = tree.getLogLike();
                double newLambda;
                while ((newLambda = oldLambda
                        + Utils.generator.nextDouble() * Utils.LAMBDA_SPAN
                        - Utils.LAMBDA_SPAN / 2.0) <= 0.0
                        || newLambda >= tree.getMu())
                    ;
                updater.updateLambda(tree, newLambda);
                double newLogLikelihood = tree.getLogLike();
                if (Utils.generator.nextDouble() < Math.exp((newLogLikelihood
                        - oldLogLikelihood - tree.getLambda() + oldLambda)
                        * tree.getHeat())
                        * (Math.min(Utils.LAMBDA_SPAN / 2.0, tree.getMu()
                        - oldLambda) + Math.min(oldLambda,
                        Utils.LAMBDA_SPAN / 2.0))
                        / (Math.min(Utils.LAMBDA_SPAN / 2.0, tree.getMu()
                        - tree.getLambda()) + Math.min(
                        tree.getLambda(), Utils.LAMBDA_SPAN / 2.0))) {
                    // accept, do nothing
                    // System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+" oldLambda: "+oldLambda+" newLambda: "+tree.hmm2.params[1]+")");
                    accepted = true;
                } else {
                    // restore
                    updater.updateLambda(tree, oldLambda);
                    // System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+" oldLambda: "+oldLambda+" newLambda: "+tree.hmm2.params[1]+")");
                }
                break;
            }
            case 2: {
                // ///////////////////////////////////////////////////////
                // System.out.print("Indel param Mu: ");
                double oldMu = tree.getMu();
                double oldLogLikelihood = tree.getLogLike();
                double newMu;
                while ((newMu = oldMu + Utils.generator.nextDouble()
                        * Utils.MU_SPAN - Utils.MU_SPAN / 2.0) <= tree.getLambda())
                    ;
                updater.updateMu(tree, newMu);
                double newLogLikelihood = tree.getLogLike();
                if (Utils.generator.nextDouble() < Math.exp((newLogLikelihood
                        - oldLogLikelihood - tree.getMu() + oldMu)
                        * tree.getHeat())
                        * (Utils.MU_SPAN / 2.0 + Math.min(oldMu
                        - tree.getLambda(), Utils.MU_SPAN / 2.0))
                        / (Utils.MU_SPAN / 2.0 + Math.min(tree.getMu()
                        - tree.getLambda(), Utils.MU_SPAN / 2.0))) {
                    // accept, do nothing
                    // System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
                    accepted = true;
                } else {
                    // restore
                    updater.updateMu(tree, oldMu);
                    // System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
                }
                break;
            }
        }
        return accepted;
    }

    protected boolean sampleEdge(Vertex vertex) {
        // System.out.print("Edge: ");
        double oldEdge = vertex.edgeLength;
        double oldLogLikelihood = tree.getLogLike();
        while ((vertex.edgeLength = oldEdge
                + Utils.generator.nextDouble() * Utils.EDGE_SPAN
                - (Utils.EDGE_SPAN / 2.0)) < 0.01)
            ;
        vertex.edgeChangeUpdate();
        // Vertex actual = tree.vertex[i];
        // while(actual != null){
        // actual.calcFelsen();
        // actual.calcOrphan();
        // actual.calcIndelLogLike();
        // actual = actual.parent;
        // }
        vertex.calcAllUp();
        double newLogLikelihood = tree.getLogLike();
        if (Utils.generator.nextDouble() < (Math.exp((newLogLikelihood
                - oldLogLikelihood - vertex.edgeLength + oldEdge)
                * tree.getHeat()) * (Math.min(oldEdge - 0.01, Utils.EDGE_SPAN / 2.0) + Utils.EDGE_SPAN / 2.0))
                / (Math.min(vertex.edgeLength - 0.01,
                Utils.EDGE_SPAN / 2.0) + Utils.EDGE_SPAN / 2.0)) {
            // acceptance, do nothing
            // System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");

            return true;
        } else {
            // reject, restore
            // System.out.print("Rejected! i: "+i+"\tOld likelihood: "+oldLogLikelihood+"\tNew likelihood: "+newLogLikelihood);
            vertex.edgeLength = oldEdge;
            vertex.edgeChangeUpdate();
            // actual = tree.vertex[i];
            // while(actual != null){
            // actual.calcFelsen();
            // actual.calcOrphan();
            // actual.calcIndelLogLike();
            // actual = actual.parent;
            // }
            vertex.calcAllUp();
            // System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");

            return false;
        }
    }

    /**
     * Samples now alignment for a given tree.
     */
    protected boolean sampleAlignment(Tree tree)
    {
        for (Vertex v : tree.vertex)
            v.selected = false;
        // System.out.print("Alignment: ");
        double oldLogLi = tree.getLogLike();
        // System.out.println("fast indel before: "+tree.root.indelLogLike);
        tree.countLeaves(); // calculates recursively how many leaves we have
        // below this node

        double[] weights = new double[tree.vertex.size()];
        for (int i = 0; i < weights.length; i++) {
            weights[i] = Math.pow(tree.vertex.get(i).leafCount, LEAFCOUNT_POWER);
        }
        int k = Utils.weightedChoose(weights, null);
        // System.out.println("Sampling from the subtree: "+tree.vertex[k].print());

        for (Vertex v : tree.vertex)
            v.checkPointers();

        tree.vertex.get(k).selectSubtree(SELTRLEVPROB, 0);
        double bpp = tree.vertex.get(k).selectAndResampleAlignment();
        double newLogLi = tree.getLogLike();

        // String[] printedAlignment = tree.printedAlignment("StatAlign");
        // for(String i: printedAlignment)
        // System.out.println(i);
        //
        // System.out.println("-----------------------------------------------------------------------------");
        // double fastFels = tree.root.orphanLogLike;
        // double fastIns = tree.root.indelLogLike;
        // report();
        // tree.root.first.seq[0] = 0.0;
        // System.out.println("Old before: "+tree.root.old.indelLogLike);
        // tree.root.calcFelsRecursivelyWithCheck();
        // tree.root.calcIndelRecursivelyWithCheck();
        // tree.root.calcIndelLikeRecursively();
        // System.out.println("Old after: "+tree.root.old.indelLogLike);
        // System.out.println("Check logli: "+tree.getLogLike()+" fastFels: "+fastFels+" slowFels: "+tree.root.orphanLogLike+
        // " fastIns: "+fastIns+" slowIns: "+tree.root.indelLogLike);
        // System.out.println("selected subtree: "+tree.vertex[k].print());
        // System.out.println("bpp: "+bpp+"old: "+oldLogLi+"new: "+newLogLi +
        // "heated diff: " + ((newLogLi - oldLogLi) * tree.heat));
        if (Math.log(Utils.generator.nextDouble()) < bpp
                + (newLogLi - oldLogLi) * tree.heat) {
            // accepted
            // System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
            return true;
        } else {
            // refused
            // String[] s = tree.printedAlignment();
            tree.vertex.get(k).alignRestore();
            // s = tree.printedAlignment();
            // System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")");
            // System.out.println("after reject fast: "+tree.root.indelLogLike);
            // tree.root.calcIndelRecursivelyWithCheck();
            // System.out.println(" slow: "+tree.root.indelLogLike);

            return false;
        }
        // tree.root.calcFelsRecursivelyWithCheck();
        // tree.root.calcIndelRecursivelyWithCheck();
    }

    @Override
    public boolean sampleSubstParameter() {
        if (tree.getSubstitutionModel().params.length == 0)
            return false;
        else {
            double oldlikelihood = tree.getOrphanLogLike();
            // Note that we should only sample from one substitution parameter,
            // since the substitution model is shared between all trees!
            double mh = tree.getSubstitutionModel().sampleParameter();
            updater.recalcSubstitutionParameters(tree);
            double newlikelihood = tree.getOrphanLogLike();
            if (Utils.generator.nextDouble() < Math.exp(mh
                    + (Math.log(tree.getSubstitutionModel().getPrior())
                    + newlikelihood - oldlikelihood))
                    * tree.getHeat()) {
                // System.out.println("Substitution parameter: accepted (old: "+oldlikelihood+" new: "+newlikelihood+")");
                return true;
            } else {
                tree.getSubstitutionModel().restoreParameter();
                updater.recalcSubstitutionParameters(tree);
                // System.out.println("Substitution parameter: rejected (old: "+oldlikelihood+" new: "+newlikelihood+")");

                return false;
            }
        }
    }


    public boolean sampleTopology(Tree tree) {
        boolean accepted = false;
        int vnum = tree.vertex.size();

        if (vnum <= 3)
            return false;

        // System.out.println("\n\n\t***\t***\t***\n\n\n");
        // System.out.print("Topology: ");
        // tree.printAllPointers();
        double oldLogLi = tree.getLogLike();
        SteinerTreeUpdater.NNIResult nni = updater.performNNI(tree);

        double newLogLi = tree.getLogLike();

        // tree.root.calcFelsRecursivelyWithCheck();
        // tree.root.calcIndelRecursivelyWithCheck();

        if (Math.log(Utils.generator.nextDouble()) < nni.bpp
                + (newLogLi - oldLogLi) * tree.heat) {
            // accepted
            // System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
            accepted = true;
        } else {
            // rejected
            if(Utils.DEBUG) {
                // Checking pointer integrity before changing back topology
                for (Vertex v : tree.vertex) {
                    if (v.left != null && v.right != null) {
                        v.checkPointers();
                        AlignColumn p;
                        // checking pointer integrity
                        for (AlignColumn c = v.left.first; c != null; c = c.next) {
                            p = v.first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + v + " "
                                                + v.print());
                        }
                        for (AlignColumn c = v.right.first; c != null; c = c.next) {
                            p = v.first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + v + " "
                                                + v.print());
                        }

                    }
                }
            }

            updater.revertNNI(tree, nni);

            if(Utils.DEBUG) {
                // Checking pointer integrity after changing back topology
                for (Vertex v : tree.vertex) {
                    if (v.left != null && v.right != null) {
                        v.checkPointers();
                        AlignColumn p;
                        // checking pointer integrity
                        for (AlignColumn c = v.left.first; c != null; c = c.next) {
                            p = v.first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + v + " "
                                                + v.print());
                        }
                        for (AlignColumn c = v.right.first; c != null; c = c.next) {
                            p = v.first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + v + " "
                                                + v.print());
                        }
                    }
                }
            }
            // uncle.swapBackUncle();
            // s = tree.root.printedMultipleAlignment();
            // System.out.println("Alignment after changing back the topology: ");
            // for(int i = 0; i < s.length; i++){
            // System.out.println(s[i]);
            // }
            // System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")");
        }

        // tree.printAllPointers();
        // System.out.println("\n\n\t***\t***\t***\n\n\n");
        if(Utils.DEBUG) {
            tree.root.calcFelsRecursivelyWithCheck();
            tree.root.calcIndelRecursivelyWithCheck();
        }

        return accepted;
    }


}
