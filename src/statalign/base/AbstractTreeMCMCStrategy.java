package statalign.base;

public abstract class AbstractTreeMCMCStrategy<T extends ITree, Updater extends ITreeUpdater<T>> implements MCMCStrategy {
    protected T tree;
    protected Updater updater;

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
            case 0:
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
            case 1:
                // ///////////////////////////////////////////////
                // System.out.print("Indel param Lambda: ");
                double oldLambda = tree.getLambda();
                oldLogLikelihood = tree.getLogLike();
                double newLambda;
                while ((newLambda = oldLambda
                        + Utils.generator.nextDouble() * Utils.LAMBDA_SPAN
                        - Utils.LAMBDA_SPAN / 2.0) <= 0.0
                        || newLambda >= tree.getMu())
                    ;
                updater.updateLambda(tree, newLambda);
                newLogLikelihood = tree.getLambda();
                if (Utils.generator.nextDouble() < Math.exp((newLogLikelihood
                        - oldLogLikelihood - tree.getLambda() + oldLambda)
                        * tree.getHeat())
                        * (Math.min(Utils.LAMBDA_SPAN / 2.0, tree.getLambda()
                        - oldLambda) + Math.min(oldLambda,
                        Utils.LAMBDA_SPAN / 2.0))
                        / (Math.min(Utils.LAMBDA_SPAN / 2.0, tree.getLambda()
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
            case 2:
                // ///////////////////////////////////////////////////////
                // System.out.print("Indel param Mu: ");
                double oldMu = tree.getMu();
                oldLogLikelihood = tree.getLogLike();
                double newMu;
                while ((newMu = oldMu + Utils.generator.nextDouble()
                        * Utils.MU_SPAN - Utils.MU_SPAN / 2.0) <= tree.getLambda())
                    ;
                updater.updateMu(tree, newMu);
                newLogLikelihood = tree.getLogLike();
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
}
