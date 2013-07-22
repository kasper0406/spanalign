package statalign.base;

public class SteinerTreeMCMCStrategy extends AbstractTreeMCMCStrategy<Tree, SteinerTreeUpdater> {
    // private Tree tree;
    private double[] weights; // for selecting internal tree node

    final static double LEAFCOUNT_POWER = 1.0;
    final static double SELTRLEVPROB[] = { 0.9, 0.6, 0.4, 0.2, 0 };

    public SteinerTreeMCMCStrategy(Tree tree) {
        super(tree, new SteinerTreeUpdater());

        this.tree = tree;
        weights = new double[tree.vertex.length];
    }

    @Override
    public boolean sampleEdge() {
        int i = Utils.generator.nextInt(tree.vertex.length - 1);
        return sampleEdge(tree.vertex[i]);
    }

    @Override
    public boolean sampleTopology() {
        boolean accepted = false;
        int vnum = tree.vertex.length;

        if (vnum <= 3)
            return false;

        // System.out.println("\n\n\t***\t***\t***\n\n\n");
        // System.out.print("Topology: ");
        // tree.printAllPointers();
        double oldLogLi = tree.getLogLike();

        int vertId, rnd = Utils.generator.nextInt(vnum - 3);
        vertId = tree.getTopVertexId(rnd);
        if (vertId != -1) {
            int lastId[] = new int[3], num = 0, newId = vertId;

            for (int i = vnum - 3; i < vnum; i++) {
                int id = tree.getTopVertexId(i);
                if (id == -1)
                    lastId[num++] = i;
                else if (id < vertId)
                    newId--;
            }
            rnd = lastId[newId];
        }
        Vertex nephew = tree.vertex[rnd];
        Vertex uncle = nephew.parent.brother();

        // for(vertId = 0; vertId < vnum; vertId++) {
        // if(tree.getTopVertexId(vertId) == -1) { // vertex eligible
        // if(rnd-- == 0)
        // break;
        // }
        // }
        // Vertex nephew = tree.vertex[vertId];

        // String[] s = tree.root.printedMultipleAlignment();
        // System.out.println("Alignment before topology changing: ");
        // for(int i = 0; i < s.length; i++){
        // System.out.println(s[i]);
        // }
        double bpp = nephew.fastSwapWithUncle();
        // double bpp = nephew.swapWithUncle();
        // s = tree.root.printedMultipleAlignment();
        // System.out.println("Alignment after topology changing: ");
        // for(int i = 0; i < s.length; i++){
        // System.out.println(s[i]);
        // }

        double newLogLi = tree.getLogLike();

        // tree.root.calcFelsRecursivelyWithCheck();
        // tree.root.calcIndelRecursivelyWithCheck();

        if (Math.log(Utils.generator.nextDouble()) < bpp
                + (newLogLi - oldLogLi) * tree.heat) {
            // accepted
            // System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
            accepted = true;
        } else {
            // rejected
            if(Utils.DEBUG) {
                // Checking pointer integrity before changing back topology
                for (int i = 0; i < tree.vertex.length; i++) {
                    if (tree.vertex[i].left != null && tree.vertex[i].right != null) {
                        tree.vertex[i].checkPointers();
                        AlignColumn p;
                        // checking pointer integrity
                        for (AlignColumn c = tree.vertex[i].left.first; c != null; c = c.next) {
                            p = tree.vertex[i].first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + tree.vertex[i] + " "
                                                + tree.vertex[i].print());
                        }
                        for (AlignColumn c = tree.vertex[i].right.first; c != null; c = c.next) {
                            p = tree.vertex[i].first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + tree.vertex[i] + " "
                                                + tree.vertex[i].print());
                        }

                    }
                }
            }

            uncle.fastSwapBackUncle();

            if(Utils.DEBUG) {
                // Checking pointer integrity after changing back topology
                for (int i = 0; i < tree.vertex.length; i++) {
                    if (tree.vertex[i].left != null && tree.vertex[i].right != null) {
                        tree.vertex[i].checkPointers();
                        AlignColumn p;
                        // checking pointer integrity
                        for (AlignColumn c = tree.vertex[i].left.first; c != null; c = c.next) {
                            p = tree.vertex[i].first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + tree.vertex[i] + " "
                                                + tree.vertex[i].print());
                        }
                        for (AlignColumn c = tree.vertex[i].right.first; c != null; c = c.next) {
                            p = tree.vertex[i].first;
                            while (c.parent != p && p != null)
                                p = p.next;
                            if (p == null)
                                throw new Error(
                                        "children does not have a parent!!!"
                                                + tree.vertex[i] + " "
                                                + tree.vertex[i].print());
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

    @Override
    public boolean sampleAlignment() {
        for (int i = 0; i < tree.vertex.length; i++) {
            tree.vertex[i].selected = false;
        }
        // System.out.print("Alignment: ");
        double oldLogLi = tree.getLogLike();
        // System.out.println("fast indel before: "+tree.root.indelLogLike);
        tree.countLeaves(); // calculates recursively how many leaves we have
        // below this node
        for (int i = 0; i < weights.length; i++) {
            weights[i] = Math.pow(tree.vertex[i].leafCount, LEAFCOUNT_POWER);
        }
        int k = Utils.weightedChoose(weights, null);
        // System.out.println("Sampling from the subtree: "+tree.vertex[k].print());
        tree.vertex[k].selectSubtree(SELTRLEVPROB, 0);
        double bpp = tree.vertex[k].selectAndResampleAlignment();
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
            tree.vertex[k].alignRestore();
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
        if (tree.substitutionModel.params.length == 0)
            return false;
        else {
            double mh = tree.substitutionModel.sampleParameter();
            double oldlikelihood = tree.root.orphanLogLike;
            for (int i = 0; i < tree.vertex.length; i++) {
                tree.vertex[i].updateTransitionMatrix();
            }
            tree.root.calcFelsRecursively();
            double newlikelihood = tree.root.orphanLogLike;
            if (Utils.generator.nextDouble() < Math.exp(mh
                    + (Math.log(tree.substitutionModel.getPrior())
                    + newlikelihood - oldlikelihood))
                    * tree.heat) {
                // System.out.println("Substitution parameter: accepted (old: "+oldlikelihood+" new: "+newlikelihood+")");
                return true;
            } else {
                tree.substitutionModel.restoreParameter();
                for (int i = 0; i < tree.vertex.length; i++) {
                    tree.vertex[i].updateTransitionMatrix();
                }
                tree.root.calcFelsRecursively();
                // System.out.println("Substitution parameter: rejected (old: "+oldlikelihood+" new: "+newlikelihood+")");

                return false;
            }
        }
    }
}
