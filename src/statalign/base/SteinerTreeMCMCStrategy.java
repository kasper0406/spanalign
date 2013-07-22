package statalign.base;

public class SteinerTreeMCMCStrategy extends AbstractTreeMCMCStrategy<Tree, SteinerTreeUpdater> {
    public SteinerTreeMCMCStrategy(Tree tree) {
        super(tree, new SteinerTreeUpdater());

        this.tree = tree;
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
        return sampleAlignment(tree);
    }
}
