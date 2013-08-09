package statalign.base;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

class Transplanter extends MCMCMove<Spannoid, TransplantResult> {
    public Transplanter(Spannoid tree) {
        super(tree);
    }

    private boolean canSample() {
        return tree.getInnerLabelledVertices().size() > 0;
    }

    public boolean sample() {
        if (!canSample())
            return false;
        return super.sample();
    }

    @Override
    protected TransplantResult jump() {
        double bpp = 0;

        // choose source vertex
        List<Vertex> sources = tree.getInnerLabelledVertices();
        Vertex source = sources.get(Utils.generator.nextInt(sources.size()));
        bpp -= -Math.log(sources.size());

        // choose destination vertex
        List<Vertex> destinations = getDestinations(source);
        Vertex destination = destinations.get(Utils.generator.nextInt(destinations.size()));
        bpp -= -Math.log(destinations.size());

        // where to move the source vertex back to in case of rejection
        Vertex prev = getArbitraryNeighbour(source);

        bpp += moveSubtree(source, destination);

        bpp += -Math.log(tree.getInnerLabelledVertices().size());
        bpp += -Math.log(getDestinations(destination).size());

        return new TransplantResult(bpp, source, destination, prev);
    }

    @Override
    protected void restore(TransplantResult result) {
        restoreSubtree(result.source, result.prev);
    }

    private Vertex getArbitraryNeighbour(Vertex vertex) {
        Set<Vertex> neighbourhood = tree.getNeighbourhood(vertex);
        for (Vertex neighbour : neighbourhood) {
            if (neighbour != vertex) {
                return neighbour;
            }
        }

        // should never reach this
        return null;
    }

    private List<Vertex> getDestinations(Vertex source) {
        List<Vertex> destinations = new ArrayList<Vertex>();

        Set<Vertex> neighbourhood = tree.getNeighbourhood(source);

        for (Vertex neighbour : neighbourhood) {
            if (neighbour == source) {
                continue;
            }

            for (Vertex v : neighbour.owner.vertex) {
                if (v != neighbour && v.name != null) {
                    destinations.add(v);
                }
            }
        }

        return destinations;
    }

    private double moveSubtree(Vertex source, Vertex dest){
        // update internal Spannoid structure
        tree.moveComponent(source, dest);

        double bpp = 0;

        source.fullWin();
        source.parent.fullWin();
        bpp += source.hmm2BackProp();

        // update vertex information
        source.seq = dest.seq;
        source.length = dest.length;
        source.name = dest.name;

        // save old alignment
        source.old.first = source.first;
        source.old.last = source.last;

        // create new alignment columns for new sequence
        source.first = new AlignColumn(source);
        source.first.seq = dest.first.seq.clone();
        source.first.parent = source.parent.last;
        AlignColumn prev = source.first;
        for (AlignColumn cur = dest.first.next; cur != dest.last; cur = cur.next) {
            AlignColumn actual = new AlignColumn(source);
            actual.seq = cur.seq.clone();
            actual.parent = source.parent.last;

            actual.prev = prev;
            prev.next = actual;

            prev = actual;
        }
        AlignColumn last = new AlignColumn(source);
        last.parent = source.parent.last;
        last.orphan = false;
        last.prev = prev;
        prev.next = last;
        source.last = last;

        // destroy parent alignment pointers
        for (AlignColumn c = source.parent.first; c != source.parent.last; c = c.next) {
            if (source.parent.left == source) {
                c.left = null;
            } else {
                c.right = null;
            }
        }

        if (source.parent.left == source) {
            source.parent.last.left = source.last;
        } else {
            source.parent.last.right = source.last;
        }

        source.fullWin();
        source.parent.fullWin();
        bpp += source.hmm2AlignWithRecalc();

        source.calcAllUp();

        return bpp;
    }

    private void restoreSubtree(Vertex source, Vertex prev){
        // update internal Spannoid structure
        tree.moveComponent(source, prev);

        // restore old vertex information
        source.seq = prev.seq;
        source.length = prev.length;
        source.name = prev.name;

        // restore old alignment
        source.first = source.old.first;
        source.last = source.old.last;

        boolean isLeft = (source.parent.left == source);

        AlignColumn c = source.first;
        AlignColumn p = source.parent.first;
        while (p != null) {
            if (c.parent != p) {
                if (isLeft) {
                    p.left = null;
                } else {
                    p.right = null;
                }
                p = p.next;
            } else if (c.orphan) {
                c = c.next;
            } else {
                if (isLeft) {
                    p.left = c;
                } else {
                    p.right = c;
                }
                p = p.next;
                c = c.next;
            }
        }

        source.calcAllUp();
    }
}

class TransplantResult extends MCMCResult {
    Vertex source, destination, prev;

    public TransplantResult(double bpp, Vertex source, Vertex destination, Vertex prev) {
        super(bpp);
        this.source = source;
        this.destination = destination;
        this.prev = prev;
    }
}