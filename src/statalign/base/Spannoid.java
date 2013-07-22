package statalign.base;

import statalign.base.hmm.HmmNonParam;
import statalign.base.hmm.HmmTkf92;
import statalign.base.thread.StoppedException;
import statalign.model.score.SubstitutionScore;
import statalign.model.subst.SubstitutionModel;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/19/13
 * Time: 2:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class Spannoid {

     ///  The chosen representation will be a tree with

    SubstitutionModel substModel;
    SubstitutionScore substScore;
    List<SpannoidComponent> components;






    /**
     * This constructor generates a tree and puts aligned sequences onto it.
     * The sequences are aligned using an iterative alignment scheme and Neighbor Joining.
     * First pairwise alignments are used to calculate distances, and these distances are used
     * to construct a tree using Neighbor Joining. As the tree is being constructed, sequences
     * are aligned together in an iterative manner.
     * @param sequences The sequences themselves.
     * @param names     The name of the sequences. The code calling this constructor is responsible to
     *                  have the same order of sequences and their names, as it is assumed that the name of
     *                  sequences[i] is names[i].
     * @param model     The time-continuous Markov model describing the substitution process of
     *                  sequences
     * @param ss        The substitution score matrix that is used to obtain pairwise distances.
     * @param filename  The name of the file that contained the sequences. This will be the name of
     *                  the tree, appearing in the graphical interface showing multiple alignments.
     */
    public Spannoid(String[] sequences, String[] names, SubstitutionModel model, SubstitutionScore ss) throws StoppedException {
        this.names = names;
        substitutionModel = model;
        hmm2 = new HmmTkf92(null);
        hmm3 = new HmmNonParam();

        //reading the sequences, transforming them into integer arrays, according to the model
        int[][][] seq = new int[sequences.length][][];
        for (int i = 0; i < sequences.length; i++) {
            //System.out.println(sequences[i]);
            //System.out.println(ss.which);
            int k = 0;
            for (int j = 0; j < sequences[i].length(); j++) {
                int sum = 0;
                char ch = sequences[i].charAt(j);
                for (int l = 0; l < ss.which[ch].length; l++) {
                    sum += ss.which[ch][l];
                }
                if (sum > 0) {
                    k++;
                }
            }
            seq[i] = new int[k][model.e.length];
            k = 0;
            for (int j = 0; j < sequences[i].length(); j++) {
                int sum = 0;
                char ch = sequences[i].charAt(j);
                for (int l = 0; l < ss.which[ch].length; l++) {
                    sum += ss.which[ch][l];
                }
                if (sum > 0) {
                    seq[i][k] = ss.which[ch];
                    k++;
                }
            }
            //			for(int j = 0; j < seq[i].length; j++){
            //	System.out.print(seq[i][j]+"\t");
            //}
            //System.out.println();
        }
        // now the pairwise distances
        int[][] dist = new int[seq.length][seq.length];
        int[][] d = null;
        int[][] p = null;
        int[][] q = null;
        int[][] charDist = ss.dist;
        for (int k = 0; k < seq.length; k++) {
            for (int l = 0; l <= k; l++) {
                stoppable();
                if (d == null || d.length < seq[k].length + 1 || d[0].length < seq[l].length + 1) {
                    d = new int[seq[k].length + 1][seq[l].length + 1];
                    p = new int[seq[k].length + 1][seq[l].length + 1];
                    q = new int[seq[k].length + 1][seq[l].length + 1];
                }
                ////////////
                d[0][0] = 0;
                for (int i = 1; i <= seq[k].length; i++) {
                    d[i][0] = q[i][0] = GAPOPEN - GAPEXT + i * GAPEXT;
                }
                for (int j = 1; j <= seq[l].length; j++) {
                    d[0][j] = p[0][j] = GAPOPEN - GAPEXT + j * GAPEXT;
                }
                for (int i = 1; i <= seq[k].length; i++) {
                    for (int j = 1; j <= seq[l].length; j++) {
                        p[i][j] = Math.min(d[i - 1][j] + GAPOPEN, p[i - 1][j] + GAPEXT);
                        q[i][j] = Math.min(d[i][j - 1] + GAPOPEN, q[i][j - 1] + GAPEXT);
                        int x, y;
                        for (x = 0; seq[k][i - 1][x] == 0; x++) {
                            ;
                        }
                        for (y = 0; seq[l][j - 1][y] == 0; y++) {
                            ;
                        }
                        d[i][j] = Math.min(d[i - 1][j - 1] + charDist[x][y],
                                Math.min(p[i][j], q[i][j]));
                    }
                }
                dist[k][l] = dist[l][k] = d[seq[k].length][seq[l].length];
                ////////////
                //System.out.println(k+"\t"+l+"\t"+dist[k][l]);
            }

        }

        //// Neighbor Joining algorithm based on the distances calculated above
        // initialization
        vertex = new Vertex[2 * seq.length - 1];
//        for (int i = 0; i < vertex.length; i++) {
//            vertex[i] = new Vertex();
//        }
        double[] sumDist = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            sumDist[i] = 0;
            for (int j = 0; j < dist.length; j++) {
                sumDist[i] += dist[i][j];
            }
        }
        int[] where = new int[dist.length];
        for (int i = 0; i < where.length; i++) {
            where[i] = i;
        }
        // the first n vertices will be the leaves
        for (int i = 0; i < seq.length; i++) {
            vertex[i] = new Vertex(this, 0.0, seq[i], names[i], sequences[i]);
        }
        // NJ main recursion
        int vnum = seq.length;
        Vertex newVert;
        for (int remN = dist.length; remN > 1; remN--) {
            stoppable();
            double minVal = BIGNUM;
            double val = 0.0;
            int i = -1;
            int j = -1;
            for (int k = 1; k < dist.length; k++) {
                for (int l = 0; l < k; l++) {
                    if (where[k] > -1 && where[l] > -1 && (val = (remN - 2) * dist[k][l] - sumDist[k] - sumDist[l]) < minVal) {
                        i = k;
                        j = l;
                        minVal = val;
                    }
                }
            }
            newVert = new Vertex(this, 0.0);    /* new vertex */
            vertex[vnum] = newVert;
            newVert.left = vertex[where[i]];
            newVert.right = vertex[where[j]];
            //System.out.println("Joining vertices "+where[i]+" and "+where[j]);
            newVert.parent = null;
            newVert.left.parent = newVert.right.parent = newVert;
            newVert.left.edgeLength = dist[i][j] / 2 - (remN > 2 ? (sumDist[i] - sumDist[j]) / (2 * remN - 4) : 0.001);
            newVert.right.edgeLength = dist[i][j] - newVert.left.edgeLength;

            val = (newVert.left.length + newVert.right.length) / 0.2;
            newVert.left.edgeLength /= val;
            newVert.right.edgeLength /= val;

            if (newVert.left.edgeLength < 0.01) {
                newVert.left.edgeLength = 0.01;
            }
            if (newVert.right.edgeLength < 0.01) {
                newVert.right.edgeLength = 0.01;
            }


            //	    newVert.left.edgeLength = 0.1;
            //newVert.right.edgeLength = 0.1;

            newVert.left.edgeChangeUpdate();
            newVert.right.edgeChangeUpdate();

            AlignColumn fake = new AlignColumn(newVert);
            newVert.first = fake;
            newVert.last = fake;
            newVert.length = 0;
            fake.left = newVert.left.last;
            fake.right = newVert.right.last;
            newVert.left.last.parent = fake;
            newVert.left.last.orphan = false;
            newVert.right.last.parent = fake;
            newVert.right.last.orphan = false;

            //			/* FAKE!!!! alignment */
            //			AlignColumn prev = new AlignColumn(newVert);
            //			newVert.first = prev;
            //			prev.seq = new double[model.e.length];
            //			AlignColumn left = newVert.left.first;
            //			AlignColumn right = newVert.right.first;
            //			left.parent = prev;
            //			left.orphan = false;
            //			right.parent = prev;
            //			right.orphan = false;
            //			left = left.next;
            //			right = right.next;
            //			while(left.next != null && right.next != null){
            //				AlignColumn actual = new AlignColumn(newVert);
            //				actual.seq = new double[model.e.length];
            //				actual.prev = prev;
            //				prev.next = actual;
            //				left.parent = actual;
            //				left.orphan = false;
            //				right.parent = actual;
            //				right.orphan = false;
            //				right = right.next;
            //				left = left.next;
            //				prev = actual;
            //			}
            //			AlignColumn actual = new AlignColumn(newVert);
            //			actual.seq = new double[model.e.length];
            //			actual.prev = prev;
            //			prev.next = actual;
            //
            //			newVert.last = actual;

            //??????

            newVert.left.fullWin();
            newVert.right.fullWin();
            newVert.fullWin();

            newVert.left.selected = true;
            newVert.right.selected = true;
            newVert.hmm3AlignWithSave();

            //System.out.println("length of the ancestral sequence: "+newVert.length);

            //String[] s = newVert.left.printedAlignment();
            //System.out.println(s[0]+"\n"+s[1]+"\n");

            //s = newVert.right.printedAlignment();
            //System.out.println(s[0]+"\n"+s[1]+"\n");

            where[i] = where[j] = -1;
            sumDist[i] = 0;
            for (int a = 0; a < dist.length; a++) {
                if (where[a] > -1) {
                    sumDist[a] -= dist[a][i] + dist[a][j];
                    dist[a][i] = dist[i][a] = (dist[a][i] + dist[a][j] - dist[i][j]) / 2;
                    sumDist[a] += dist[a][i];
                    sumDist[i] += dist[a][i];
                }
            }

            where[i] = vnum;
            vnum++;
        }
        root = vertex[vnum - 1];
        root.calcOrphan();
        root.calcFelsRecursively();
        /////////////

    }





}
