package statalign.base;

import statalign.model.score.SubstitutionScore;
import statalign.postprocess.plugins.TreeNode;

public class NJTree {
    private static final int GAPOPEN = 9;
    private static final int GAPEXT = 2;

    public static String getNJTree(int[][][] seq, String[] names, SubstitutionScore ss) {
        // now the pairwise distances
        int[][] dist = new int[seq.length][seq.length];
        int[][] d = null;
        int[][] p = null;
        int[][] q = null;
        int[][] charDist = ss.dist;
        for (int k = 0; k < seq.length; k++) {
            for (int l = 0; l <= k; l++) {
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
            }

        }

        //// Neighbor Joining algorithm based on the distances calculated above
        // initialization
        TreeNode[] nodes = new TreeNode[2 * seq.length - 1];

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
            // vertex[i] = new Vertex(this, 0.0, seq[i], names[i], sequences[i]);
            nodes[i] = new TreeNode(names[i], 0.0);
        }
        // NJ main recursion
        int vnum = seq.length;
        TreeNode newNode;
        for (int remN = dist.length; remN > 1; remN--) {
            double minVal = Double.MAX_VALUE;
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

            newNode = new TreeNode(null, 0.0);
            nodes[vnum] = newNode;
            newNode.addChild(nodes[where[i]]);
            newNode.addChild(nodes[where[j]]);

            TreeNode left = newNode.children.get(0);
            TreeNode right = newNode.children.get(1);
            left.edgeLength = dist[i][j] / 2 - (remN > 2 ? (sumDist[i] - sumDist[j]) / (2 * remN - 4) : 0.001);
            right.edgeLength = dist[i][j] - left.edgeLength;

            double length = (seq[i].length + seq[j].length) / 2;
            double scale = 0.1 / length;
            left.edgeLength *= scale;
            right.edgeLength *= scale;

            left.edgeLength = Math.max(left.edgeLength, 0.01);
            right.edgeLength = Math.max(right.edgeLength, 0.01);

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
        TreeNode root = nodes[vnum - 1];
        return root.toString();
    }
}
