package statalign.base;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;
import java.util.List;

/**
 * This is a vertex of the tree.
 * The 'hardcore' functions are implemented in this class, developers are suggested
 * not change functions in it. The implemented functions are quite unreadable, since we
 * opted for efficiency and not readability. You should be able to develop novel
 * functionality of the software package (postprocessing, substitution models, etc.)
 * without touching this class.
 * @author miklos, novak
 */
public class Vertex extends TopologyNode{

    public Vertex parent;
    public Vertex left;
    public Vertex right;
    public Vertex old;


    Vertex(){
        
    }

    Vertex(Tree owner, double edgeLength) {
        this.owner = owner;
        this.edgeLength = edgeLength;
        updateTransitionMatrix();
        indelLogLike = 0.0;
        old = new Vertex();
    }

    Vertex(Tree owner, double edgeLength, int[][] seq, String name, String origSeq) {
        this.owner = owner;
        this.edgeLength = edgeLength;
        this.name = new String(name);
        this.seq = origSeq;
        // specify if the Vertex is labeled or unlabeled
        if (origSeq != null){
            this.isLabeled = true;
        }
        else{
            this.isLabeled = false;

        }
        old = new Vertex();
        int size = owner.substitutionModel.e.length;
        first = new AlignColumn(this);
        first.seq = new double[size];
        for (int i = 0; i < first.seq.length; i++) {
            first.seq[i] = seq[0][i];
        }
        first.prev = null;
        AlignColumn prev = first;
        for (int i = 1; i < seq.length; i++) {
            AlignColumn actual = new AlignColumn(this);
            prev.next = actual;
            actual.prev = prev;
            actual.seq = new double[size];
            for (int j = 0; j < first.seq.length; j++) {
                actual.seq[j] = seq[i][j];
            }
            prev = actual;
        }
        last = new AlignColumn(this);
        last.prev = prev;
        prev.next = last;
        length = seq.length;
        edgeChangeUpdate();
        indelLogLike = 0.0;
    }

    Vertex(String descriptor, Vertex parent) {

        /* descriptor is a string describing the Vertexs below this Vertex */

        //System.out.println(descriptor);

        this.parent = parent;
        owner = parent.owner;
        old = new Vertex();


        if (descriptor.charAt(0) == '(') {
            String leftDescriptor = "";
            int counter = 0;
            int j;
            for (j = 1; counter != 0 || descriptor.charAt(j) != ','; j++) {
                leftDescriptor += descriptor.charAt(j);
                if (descriptor.charAt(j) == '(') {
                    counter++;
                }
                if (descriptor.charAt(j) == ')') {
                    counter--;
                }
            }
            Vertex leftChild = new Vertex(leftDescriptor, this);
            left = leftChild;

            String rightDescriptor = "";
            for (j += 1; counter != -1; j++) {
                rightDescriptor += descriptor.charAt(j);
                if (descriptor.charAt(j) == '(') {
                    counter++;
                }
                if (descriptor.charAt(j) == ')') {
                    counter--;
                }
            }
            rightDescriptor = rightDescriptor.substring(0, rightDescriptor.length() - 1);
            Vertex rightChild = new Vertex(rightDescriptor, this);
            right = rightChild;
            String tempStringlength = "";
            for (j += 1; j < descriptor.length(); j++) {
                tempStringlength += descriptor.charAt(j);
            }
            //System.out.println(tempStringlength);
            if (tempStringlength.length() > 0) {
                edgeLength = Math.max(Double.parseDouble(tempStringlength), 0.01);
            } else {
                edgeLength = 0.01;
            }
            //  edgeLength = 0.01;
            //else we have the root Vertex, which does not need an edgeLength
        } else {
            left = null;
            right = null;
            name = "";
            int i = 0;
            while (descriptor.charAt(i) != ':') {
                name += descriptor.charAt(i);
                i++;
            }
            String tempString = "";
            for (i++; i < descriptor.length(); i++) {
                tempString += descriptor.charAt(i);
            }
            //System.out.println(tempString);
            edgeLength = Math.max(Double.parseDouble(tempString), 0.01);

        }

    }

    Vertex brother() {
        return parent.left == this ? parent.right : parent.left;
    }
}