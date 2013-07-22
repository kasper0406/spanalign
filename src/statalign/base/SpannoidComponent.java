package statalign.base;

import statalign.base.thread.StoppedException;
import statalign.model.score.SubstitutionScore;
import statalign.model.subst.SubstitutionModel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Aldo Pacchiano
 * Date: 7/19/13
 * Time: 2:33 PM
 * To change this template use File | Settings | File Templates.
 */

///// This object represents a component of a spannoid.
public class SpannoidComponent {
    List<Vertex> vertices = new ArrayList<Vertex>(); /// List of vertices within this component
    Tree tree = new Tree();   /// Tree within this component
    /// We need a map between this component and the ones to which it is connected.
    HashMap<Vertex, HashSet<SpannoidComponent>> connections = new HashMap<Vertex, HashSet<SpannoidComponent>>();
    int k = 0; /// k is the number of unlabeled nodes within the component.



    SpannoidComponent(String[] sequences, String[] names, SubstitutionModel model, SubstitutionScore ss) throws StoppedException{
      tree = new Tree(sequences, names, model, ss);




    }

    /// Return a random vertex of the component
    public Vertex getRandomVertex(){
        Vertex v = new Vertex();
        return v;
    }

    /// Removes the specified vertex from the component
    public void removeVertex(Vertex v){
        vertices.remove(v);
        connections.remove(v);
        /// remove the vertex from the tree IMPORTANT


    }


    /// Adds the specified vertex into the component.
    public void addVertex(Vertex v, HashSet<SpannoidComponent> neighbors){
            vertices.add(v);
            connections.put(v, neighbors);
            /// remove the vertex from the tree IMPORTANT
    }







}
