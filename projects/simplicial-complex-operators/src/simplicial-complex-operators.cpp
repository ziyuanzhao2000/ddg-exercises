// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <set>

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        geometry->vertexIndices[v] = idx;
        ++idx;
    }

    idx = 0;
    for (Edge e : mesh->edges()) {
        geometry->edgeIndices[e] = idx;
        ++idx;
    }

    idx = 0;
    for (Face f : mesh->faces()) {
        geometry->faceIndices[f] = idx;
        ++idx;
    }
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> triplets;
    for (Edge e : mesh->edges()) {
        triplets.push_back(T(e.getIndex(), e.firstVertex().getIndex(), 1));
        triplets.push_back(T(e.getIndex(), e.secondVertex().getIndex(), 1));
    }
    Eigen::SparseMatrix<size_t> mat(mesh->nEdges(), mesh->nVertices());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> triplets;
    for (Face f : mesh->faces()) {
        for (Edge e: f.adjacentEdges()) {
            triplets.push_back(T(f.getIndex(),  e.getIndex(), 1));
        }
    }
    Eigen::SparseMatrix<size_t> mat(mesh->nFaces(), mesh->nEdges());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t i=0; i < mesh->nVertices(); i++){
        if (subset.vertices.count(i)){
            vec[i] = 1;
        }
    }
    return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t i=0; i < mesh->nEdges(); i++){
        if (subset.edges.count(i)){
            vec[i] = 1;
        }
    }
    return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t i=0; i < mesh->nFaces(); i++){
        if (subset.faces.count(i)){
            vec[i] = 1;
        }
    }
    return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO

    // Half-edge implementation, but not allowed!
    // std::set<size_t> V;
    // std::set<size_t> E;
    // std::set<size_t> F;

    // for (size_t vid : subset.vertices){
    //     V.insert(vid);
    //     Vertex v = mesh->vertex(vid);
    //     for (Edge e : v.adjacentEdges()){
    //         E.insert(e.getIndex());
    //     }
    //     for (Face f : v.adjacentFaces()){
    //         F.insert(f.getIndex());
    //     }
    // }

    // for (size_t eid : subset.edges){
    //     E.insert(eid);
    //     Edge e = mesh->edge(eid);
    //     for (Face f : e.adjacentFaces()){
    //         F.insert(f.getIndex());
    //     }
    // }

    // for (size_t fid : subset.faces){
    //     F.insert(fid);
    // }

    // return MeshSubset(V, E, F);

    // Adjacency matrix implementation, as required
    MeshSubset new_subset = subset.deepCopy();
    Vector<size_t> v = buildVertexVector(subset);
    Vector<size_t> e = A0 * v;
    Vector<size_t> f = A1 * e;
    for (size_t i=0; i<mesh->nEdges(); i++) {
        if (e[i]) {
            new_subset.addEdge(i);
        }
    }
    for (size_t i=0; i<mesh->nFaces(); i++) {
        if (f[i]) {
            new_subset.addFace(i);
        }
    }

    e = buildEdgeVector(subset);
    f = A1 * e;
    for (size_t i=0; i<mesh->nFaces(); i++) {
        if (f[i]) {
            new_subset.addFace(i);
        }
    }
    return new_subset;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset new_subset = subset.deepCopy();
    Vector<size_t> f = buildFaceVector(subset);
    Vector<size_t> e = A1.transpose() * f;
    Vector<size_t> v = A0.transpose() * e;
    for (size_t i=0; i<mesh->nEdges(); i++) {
        if (e[i]) {
            new_subset.addEdge(i);
        }
    }
    for (size_t i=0; i<mesh->nVertices(); i++) {
        if (v[i]) {
            new_subset.addVertex(i);
        }
    }

    e = buildEdgeVector(subset);
    v = A0.transpose() * e;
    for (size_t i=0; i<mesh->nVertices(); i++) {
        if (v[i]) {
            new_subset.addVertex(i);
        }
    }
    return new_subset;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    // Compute Cl * St - St * Cl
    MeshSubset set1 = closure(star(subset));
    MeshSubset set2 = star(closure(subset));
    set1.deleteSubset(set2);
    return set1;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    // The idea is to check if the subset of simplices is equal to its closure
    // if not, then some face does not include all its subsets so the subset
    // doesn't correspond to a simplicial complex.
    MeshSubset cl = closure(subset);
    return cl.equals(subset);
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    if (isComplex(subset) == false) {
        return -1;
    }

    Vector<size_t> e = buildEdgeVector(subset);
    Vector<size_t> f = buildFaceVector(subset);

    if (e.sum()==0 && f.sum()==0) {
        return 0; // only have vertices 
    } 
    
    // To check that it's pure, make sure every vertex in the subset 
    // is in some edge in the subset (deg >= 1)
    // and every edge is in some face (deg >= 2)
    Vector<size_t> v_all = buildVertexVector(subset);
    Vector<size_t> v_from_e = A0.transpose() * e;
    for (auto idx=0; idx<v_all.size(); idx++){
        if (v_from_e[idx] < v_all[idx]) {   // not ideal, how to cast into matrix operation?
            return -1;
        }
    }

    if (f.sum()==0) {
        return 1; // only have vertices and edges 
    } 

    Vector<size_t> e_all = buildEdgeVector(subset);
    Vector<size_t> e_from_f = A1.transpose() * f;
    for (auto idx=0; idx<e_all.size(); idx++){
        if (e_from_f[idx] < e_all[idx]) {
            return -1;
        }
    } 

    return 2;

}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    int degree = isPureComplex(subset);

    MeshSubset pre_bd;

    // we do not have a straightforward definition for the boundary of non-pure simplicial complex
    // Boundary of vertices is empty

    // If we have edges, the boundary will consist of only vertices
    // If we have faces, the boundary will consist of only edges
    if (degree == 1) { 
        Vector<size_t> v_from_e = A0.transpose() * buildEdgeVector(subset);
        for (size_t idx=0; idx<mesh->nVertices(); idx++){
            if (v_from_e[idx]==1) {
                pre_bd.addVertex(idx);
            }
        }
    } else if (degree == 2) {
        Vector<size_t> e_from_f = A1.transpose() * buildFaceVector(subset);
        for (size_t idx=0; idx<mesh->nEdges(); idx++){
            if (e_from_f[idx]==1) {
                pre_bd.addEdge(idx);
            }
        }
    }

    return closure(pre_bd);
}