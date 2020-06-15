//
// Created by Adam Christa on 2/26/20.
//

#ifndef GRAPHPRACTICE_EDGE_H
#define GRAPHPRACTICE_EDGE_H

#include "Node.h"
#include "Graph.h"

/*
 *   This is the Edge class. The edge class represents a an edge within a graph class. The edge connects two
 *   adjacent nodes together. The edge has a begginning node name and an ending Node. The ending Node is an actual
 *   node data structure
 */

template <class T> class Node;
template <class T> class Graph;

template<class T>
class Edge {


    int beginningNodeName;
    int distance;
    Node<T>* end;
    int numPaths;
    double betweennessValue;


    explicit Edge(Node<T>*);


public:

    int getEndName();
    Node<T>* getEndNode();
    void setBeginningNodeName(int);
    int getBeginningNodeName();
    double getBetweennessValue() {return betweennessValue;}
    void setBetweennessValue(double temp) {betweennessValue = temp;}

    void setNumPaths(int x) {numPaths = x;}
    int getNumPaths(){return numPaths;}

    friend class Graph<T>;
    friend class Node<T>;


    int getDistance();
    void setDistance(int);
};

template<class T>
Node<T> *Edge<T>::getEndNode() {
    return end;
}

template<class T>
int Edge<T>::getEndName() {
    return end->getNodeName();
}


template <class T>
Edge<T>::Edge(Node<T> *e)
{
    end = e;
    beginningNodeName = 0;
    distance = 0;
    betweennessValue = 0;
    numPaths = 0;
}

template<class T>
void Edge<T>::setBeginningNodeName(int temp) {
    beginningNodeName = temp;
}

template<class T>
int Edge<T>::getBeginningNodeName() {
    return beginningNodeName;

}

template<class T>
int Edge<T>::getDistance() {
    return distance;
}

template<class T>
void Edge<T>::setDistance(int temp) {
    distance = temp;
}


#endif //GRAPHPRACTICE_EDGE_H
