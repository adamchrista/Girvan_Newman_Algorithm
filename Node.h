//
// Created by Adam Christa on 2/26/20.
//

#ifndef GRAPHPRACTICE_NODE_H
#define GRAPHPRACTICE_NODE_H

#include <vector>
#include "Graph.h"

/*
 *  This is the Node class. The node class represents a node in an adjacency list that contains a vector of edges.
 *  It is a large part of the adjacency list data structure. The Node has a name, an iterator number, and most
 *  importantly a vector of Edges. The node therefore can be used to find all edges that are adjacent to it
 *
 */

using namespace std;

template <class T> class Edge;
template<class T>  class Graph;

template<class T>
class Node {

public:

private:

    int nodeName;
    int iteratorNumber;
    int distanceFromTop;

    T data;

    vector<Edge<T>*> adj;
    bool visited;

    Node<T>* path;

public:

    int getNodeName();
    int getIteratorNumber();
    void incrementDistanceFromTop();
    void setDistanceFromTop(int);
    int getDistanceFromTop();
    void setIteratorNumber(int);
    int getAdjSize();
    bool contains(int);

    explicit Node(int);
    //explicit Node(T in);
    Node(const Node<T> &in);
    Node();
    void addEdge(Node<T>*);
    Edge<T>* getEdge(int);
    void printEdges();

    friend class Graph<T>;
};


template <class T>
Node<T>::Node(): visited(false), path(NULL), iteratorNumber(0), distanceFromTop(0), nodeName(0){}


template <class T>
Node<T>::Node(int name): nodeName(name), visited(false), path(NULL), iteratorNumber(0), distanceFromTop(0){}


template <class T>
Node<T>::Node(const Node<T> & in): data(in.info), visited(in.visited), path(NULL), iteratorNumber(0), distanceFromTop(0), nodeName(in.nodeName){}


template<class T>
void Node<T>::addEdge(Node<T> *end) {

    auto* edgeD = new Edge<T>(end);

    adj.push_back(edgeD);

}

template<class T>
int Node<T>::getNodeName() {
    return nodeName;
}

template<class T>
bool Node<T>::contains(int nameNumber) {

    int counter = 0;
    while (counter < adj.size())
    {
        if (adj[counter]->getEndName() == nameNumber)
        {
            return true;
        }
        counter++;
    }

    return false;

}

template<class T>
void Node<T>::printEdges() {

    for (int i = 0; i < adj.size(); i++)
    {
        cout << "Node " << adj[i]->getEndName() << " -> ";
    }

}

template<class T>
int Node<T>::getIteratorNumber() {
    return iteratorNumber;
}

template<class T>
void Node<T>::setIteratorNumber(int temp) {
    iteratorNumber = temp;
}

template<class T>
int Node<T>::getAdjSize() {
    return adj.size();
}

template<class T>
Edge<T>* Node<T>::getEdge(int number) {
    return adj[number];
}

template<class T>
void Node<T>::incrementDistanceFromTop() {
    distanceFromTop++;
}

template<class T>
void Node<T>::setDistanceFromTop(int temp) {
    distanceFromTop = temp;
}

template<class T>
int Node<T>::getDistanceFromTop() {
    return distanceFromTop;
}



#endif //GRAPHPRACTICE_NODE_H
