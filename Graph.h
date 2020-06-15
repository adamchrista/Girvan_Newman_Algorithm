//
// Created by Adam Christa on 2/25/20.
//

#ifndef GRAPHPRACTICE_GRAPH_H
#define GRAPHPRACTICE_GRAPH_H

#include <utility>
#include <string>
#include <iostream>
#include <queue>
#include <vector>
#include <stack>
#include "Node.h"
#include "Edge.h"

/*
 * This is the graph class. The main data structure behind the Graph class is an adjacency list
 * that is made up of a vector of Node objects
 *
 */

using namespace std;

template<class T>
class Graph {

    friend class Node<T>;

private:
    vector<Node<T>*> adjListOfNodes;
    void createAndPushNode(int, int);
    void dfsInternal(Node<T>*, vector<T>&) const;
    bool isThereAZeroVertex();
    void setNodeIteratorNumber(int, int);
    void setNodeDistanceFromTop(int, int);
    vector<int> TarjansAlgorithmInner(int, int[], int[], stack<int>*, bool[]);
    vector<Edge<T>*> calculateBetweenness(vector<Edge<T>*>, int[]);


public:

    explicit Graph();

    void addRelationship(int, int);
    void printGraph();
    Node<T>* getNode(int tempNodeName);
    int getNumberOfNodes() {return adjListOfNodes.size();}

    void makeAConnection(int, int, map<int, string>, ofstream&);
    vector<Edge<T>*> dfs(int) ;

    vector<vector<int>> TargansAlgorithm();

    vector<Edge<T>*> bfs(int, int&);

    void printBFSAndDFS(vector<Edge<T>*>);
    void clear() {adjListOfNodes.clear();}

    void resetIteratorNumbers();



};


template<class T>
Graph<T>::Graph() {



}

/*
 * This is the resetIteratorNumbers function. It resets iterator numbers back to 0 in order to perform
 * a breadth first search on a graph more than once.
 *
 */
template<class T>
void Graph<T>::resetIteratorNumbers() {

    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        adjListOfNodes[i]->setIteratorNumber(0);
        for (int j = 0; j < adjListOfNodes[i]->getAdjSize(); j++)
        {
            adjListOfNodes[i]->getEdge(j)->getEndNode()->setIteratorNumber(0);
        }
    }

}

/*
 * This is the create and push node function. It creates a node with the proper edge
 * and pushes it onto the adjacency list
 */
template <class T>
void Graph<T>::createAndPushNode(int name1, int name2) {

    auto *v1 = new Node<T>(name1);
    auto *v2 = new Node<T>(name2);
    v1->addEdge(v2);

    adjListOfNodes.push_back(v1);

}

/*
 *  This is the make a connection function. It performs the mc command that is required for the project. It
 *  uses a breadth first search in order to find the shortest path between two nodes
 */
template<class T>
void Graph<T>::makeAConnection(int firstNode, int secondNode, map<int, string> reversedMap, ofstream& outputFile) {

    int pred[adjListOfNodes.size()];
    Node<T>* tempNode = getNode(firstNode);
    vector<Edge<T>*> vectorOfEdges;
    queue<Node<T>*> theQueue;
    int i = 1;

    int dist[adjListOfNodes.size()];
    int paths[adjListOfNodes.size()];
    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        dist[i] = 100000;
        pred[i] = -1;
    }
    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        paths[i] = 0;
    }


    dist[tempNode->getNodeName() - 1] = 0;
    paths[tempNode->getNodeName() - 1] = 1;

    int stopEverything = 0;

    int counter = 0;

    while (isThereAZeroVertex() && counter < 1000 && stopEverything == 0)
    {
        setNodeIteratorNumber(tempNode->getNodeName(), i++);
        theQueue.push(tempNode);

        counter++;

        while (!theQueue.empty() && stopEverything == 0)
        {
            tempNode = theQueue.front();
            Node<T>* tempNode2 = getNode(tempNode->getNodeName());
            theQueue.pop();
            vector<int> distances;

            for (int j = 0; j < tempNode2->getAdjSize(); j++)
            {

                if (tempNode2->getEdge(j)->getEndNode()->getNodeName() == secondNode)
                {
                    stopEverything = 1;
                }
                if (tempNode2->getEdge(j)->getEndNode()->getIteratorNumber() == 0)
                {
                    setNodeIteratorNumber(tempNode2->getEdge(j)->getEndNode()->getNodeName(), i++);
                    theQueue.push(tempNode2->getEdge(j)->getEndNode());
                    tempNode2->getEdge(j)->setBeginningNodeName((tempNode->getNodeName()));
                    vectorOfEdges.push_back(tempNode2->getEdge(j));
                    tempNode2->getEdge(j)->getEndNode()->setDistanceFromTop(tempNode2->getDistanceFromTop());
                    tempNode2->getEdge(j)->getEndNode()->incrementDistanceFromTop();
                    setNodeDistanceFromTop(tempNode2->getEdge(j)->getEndNode()->getNodeName(), tempNode2->getEdge(j)->getEndNode()->getDistanceFromTop());
                    pred[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] = tempNode2->getNodeName();
                }


                if (dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] > dist[tempNode2->getNodeName() - 1] + 1)
                {
                    dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] = dist[tempNode2->getNodeName() - 1] + 1;
                    paths[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] = paths[tempNode2->getNodeName() - 1];

                }
                else if (dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] == dist[tempNode2->getNodeName() - 1] + 1)
                {
                    paths[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] += paths[tempNode2->getNodeName() - 1];
                    tempNode2->getEdge(j)->setBeginningNodeName((tempNode->getNodeName()));
                    vectorOfEdges.push_back(tempNode2->getEdge(j));
                }



            }

        }

    }


    vector<int> path;
    int crawl = secondNode;
    path.push_back(crawl);
    while (pred[crawl - 1] != -1)
    {
        path.push_back(pred[crawl - 1]);
        crawl = pred[crawl -1];
    }

    cout << "{";
    outputFile << "{";

    for (int i = path.size() - 1; i > 0; i--)
    {
        if (i > 1)
        {
            cout << "(" << reversedMap[path[i]] << " - " << reversedMap[path[i-1]] << "), ";
            outputFile << "(" << reversedMap[path[i]] << " - " << reversedMap[path[i-1]] << "), ";
        }

        if (i == 1)
        {
            cout << "(" << reversedMap[path[i]] << " - " << reversedMap[path[i-1]] << ")";
            outputFile << "(" << reversedMap[path[i]] << " - " << reversedMap[path[i-1]] << ")";
        }
    }

    cout << "}" << endl;
    outputFile << "}" << endl;

}

/*
 *  This is the addRelationship function. It does the heavy lifting when a new graph is being created by accepting
 *  two node names and putting those node names into the adjacency list with the proper relationship.
 */
template<class T>
void Graph<T>::addRelationship(int nameNumber1, int nameNumber2) {


    bool doesNodeOneNeedToBeAdded = true;
    bool doesNodeTwoNeedToBeAdded = true;


    if (adjListOfNodes.size() != 0)
    {
        int counter = 0;
        while (counter < adjListOfNodes.size())
        {
            int temp = adjListOfNodes[counter]->getNodeName();
            if ((temp == nameNumber1) && (adjListOfNodes[counter]->contains(nameNumber2) == false))
            {
                auto *v1 = new Node<T>(nameNumber2);
                adjListOfNodes[counter]->addEdge(v1);
                doesNodeOneNeedToBeAdded = false;
            }
            if ((temp == nameNumber2) && (adjListOfNodes[counter]->contains(nameNumber1) == false))
            {
                auto *v1 = new Node<T>(nameNumber1);
                adjListOfNodes[counter]->addEdge(v1);
                doesNodeTwoNeedToBeAdded = false;
            }

            if ((temp == nameNumber1) && (adjListOfNodes[counter]->contains(nameNumber2) == true))
            {
                doesNodeOneNeedToBeAdded = false;
            }

            if ((temp == nameNumber2) && (adjListOfNodes[counter]->contains(nameNumber1) == true))
            {
                doesNodeTwoNeedToBeAdded = false;
            }

            counter++;
        }

    }

    if (adjListOfNodes.size() == 0)
    {
        createAndPushNode(nameNumber1, nameNumber2);
        createAndPushNode(nameNumber2, nameNumber1);

        doesNodeOneNeedToBeAdded = false;
        doesNodeTwoNeedToBeAdded = false;
    }

    if (doesNodeOneNeedToBeAdded)
    {
        createAndPushNode(nameNumber1, nameNumber2);
    }

    if (doesNodeTwoNeedToBeAdded)
    {
        createAndPushNode(nameNumber2, nameNumber1);
    }

}

/*
 * Does a depth first search internally. Uses a stack instead of a queue in order to pop off nodes in the right
 * order.
 */
template <class T>
void Graph<T>::dfsInternal(Node<T> *v,vector<T> &res) const {
    v->visited = true;
    res.push_back(v->getNodeName());
    typename vector<Edge<T> >::iterator it = (v->adj).begin();
    typename vector<Edge<T> >::iterator ite = (v->adj).end();
    for (; it != ite; it++)
        if (!it->end->visited)
        {
            dfsInternal(it->end, res);
        }
}


template<class T>
void Graph<T>::printGraph() {

    int counter = 0;
    while (counter < adjListOfNodes.size())
    {
        int nodeName = adjListOfNodes[counter]->getNodeName();
        cout << "Node " << nodeName << " -> ";
        adjListOfNodes[counter]->printEdges();
        cout << endl;
        counter++;
    }



}

/*
 *  A depth first search. The function accepts a node name and outputs all of the edges that the
 *  are found from the results of the depth first search.
 */
template<class T>
vector<Edge<T>*> Graph<T>::dfs(int tempNodeName) {

    Node<T>* tempNode = getNode(tempNodeName);
    vector<Edge<T>*> vectorOfEdges;
    stack<Node<T>*> theStack;
    int i = 1;
    int dist[adjListOfNodes.size()];
    int paths[adjListOfNodes.size()];
    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        dist[i] = 100000;
    }
    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        paths[i] = 0;
    }


    dist[tempNode->getNodeName() - 1] = 0;
    paths[tempNode->getNodeName() - 1] = 1;

    while (isThereAZeroVertex())
    {
        setNodeIteratorNumber(tempNode->getNodeName(), i++);
        theStack.push(tempNode);



        while (!theStack.empty())
        {
            tempNode = theStack.top();
            Node<T>* tempNode2 = getNode(tempNode->getNodeName());
            theStack.pop();
            vector<int> distances;

            for (int j = 0; j < tempNode2->getAdjSize(); j++)
            {
                if (tempNode2->getEdge(j)->getEndNode()->getIteratorNumber() == 0)
                {
                    setNodeIteratorNumber(tempNode2->getEdge(j)->getEndNode()->getNodeName(), i++);
                    theStack.push(tempNode2->getEdge(j)->getEndNode());
                    tempNode2->getEdge(j)->setBeginningNodeName((tempNode->getNodeName()));
                    vectorOfEdges.push_back(tempNode2->getEdge(j));
                    tempNode2->getEdge(j)->getEndNode()->setDistanceFromTop(tempNode2->getDistanceFromTop());
                    tempNode2->getEdge(j)->getEndNode()->incrementDistanceFromTop();
                    setNodeDistanceFromTop(tempNode2->getEdge(j)->getEndNode()->getNodeName(), tempNode2->getEdge(j)->getEndNode()->getDistanceFromTop());
                }

                if (dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] > dist[tempNode2->getNodeName() - 1] + 1)
                {
                    dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] = dist[tempNode2->getNodeName() - 1] + 1;
                    paths[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] = paths[tempNode2->getNodeName() - 1];
                }
                else if (dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] == dist[tempNode2->getNodeName() - 1] + 1)
                {
                    paths[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] += paths[tempNode2->getNodeName() - 1];
                    tempNode2->getEdge(j)->setBeginningNodeName((tempNode->getNodeName()));
                    vectorOfEdges.push_back(tempNode2->getEdge(j));
                }

            }

        }
    }


    printBFSAndDFS(vectorOfEdges);
    return vectorOfEdges;
}

/*
 *  A typical breadth first search. It will output a list of edges between nodes that were found
 */
template<class T>
vector<Edge<T>*> Graph<T>::bfs(int tempNodeName, int& stop) {

    Node<T>* tempNode = getNode(tempNodeName);
    vector<Edge<T>*> vectorOfEdges;
    queue<Node<T>*> theQueue;
    int i = 1;
    int dist[adjListOfNodes.size()];
    int paths[adjListOfNodes.size()];
    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        dist[i] = 100000;
    }
    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        paths[i] = 0;
    }


    dist[tempNode->getNodeName() - 1] = 0;
    paths[tempNode->getNodeName() - 1] = 1;

    int counter = 0;

    while (isThereAZeroVertex() && counter < 1000)
    {
        setNodeIteratorNumber(tempNode->getNodeName(), i++);
        theQueue.push(tempNode);

        if (counter == 990)
        {
         //   cout << "Hi";
        }

        counter++;

        while (!theQueue.empty())
        {
            tempNode = theQueue.front();
            Node<T>* tempNode2 = getNode(tempNode->getNodeName());
            theQueue.pop();
            vector<int> distances;

            for (int j = 0; j < tempNode2->getAdjSize(); j++)
            {
                if (tempNode2->getEdge(j)->getEndNode()->getIteratorNumber() == 0)
                {
                    setNodeIteratorNumber(tempNode2->getEdge(j)->getEndNode()->getNodeName(), i++);
                    theQueue.push(tempNode2->getEdge(j)->getEndNode());
                    tempNode2->getEdge(j)->setBeginningNodeName((tempNode->getNodeName()));
                    vectorOfEdges.push_back(tempNode2->getEdge(j));
                    tempNode2->getEdge(j)->getEndNode()->setDistanceFromTop(tempNode2->getDistanceFromTop());
                    tempNode2->getEdge(j)->getEndNode()->incrementDistanceFromTop();
                    setNodeDistanceFromTop(tempNode2->getEdge(j)->getEndNode()->getNodeName(), tempNode2->getEdge(j)->getEndNode()->getDistanceFromTop());
                }

                if (dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] > dist[tempNode2->getNodeName() - 1] + 1)
                {
                    dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] = dist[tempNode2->getNodeName() - 1] + 1;
                    paths[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] = paths[tempNode2->getNodeName() - 1];
                }
                else if (dist[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] == dist[tempNode2->getNodeName() - 1] + 1)
                {
                    paths[tempNode2->getEdge(j)->getEndNode()->getNodeName() - 1] += paths[tempNode2->getNodeName() - 1];
                    tempNode2->getEdge(j)->setBeginningNodeName((tempNode->getNodeName()));
                    vectorOfEdges.push_back(tempNode2->getEdge(j));
                }
            }
        }
    }

    if (counter > 990)
    {
        stop = 1;

    }
    vector<Edge<T>*> newVectorOfEdges = calculateBetweenness(vectorOfEdges, paths);
    return newVectorOfEdges;
}

template<class T>
Node<T> *Graph<T>::getNode(int tempNodeName) {

    int i = 0;
    while (i < adjListOfNodes.size())
    {
        if (adjListOfNodes[i]->getNodeName() == tempNodeName)
        {
            return adjListOfNodes[i];
        }
        i++;
    }

    return nullptr;

}

/*
 *  Iterates through the adjacency list to determine whether or not any nodes haven't been visited
 *  yet. It is essential to breadth first search
 */
template<class T>
bool Graph<T>::isThereAZeroVertex() {

    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        if (adjListOfNodes[i]->getIteratorNumber() == 0)
        {
            return true;
        }
    }

    return false;
}

/*
 *  Iterates through the graph in order to give all nodes and edges containing node name the same
 *  iterator number so that they are not zero
 */
template<class T>
void Graph<T>::setNodeIteratorNumber(int nodeName, int iteratorNumber) {

    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        if (adjListOfNodes[i]->getNodeName() == nodeName)
        {
            adjListOfNodes[i]->setIteratorNumber(iteratorNumber);
        }
        else
        {
            for (int j = 0; j < adjListOfNodes[i]->getAdjSize(); j++)
            {
                if (adjListOfNodes[i]->getEdge(j)->getEndNode()->getNodeName() == nodeName)
                {
                    adjListOfNodes[i]->getEdge(j)->getEndNode()->setIteratorNumber(iteratorNumber);
                }
            }
        }
    }

}

/*
 *  Sets all nodes of same node name to be the same distance from top of searching node
 */
template<class T>
void Graph<T>::setNodeDistanceFromTop(int nodeName, int theDistanceFromTop) {

    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        if (adjListOfNodes[i]->getNodeName() == nodeName)
        {
            adjListOfNodes[i]->setDistanceFromTop(theDistanceFromTop);
        }
        else
        {
            for (int j = 0; j < adjListOfNodes[i]->getAdjSize(); j++)
            {
                if (adjListOfNodes[i]->getEdge(j)->getEndNode()->getNodeName() == nodeName)
                {
                    adjListOfNodes[i]->getEdge(j)->getEndNode()->setDistanceFromTop(theDistanceFromTop);
                }
            }
        }
    }

}

/*
 * Prints the results of breadth first search and depth first search
 */
template<class T>
void Graph<T>::printBFSAndDFS(vector<Edge<T>*> results) {

    for (int i = 0; i < results.size(); i++)
    {
        cout << results[i]->getBeginningNodeName() << " -- " << results[i]->getEndNode()->getNodeName() << endl;
    }


}

/*
 * Iterates through all the results of the edges in order to determine the betweenness for all of the edges
 * gathered from breadth first search. Used heavily in Girvan Newman Algorithm
 */
template<class T>
vector<Edge<T>*> Graph<T>::calculateBetweenness(vector<Edge<T>*> vectorOfEdges, int paths[]) {

    vector<Edge<T>*> finalEdgesWithBetweenness;

    queue<int> theQueue;
    stack<int> theStack;

    theQueue.push(vectorOfEdges[0]->getBeginningNodeName());

    while (!theQueue.empty()) {
        int tempNodeName = theQueue.front();
        theQueue.pop();

        for (int i = 0; i < vectorOfEdges.size(); i++) {
            if (vectorOfEdges[i]->getBeginningNodeName() == tempNodeName) {
                theQueue.push(vectorOfEdges[i]->getEndNode()->getNodeName());
            }
        }

        theStack.push(tempNodeName);
    }

    while (!theStack.empty()) {
        int temp = theStack.top();
        theStack.pop();

        if (!theStack.empty()) {
            while (theStack.top() == temp) {
                theStack.pop();
            }
        }


        vector<int> children;

        for (int i = 0; i < vectorOfEdges.size(); i++) {
            if (vectorOfEdges[i]->getBeginningNodeName() == temp) {
                children.push_back(vectorOfEdges[i]->getEndNode()->getNodeName());
            }
        }

        if (children.empty()) {
            for (int i = 0; i < vectorOfEdges.size(); i++) {
                if (vectorOfEdges[i]->getEndNode()->getNodeName() == temp) {
                    auto node = new Node<T>(vectorOfEdges[i]->getBeginningNodeName());
                    auto edge = new Edge<T>(node);
                    edge->setBeginningNodeName(temp);
                    double numPaths = paths[temp - 1];
                    double betweennessValue = 1 / numPaths;
                    edge->setBetweennessValue(betweennessValue);
                    edge->setNumPaths(numPaths);
                    finalEdgesWithBetweenness.push_back(edge);
                }
            }
        }


        if (!children.empty()) {
            double betweenness = 0;

            for (int j = 0; j < children.size(); j++) {

                for (int k = 0; k < finalEdgesWithBetweenness.size(); k++) {
                    if (finalEdgesWithBetweenness[k]->getBeginningNodeName() == children[j]) {
                        betweenness = betweenness + finalEdgesWithBetweenness[k]->getBetweennessValue();
                        break;
                    }

                }
            }

            for (int l = 0; l < vectorOfEdges.size(); l++) {
                if (vectorOfEdges[l]->getEndNode()->getNodeName() == temp) {
                    auto node = new Node<T>(vectorOfEdges[l]->getBeginningNodeName());
                    auto edge = new Edge<T>(node);
                    edge->setBeginningNodeName(temp);
                    double numPaths = paths[temp - 1];
                    double betweennessValue = betweenness / numPaths;
                    edge->setBetweennessValue(betweennessValue + 1);
                    edge->setNumPaths(numPaths);
                    finalEdgesWithBetweenness.push_back(edge);
                }


            }

        }
    }

    return finalEdgesWithBetweenness;
}

/*
 *  This is Targans Algorithm. It used once the necessary amount of edges have been removed. The algorithm
 *  determines all of the strongly connected communities in a graph. It is crucial to the GirvanNewman Algorithm.
 *  The Targans Algorithm uses a depth first search in order to find all of the strongly connected communities
 */
template<class T>
vector<vector<int>> Graph<T>::TargansAlgorithm() {

    int *disc = new int[adjListOfNodes.size()];
    int *low = new int[adjListOfNodes.size()];
    bool *stackMember = new bool[adjListOfNodes.size()];
    auto st = new stack<int>();
    vector<vector<int>> finalCommunities;

    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        disc[i] = -1;
        low[i] = -1;
        stackMember[i] = false;
    }

    for (int i = 0; i < adjListOfNodes.size(); i++)
    {
        if (disc[i] == -1)
        {
            vector<int> temp;
            temp = TarjansAlgorithmInner(i+1, disc, low, st, stackMember);
            finalCommunities.push_back(temp);
        }
    }

    return finalCommunities;



}

/*
 *  This is Targans Algorithm Inner. It used once the necessary amount of edges have been removed. The algorithm
 *  determines all of the strongly connected communities in a graph. It is crucial to the GirvanNewman Algorithm.
 *  The Targans Algorithm uses a depth first search in order to find all of the strongly connected communities
 */
template<class T>
vector<int> Graph<T>::TarjansAlgorithmInner(int u, int disc[], int low[], stack<int>* st, bool stackMember[]) {

    static int time = 0;

    disc[u - 1] = low[u - 1] = ++time;
    st->push(u);
    stackMember[u - 1] = true;

    Node<int>* currNode = getNode(u);

    for (int i = 0; i < currNode->getAdjSize(); i++) {
        Node<int> *currEdgeNode = getNode(currNode->getEdge(i)->getEndName());

        if (disc[currEdgeNode->getNodeName() - 1] == -1) {
            TarjansAlgorithmInner(currEdgeNode->getNodeName(), disc, low, st, stackMember);
            low[u - 1] = min(low[u - 1], low[currEdgeNode->getNodeName() - 1]);
        } else {
            if (!stackMember[currEdgeNode->getNodeName() - 1]) continue;
            low[u - 1] = min(low[u - 1], disc[currEdgeNode->getNodeName() - 1]);
        }

    }

        Node<int>* w = nullptr;

        vector<int> community;

        if (low[u - 1] == disc[u - 1])
        {
            while (st->top() != u)
            {
                w = getNode(st->top());
                cout << w->getNodeName() << " ";
                community.push_back(w->getNodeName());
                stackMember[w->getNodeName() - 1] = false;
                st->pop();
            }

            w = getNode(st->top());
            cout << w->getNodeName() << endl;
            community.push_back(w->getNodeName());
            stackMember[w->getNodeName() - 1] = false;
            st->pop();

        }

        return community;


}


#endif //GRAPHPRACTICE_GRAPH_H
