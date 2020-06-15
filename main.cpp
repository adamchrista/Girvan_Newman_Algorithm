#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include "Graph.h"
#include <algorithm>
#include <iterator>


using namespace std;

/*
 *  This is the loadFileIntoGraph function. This function takes in a file to load into a graph and then creates
 *  all of the necessary relationships between nodes and edges for the graph.
 */
void loadFileIntoGraph(Graph<int>& theGraph, map<string,int>& normalMap, map<int, string>& reversedMap, string const& fileToLoad)
{

    fstream myFile;
    myFile.open(fileToLoad);
    int count = 1;

    char unnecessary[500];
    //string str;
    myFile.getline(unnecessary, 500);

    int index = 1;
    char c;
    string theString;
    c = unnecessary[index];
    while (c != ']')
    {
        if (index == 1)
        {
            theString += c;

        }

        if (index != 1)
        {
            c = unnecessary[index];
            theString += c;
        }
        index++;
        c = unnecessary[index];
    }

    int numOfRows = stoi(theString);

    char temp[500];

    for (int i = 0; i < numOfRows; i++)
    {
        myFile.getline(temp, 500);
        string stringTemp = temp;
        normalMap[stringTemp] = count;
        reversedMap[count] = stringTemp;
        count++;
    }

    myFile.getline(temp, 500);

    while (myFile.getline(temp, 500))
    {
        string AString = temp;
        vector<string> results;
        istringstream iss(AString);
        copy(istream_iterator<string>(iss),
                  istream_iterator<string>(),
                  back_inserter(results));

        int nodeName1 = normalMap[results[0]];
        int nodeName2 = normalMap[results[2]];
        theGraph.addRelationship(nodeName1, nodeName2);

    }

    myFile.close();
}


/*
 *  This is the removeHighestBetweennessAndCreateNewGraph function. This function takes in an old graph
 *  object (with betweenness already calculated for each edge) and removes the edge with the highest betweeness.
 *  Using the edges with a removed edge, a new graph is constructed and returned.
 *
 */
Graph<int> removeHighestBetweennessAndCreateNewGraph(vector<Edge<int>*> vectorOfEdges, Graph<int> oldGraph,
        ofstream& outputFile, map<int, string> reversedMap, int& finish)
{
double highestBetweennes = 0;
int hasRemovedSingleHighestBetweenness = 0;

Graph<int> newGraph;

    for (auto & vectorOfEdge : vectorOfEdges)
    {
        if (vectorOfEdge->getBetweennessValue() > highestBetweennes)
        {

             highestBetweennes = vectorOfEdge->getBetweennessValue();
        }
    }

     cout << highestBetweennes << endl;


    for (auto & vectorOfEdge : vectorOfEdges)
    {
         if (vectorOfEdge->getBetweennessValue() == highestBetweennes && hasRemovedSingleHighestBetweenness == 0)
            {
                 hasRemovedSingleHighestBetweenness = 1;
      }       else
            {
            newGraph.addRelationship(vectorOfEdge->getBeginningNodeName(), vectorOfEdge->getEndName());
            }

         }

   if (newGraph.getNumberOfNodes() < oldGraph.getNumberOfNodes())
   {
       finish = 1;
       vector<vector<int>> finalCommunnities;
       finalCommunnities = oldGraph.TargansAlgorithm();

       for (int i = 0; i < finalCommunnities.size(); i++)
       {
           cout << "Community " << i+1 << ": " << endl;
           outputFile << "Community " << i+1 << ": " << endl;
           vector<string> tempStringVectorToSort;
           for (int j : finalCommunnities[i])
           {
               tempStringVectorToSort.push_back(reversedMap[j]);
           }

           for(int k = 0; k < tempStringVectorToSort.size(); k++)
           {
               sort(tempStringVectorToSort.begin(),tempStringVectorToSort.end());
           }

           for (int l = 0; l < finalCommunnities[i].size(); l++)
           {
               cout << tempStringVectorToSort[l] << endl;
               outputFile << tempStringVectorToSort[l] << endl;
           }

           cout << endl;
           outputFile << endl;
       }
   }



return newGraph;


}

/*
 *  This is the Girvan Newman Algorithm Function. The function recieves a graph, calculates betweenness for all
 *  edges, removes edge with highest betweenness. Once a node does not have any more edges, the algorithm uses
 *  Tarjans Strong Communities algorithm to produce the communities.
 *
 */
void GirvanNewmanAlgorithm(Graph<int> graph, int& stop, const map<int, string>& reversedMap, ofstream& outputFile)
{

    vector<vector<int>> finalCommunnities;

    vector<Edge<int>*> finalEdges;

    for (int i = 0; i < graph.getNumberOfNodes(); i++)
    {
        graph.resetIteratorNumbers();


        vector<Edge<int>*> tempEdgesFromBFS = graph.bfs(i + 1, stop);

        vector<Edge<int>*> tempEdgesFromBFS2;

        for (auto & a : tempEdgesFromBFS)
        {
            tempEdgesFromBFS2.push_back(a);
            for (int b = 0; b < tempEdgesFromBFS2.size() - 1; b++)
            {
                if (((a->getBeginningNodeName() == tempEdgesFromBFS2[b]->getBeginningNodeName()) &&
                     (a->getEndName() == tempEdgesFromBFS2[b]->getEndName())) && (a->getBetweennessValue()
                                                                                                    == tempEdgesFromBFS2[b]->getBetweennessValue()))
                {
                    tempEdgesFromBFS2.pop_back();
                }
            }
        }


        for (int j = 0; j < tempEdgesFromBFS2.size(); j++)
        {
            int found = 0;
            for (auto & finalEdge : finalEdges)
            {
                int begNodeNameFromBFS = tempEdgesFromBFS2[j]->getBeginningNodeName();
                int endNodeNameFromBFS = tempEdgesFromBFS2[j]->getEndNode()->getNodeName();
                int begNodeNameOnFinalEdge = finalEdge->getBeginningNodeName();
                int endNodeNameOnFinalEdge = finalEdge->getEndNode()->getNodeName();
                if (((begNodeNameFromBFS == begNodeNameOnFinalEdge) || (begNodeNameFromBFS == endNodeNameOnFinalEdge))
                    && ((endNodeNameFromBFS == begNodeNameOnFinalEdge) || (endNodeNameFromBFS == endNodeNameOnFinalEdge)))
                {
                    double newBetweenness = finalEdge->getBetweennessValue() + tempEdgesFromBFS2[j]->getBetweennessValue();

                    int hasItBeenSeen = 0;
                    for (int l = 0; l < j; l++)
                    {
                        int begNodeNameFromBFS2 = tempEdgesFromBFS2[l]->getBeginningNodeName();
                        int endNodeNameFromBFS2 = tempEdgesFromBFS2[l]->getEndNode()->getNodeName();
                        if (((begNodeNameFromBFS2 == begNodeNameOnFinalEdge) || (begNodeNameFromBFS2 == endNodeNameOnFinalEdge))
                            && ((endNodeNameFromBFS2 == begNodeNameOnFinalEdge) || (endNodeNameFromBFS2 == endNodeNameOnFinalEdge)))
                        {
                            hasItBeenSeen = 1;
                        }
                    }

                    if (hasItBeenSeen == 0)
                    {
                        finalEdge->setBetweennessValue(newBetweenness);
                    }

                    found = 1;
                    break;

                }

            }

            if (found == 0 && i != 0)
            {
                finalEdges.push_back(tempEdgesFromBFS2[j]);
            }

        }

        if (i == 0)
        {
            for (auto b : tempEdgesFromBFS2)
            {
                finalEdges.push_back(b);
            }
        }

    }

    if (finalEdges.size() == 1)
    {
        finalCommunnities = graph.TargansAlgorithm();
    }

    int finish = 0;

    Graph<int> newGraph = removeHighestBetweennessAndCreateNewGraph(finalEdges, graph, outputFile,
            reversedMap, finish);


    if (finish == 0)
    {
        GirvanNewmanAlgorithm(newGraph, stop, reversedMap, outputFile);
    }

}

/*
 * This is the main function for my Algorithms Project 2. It contains the logic
 *  for reading in the input control file in order to perform the required commands.
 */
int main(int argc, char* argv[]) {

    string currentFile;
    cout << "ACTUAL RESULTS ARE IN OUTPUT FILE" << endl;
    fstream inputFile;
    ofstream outputFile;
    inputFile.open(argv[1]);

    Graph<int> theGraph;
    map<string, int> normalMap;
    map<int, string> reversedMap;

    char text[500];

    while(inputFile.getline(text, 500))
    {
        string AString = text;

        if (text == "")
        {
            break;
        }

        vector<string> results;
        istringstream iss(AString);
        copy(istream_iterator<string>(iss),
                  istream_iterator<string>(),
                  back_inserter(results));


        if (results[0] == "mc")
        {
            theGraph.clear();
            normalMap.clear();
            reversedMap.clear();

            loadFileIntoGraph(theGraph, normalMap, reversedMap, currentFile);

            cout << "Shortest Path Connection: ";
            outputFile << "Shortest Path Connection: ";
            theGraph.makeAConnection(normalMap[results[1]],normalMap[results[2]], reversedMap, outputFile);
        }

        if (results[0] == "or")
        {
            normalMap.clear();
            reversedMap.clear();
            theGraph.clear();
            loadFileIntoGraph(theGraph, normalMap, reversedMap, results[1]);
            currentFile = results[1];
        }

        if (results[0] == "dc")
        {
            int stop = 0;
            GirvanNewmanAlgorithm(theGraph, stop, reversedMap, outputFile);
        }

        if (results[0] == "ow")
        {
            outputFile.open(results[1]);
        }

        if (results[0] == "dfs")
        {
            cout << "DFS RESULTS: " << endl << endl;
            outputFile << "DFS RESULTS: " << endl << endl;
            theGraph.clear();
            normalMap.clear();
            reversedMap.clear();

            cout << "{";
            outputFile << "{";

            loadFileIntoGraph(theGraph, normalMap, reversedMap, currentFile);

            vector<Edge<int>*> tempEdges = theGraph.dfs(normalMap[results[1]]);
            for (auto & tempEdge : tempEdges)
            {
                cout << "(" << reversedMap[tempEdge->getBeginningNodeName()] << " -- " << reversedMap[tempEdge->getEndNode()->getNodeName()] << "), ";
                outputFile << "(" << reversedMap[tempEdge->getBeginningNodeName()] << " -- " << reversedMap[tempEdge->getEndNode()->getNodeName()] << "), ";
            }

            cout << "}" << endl;
            outputFile << "}" << endl;
        }

        if (results[0] == "bfs")
        {
            cout << "BFS RESULTS: " << endl << endl;
            outputFile << "BFS RESULTS: " << endl << endl;
            cout << "{";
            outputFile << "{";
            theGraph.clear();
            normalMap.clear();
            reversedMap.clear();

            loadFileIntoGraph(theGraph, normalMap, reversedMap, currentFile);
            int stop = 0;
            vector<Edge<int>*> tempEdges = theGraph.bfs(normalMap[results[1]], stop);
            for (auto & tempEdge : tempEdges)
            {
                cout << "(" << reversedMap[tempEdge->getBeginningNodeName()] << " -- " << reversedMap[tempEdge->getEndNode()->getNodeName()] << "), ";
                outputFile << "(" << reversedMap[tempEdge->getBeginningNodeName()] << " -- " << reversedMap[tempEdge->getEndNode()->getNodeName()] << "), ";
            }

            cout << "}" << endl;
            outputFile << "}" << endl;

        }


    }



    return 0;
}

