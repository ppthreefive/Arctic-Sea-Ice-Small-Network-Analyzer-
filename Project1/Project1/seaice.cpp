/*	Name: Phillip Pham
	Course: CSE310, Section: 84794
	Instructor: Dr. Violet R. Syrotiuk

	Program Title: seaice
	Program Description: This program is designed to represent climatological data of the beaufort sea region as a graph, 
		and then we compute characteristics of it to determine if it is a small-world graph. This program specifically
		uses the data set gathered from the SMMR-SSM/I passive microwave data set and finds anomalies, so that it is
		more amenable to statistical analysis.

	The binary files used for this data calculation was provided by the instructor on Canvas.

	Extra resources used for building this application:
		http://www.cplusplus.com/forum/general/833/
		https://www.geeksforgeeks.org/connected-components-in-an-undirected-graph/
		https://www.softwaretestinghelp.com/graph-implementation-cpp/
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>

using namespace std;

vector<float> avgDegrees;

// Node structure for the adjacency lists
struct node 
{
	int vertex;
	struct node * next;
};

/* This is the adjacency linked lists that are going to be stored at every spot in the graph */
class adjacencyList 
{
	public:

		struct node * head, *tail;
		int size;

		adjacencyList()
		{
			head = NULL;
			tail = NULL;
			size = 0;
		}

		void push_back(int data) 
		{
			node *temp = new node;
			temp->vertex = data;
			temp->next = NULL;

			if (head == NULL)
			{
				head = temp;
				tail = temp;
				temp = NULL;
			}
			else
			{
				tail->next = temp;
				tail = temp;
			} 
			
			size++;
		}

		void push_front(int data)
		{
			node *temp = new node;

			temp->vertex = data;

			if (head != NULL)
			{
				temp->next = head;
				head = temp;
			}
			else
			{
				head = temp;
				head->next = NULL;
			}

			size++;
		}

		void pop() 
		{
			if (head != NULL) 
			{
				node * temp = new node;
				temp = head;
				head = head->next;
				delete temp;
				size--;
			}
		}

		int top() 
		{
			if (head != NULL) 
			{
				return head->vertex;
			}
			else 
			{
				return -100;
			}
		}

		bool empty() 
		{
			if (size > 0) 
			{
				return false;
			}
			else 
			{
				return true;
			}
		}
};

/* This class stores all of our edges and vertices, and has functions to create a degree distribution, and to figure out how many connected components there are. */
class Graph 
{
	public:
		
		int numVertices;
		adjacencyList *list;
		int ** adjMatrix;
	
		Graph(int vertices)
		{
			// Initialize our adjacency lists
			numVertices = vertices;
			list = new adjacencyList[vertices];

			// Initialize our adjacency matrix as well
			adjMatrix = new int*[vertices];
			for (int i = 0; i < vertices; i++)
			{
				adjMatrix[i] = new int[numVertices];
				adjMatrix[i][i] = 0;

				for (int j = 0; j < numVertices; j++)
				{
					if (j != i) 
					{
						adjMatrix[i][j] = INT_MAX;
					}
				}				
			}
		}

		void addEdge(int u, int v) 
		{
			list[u].push_back(v);
			adjMatrix[u][v] = 1;
			adjMatrix[v][u] = 1;
		}

		bool isEdge(int u, int v) 
		{
			return (adjMatrix[u][v] == 1);
		}

		int getIndexDegree(int index)
		{
			node *p = list[index].head;

			if (p != NULL && p->next == NULL)
			{
				return 1;
			}
			else if (p == NULL)
			{
				return 0;
			}
			else
			{
				int count = 1;

				while (p->next != NULL)
				{
					count++;
					p = p->next;
				}

				return count;
			}
		}

		int getMaxChainCount()
		{
			int maxChainSize = 0;
			int index = 0;

			for (int i = 0; i < numVertices; i++)
			{
				int currentMax = getIndexDegree(i);

				if (currentMax > maxChainSize)
				{
					maxChainSize = currentMax;
				}
			}

			return maxChainSize;
		}

		int getComponentNumbers()
		{
			int count = 0;

			bool * visited = new bool[numVertices];

			for (int v = 0; v < numVertices; v++)
			{
				visited[v] = false;
			}

			for (int v = 0; v < numVertices; v++)
			{
				int compCount = 0;

				if (visited[v] == false)
				{
					depthFirstSearch(v, visited, compCount);
					count++;
					cout << "Component " << to_string(count) << " size: " << to_string(compCount) << endl;
				}				
			}

			return count;
		}

		void depthFirstSearch(int v, bool visited[], int &count)
		{
			visited[v] = true;
			count++;

			node * p = list[v].head;

			while (p != NULL)
			{
				v = p->vertex;

				if (!visited[v])
				{
					depthFirstSearch(v, visited, count);
				}

				p = p->next;
			}
		}

		// Computes the characteristic path length (average path length) of a graph
		float computeCharacteristicPathLength() 
		{
			float result = 0.0f;
			float count = 0.0f;
			float sum = 0.0f;

			floydWarshall();

			for (int i = 0; i < numVertices; i++)
			{
				for (int j = 0; j < numVertices; j++)
				{
					if (i != j)
					{
						if (adjMatrix[i][j] != INT_MAX)
						{
							sum += (float)adjMatrix[i][j];
							count++;
						}
					}
				}
			}

			if (count > 0) 
			{
				result = sum / count;
			}

			return result;
		}

		// Shortened version of Floyd-Warshall algorithm
		void floydWarshall() 
		{
			for (int k = 0; k < 3186; k++)
			{
				for (int i = 0; i < 3186; i++) 
				{
					for (int j = i + 1; j < 3186; j++)
					{
						// IMPORTANT TO CHECK THAT WE ARE NOT ADDING ANYTHING TO INT_MAX, OR IT WILL OVERFLOW INTO NEGATIVE INTS
						if ((!(adjMatrix[i][k] == INT_MAX) && !(adjMatrix[k][j] == INT_MAX))
							&& (!(adjMatrix[k][i] == INT_MAX) && !(adjMatrix[j][k] == INT_MAX)))
						{
							if (adjMatrix[i][j] > adjMatrix[i][k] + adjMatrix[k][j])
							{
								adjMatrix[i][j] = adjMatrix[i][k] + adjMatrix[k][j];
							}

							if (adjMatrix[j][i] > adjMatrix[k][i] + adjMatrix[j][k])
							{
								adjMatrix[j][i] = adjMatrix[k][i] + adjMatrix[j][k];
							}
						}						
					}
				}
			}
		}

		float computeMeanClusteringCoefficient()
		{
			// Initializing local variables
			float sum = 0.0;
			float result = 0.0;

			for (int i = 0; i < numVertices; i++)
			{
				// Number of vertices in the neighborhood
				float k = (float)list[i].size;

				// Number of edges initialized, we havent counted yet
				float e = 0.0;

				// Set a pointer to the beginning of the i vertice's adjacency list
				node * p = list[i].head;

				if (k > 0)
				{
					// Check if the individual vertexes in vertex i's adjacency list have edges between them
					// If they do, increment e (We are counting edges in the neighborhood of vertex i)
					while (p != NULL)
					{
						int vertex = p->vertex;

						node * q = list[vertex].head;

						while (q != NULL)
						{
							node * r = list[i].head;
							int vertexTwo = q->vertex;
							int vertexThree = r->vertex;

							while (r != NULL)
							{
								if (vertexTwo != vertexThree)
								{
									r = r->next;
								}
								else
								{
									e++;

									r = r->next;
								}
							}

							q = q->next;
						}

						p = p->next;
					}

					if (e > 0)
					{
						// Clustering Coefficient formula
						sum = sum + ((2 * e) / (k * (k - 1)));
					}
				}
			}

			// Get the mean of clustering coefficients
			result = (sum / (float)numVertices);

			return result;
		}

		void computeAverageVertexDegree(int maxDegree, int * degrees) 
		{
			float sum = 0.0f;
			float count = 0.0f;
			float result = 0.0f;

			for (int i = 0; i < maxDegree; i++)
			{
				for (int j = 0; j < degrees[i]; j++)
				{
					sum += i;
					count++;
				}
			}

			result = sum / count;

			avgDegrees.push_back(result);
		}

		string print() 
		{
			string output;
			int maxChain = getMaxChainCount() + 1;
			int *chains = new int[maxChain];
			int chainsNotZero = 0;

			// Initialize all elements in the array to zero first
			for (int i = 0; i < maxChain; i++)
			{
				chains[i] = 0;
			}

			int n = 0;
			int index = 0;

			while (n <= maxChain)
			{
				index = 0;

				while (index < numVertices)
				{
					int count = getIndexDegree(index);

					if (count == n)
					{
						chains[n] = chains[n] + 1;
					}

					index++;
				}

				n++;
			}

			computeAverageVertexDegree(maxChain, chains);

			// Print the number of chains for the respective number of chain up to the max chain
			for (int i = 0; i < maxChain; i++)
			{
				output += "Degree of " + to_string(i) + "\t|";

				for (int j = 0; j < (chains[i] / 10); j++)
				{
					output += "*";
				}

				output += " (" + to_string(chains[i]) + ")\n";
			}

			return output;
		}
};

/* Reads all the cells in the binary file given, then puts the data into the correct place in our table. Each cell in our table is a vector. */
void readBinary(string fileName, vector<vector<vector<float>>> &array)
{
	ifstream inputFile(fileName, ios::in | ios::binary);

	float dataIn = 0;

	for (int i = 0; i < 63; i++)
	{
		for (int j = 0; j < 63; j++)
		{
			inputFile.read((char*)&dataIn, 4); // Read 4 bytes of data

			array[i][j].push_back(dataIn);
		}
	}

	inputFile.close();
}

/* This function calculates all of the means for every (x,y) cell and stores it into a 1D vector */
void calculateMean(vector<float> &means, vector<vector<float>> &array)
{
	for (int i = 0; i < 3186; i++)
	{
		float addition = 0.0f;
		float result = 0.0f;
		int count = 0;

		for (int j = 0; j < array[i].size(); j++)
		{
			addition += array[i][j];
			count++;
		}

		// Calculated mean of the cell
		result = addition / (float)count;

		//cout << "The mean for cell " << i << ": " << to_string(result) << endl;

		means.at(i) = result;
	}
}

/* This function calculates the Sxx values for every (x,y) cell and stores it into a 1D vector */
void calculateTimeSeries(vector<float> &timeSeries, vector<float> &means, vector<vector<float>> &array)
{
	for (int i = 0; i < 3186; i++)
	{
		float result = 0.0f;

		for (int j = 0; j < array[i].size(); j++)
		{
			float num = array[i][j];
			float mean = means[i];

			float numTwo = array[i][j] - mean;

			num = num - mean;
			num = num * num;

			result += num;
		}

		timeSeries.at(i) = result;
	}
}

/* This function simultaneously calculates all of the pearson coefficients and builds the graphs */
void calculateCoefficientsAndGraph(vector<float> &means, vector<vector<float>> &array, vector<float> &timeSeries, float (&thresh)[3], Graph (&graphs)[3])
{
	// In case we want to count how many edges we made, I initialized these.
	int edges1 = 0;
	int edges2 = 0;
	int edges3 = 0;

	for (int i = 0; i < 3186; i++)
	{
		for (int j = 0; j < 3186; j++)
		{
			float coef = 0.0f;

			if (j != i) 
			{
				float Sxy = 0.0f;

				// Here we calculate our Sxy value, and then the following pearson correlation coefficient
				for (int k = 0; k < 832; k++)
				{
					Sxy += (array[i][k] - means[i]) * (array[j][k] - means[j]);
				}

				coef = fabs(Sxy / sqrtf((timeSeries.at(i)*timeSeries.at(j))));

				//cout << "The coefficient for " << i << " and " << j << ": " << to_string(coef) << endl;

				// We check where we need to place an edge. We can choose up to 3 different graphs in the same loop.
				if (coef >= thresh[0])
				{
					graphs[0].addEdge(i, j);
					edges1++;
				}
				if (coef >= thresh[1])
				{
					graphs[1].addEdge(i, j);
					edges2++;
				}
				if (coef >= thresh[2])
				{
					graphs[2].addEdge(i, j);
					edges3++;
				}
			}
		}
	}

	/*cout << "Number of edges for graph 1: " << to_string(edges1) << endl;
	cout << "Number of edges for graph 2: " << to_string(edges2) << endl;
	cout << "Number of edges for graph 3: " << to_string(edges3) << endl;*/
}

void cullLandCells(vector<vector<vector<float>>> &array, vector<vector<float>> &newArray)
{
	for (int i = 0; i < 63; i++)
	{
		for (int j = 0; j < 63; j++) 
		{
			if (array[i][j][0] != 168.0f) 
			{
				newArray.push_back(array[i][j]);
			}
		}
	}

	cout << "The number of non-land cells after culling: " << newArray.size() << endl;
}

float computeRandomClusteringCoefficient(float degreeAvg) 
{
	float result = 0.0f;

	result = degreeAvg / 3186.0f;

	return result;
}

float computeRandomCharacteristicPathLength(float degreeAvg)
{
	float result = 0.0f;

	result = logf(3186.0f) / logf(degreeAvg);

	return result;
}

int main() 
{
	// Initialize a 3D array structure, in this case we are using vector for dynamic allocation
	vector<vector<vector<float>>> array; // This will store a 3D array of 3969 values before we cull it to 3186

	// This is a 2D array or table containing only our sea/ice cells
	vector<vector<float>> nonLand;

	// Initialize 1D array to store the means of all (x,y) cells
	vector<float> means;
	//means.resize(3969);
	means.resize(3186); // Resize it to the amount of sea/ice values only

	// Initialize 1D array to store the Sxx values of all (x,y) cells
	vector<float> timeSeries;
	//timeSeries.resize(3969);
	timeSeries.resize(3186); // Resize it to the amount of sea/ice values only

	// Resizes our 3D vector array properly
	int height = 63;
	int width = 63;
	//int depth = 832;	
	array.resize(height);
	for (int i = 0; i < height; ++i)
	{
		array[i].resize(width);
	}

	// Thresholds that we use for our pearson correlation coefficient calculation
	float thresh[3] = { .95f, .925f, .9f };

	// Initilize the year numbers
	int yearNumbers[16] = {1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005};

	cout << "Currently reading binary files..." << endl;
	int numFiles = 0;

	// Calculates all of the filenames that we are including in our calculations and then reads the files
	for (int i = 0; i < 16; i++)
	{
		for (int j = 1; j < 53; j++) 
		{
			string fileName = "";
			string weekNumber = to_string(j);
			string yearNumber = to_string(yearNumbers[i]);

			if (weekNumber.length() == 1) 
			{
				weekNumber.insert(0, "0");
			}

			fileName = "CS310_project_subregion/" + yearNumber + "/Beaufort_Sea_diffw" + weekNumber + "y" + yearNumber + "+landmask";

			readBinary(fileName, array);

			numFiles++;
		}
	}

	cout << "Files read: " << numFiles << ". Done."<< endl;

	// Here we will cull our 3D array into a 2D array for faster calculation, also removing land cells
	cullLandCells(array, nonLand);
	
	// Initialize the graphs
	Graph graphs[3] = { Graph(3186), Graph(3186), Graph(3186) };

	cout << "Calculating means of all indexes...";

	// After reading all of the files, calculate the mean for each spot in the table
	calculateMean(means, nonLand);

	cout << " Done.\nCalculating time series of all indexes...";
	// After reading all of the files, calculate the Sxx value for each spot in the table
	calculateTimeSeries(timeSeries, means, nonLand);

	cout << " Done.\nCalculating coefficients of all indexes...";
	// Now we calculate the coefficients
	calculateCoefficientsAndGraph(means, nonLand, timeSeries, thresh, graphs);

	cout << " Done.\nCreating graphs..." << endl;

	cout << endl;

	// Print out graph data here
	cout << "----------Threshold 0.950 Histogram----------" << endl;
	cout << graphs[0].print() << endl;

	cout << "Summary of Connected Components for threshold 0.950" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << "Number of connected components: " << to_string(graphs[0].getComponentNumbers()) << endl << endl;

	cout << "----------Threshold 0.925 Histogram----------" << endl;
	cout << graphs[1].print() << endl;

	cout << "Summary of Connected Components for threshold 0.925" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << "Number of connected components: " << to_string(graphs[1].getComponentNumbers()) << endl << endl;

	cout << "----------Threshold 0.900 Histogram----------" << endl;
	cout << graphs[2].print() << endl;

	cout << "Summary of Connected Components for threshold 0.900" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << "Number of connected components: " << to_string(graphs[2].getComponentNumbers()) << endl << endl;

	cout << "-----------------------------------------------------------------" << endl;

	cout << "The mean of clustering coefficients for Graph of Threshold 0.950: " << to_string(graphs[0].computeMeanClusteringCoefficient()) << endl;
	cout << "The mean of clustering coefficients for Graph of Threshold 0.925: " << to_string(graphs[1].computeMeanClusteringCoefficient()) << endl;
	cout << "The mean of clustering coefficients for Graph of Threshold 0.900: " << to_string(graphs[2].computeMeanClusteringCoefficient()) << endl;

	cout << "-----------------------------------------------------------------" << endl;

	cout << "The characteristic path length for Graph of Threshold 0.950: " << to_string(graphs[0].computeCharacteristicPathLength()) << endl;
	cout << "The characteristic path length for Graph of Threshold 0.925: " << to_string(graphs[1].computeCharacteristicPathLength()) << endl;
	cout << "The characteristic path length for Graph of Threshold 0.900: " << to_string(graphs[2].computeCharacteristicPathLength()) << endl;

	cout << "-----------------------------------------------------------------" << endl;

	cout << "Random Graph 1 clustering coefficient: " << to_string(computeRandomClusteringCoefficient(avgDegrees[0])) << endl;
	cout << "Random Graph 1 clustering coefficient: " << to_string(computeRandomClusteringCoefficient(avgDegrees[1])) << endl;
	cout << "Random Graph 1 clustering coefficient: " << to_string(computeRandomClusteringCoefficient(avgDegrees[2])) << endl;

	cout << "-----------------------------------------------------------------" << endl;

	cout << "Random Graph 1 characteristic path length: " << to_string(computeRandomCharacteristicPathLength(avgDegrees[0])) << endl;
	cout << "Random Graph 1 characteristic path length: " << to_string(computeRandomCharacteristicPathLength(avgDegrees[1])) << endl;
	cout << "Random Graph 1 characteristic path length: " << to_string(computeRandomCharacteristicPathLength(avgDegrees[2])) << endl;
	
	return 0;
}