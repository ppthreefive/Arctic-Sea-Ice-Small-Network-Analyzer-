#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int landCells = 0;

struct node 
{
	int vertex;
	struct node * next;
};

class adjacencyList 
{
	public:

		struct node * head;
		int size;

		adjacencyList()
		{
			head = NULL;
			size = 0;
		}

		void insertBack(int data)
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
};

class Graph 
{
	public:
		
		int numVertices;
		int time;
		adjacencyList *list;
	
		Graph(int vertices)
		{
			numVertices = vertices;
			list = new adjacencyList[vertices];
		}

		void addEdge(int u, int v) 
		{
			list[u].insertBack(v);
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
				if (visited[v] == false) 
				{
					depthFirstSearch(v, visited);
					count++;
				}
			}

			return count;
		}

		void depthFirstSearch(int v, bool visited[]) 
		{
			visited[v] = true;

			node * p = list[v].head;

			while (p != NULL) 
			{
				v = p->vertex;

				if (!visited[v]) 
				{
					depthFirstSearch(v, visited);
				}

				p = p->next;
			}
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

			chains[0] -= landCells;

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

/* Reads all the cells in the binary file given, then puts the data into the correct place in our table. Each cell in our table is a doubly linked list. */
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

void calculateMean(vector<float> &means, vector<vector<vector<float>>> &array)
{
	for (int i = 0; i < 63; i++)
	{
		for (int j = 0; j < 63; j++)
		{
			float addition = 0.0f;
			float result = 0.0f;
			int count = 0;

			if (array[i][j][0] != 168.0f) 
			{
				for (int k = 0; k < 832; k++)
				{
					if (array[i][j][k] != 157.0f)
					{
						addition += array[i][j][k];
						count++;
					}
				}

				result = addition / (float)count;
			}
			else 
			{
				landCells++;
			}

			means.at((i * 63) + j) = result;
		}
	}
}

void calculateTimeSeries(vector<float> &timeSeries, vector<float> &means, vector<vector<vector<float>>> &array)
{
	for (int i = 0; i < 63; i++)
	{
		for (int j = 0; j < 63; j++)
		{
			float result = 0.0f;
			float notSq = 0.0f;

			if (array[i][j][0] != 168.0f) 
			{
				for (int k = 0; k < 832; k++)
				{
					if (array[i][j][k] != 157.0f)
					{

						float num = array[i][j][k];
						float mean = means[(i * 63) + j];

						float numTwo = array[i][j][k] - mean;

						num = num - mean;
						num = num * num;

						result += num;
					}
				}
			}

			timeSeries.at((i * 63) + j) = result;
		}
	}
}

void calculateCoefficientsAndGraph(vector<float> &means, vector<vector<vector<float>>> &array, vector<float> &timeSeries, float (&thresh)[3], Graph (&graphs)[3])
{
	int edges1 = 0;
	int edges2 = 0;
	int edges3 = 0;

	for (int i = 0; i < 3969; i++)
	{
		if (means[i] != 0.0f) 
		{
			for (int j = 0; j < 3969; j++)
			{
				float coef = 0.0f;

				if (j != i)
				{
					int index_X = 0;
					int index_Y = 0;

					int index_X2 = 0;
					int index_Y2 = 0;

					float Sxy = 0.0f;

					if (i > 62)
					{
						index_X = i / 63;
						index_Y = i % 63;
					}
					else
					{
						index_X = 0;
						index_Y = i;
					}

					if (j > 62)
					{
						index_X2 = j / 63;
						index_Y2 = j % 63;
					}
					else
					{
						index_X2 = 0;
						index_Y2 = j;
					}

					if (means[i] != 0.0f) 
					{
						if (means[j] != 0.0f)
						{
							for (int k = 0; k < 832; k++)
							{
								Sxy += (array[index_X][index_Y][k] - means[i]) * (array[index_X2][index_Y2][k] - means[j]);
							}

							coef = fabs(Sxy / sqrtf((timeSeries.at(i)*timeSeries.at(j))));
						}
					}

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
	}
}

int main() 
{
	vector<vector<vector<float>>> array;
	vector<float> means;
	means.resize(3969);

	vector<float> timeSeries;
	timeSeries.resize(3969);

	vector<float> notSquared;
	notSquared.resize(3969);

	int height = 63;
	int width = 63;
	//int depth = 832;

	// Resizes our 3D vector array properly
	array.resize(height);
	for (int i = 0; i < height; ++i)
	{
		array[i].resize(width);
	}

	// Thresholds
	float thresh[3] = { .95f, .925f, .9f };

	// Initilize the year numbers
	int yearNumbers[16] = {1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005};

	cout << "Currently reading binary files..." << endl;
	int numFiles = 0;

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

	cout << "Files read: " << numFiles << endl;

	// Initialize the graphs
	Graph graphs[3] = { Graph(3969), Graph(3969), Graph(3969) };

	cout << "Calculating means of all indexes..." << endl;

	// After reading all of the files, calculate the mean for each spot in the table
	calculateMean(means, array);

	cout << "Calculating time series of all indexes..." << endl;
	// After reading all of the files, calculate the Sxx value for each spot in the table
	calculateTimeSeries(timeSeries, means, array);

	cout << "Calculating coefficients of all indexes..." << endl;
	// Now we calculate the coefficients
	calculateCoefficientsAndGraph(means, array, timeSeries, thresh, graphs);

	cout << endl;

	cout << "----------Threshold 0.950 Histogram----------" << endl;
	cout << graphs[0].print() << endl;

	cout << "Summary of Connected Components for threshold 0.950" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << "Number of Connected Nodes: " << to_string(graphs[0].getComponentNumbers()) << endl << endl;

	cout << "----------Threshold 0.925 Histogram----------" << endl;
	cout << graphs[1].print() << endl;

	cout << "Summary of Connected Components for threshold 0.925" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << "Number of Connected Nodes: " << to_string(graphs[1].getComponentNumbers()) << endl << endl;

	cout << "----------Threshold 0.900 Histogram----------" << endl;
	cout << graphs[2].print() << endl;

	cout << "Summary of Connected Components for threshold 0.900" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << "Number of Connected Nodes: " << to_string(graphs[2].getComponentNumbers()) << endl << endl;

	return 0;
}