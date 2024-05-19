// 21I-0523 Sheraz Sadiq
// 21I-0586 Huzaifa Tahir
// 21I-0737 Ali Hassan
// Section G

// Paralled and Distributed Computing Project

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>


#define MAX_COLS 4
#define MAX_CITIES 1000
#define INFINITY INT_MAX

int numCities = 0;
int globalK;
int * globalKPath;
int ** memoizationTable;
// Structure to store each row of data
struct Edge {
    char source[100]; 
    char target[100];
    int weight;
};

// Structure to store the Unique Cities
struct CityNode {
    int index;
    char cityName[100];
};

// Readign File From the csv file
int readNoLines(const char *filename){
    FILE *file = fopen(filename, "r");
    if (file == NULL) {             // If Fails to open File
        printf("Error opening file.\n");
        return -1;
    }

    char line[256];
    int numLines = 0;
    fgets(line, sizeof(line), file);
    while (fgets(line, sizeof(line), file)) {
        numLines++;     // Calculating no of lines of the csv File
    }

    return numLines;
}

// Function to read data from CSV file and transform it
void read_csv_and_transform(const char *filename, struct Edge *edges) {
    FILE *file = fopen(filename, "r");

    char line[256];
    int row = 0;
    while (fgets(line, sizeof(line), file)) {
        // Skip the header row
        if (row == 0) {
            row++;
            continue;
        }

        // Parse the line
        char *token = strtok(line, ",");
        int col = 0;
        while (token != NULL) {
            // Store source, target, and weight in the Edge structure
            if (col == 0) {
                strcpy(edges[row - 1].source, token);
            } else if (col == 1) {
                strcpy(edges[row - 1].target, token);
            } else if (col == 2) {
                edges[row - 1].weight = atoi(token);
            }
            token = strtok(NULL, ",");
            col++;
        }
        row++;
    }

    fclose(file);
}



// Dusplaying CSV Data
void displayCSVData(const struct Edge *edges, int size){
    // Displaying the edges in the desired format
    printf("{");
    for (int i = 0; i < size; i++) {
        printf("{ \"%s\", \"%s\", %d }, ", edges[i].source, edges[i].target, edges[i].weight);
    }
    printf("}\n");
}

// Displaying Uniquue Cities as Nodes
void displayCityNodes(const struct CityNode *cities, int size){
    // Displaying the city nodes
    printf("{\n");
    for (int i = 0; i < size; i++) {
        printf("  { \"%d\", \"%s\" },\n", cities[i].index, cities[i].cityName);
    }
    printf("}\n");
}

// Transforming the data to the Adjency Matrix
int **dataTODistanceMatrix(const struct Edge* edges, struct CityNode *cities, const int numLines) {

    // Initialize cities array
    for (int i = 0; i < MAX_CITIES; i++) {
        cities[i].index = -1;
        cities[i].cityName[0] = '\0';
    }

    // Store unique source cities
    for (int i = 0; i < numLines; i++) {
        bool cityExists = false;
        for (int j = 0; j < numCities; j++) {
            if (strcmp(edges[i].source, cities[j].cityName) == 0) {
                cityExists = true;      // Checking if the city Exists or not
                break;
            }
        }

        if (!cityExists) {      // If not exists then store it
            strcpy(cities[numCities].cityName, edges[i].source);
            cities[numCities].index = numCities;
            numCities++;
        }
    }

    // Store unique target cities
    for (int i = 0; i < numLines; i++) {
        bool cityExists = false;
        for (int j = 0; j < numCities; j++) {
            if (strcmp(edges[i].target, cities[j].cityName) == 0) {
                cityExists = true;      // Checking if the city Exists or not
                break;
            }
        }

        if (!cityExists) {      // If not exists then store it
            strcpy(cities[numCities].cityName, edges[i].target);
            cities[numCities].index = numCities;
            numCities++;
        }
    }

    // Dynamically allocate memory for distance matrix
    int **distCities = (int **)malloc(numCities * sizeof(int *));
    for (int i = 0; i < numCities; i++) {
        distCities[i] = (int *)malloc(numCities * sizeof(int));
    }

    // Initialize distance matrix
    for (int i = 0; i < numCities; i++) {
        for (int j = 0; j < numCities; j++) {
            distCities[i][j] = INFINITY; // Initialize with default value (0 on Diagonal and Other infinity)
            if(i == j)
                distCities[i][j] = 0;
        }
    }

    // Fill distance matrix
    for (int i = 0; i < numCities; i++) {
        for (int j = 0; j < numLines; j++) {
            if (strcmp(cities[i].cityName, edges[j].source) == 0) {
                int targetIndex = -1; // Initialize target index
                for (int k = 0; k < numCities; k++) {
                    if (strcmp(cities[k].cityName, edges[j].target) == 0) {
                        targetIndex = cities[k].index;
                        break;
                    }
                }
                if (targetIndex != -1) {
                    distCities[cities[i].index][targetIndex] = edges[j].weight;     // Storing the weight
                    distCities[targetIndex][cities[i].index] = edges[j].weight;
                }
            }
        }
    }
    return distCities;
}




// Calculating the Shortest Path from one node to other
int shortestPath(int src, int dest, int **adj, int V) {
    int dist[numCities]; // Initialize distances to infinity
    int parent[numCities];
    bool visited[numCities];
    for (int i = 0; i < V; ++i) {
        visited[i] = false;
        parent[i] = -1; // Initialize parent array with -1 (indicating no parent)
        dist[i] = INFINITY; // Initialize distances to infinity
    }
    int queue[numCities];
    int front = 0, rear = 0;

    dist[src] = 0; // Distance from source to source is 0
    visited[src] = true;
    queue[rear++] = src;

    while (front != rear) {
        int u = queue[front++];

        // Check all adjacent nodes of u
        #pragma omp parallel for
        for (int v = 0; v < V; ++v) {
            if (adj[u][v] && adj[u][v] != INFINITY && adj[u][v]!=0 ) {
                int newDist = dist[u] + adj[u][v];
                if (newDist < dist[v]) { // Update distance only if new distance is shorter
                    dist[v] = newDist;
                    parent[v] = u; // Set parent of v to u
                    if (!visited[v]) {
                        visited[v] = true;
                        queue[rear++] = v; // Enqueue v only if it hasn't been visited yet
                    }
                }
            }
        }


    }

    /*
    // Print the shortest path from source to destination
    printf("Shortest path from node %d to node %d: ", src, dest);
    int u = dest;
    while (u != -1) {
        printf("%d ", u);
        u = parent[u]; // Move to the parent node
       
    }
    printf("\n");

    */


    return dist[dest];
}



// Initialize memoization table with INFINITY
void initializeMemoizationTable() {
    #pragma omp parallel for
    for (int i = 0; i < numCities; i++) {
        #pragma omp parallel for
        for (int j = 0; j < numCities; j++) {
            memoizationTable[i][j] = INFINITY;
        }
    }
}

// Function to get shortest path length using memoization
int shortestPathMemoized(int src, int dest, int **adj, int V) {
    #pragma omp task
    {
        if (memoizationTable[src][dest] != INFINITY) {
            return memoizationTable[src][dest]; // If already computed, return memoized value
        }
    }
    #pragma omp task
    {
        if (memoizationTable[src][dest] == INFINITY) {
            int pathLength = shortestPath(src, dest, adj, V); // Compute shortest path
            memoizationTable[src][dest] = pathLength; // Memoize the result
             memoizationTable[dest][src] = pathLength; // Memoize the result
            return pathLength;
        }
    }
}



// Calculating the kth shortest paths from one node to other
/*
    For Example Source = 0, Destination = 23
    Creating Array of K length initializing with INFINITY
    Calculating as
    {0 -> Directly Connected Node} + ..... + {Connected Node -> destination}

*/
void findKthPath(int **adjacencyMatrix, int numNodes, int source, int destination, int k) {
   
    int *dist = (int *)malloc(numNodes * sizeof(int));
    int *prev = (int *)malloc(numNodes * sizeof(int));

    int mainSource = source;
    #pragma omp parallel for
    for (int j = 0; j < numCities; j++) {
        if (adjacencyMatrix[source][j] != INFINITY && j != destination ) { 
            int local_source = j;
            #pragma omp parallel for
            for (int i = 0; i < numCities; i++) {
                int cost0 = 0;
                int cost1 = 0;
                int cost2 = 0;
                int totalCost = 0;

                if (adjacencyMatrix[local_source][i] != 0 && adjacencyMatrix[local_source][i] != INFINITY && i != mainSource) {
                    int temp = destination;
                    destination = i;

                    // Compute the cost of the first segment of the path
                    #pragma omp task
                    {
                        cost1 = shortestPathMemoized(local_source, destination, adjacencyMatrix, numCities);
                    }

                    destination = temp;
                    temp = local_source;

                    // Compute the cost of the second segment of the path
                    #pragma omp task
                    {
                        if (local_source != mainSource) {
                            int temp_destination = destination;
                            destination = local_source;
                            int temp_source = local_source;
                            local_source = mainSource;
                            cost0 = shortestPathMemoized(local_source, destination, adjacencyMatrix, numCities);
                            destination = temp_destination;
                            local_source = temp_source;
                        }
                    }

                    local_source = i;

                    // Compute the cost of the third segment of the path
                    #pragma omp task
                    {
                        if (i != destination) {
                            cost2 = shortestPathMemoized(local_source, destination, adjacencyMatrix, numCities);
                        }
                    }

                    // Wait for all tasks to complete
                    #pragma omp taskwait

                    // Compute the total cost
                    totalCost = cost0 + cost1 + cost2;

                    // Update the global Kth shortest paths
                    #pragma omp critical
                    {
                        int maxIdx = 0;
                        for (int l = 0; l < k; l++) {
                            if (globalKPath[l] > globalKPath[maxIdx]) {
                                maxIdx = l;
                            }
                        }

                        if (globalKPath[maxIdx] > totalCost) {
                            globalKPath[maxIdx] = totalCost;
                        }
                    }
                }
            }
        }
    }

    free(dist);
    free(prev);
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(NULL) + rank); 

    int numLines = readNoLines("doctorwho.csv");
    if (numLines == -1) {
        MPI_Finalize();
        return 1; // Error reading file
    }

    struct Edge edges[numLines];
    read_csv_and_transform("doctorwho.csv", edges);
    int k = 3; // Find K shortest paths

    struct CityNode cities[MAX_CITIES];
    int **distanceMatrix = NULL;


    if (rank == 0) {    // Root process initializes the distance matrix and memoization table
        distanceMatrix = dataTODistanceMatrix(edges, cities, numLines);
        memoizationTable = (int **)malloc(numCities * sizeof(int *));
        for (int i = 0; i < numCities; i++) {
            memoizationTable[i] = (int *)malloc(numCities * sizeof(int));
        }
        initializeMemoizationTable(); // Initialize memoization table
    }

    // Broadcast distance matrix and memoization table from root process to all other processes
    MPI_Bcast(&numCities, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        distanceMatrix = (int **)malloc(numCities * sizeof(int *));
        for (int i = 0; i < numCities; i++) {
            distanceMatrix[i] = (int *)malloc(numCities * sizeof(int));
        }
        memoizationTable = (int **)malloc(numCities * sizeof(int *));
        for (int i = 0; i < numCities; i++) {
            memoizationTable[i] = (int *)malloc(numCities * sizeof(int));
        }
    }
    for (int i = 0; i < numCities; i++) {
        MPI_Bcast(distanceMatrix[i], numCities, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(memoizationTable[i], numCities, MPI_INT, 0, MPI_COMM_WORLD);
    }

    globalK = k;
    globalKPath = (int *)malloc(k * sizeof(int));


    for(int i = 0; i < k; i++){
        globalKPath[i] = INFINITY;
    }
    
    clock_t start_time, end_time;
    double execution_time;

    start_time = clock();

    // Determine the number of iterations each process will handle
    int numIterationsPerProcess = 5 / numProcesses;
    int remainder = 5 % numProcesses;

    // Calculate the starting and ending indices for the loop for this process
    int startIndex = rank * numIterationsPerProcess + (rank < remainder ? rank : remainder);
    int endIndex = startIndex + numIterationsPerProcess + (rank < remainder ? 1 : 0);

    int localSource = (rand() % numCities) + 1; // Generate random source for this iteration
    int localDestination = (rand() % numCities) + 1; // Generate random destination for this iteration
    while (localDestination == localSource) {
        localDestination = (rand() % numCities) + 1; // Ensure destination is different from source
    }
    
    printf("{%d ,%d}\n",localSource ,localDestination);
    // Compute Kth path for this source and destination pair
    findKthPath(distanceMatrix, numCities, localSource, localDestination, k);

    
    printf("Shortest K(%d) Paths: ", k);
    for(int j = 0; j < k; j++){
        printf("%d ", globalKPath[j]);
    }

    printf("\n");

    end_time = clock();

    execution_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;       // CalCulating time of Execution
    printf("Execution time: %f seconds\n\n", execution_time);

    MPI_Finalize();

    return 0;
}
