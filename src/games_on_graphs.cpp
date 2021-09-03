#include <iostream>
#include <csignal>
#include <igraph/igraph.h>
#include <sstream>
#include <string>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/filereadstream.h"
#include <fstream>
#include <stdio.h>
#include <string.h>
#include "simulation.h"


using namespace std;
using namespace rapidjson;



void printAdjacencyMatrix (igraph_t *graph)
{
	igraph_matrix_t A;
	igraph_matrix_init(&A, 0, 0);
	igraph_get_adjacency(graph, &A, IGRAPH_GET_ADJACENCY_BOTH, false);
	int N = igraph_vcount(graph);

	std::cout << "\n";
	for (int i=0; i<N; i++)
	{
		std::cout << "| ";
		for (int j=0; j<N; j++)
		{
			std::cout << MATRIX(A, i, j) << " ";
		}
		std::cout << "|\n";
	}
	std::cout << "\n";
}


int main(int argc, char* argv[])
{

	Document document;
	ostringstream strCout;
	streambuf* oldCoutStreamBuf = cout.rdbuf();

	if (argc == 1) // Configuration via stdin
	{
		stringstream config;
		for (string line; getline(std::cin, line);)
		{
			config << line << endl;
		}

		cout.rdbuf(strCout.rdbuf());				// Redirect stdout to string buffer (avoid unwanted prints in response)
		document.Parse(config.str().c_str());		// Parse JSON configuration
	}
	else if (argc == 3 && strcmp(argv[1], "-c") == 0)	// Configuration file
	{
		char *configFile = argv[2];
		if (FILE *file = fopen(configFile, "r"))
		{
			fclose(file);
			ifstream ifs(configFile);
			string content( (istreambuf_iterator<char>(ifs) ),(istreambuf_iterator<char>()));
			document.Parse(content.c_str());
		}
		else
		{
			cerr << "Configuration file '" << configFile << "' not found.\n";
			return 1;
		}
	}
	else
	{
		cerr << "Invalid command line parameters supplied." << "\n";
		exit(1);
	}


	if (!document.IsObject())
	{
		cout << "Failed to load configuration\n";
		exit(1);
	}

	bool errorOccured = false;

	// Validate parameters
	if (!document.HasMember("NETWORK") || !document["NETWORK"].IsString()) errorOccured = true;
	if (!document.HasMember("STRATEGIES") || !document["STRATEGIES"].IsString()) errorOccured = true;
	if (!document.HasMember("STRATEGIES_ENFORCED") || !document["STRATEGIES_ENFORCED"].IsString()) errorOccured = true;
	if (!document.HasMember("B") || !document["B"].IsDouble()) errorOccured = true;
	if (!document.HasMember("T") || !document["T"].IsInt()) errorOccured = true;
	if (!document.HasMember("SAMPLE_TIME") || !document["SAMPLE_TIME"].IsInt()) errorOccured = true;
	if (!document.HasMember("RANDOM_SEED") || !document["RANDOM_SEED"].IsInt()) errorOccured = true;

	if (errorOccured)
	{
		cerr << "Invalid parameters given" << "\n";
		exit(1);
	}

	unsigned T, sampleTime;
	int randomSeed;
	float b;
	string network;
	string strategies;
	string strategiesEnforced;

	T 			= document["T"].GetInt();
	b 			= document["B"].GetFloat();
	strategies 	= document["STRATEGIES"].GetString();
	sampleTime 	= document["SAMPLE_TIME"].GetInt();
	randomSeed 	= document["RANDOM_SEED"].GetInt();
	network 	= document["NETWORK"].GetString();
	strategiesEnforced = document["STRATEGIES_ENFORCED"].GetString();



	srand (randomSeed);
	char * strategies_c = new char [strategies.length()+1];
	strcpy (strategies_c, strategies.c_str());

	cout << "T: " 			<< T << "\n";
	cout << "b: " 			<< b << "\n";
	cout << "Random seed: " << randomSeed << "\n";


	// Load the graph from the given GML definition
	char * network_c = new char [network.length()+1];
	strcpy (network_c, network.c_str());
	igraph_t graph;
	FILE *file = fmemopen(network_c, strlen(network_c), "r");
	if (file == NULL)
	{
		cerr << "Failed to load network\n";
		exit(1);
	}

	igraph_read_graph_gml(&graph, file);
	//printAdjacencyMatrix(&graph);

	// Run Simulation
	Simulation* simulation;

	try {
		simulation = new Simulation(graph, strategies_c, T, b, randomSeed, sampleTime);

		// Activate strategy enforcement in case the strategies to enforce are given as a string
		if (strategiesEnforced.length()==strategies.length()) {
			char * strategiesEnforced_c = new char [strategiesEnforced.length()+1];
			strcpy (strategiesEnforced_c, strategiesEnforced.c_str());
			simulation->setEnforceStrategies(strategiesEnforced_c);
		}


		simulation->run();

		// Return results
		cout.rdbuf( oldCoutStreamBuf ); // Recover old stdout

		cout << "TIME,NUMBER_OF_COOPERATORS,NUMBER_OF_DEFECTORS,STRATEGIES";

		if (simulation->doEnforceStrategies) {
			cout << ",NUMBER_STRATEGY_ENFORCEMENTS";
		}

		for (unsigned int i=0; i<simulation->N; i++)
		{
			cout << ",P_" << i;
		}
		cout << "\n";
		for (unsigned int t=T-sampleTime; t<T; t++)
		{
			if (simulation->fixationReached && t > simulation->t) {
				// If fixation was reached for t < T -> just output last state
				cout << t << "," << simulation->numCooperators_t[simulation->t] << "," << simulation->numDefectors_t[simulation->t] << "," << simulation->strategies_t[simulation->t];

				if (simulation->doEnforceStrategies) {
					cout << ",0";
				}

				for (unsigned int i=0; i<simulation->N; i++)
				{
					cout << "," << simulation->payoffs_t[simulation->t][i];
				}

			} else {
				cout << t << "," << simulation->numCooperators_t[t] << "," << simulation->numDefectors_t[t] << "," << simulation->strategies_t[t];

				// Print number of times of strategy enforcement
				if (simulation->doEnforceStrategies) {
					cout << "," << simulation->numStrategyEnforcements_t[t];
				}

				for (unsigned int i=0; i<simulation->N; i++)
				{
					cout << "," << simulation->payoffs_t[t][i];
				}
			}
			cout << "\n";
		}
	}
	catch (std::exception &e)
	{
		return 1;
	}
	return 0;

}
