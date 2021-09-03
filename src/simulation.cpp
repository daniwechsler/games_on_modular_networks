#include "simulation.h"
#include <iostream>
#include <igraph/igraph.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <unistd.h>
#include <string.h>
#include "gsl/gsl_rng.h"


Simulation::Simulation(igraph_t graph, char *strategies, unsigned int T, double b, int randomSeed, unsigned int sampleTime)
{
	this->graph 				= graph;
	this->strategies 			= std::string(strategies);
	this->T 					= T;
	this->b 					= b;
	this->randomSeed 			= randomSeed;
	this->sampleTime 			= sampleTime;
	this->doEnforceStrategies 	= false;
	this->N 					= igraph_vcount(&graph);;
	this->t 					= 0;
	this->payoffs.resize(this->N);
	this->numCooperators_t.resize(this->T);
	this->numDefectors_t.resize(this->T);
	this->strategies_t.resize(this->T);
	this->payoffs_t.resize(this->T);
	this->fixationReached = false;
	this->numStrategyEnforcements_t.resize(this->T);

	int size = sizeof(std::vector<std::string>) + (sizeof(std::string) * strategies_t.size());
	std::cout << "Size: " << size << "\n";


	// Do some verification
	if (this->N != this->strategies.length()) {
		this->abortWithError("Error: Number of nodes not equal to number of strategies.");
	}

	if (this->sampleTime > this->T) {
			this->abortWithError("Error: Sample time larger than simulation time.");
	}

	for (unsigned int i = 0; i < this->N; i++)
	{
		this->payoffs[i] = 0.0;
	}
	// Fill list of node degrees
	igraph_vs_t v;
	igraph_vector_init(&this->k, 0);
	v = igraph_vss_all();
	igraph_degree(&this->graph, &this->k, v, IGRAPH_ALL, true);

	// Initialize GSL random number generator
	gsl_rng_env_setup();
	this->randomNumberGenerator = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set (this->randomNumberGenerator, this->randomSeed);
}



void Simulation::destroy(void)
{
	igraph_destroy(&this->graph);
	igraph_vector_destroy(&this->k);
	delete &this->strategies;
	delete &this->payoffs;
	delete &this->numCooperators_t;
	delete &this->numDefectors_t;
	delete &this->strategies_t;
	delete &this->payoffs_t;
}

Simulation::~Simulation(void)
{
	this->destroy();
}

void Simulation::run(void)
{

	while (this->t < this->T)
	{
		//std::cout << "Time: " << t <<  "\t " << this->strategies.length() << " \t " << this->strategies << "\n";

		// Enforce strategies in case a strategy enforecement is defined
		if (this->doEnforceStrategies)
		{
			int numberOfStrategyEnforcements = 0;
			for (unsigned int i = 0; i < this->N; i++)
			{
				if (this->strategiesEnforced[i] == 'C' && this->strategies[i] == 'D') {
					this->strategies[i] = 'C';
					numberOfStrategyEnforcements += 1;
					continue;
				}
				if (this->strategiesEnforced[i] == 'D' && this->strategies[i] == 'C') {
					this->strategies[i] = 'D';
					numberOfStrategyEnforcements += 1;
					continue;
				}
			}
			this->numStrategyEnforcements_t[this->t] = numberOfStrategyEnforcements;
		}

		this->computeStats();
		// Abort if fixation is reached
		if (this->numCooperators_t[this->t] == 0 || this->numCooperators_t[this->t] == this->N)
		{
			this->fixationReached = true;
			this->playGame();	// Compute the payoff for the current generation
			return;
		}
		this->iterate();
		this->t++;
	}
}

void Simulation::iterate(void)
{
	this->playGame();
	this->updateStrategies();
}

void Simulation::playGame(void)
{
	igraph_vector_t neighbors;
	double P_i, P_i_j;

	for (unsigned int i = 0; i < this->N; i++)
	{
		igraph_vector_init(&neighbors, 0);
		igraph_neighbors(&this->graph, &neighbors, i, IGRAPH_ALL);

		P_i = 0.0;
		int j;
		for (int j_index = 0; j_index < igraph_vector_size(&neighbors); j_index++)
		{
			j = VECTOR(neighbors)[j_index]; // Get neighbor vertex
			P_i_j = 0.0;

			if (this->strategies[i] == 'C' && this->strategies[j] == 'C') P_i_j += R;
			if (this->strategies[i] == 'D' && this->strategies[j] == 'C') P_i_j += this->b;
			if (this->strategies[i] == 'C' && this->strategies[j] == 'D') P_i_j += S;
			if (this->strategies[i] == 'D' && this->strategies[j] == 'D') P_i_j += P;

			P_i += P_i_j;
		}
		this->payoffs[i] = P_i;
		igraph_vector_destroy(&neighbors);
	}
	this->payoffs_t[this->t] = this->payoffs;
}

void Simulation::updateStrategies()
{

	std::string strategiesNext = std::string(this->N, 'a');
	igraph_vector_t neighbors;
	igraph_vector_init(&neighbors, 0);


	for (unsigned int i = 0; i < this->N; i++)
	{
		int k_i = VECTOR(this->k)[i];
		if (k_i == 0)
			this->abortWithError("Error: There is a disconnected node in the network.");

		igraph_neighbors(&this->graph, &neighbors, i, IGRAPH_ALL);

		if (k_i != igraph_vector_size(&neighbors))
			this->abortWithError("Error: Node degree k is different from number of neighbors.");

		int j;
		if (USE_GSL_RANDOM) { // Select a random neighbor
			j = VECTOR(neighbors)[gsl_rng_uniform_int(this->randomNumberGenerator, k_i)];
		} else {
			j = VECTOR(neighbors)[rand() % k_i];
		}

		int k_j 		= VECTOR(this->k)[j];				// Get degree of the neighbor
		double P_i 		= this->payoffs[i];
		double P_j 		= this->payoffs[j];
		double k_max;

		if (k_j > k_i)
		{
			k_max = double(k_j);
		}
		else
		{
			k_max = double(k_i);
		}

		double p = (P_j - P_i) / (k_max * (this->b - S));

		if (p <= 0.0)
		{
			strategiesNext[i] = this->strategies[i];
			continue;
		}

		double r;

		if (USE_GSL_RANDOM) { // Depending on which random generator to use
			r = gsl_rng_uniform (this->randomNumberGenerator);
		} else {
			long randMax = (long) RAND_MAX;
			r = rand() / ((double) (RAND_MAX));
		}

		if (r <= p)
		{
			strategiesNext[i] = this->strategies[j];
		}
		else
		{
			strategiesNext[i] = this->strategies[i];
		}

	}

	this->strategies = strategiesNext;
	igraph_vector_destroy(&neighbors);
}

void Simulation::computeStats(void) {

	this->numCooperators_t[this->t] = 0;
	this->numDefectors_t[this->t] = 0;

	for (unsigned int i = 0; i < this->N; i++)
	{
		if (this->strategies[i] == 'C')
		{
			this->numCooperators_t[this->t]++;
		}
		else
		{
			this->numDefectors_t[this->t]++;
		}
	}
	this->strategies_t[this->t] = this->strategies;

	return;
}

void Simulation::setEnforceStrategies(char *strategiesEnforced) {

	this->strategiesEnforced = std::string(strategiesEnforced);
	this->doEnforceStrategies = true;
}


void Simulation::print(bool printAdjacency)
{
	std::cout << "N: " << this->N << "\n";
	std::cout << "T: " << this->T << "\n";
	std::cout << "b: " << this->b << "\n";

	int k_max = 0;
	int i_k_max = 0;
	for (unsigned int i = 0; i < this->N; i++)
	{
		if (VECTOR(this->k)[i] > k_max)
		{
			k_max = VECTOR(this->k)[i];
			i_k_max = i;
		}
	}
	std::cout << "Highest degree node k=" << k_max << " S="
			<< this->strategies[i_k_max] << "\n";

	if (!printAdjacency)
		return;

	igraph_matrix_t A;
	igraph_matrix_init(&A, 0, 0);
	igraph_get_adjacency(&this->graph, &A, IGRAPH_GET_ADJACENCY_BOTH, false);
	std::cout << "\n    ";

	for (unsigned int i = 0; i < this->N; i++)
	{
		std::cout << this->strategies[i] << " ";
	}
	std::cout << "\n";
	for (unsigned int i = 0; i < this->N; i++)
	{
		std::cout << this->strategies[i] << " ";
		std::cout << "| ";
		for (unsigned int j = 0; j < this->N; j++)
		{
			std::cout << MATRIX(A, i, j) << " ";
		}
		std::cout << "|\n";
	}
	std::cout << "\n";

}

void Simulation::printStrategies(void)
{
	for (unsigned int i = 0; i < this->N; i++)
	{
		std::cout << this->strategies[i] << " ";
	}
	std::cout << "\n";
}

void Simulation::writeGraphToFile(void)
{
	std::cout << "Write to file";

	FILE * outFile;
	outFile = fopen("graph.txt", "w");

	igraph_write_graph_edgelist(&this->graph, outFile);
}

void Simulation::abortWithError (std::string msg)
{
	std::cerr << msg << "\n";
	throw std::exception();
}

