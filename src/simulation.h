
#ifndef INC_SIMULATION_H
#define INC_SIMULATION_H

#include <igraph/igraph.h>
#include <iostream>
#include <vector>
#include "gsl/gsl_rng.h"

class Simulation {

public:

	const static bool USE_GSL_RANDOM = false;

	constexpr static double R = 1.0;
	constexpr static double P = 0.0;
	constexpr static double S = 0.0;

	Simulation(igraph_t graph, char *strategies, unsigned int T, double b, int randomSeed, unsigned int sampleTime);
	~Simulation();

	void destroy(void);
	void run(void);
	void iterate(void);
	void playGame(void);
	void updateStrategies(void);
	void computeStats(void);
	void setEnforceStrategies(char *strategiesEnforced);
	void print(bool printAdjacency);
	void printStrategies(void);
	void writeGraphToFile(void);
	void abortWithError(std::string msg);

	igraph_t graph;							// NOC
	std::string strategies;					// The strategies of the player [0, N-1]
	unsigned int T;							// Number of time steps to simulate
	double b;								// Defector advantage (T in game matrix)
	int randomSeed;							// The random seed used for the simulation
	unsigned int sampleTime;				// Number of time steps for which the results should be returned
	bool doEnforceStrategies;				// Whether to enforce strategies (see strategiesEnforced)
	std::string strategiesEnforced;			// The strategies to enforce. After every iteration they overwrite
											// the current strategies. It is a string consisting of C|D|X
											// C or D at position i overwrite strategy[i] after update, while X (and in fact any
											// other character) will not cause strategy[i] to be overwritten.

	unsigned int N;							// Number of nodes in the NOC
	unsigned int t;							// Current time
	std::vector<double> payoffs;			// The payoffs of the players [0, N-1]
	igraph_vector_t k;						// The degrees of the nodes [0, N-1]
	std::vector<int> numCooperators_t; 		// For each time step the number of cooperators
	std::vector<int> numDefectors_t;		// For each time step the number of defectors
	std::vector<std::string> strategies_t;	// For each time step the strategy vector
	std::vector< std::vector<double> > payoffs_t;		// For each time step the payoffs of each node
	bool fixationReached;					// If true, all nodes are occupied by either coop. or defectors.
	std::vector<int> numStrategyEnforcements_t;		// For each time steps the number of individuals for which the strategy was enforced

	gsl_rng * randomNumberGenerator;

};

#endif  /* INC_SIMULATION_H */
