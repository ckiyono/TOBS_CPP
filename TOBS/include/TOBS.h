#ifndef _TOBS_H
#define _TOBS_H

#include <vector>

// Magic tricks to have CPLEX behave well:
#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// End magic tricks

class TOBS{
private:

	// TOBS parameters
	double epsilons_ = 0.001; // constraints relaxation 
	double flipLimits_ = 0.05; // how many variables are going to be changed at each iteration

	// Design variables {0,1}
	std::vector<double>* designVariables_;

	// Objective and sensitivity
	double objValue_;
	std::vector<double> objSensitivity_;

	// Constraints and sensitivities
	std::vector<double> constraintValues_;
	std::vector<double> constraintTargets_;
	std::vector<std::vector<double>> contraintSensitivities_;

	
public:

	TOBS();
	TOBS(double epsilons, double flipLimits);
	TOBS(std::vector<double>* designVariables, double epsilons, double flipLimits);

	void setDesignVariables(std::vector<double>* designVariables);
	void createDesignVariables();
	std::vector<double>* getDesignVariablesPtr();
	std::vector<double> getDesignVariables();

	void setObjective(double value, std::vector<double> sens);

	void addConstraint(double value, double target, std::vector<double> sens);

};

#endif