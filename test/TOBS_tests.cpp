#include "gtest/gtest.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <iostream>

#include "TOBS.h"

#include <vector>

TEST(TOBS_Test, Definition_test)
{
	// Integer Programming test with 4 constraints
	// Maximixe     8*x[0] + 11*x[1] + 6*x[2] + 4*x[3]
	// Subject to   5*x[0] +  7*x[1] + 4*x[2] + 3*x[3] <= 14
	//                x[0] +    x[1] +   x[2] +   x[3] <=  2
	//					        x[1]          -   x[3] <=  0
	//                x[0]           +   x[2]          <=  1
	// Bounds       0 <= x[i] <= 1

	int nVars = 4;
	int nConst = 4;

	std::vector<double>* designVariables = new std::vector<double>;
	designVariables->push_back(1);
	designVariables->push_back(1);
	designVariables->push_back(1);
	designVariables->push_back(1);

	double objValue = 0;

	std::vector<double> objSens;
	objSens.push_back(8);
	objSens.push_back(11);
	objSens.push_back(6);
	objSens.push_back(4);

	std::vector<double> constraintValues;
	constraintValues.push_back(0);
	constraintValues.push_back(0);
	constraintValues.push_back(0);
	constraintValues.push_back(0);

	std::vector<double> constraintTargets;
	constraintTargets.push_back(14);
	constraintTargets.push_back(2);
	constraintTargets.push_back(0);
	constraintTargets.push_back(1);

	std::vector<std::vector<double>> constraintSens;
	std::vector<double> cSens;

	// constraint 1
	cSens.push_back(5);
	cSens.push_back(7);
	cSens.push_back(4);
	cSens.push_back(3);
	constraintSens.push_back(cSens);
	cSens.clear();
	
	// constraint 2
	cSens.push_back(1);
	cSens.push_back(1);
	cSens.push_back(1);
	cSens.push_back(1);
	constraintSens.push_back(cSens);
	cSens.clear();

	// constraint 3
	cSens.push_back(0);
	cSens.push_back(1);
	cSens.push_back(0);
	cSens.push_back(-1);
	constraintSens.push_back(cSens);
	cSens.clear();

	// constraint 4
	cSens.push_back(1);
	cSens.push_back(0);
	cSens.push_back(1);
	cSens.push_back(0);
	constraintSens.push_back(cSens);
	cSens.clear();

	TOBS tobs(designVariables, 0.01, 0.25);
	
	for (int i = 0; i < nVars; i++)
	{
		objValue += (*designVariables)[i] * objSens[i];
	}
	tobs.setObjective(&objValue, &objSens, true);
	

	for (int j = 0; j < nConst; j++)
	{
		for (int i = 0; i < nVars; i++)
		{
			constraintValues[j] += (*designVariables)[i] * constraintSens[j][i];
		}
		tobs.addConstraint(&constraintValues[j], constraintTargets[j], &constraintSens[j]);
	}

	double diff = 1;
	double objValueOld = objValue;

	while (diff > 1e-8)
	{
		std::cout << std::endl;
		std::cout << "Solution value = " << objValue << std::endl;
		std::cout << "Values = [";
		for (int i = 0; i < nVars; i++)
		{
			std::cout << (*designVariables)[i];
			if (i != nVars - 1)
				std::cout << ", ";
		}
		std::cout << "]" << std::endl;
		std::cout << std::endl;

		tobs.solve();

		objValue = 0.0;

		for (int i = 0; i < nVars; i++)
		{
			objValue += (*designVariables)[i] * objSens[i];
		}

		for (int j = 0; j < nConst; j++)
		{
			constraintValues[j] = 0.0;

			for (int i = 0; i < nVars; i++)
			{
				constraintValues[j] += (*designVariables)[i] * constraintSens[j][i];
			}
		}

		diff = abs(objValue - objValueOld) / abs(objValueOld);

		objValueOld = objValue;

	}
}