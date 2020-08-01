#include "gtest/gtest.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <iostream>

// Magic tricks to have CPLEX behave well:
#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// End magic tricks

TEST(Cplex_Test, LP_test)
{
	// Linear Programming test
	// Maximixe     x[0] + 2*x[1] + 3*x[2]
	// Subject to  -x[0] +   x[1] +   x[2] <= 20
	//              x[0] - 3*x[1] +   x[2] <= 30
	// Bounds       0 <= x[0] <= 40 
	//              0 <= x[1] <= infinity
	//              0 <= x[2] <= infinity

	std::cout << std::endl;

	IloEnv env;

	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env, 0.0, 40.0));
		vars.add(IloNumVar(env));
		vars.add(IloNumVar(env));
		
		IloRange constraint1 = IloRange(env, -IloInfinity, 20);
		IloNumArray con1Coefs(env, 3);
		con1Coefs[0] = -1;
		con1Coefs[1] = 1;
		con1Coefs[2] = 1;
		constraint1.setLinearCoefs(vars, con1Coefs);
		IloRange constraint2 = IloRange(env, -IloInfinity, 30);
		IloNumArray con2Coefs(env, 3);
		con2Coefs[0] = 1;
		con2Coefs[1] = -3;
		con2Coefs[2] = 1;
		constraint2.setLinearCoefs(vars, con2Coefs);

		IloObjective obj = IloMaximize(env, 0);
		IloNumArray objCoefs(env,3);
		objCoefs[0] = 1;
		objCoefs[1] = 2;
		objCoefs[2] = 3;
		obj.setLinearCoefs(vars, objCoefs);

		model.add(obj);
		model.add(constraint1);
		model.add(constraint2);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << std::endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution value = " << cplex.getObjValue() << std::endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << std::endl << std::endl;

	}
	catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}

	env.end();
}

TEST(Cplex_Test, ILP_test_1_Constraint)
{
	// Integer Programming test with 1 constraint
	// Maximixe     8*x[0] + 11*x[1] + 6*x[2] + 4*x[3]
	// Subject to   5*x[0] +  7*x[1] + 4*x[2] + 3*x[3] <= 14
	// Bounds       0 <= x[i] <= 1

	std::cout << std::endl;

	IloEnv env;

	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));

		IloObjective obj = IloMaximize(env, 0);
		IloNumArray objCoefs(env, 4);
		objCoefs[0] = 8;
		objCoefs[1] = 11;
		objCoefs[2] = 6;
		objCoefs[3] = 4;
		obj.setLinearCoefs(vars, objCoefs);

		IloRange constraint1 = IloRange(env, -IloInfinity, 14);
		IloNumArray con1Coefs(env, 4);
		con1Coefs[0] = 5;
		con1Coefs[1] = 7;
		con1Coefs[2] = 4;
		con1Coefs[3] = 3;
		constraint1.setLinearCoefs(vars, con1Coefs);
		
		model.add(obj);
		model.add(constraint1);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << std::endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution value = " << cplex.getObjValue() << std::endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << std::endl << std::endl;

	}
	catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}

	env.end();
}

TEST(Cplex_Test, ILP_test_2_Constraints)
{
	// Integer Programming test with 2 constraints
	// Maximixe     8*x[0] + 11*x[1] + 6*x[2] + 4*x[3]
	// Subject to   5*x[0] +  7*x[1] + 4*x[2] + 3*x[3] <= 14
	//                x[0] +    x[1] +   x[2] +   x[3] <=  2
	// Bounds       0 <= x[i] <= 1

	std::cout << std::endl;

	IloEnv env;

	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));

		IloObjective obj = IloMaximize(env, 0);
		IloNumArray objCoefs(env, 4);
		objCoefs[0] = 8;
		objCoefs[1] = 11;
		objCoefs[2] = 6;
		objCoefs[3] = 4;
		obj.setLinearCoefs(vars, objCoefs);

		IloRange constraint1 = IloRange(env, -IloInfinity, 14);
		IloNumArray con1Coefs(env, 4);
		con1Coefs[0] = 5;
		con1Coefs[1] = 7;
		con1Coefs[2] = 4;
		con1Coefs[3] = 3;
		constraint1.setLinearCoefs(vars, con1Coefs);

		IloRange constraint2 = IloRange(env, -IloInfinity, 2);
		IloNumArray con2Coefs(env, 4);
		con2Coefs[0] = 1;
		con2Coefs[1] = 1;
		con2Coefs[2] = 1;
		con2Coefs[3] = 1;
		constraint2.setLinearCoefs(vars, con2Coefs);

		model.add(obj);
		model.add(constraint1);
		model.add(constraint2);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << std::endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution value = " << cplex.getObjValue() << std::endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << std::endl << std::endl;

	}
	catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}

	env.end();
}

TEST(Cplex_Test, ILP_test_3_Constraints)
{
	// Integer Programming test with 3 constraints
	// Maximixe     8*x[0] + 11*x[1] + 6*x[2] + 4*x[3]
	// Subject to   5*x[0] +  7*x[1] + 4*x[2] + 3*x[3] <= 14
	//                x[0] +    x[1] +   x[2] +   x[3] <=  2
	//					        x[1]          -   x[3] <=  0
	// Bounds       0 <= x[i] <= 1

	std::cout << std::endl;

	IloEnv env;

	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));

		IloObjective obj = IloMaximize(env, 0);
		IloNumArray objCoefs(env, 4);
		objCoefs[0] = 8;
		objCoefs[1] = 11;
		objCoefs[2] = 6;
		objCoefs[3] = 4;
		obj.setLinearCoefs(vars, objCoefs);

		IloRange constraint1 = IloRange(env, -IloInfinity, 14);
		IloNumArray con1Coefs(env, 4);
		con1Coefs[0] = 5;
		con1Coefs[1] = 7;
		con1Coefs[2] = 4;
		con1Coefs[3] = 3;
		constraint1.setLinearCoefs(vars, con1Coefs);

		IloRange constraint2 = IloRange(env, -IloInfinity, 2);
		IloNumArray con2Coefs(env, 4);
		con2Coefs[0] = 1;
		con2Coefs[1] = 1;
		con2Coefs[2] = 1;
		con2Coefs[3] = 1;
		constraint2.setLinearCoefs(vars, con2Coefs);

		IloRange constraint3 = IloRange(env, -IloInfinity, 0);
		IloNumArray con3Coefs(env, 4);
		con3Coefs[0] =  0;
		con3Coefs[1] =  1;
		con3Coefs[2] =  0;
		con3Coefs[3] = -1;
		constraint3.setLinearCoefs(vars, con3Coefs);

		model.add(obj);
		model.add(constraint1);
		model.add(constraint2);
		model.add(constraint3);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << std::endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution value = " << cplex.getObjValue() << std::endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << std::endl << std::endl;

	}
	catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}

	env.end();
}

TEST(Cplex_Test, ILP_test_4_Constraints)
{
	// Integer Programming test with 4 constraints
	// Maximixe     8*x[0] + 11*x[1] + 6*x[2] + 4*x[3]
	// Subject to   5*x[0] +  7*x[1] + 4*x[2] + 3*x[3] <= 14
	//                x[0] +    x[1] +   x[2] +   x[3] <=  2
	//					        x[1]          -   x[3] <=  0
	//                x[0]           +   x[2]          <=  1
	// Bounds       0 <= x[i] <= 1

	std::cout << std::endl;

	IloEnv env;

	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));
		vars.add(IloNumVar(env, 0, 1, ILOINT));

		IloObjective obj = IloMaximize(env, 0);
		IloNumArray objCoefs(env, 4);
		objCoefs[0] = 8;
		objCoefs[1] = 11;
		objCoefs[2] = 6;
		objCoefs[3] = 4;
		obj.setLinearCoefs(vars, objCoefs);

		IloRange constraint1 = IloRange(env, -IloInfinity, 14);
		IloNumArray con1Coefs(env, 4);
		con1Coefs[0] = 5;
		con1Coefs[1] = 7;
		con1Coefs[2] = 4;
		con1Coefs[3] = 3;
		constraint1.setLinearCoefs(vars, con1Coefs);

		IloRange constraint2 = IloRange(env, -IloInfinity, 2);
		IloNumArray con2Coefs(env, 4);
		con2Coefs[0] = 1;
		con2Coefs[1] = 1;
		con2Coefs[2] = 1;
		con2Coefs[3] = 1;
		constraint2.setLinearCoefs(vars, con2Coefs);

		IloRange constraint3 = IloRange(env, -IloInfinity, 0);
		IloNumArray con3Coefs(env, 4);
		con3Coefs[0] = 0;
		con3Coefs[1] = 1;
		con3Coefs[2] = 0;
		con3Coefs[3] = -1;
		constraint3.setLinearCoefs(vars, con3Coefs);

		IloRange constraint4 = IloRange(env, -IloInfinity, 1);
		IloNumArray con4Coefs(env, 4);
		con4Coefs[0] = 0;
		con4Coefs[1] = 1;
		con4Coefs[2] = 0;
		con4Coefs[3] = -1;
		constraint4.setLinearCoefs(vars, con4Coefs);

		model.add(obj);
		model.add(constraint1);
		model.add(constraint2);
		model.add(constraint3);
		model.add(constraint4);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << std::endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution value = " << cplex.getObjValue() << std::endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << std::endl << std::endl;

	}
	catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}

	env.end();
}