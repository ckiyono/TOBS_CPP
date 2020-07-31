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

int main()
{
	std::cout << "Hello TOBS" << std::endl;
	std::cin.ignore();

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
		env.out() << "Values = " << vals << std::endl;

	}
	catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}

	env.end();

	std::cin.ignore();

	return 0;
}