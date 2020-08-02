#include "TOBS.h"
#include <algorithm>

// private methods

void TOBS::calculateVariableLimits()
{
	std::fill(lowerLimits_.begin(), lowerLimits_.end(), 0);
	std::fill(upperLimits_.begin(), upperLimits_.end(), 0);

	size_t nVars = designVariables_->size();

	for (int i = 0; i < nVars; i++)
	{
		if (designVariables_->at(i) < 0.5)
			upperLimits_[i] = 1;
		else
			lowerLimits_[i] = -1;
	}
}

void TOBS::constraintRelaxation(std::vector<double>& g,
	std::vector<double>& gbar)
{
	size_t nConst = constraintValues_.size();

	for (int j = 0; j < nConst; j++)
	{
		if (gbar[j] < (1 - epsilons_)*g[j])
		{
			gbar[j] = -epsilons_*g[j];
		}
		else if (gbar[j] > (1 + epsilons_)*g[j])
		{
			gbar[j] = epsilons_*g[j];
		}
		else
		{
			gbar[j] -= g[j];
		}
	}

}

void TOBS::constraintNormalization(std::vector<double>& g,
	std::vector<double>& gbar,
	std::vector<std::vector<double>>& dgdx)
{
	size_t nVars = designVariables_->size();
	size_t nConst = constraintValues_.size();

	double norm = 1e-22;

	// getting norm value
	for (int j = 0; j < nConst; j++)
	{
		for (int i = 0; i < nVars; i++)
		{
			norm = max(norm, abs(dgdx[j][i]));
		}
	}

	// normalizing constraints
	for (int j = 0; j < nConst; j++)
	{
		g[j] /= norm;
		gbar[j] /= norm;

		for (int i = 0; i < nVars; i++)
		{
			dgdx[j][i] /= norm;
		}
	}
}

void TOBS::addFlipLimitConstraint(std::vector<double>& gbar,
	std::vector<std::vector<double>>& dgdx)
{
	size_t nVars = designVariables_->size();
	gbar.push_back(flipLimits_*nVars);

	std::vector<double> truncation;

	for (int i = 0; i < nVars; i++)
	{
		if (designVariables_->at(i) < 0.5)
			truncation.push_back(1.0);
		else
			truncation.push_back(-1.0);
	}

	dgdx.push_back(truncation);
}

void TOBS::cplexSolve(const std::vector<double> dfdx, 
	const std::vector<double> gbar,
	const std::vector<std::vector<double>> dgdx)
{
	IloEnv env;

	try {
		IloModel model(env);
		IloNumVarArray vars(env);

		size_t nVars = designVariables_->size();
		size_t nConst = gbar.size();

		IloObjective obj;
		IloNumArray objCoefs(env, nVars);

		if (isMaximization_)
		{
			obj = IloMaximize(env, 0);
		}
		else
		{
			obj = IloMinimize(env, 0);
		}

		for (int i = 0; i < nVars; i++)
		{
			vars.add(IloNumVar(env, lowerLimits_[i], upperLimits_[i], ILOINT));
			objCoefs[i] = dfdx[i];
		}
		obj.setLinearCoefs(vars, objCoefs);

		model.add(obj);

		for (int j = 0; j < nConst; j++)
		{
			IloRange constraint = IloRange(env, -IloInfinity, gbar[j]);
			IloNumArray conCoefs(env, nVars);

			for (int i = 0; i < nVars; i++)
			{
				conCoefs[i] = dgdx[j][i];
			}

			constraint.setLinearCoefs(vars, conCoefs);

			model.add(constraint);
		}

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		if (!cplex.solve()) {
			env.error() << "Failed to optimize Cplex." << std::endl;
			throw(-1);
		}

		IloNumArray vals(env);
		cplex.getValues(vals, vars);

		for (int i = 0; i < nVars; i++)
		{
			designVariables_->at(i) = round(designVariables_->at(i) + vals[i]);
		}
	}
	catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}

	env.end();
}

void TOBS::ILP()
{
	size_t nVars = designVariables_->size();
	size_t nConst = constraintValues_.size();

	std::vector<double> dfdx = *objSensitivity_;

	std::vector<double> g;
	std::vector<double> gbar;
	std::vector<std::vector<double>> dgdx;

	for (int j = 0; j < nConst; j++)
	{
		g.push_back(*constraintValues_[j]);
		gbar.push_back(constraintTargets_[j]);
		dgdx.push_back(*contraintSensitivities_[j]);
	}

	constraintNormalization(g, gbar, dgdx);

	constraintRelaxation(g, gbar);

	addFlipLimitConstraint(gbar, dgdx);

	cplexSolve(dfdx, gbar, dgdx);

}



// public methods

TOBS::TOBS()
{

}

TOBS::TOBS(double epsilons, double flipLimits)
{
	epsilons_ = epsilons;
	flipLimits_ = flipLimits;
}

TOBS::TOBS(std::vector<double>* designVariables, double epsilons, double flipLimits)
{
	designVariables_ = designVariables;
	epsilons_ = epsilons;
	flipLimits_ = flipLimits;

	size_t nVars = designVariables_->size();
	lowerLimits_.resize(nVars, 0.0);
	upperLimits_.resize(nVars, 0.0);
}

void TOBS::setDesignVariables(std::vector<double>* designVariables)
{
	designVariables_ = designVariables;

	size_t nVars = designVariables_->size();
	lowerLimits_.resize(nVars, 0.0);
	upperLimits_.resize(nVars, 0.0);
}

void TOBS::createDesignVariables()
{
	if (designVariables_ == nullptr)
		designVariables_ = new std::vector<double>;

	size_t nVars = designVariables_->size();
	lowerLimits_.resize(nVars, 0.0);
	upperLimits_.resize(nVars, 0.0);
}

std::vector<double>* TOBS::getDesignVariablesPtr()
{
	return designVariables_;
}

std::vector<double> TOBS::getDesignVariables()
{
	return *designVariables_;
}

void TOBS::setObjective(double* value, std::vector<double>* sens, bool isMaximization)
{
	objValue_ = value;
	objSensitivity_ = sens;
	isMaximization_ = isMaximization;
}

void TOBS::addConstraint(double* value, double target, std::vector<double>* sens)
{
	constraintValues_.push_back(value);
	constraintTargets_.push_back(target);
	contraintSensitivities_.push_back(sens);
}

void TOBS::solve()
{
	calculateVariableLimits();

	ILP();
}