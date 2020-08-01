#include "TOBS.h"

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
}

void TOBS::setDesignVariables(std::vector<double>* designVariables)
{
	designVariables_ = designVariables;
}

void TOBS::createDesignVariables()
{
	if (designVariables_ == nullptr)
		designVariables_ = new std::vector<double>;
}

std::vector<double>* TOBS::getDesignVariablesPtr()
{
	return designVariables_;
}

std::vector<double> TOBS::getDesignVariables()
{
	return *designVariables_;
}

void TOBS::setObjective(double value, std::vector<double> sens)
{
	objValue_ = value;
	objSensitivity_ = sens;
}

void TOBS::addConstraint(double value, double target, std::vector<double> sens)
{
	constraintValues_.push_back(value);
	constraintTargets_.push_back(target);
	contraintSensitivities_.push_back(sens);
}