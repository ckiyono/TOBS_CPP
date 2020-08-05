#include "gtest/gtest.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <dirent.h>

#include "TOBS.h"

#include <vector>
#include <chrono>

// Eigen libs
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

typedef Eigen::Triplet<double> T;
namespace chronos = std::chrono;

void checkFolder(std::string folder)
{
	DIR* dir = opendir(folder.c_str());
	if (dir)
	{
		struct dirent *next_file;
		char filepath[256];
		while ((next_file = readdir(dir)) != NULL)
		{
			sprintf(filepath, "%s/%s", folder.c_str(), next_file->d_name);
			remove(filepath);
		}
		closedir(dir);
	}
	else
	{
		_mkdir(folder.c_str());
	}
}

void writeResult(Eigen::MatrixXd noCoord, 
	Eigen::MatrixXi elConnect, 
	Eigen::VectorXd values,
	std::string fileName)
{
	std::ofstream file;
	file.open(fileName.c_str());

	file << "# vtk DataFile Version 1.0\n";
	file << "Unstructured Grid Example\n";
	file << "ASCII\n\n";

	file << "DATASET UNSTRUCTURED_GRID\n";
	file << "POINTS " << noCoord.rows() << " float\n";
	for (int row = 0; row < noCoord.rows(); row++)
	{
		file << noCoord(row, 0) << " " << noCoord(row, 1) << " 0.0\n";
	}
	file << "\n";

	file << "CELLS " << elConnect.rows() << " " << elConnect.rows()*(elConnect.cols()+1) << "\n";
	for (int row = 0; row < elConnect.rows(); row++)
	{
		file << "4 " << elConnect(row, 0) << " " << 
			elConnect(row, 1) << " " <<
			elConnect(row, 2) << " " << 
			elConnect(row, 3) << "\n";
	}
	file << "\n";

	file << "CELL_TYPES " << elConnect.rows() << "\n";
	for (int row = 0; row < elConnect.rows(); row++)
	{
		file << "9\n";
	}
	file << "\n";

	file << "CELL_DATA " << elConnect.rows() << "\n";
	file << "SCALARS rho float\n";
	file << "LOOKUP_TABLE xray\n";
	for (int row = 0; row < elConnect.rows(); row++)
	{
		file << values(row) << "\n";
	}
	file << "\n";


	file.close();
}

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

TEST(TOBS_Test, MBB_Test)
{
	chronos::steady_clock::time_point begin_total = chronos::steady_clock::now();
	chronos::steady_clock::time_point begin;
	double timeEl;

	// Input parameters
	int nelx = 120; // number of elements in x
	int nely = 40; // number of elements in y
	double gbar = 0.5; // volume cosntraint target
	double epsilons = 0.01; // constraint relaxation
	double beta = 0.05; // variable flip limits
	double rmin = 4; // filter radius in elements
	int maxLoop = 150;
	std::string workingDir = "E:\\OneDrive\\Trabalho\\Post-Doc2\\TOBS\\Result\\";
	std::string fileName = "topology";
	std::string extension = ".vtk";
	checkFolder(workingDir);

	// Material properties
	double E0 = 1;
	double Emin = 1e-9;
	double nu = 0.3;
	double penal = 3;

	// Preparing Finite Element Analysis
	Eigen::MatrixXd A11(4, 4), A12(4, 4), A(8, 8), B11(4, 4), B12(4, 4), B(8, 8);
	A11 << 12, 3, -6, -3, 3, 12, 3, 0, -6, 3, 12, -3, -3, 0, -3, 12;
	A12 << -6, -3, 0, 3, -3, -6, -3, -6, 0, -3, -6, 3, 3, -6, 3, -6;
	B11 << -4, 3, -2, 9, 3, -4, -9, 4, -2, -9, -4, -3, 9, 4, -3, -4;
	B12 << 2, -3, 4, -9, -3, 2, 9, -2, 4, 9, 2, 3, -9, -2, 3, 2;
	A << A11, A12, A12.transpose(), A11;
	B << B11, B12, B12.transpose(), B11;

	Eigen::MatrixXd KE = 1 / (1 - nu*nu) / 24 * (A + nu*B);

	Eigen::VectorXi nodes = Eigen::VectorXi::LinSpaced((nelx + 1)*(nely + 1), 0, (nelx + 1)*(nely + 1) - 1);
	Eigen::MatrixXi dofs(nodes.size(), 2);
	dofs << nodes, nodes + Eigen::VectorXi::Ones(nodes.size());

	Eigen::VectorXd coordX = Eigen::VectorXd::LinSpaced((nelx + 1), 0, 120);
	Eigen::VectorXd coordY = Eigen::VectorXd::LinSpaced((nely + 1), 40, 0);

	Eigen::MatrixXd coordXMat = coordX.transpose().replicate(nely + 1, 1);
	Eigen::MatrixXd coordYMat = coordY.replicate(1, nelx + 1);

	coordX = Eigen::Map<Eigen::VectorXd>(coordXMat.data(), (nelx + 1)*(nely + 1));
	coordY = Eigen::Map<Eigen::VectorXd>(coordYMat.data(), (nelx + 1)*(nely + 1));

	Eigen::MatrixXd coordMat((nelx + 1)*(nely + 1), 2);
	coordMat << coordX, coordY;

	//std::cout << "coordMat = " << std::endl << coordMat << std::endl;
	
	Eigen::MatrixXi elConnect(nelx*nely, 4);
	int el = 0;
	for (int elx = 0; elx < nelx; elx++)
	{
		for (int ely = 0; ely < nely; ely++)
		{
			elConnect(el, 0) = elx*(nely + 1) + 1 + ely;
			elConnect(el, 1) = (elx+1)*(nely + 1) + 1 + ely;
			elConnect(el, 2) = (elx+1)*(nely + 1) + ely;
			elConnect(el, 3) = elx*(nely + 1) + ely;
			el++;
		}
	}

	//std::cout << "elConnect = " << std::endl << elConnect << std::endl;

	Eigen::MatrixXi nodesM = Eigen::Map<Eigen::MatrixXi>(nodes.data(), 1 + nely, 1 + nelx);
	Eigen::MatrixXi nodesM2 = 2 * nodesM.block(0, 0, nely, nelx);
	
	Eigen::VectorXi edofVec = Eigen::Map<Eigen::VectorXi>(nodesM2.data(), nelx*nely, 1);
	Eigen::VectorXi edofVec1(8);
	edofVec1 << 2, 3, 2 * nely + 4, 2 * nely + 5, 2 * nely + 2, 2 * nely + 3, 0, 1;
	Eigen::MatrixXi edofMat = edofVec.replicate(1, 8) + edofVec1.transpose().replicate(nelx*nely, 1);
	
	Eigen::MatrixXi kronI = Eigen::kroneckerProduct(edofMat, Eigen::MatrixXi::Ones(8, 1)).eval().transpose();
	Eigen::MatrixXi kronJ = Eigen::kroneckerProduct(edofMat, Eigen::MatrixXi::Ones(1, 8)).eval().transpose();

	Eigen::VectorXi iK = Eigen::Map<Eigen::VectorXi>(kronI.data(), 64 * nelx*nely, 1);
	Eigen::VectorXi jK = Eigen::Map<Eigen::VectorXi>(kronJ.data(), 64 * nelx*nely, 1);

	Eigen::SparseVector<double> F(2 * (nely + 1)*(nelx + 1));
	F.coeffRef(2 * nely + 1) = -1;
	Eigen::VectorXd U = Eigen::VectorXd::Zero(2 * (nely + 1)*(nelx + 1));

	Eigen::VectorXi fixeddofs(nely + 1 + 1);
	fixeddofs << Eigen::VectorXi::LinSpaced(nely + 1, 0, 2 * nely), 2 * (nelx + 1) *(nely + 1) - 1;

	Eigen::VectorXi allddofs = Eigen::VectorXi::LinSpaced(2 * (nelx + 1)*(nely + 1), 0, 2 * (nelx + 1)*(nely + 1) - 1);

	Eigen::VectorXi freeddofs(allddofs.size());
	auto it = std::set_difference(allddofs.data(), allddofs.data() + allddofs.size(),
		fixeddofs.data(), fixeddofs.data() + fixeddofs.size(),
		freeddofs.data());
	freeddofs.conservativeResize(std::distance(freeddofs.data(), it));

	// Prepare Filter
	std::vector<T> tripletList;
	for (int i1 = 0; i1 < nelx; i1++)
	{
		for (int j1 = 0; j1 < nely; j1++)
		{
			int e1 = i1*nely + j1;
			for (int i2 = max(i1 - (int(ceil(rmin)) - 1), 0); i2 <= min(i1 + (int(ceil(rmin)) - 1), nelx - 1); i2++)
			{
				for (int j2 = max(j1 - (int(ceil(rmin)) - 1), 0); j2 <= min(j1 + (int(ceil(rmin)) - 1), nely - 1); j2++)
				{
					int e2 = i2*nely + j2;
					double value = max(0.0, rmin - sqrt(pow(i1 - i2, 2) + pow(j1 - j2, 2)));
					tripletList.push_back(T(e1, e2, value));
				}
			}
		}
	}
	Eigen::SparseMatrix<double> H(nelx*nely, nelx*nely);
	H.setFromTriplets(tripletList.begin(), tripletList.end());
	Eigen::VectorXd Hs = H*Eigen::VectorXd::Ones(nelx*nely);

	// Initialize Iteration
	Eigen::VectorXd x(nelx*nely);
	int loop = 0;
	double change = 1;
	double obj = 0;
	double g = 0;
	Eigen::VectorXd objVector(maxLoop);
	Eigen::VectorXd objSens(nelx*nely);
	Eigen::VectorXd oldObjSens = Eigen::VectorXd::Zero(nelx*nely,1);
	Eigen::VectorXd gSens(nelx*nely);

	// TOBS initialization
	std::vector<double>* designVariables = new std::vector<double>(nelx*nely);
	std::fill(designVariables->begin(), designVariables->end(), 1);
	std::vector<double> dobj(nelx*nely), dcon(nelx*nely);
	TOBS tobs(designVariables, epsilons, beta);
	tobs.setObjective(&obj, &dobj, false);
	tobs.addConstraint(&g, gbar, &dcon);

	while (change > 1e-4)
	{
		// copy std::vector to Eigen::VectorXd
		x = Eigen::Map<Eigen::VectorXd>(designVariables->data(), designVariables->size());

		// Assemble K
		std::cout << std::endl << "---Assembling K" << std::endl;
		begin = chronos::steady_clock::now();
		Eigen::VectorXd Ke = Eigen::Map<Eigen::VectorXd>(KE.data(), 64, 1);
		Eigen::VectorXd E = Emin*Eigen::VectorXd::Ones(x.size()) + x.cwiseProduct((E0 - Emin)*Eigen::VectorXd::Ones(x.size()));
		Eigen::MatrixXd Ke_E = Ke* E.transpose();
		Eigen::VectorXd sK = Eigen::Map<Eigen::VectorXd>(Ke_E.data(), 64 * nelx*nely, 1);

		tripletList.clear();
		for (int ind = 0; ind < sK.size(); ind++)
		{
			tripletList.push_back(T(iK(ind), jK(ind), sK(ind)));
		}
		Eigen::SparseMatrix<double> K(2 *(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
		K.setFromTriplets(tripletList.begin(), tripletList.end());
		K = 0.5*(K + Eigen::SparseMatrix<double>(K.transpose()));
		Eigen::SparseMatrix<double> Kdiag(2 * (nelx + 1)*(nely + 1), 2 * (nelx + 1)*(nely + 1));

		// Direchlet
		for (int ind = 0; ind < fixeddofs.size(); ind++)
		{
			K.row(fixeddofs(ind)) *= 0;
			K.col(fixeddofs(ind)) *= 0;
			Kdiag.insert(fixeddofs(ind), fixeddofs(ind)) = 1;
		}
		K.makeCompressed();
		Kdiag.makeCompressed();
		K = K + Kdiag;

		timeEl = chronos::duration_cast<chronos::microseconds>(chronos::steady_clock::now() - begin).count() / 1e6;
		std::cout << "---Finished. TimeElapsed: " << timeEl  << "s" << std::endl << std::endl;

		// Solve FE Analysis
		std::cout << std::endl << "---Solving FE" << std::endl;
		begin = chronos::steady_clock::now();
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
		U = solver.compute(K).solve(F);
		timeEl = chronos::duration_cast<chronos::microseconds>(chronos::steady_clock::now() - begin).count() / 1e6;
		std::cout << "---Finished. TimeElapsed: " << timeEl << "s" << std::endl << std::endl;
		
		std::cout << std::endl << "---Computing Objective, Constraints, and sensitivities" << std::endl;
		begin = chronos::steady_clock::now();
		// Objective function
		obj = F.transpose()*U;
		objVector(loop) = obj;

		// Objective sensitivity
		for (int e = 0; e < nelx*nely; e++)
		{
			Eigen::VectorXi edofs = edofMat.row(e);
			Eigen::VectorXd Ue(edofs.size());

			for (int i = 0; i < edofs.size(); i++)
			{
				Ue(i) = U(edofs(i));
			}
			
			objSens(e) = -penal*(E0 - Emin)*pow(x(e),(penal - 1))*Ue.transpose()*KE*Ue;
		}

		// Constraint
		g = x.mean();

		// Constraint sensitivity
		gSens = (1.0/(double(nelx*nely)))*Eigen::VectorXd::Ones(nelx*nely, 1);

		// Filtering/Stabilization
		objSens = H*objSens.cwiseQuotient(Hs);
		if (loop >= 1)
		{
			objSens += oldObjSens;
			objSens /= 2;
			oldObjSens = objSens;
		}
		else
		{
			oldObjSens = objSens;
		}
		timeEl = chronos::duration_cast<chronos::microseconds>(chronos::steady_clock::now() - begin).count() / 1e6;
		std::cout << "---Finished. TimeElapsed: " << timeEl << "s" << std::endl << std::endl;

		dobj = std::vector<double>(&objSens[0], objSens.data() + objSens.size());
		dcon = std::vector<double>(&gSens[0], gSens.data() + gSens.size());

		std::cout << std::endl << "---Solving TOBS" << std::endl;
		begin = chronos::steady_clock::now();
		tobs.solve();
		timeEl = chronos::duration_cast<chronos::microseconds>(chronos::steady_clock::now() - begin).count() / 1e6;
		std::cout << "---Finished. TimeElapsed: " << timeEl << "s" << std::endl << std::endl;

		if (loop > 8)
		{
			change = abs((objVector.segment(loop-9,5).sum()-objVector.segment(loop-4,5).sum())/objVector.segment(loop-4,5).sum());
		}

		loop++;

		std::string tempFileName = workingDir + fileName + std::to_string(loop) + extension;
		writeResult(coordMat, elConnect, x, tempFileName);

		std::cout << std::endl;
		std::cout << "It.: " << loop << std::endl;
		std::cout << "Obj value = " << obj << std::endl;
		std::cout << "Volume = " << g << std::endl;
		std::cout << "Change = " << change << std::endl;

		if (loop == maxLoop)
		{
			break;
		}
	}

	timeEl = chronos::duration_cast<chronos::microseconds>(chronos::steady_clock::now() - begin_total).count() / 1e6;
	std::cout << "---Optimization Finished. TimeElapsed: " << timeEl << "s" << std::endl << std::endl;

}