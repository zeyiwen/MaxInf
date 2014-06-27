#include "BIP/BIPSolver.h"
#include "EEP/EEPSolver.h"
#include "EEP/NASolver.h"

void main()
{
	string CustomerFile = FILE_M;
	string FacilityFile = FILE_F;
	string CandidateFile = FILE_C;
	int nNumofCustomer = FILE_M_SIZE;
	int nNumofFacility = FILE_F_SIZE;
	int nNumofCandidate = FILE_C_SIZE;

	BIPSolver * bipSolver;
	EEPSolver * eepSolver;
	NASolver * naSolver;

	bipSolver = new BIPSolver(CustomerFile, FacilityFile, CandidateFile, nNumofCustomer, nNumofFacility,nNumofCandidate);
	bipSolver->Solve();
	bipSolver->Print();
	bipSolver->PrintDuration();
	bipSolver->PrintComputation();
	cout<< "BIP end \n" << endl;

	eepSolver = new EEPSolver(CustomerFile, FacilityFile, CandidateFile, nNumofCustomer, nNumofFacility,nNumofCandidate);
	eepSolver->Solve();
	eepSolver->PrintAnswerSet();
	eepSolver->PrintDuration();
	eepSolver->PrintDistComputation();
	cout <<  "EEP end \n" << endl;

	naSolver = new NASolver(CustomerFile, FacilityFile, CandidateFile, nNumofCustomer, nNumofFacility,nNumofCandidate);
	naSolver->Solve();
	naSolver->PrintAnswerSet();
	naSolver->PrintDuration();
	naSolver->PrintComputation();
	cout << "NA end \n" << endl;

}