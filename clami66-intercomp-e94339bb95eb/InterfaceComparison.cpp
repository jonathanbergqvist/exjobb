/*
 * InterfaceComparison2.cpp
 *
 *  Created on: Feb 10, 2016
 *      Author: claudio
 */

#include "InterfaceComparison.h"
#include <fstream>
#include <sstream>
#include <vector>

#define FASTRAND_MAX 32767.00

using namespace std;

bool verbose = false;
bool superimpose = false;
bool apply = false;
bool map = false;

double best = -1;

int rounds = 1;
double chainMult = 10.0;

double T0 = 0.2;
double d0 = 0.5;
double epsilon = 0.5;
double dW = 0.5; //DEFAULT VALUES
long double pvalues[30][1001];
char blosum[200] = "BLOSUM62";

//srand (time(NULL));

int nulls;
RunningStat rs;

Molecule* molA;
Molecule* molB;
InterfaceComparison* compare;

static unsigned int g_seed;
int seed = 0;
//Used to seed the generator.

inline void fast_srand(int seed) {

	g_seed = seed;

}

//fastrand routine returns one integer, similar output value range as C lib.

inline int fastrand() {

	g_seed = (214013 * g_seed + 2531011);

	return (g_seed >> 16) & 0x7FFF;

}

void p_values() {
	std::ifstream file("p_values_intsize_shellsize_step5_0.5_0.5_0.5");

	int intl = 0;
	int shel = 0;

	for (int row = 0; row < 30; row++) {
		std::string line;
		std::getline(file, line);
		if (!file.good())
			break;

		std::stringstream iss(line);

		for (int col = 0; col < 1002; ++col) {
			std::string val = "";
			std::getline(iss, val, ' ');
			if (!iss.good())
				break;
			//cout << val << " " << flush;
			std::stringstream convertor(val);

			if (col == 0) {
				convertor >> intl;
			} else {
				convertor >> pvalues[int(intl / 5)][col - 1];
				//cout << pvalues[int(intl/5)][col-1] << " ";
			}
		}
		//cout << "\n" << flush;
	}
	return;
}

double anneal() {

	//shuffle B, reset null correspondances
	molB->shuffle(1000);
	molA->resetNull();

	//calculate the difference between the two (only where A and B overlap and A has no null corrispondences)

	//rs.Clear();

	//cout << "Initial energy: " << score << "\n" << flush;

	int r;
	int c;

	int loop = 1;

	double p = 0;
	double de = 0;


	double T = T0;

	double dice = 0;

	double mean = 0;
	double stdev = 0;

	double accRate = 1;
	int accepted = 1;
	int rejected = 0;

	int maxChainLength = max(400,
			int(chainMult * molB->getLength() * (molB->getLength() - 1))); //2*lengthMin*(lengthMin-1);
	//cout << maxChainLength << "\n";
	int maxAccept = molB->getLength() * (molB->getLength() - 1);

	int toNull = -1;

	for (int i = 0; i < nulls; i++) {

		while (toNull < 0 || molA->getNull()[toNull]) {
			toNull = fastrand() % molA->getLength();
		}
		//if(verbose)
		//cout << "nullA["<< toNull << "] " << "true" << "\n";

		molA->setNull(toNull);
	}
	if (verbose)
		cout << "Caching difference..." << "\n" << flush;
	compare->loadBLOSUM62(blosum);
	compare->cacheDifference();

	if (verbose)
		cout << "Difference cached" << "\n" << flush;

	//a or b selection parameters
	//double ab = 0;

	char AorB = 'A';

	//compare->printBLOSUM62();
	while (1) {

		if (nulls == 0 || loop % 3) {

			//swapping on B
			AorB = 'B';

			r = 0;
			c = 0;

			/*illegal moves are:				(for c, r in [0; lengthB), c > r)
			 1) swap two residues aligned to null correspondences in A (nullA[c] and nullA[r])
			 2) swap a residue in a null correspondence in A with a residue in the outer subset in B (nullA[r] and c >= lengthA)
			 3) swap two residues in the outer subset of B (r >= lengthA, c will be too then)
			 */
			while (r >= c || (molA->getNull()[c] and molA->getNull()[r])
					|| (molA->getNull()[r] and c >= molA->getLength())
					|| r >= molA->getLength()) {
				//while(r >= c ||  molA->getNull()[r] || r >= molA->getLength()){
				r = fastrand() % molB->getLength();
				c = fastrand() % molB->getLength();

			}

			/*if(verbose){
			 cout << "B: " << r << " " << c << " " << molA->getLength() << " " << molA->getNull()[r] << " " << molA->getNull()[c] << "\n" << flush;
			 }*/

		} else {
			//null correspondences on A

			AorB = 'A';

			r = 0;
			c = 0;

			/*illegal moves are:				(for c, r in [0; lengthA), c > r)
			 1) swap two residues aligned to null correspondences in A (nullA[c] and nullA[r])
			 2) swap two non null residues
			 */
			while (r >= c || (molA->getNull()[c] and molA->getNull()[r])
					|| (!(molA->getNull()[c]) and !(molA->getNull()[r]))) {
				r = fastrand() % molA->getLength();
				c = fastrand() % molA->getLength();
			}

		}

		if (AorB == 'A') {

			compare->deltaA(r, c);
			compare->deltaSeqA(r, c);
		}

		if (AorB == 'B') {

			compare->deltaB(r, c);
			compare->deltaSeqB(r, c);

		}

		de = epsilon * -compare->getDeltaN()
								+ dW * (double) compare->getDeltaS();// / compare->getDeltaD();
		rs.Push(compare->getScore() + de);

		if (accepted > maxAccept || accepted + rejected > maxChainLength) {

			//starting a new markov chain
			mean = rs.Mean();
			stdev = rs.StandardDeviation();

			T = compare->temperatureSD(mean, stdev, T);

			accepted = 1;
			rejected = 0;

			rs.Clear();
		}

		if (de < 0) {

			p = 1;

		} else {

			p = exp(-de / T);

		}

		dice = (((double) fastrand() / (FASTRAND_MAX)));

		if (p > dice) {	//probability

			/*if(verbose){
			 cout << "Accepted with p: " << p << "\n";

			 cout << "Delta: " << de <<"\n" << flush;
			 cout << "DeltaS: " <<  -compare->getDeltaS() <<"\n" << flush;
			 cout << "Temp: " << T << "\n" << flush;
			 }*/
			accepted++;

			if (AorB == 'A') {

				molA->swapA(r, c);
				compare->updateDifferenceA2(r, c);

			}

			if (AorB == 'B') {

				molB->swapB(r, c);
				compare->updateDifferenceB2(r, c);
			}

			if (compare->getScore() > best) {

				best = compare->getScore();

				if (verbose)
					cout << "New best: " << best << " deltaN: "
					<< -compare->getDeltaN() << " deltaS: "
					<< compare->getDeltaS() << "\n";

			}

		} else {

			rejected++;

		}

		accRate = (double) accepted / (double) rejected;

		if (accRate < 0.003) {

			if (verbose) {
				molA->printSeq(molA->getLength());
				molB->printSeq(molB->getLength());
			}

			break;

		}
		loop++;
	}

	return best;
}

void run_pdb(int argc, char* argv[]) {

	p_values(); //load the pvalues from file into array

	char to_apply[200];

	double percent = 0;

	for (int arg = 4; arg < argc; arg++) {
		if (!strcmp(argv[arg], "-nullP")) {

			percent = atoi(argv[arg + 1]);


		} else if (!strcmp(argv[arg], "-anneal")) {

			rounds = atoi(argv[arg + 1]);
			if (verbose)
				cout << "Annealing rounds: " << rounds << "\n";

		} else if (!strcmp(argv[arg], "-d0")) {

			if (!strcmp(argv[arg + 1], "auto")) {

				d0 = 0.5;

			} else {

				d0 = atof(argv[arg + 1]);

			}
			if (verbose)
				cout << "d0: " << d0 << "\n";

		} else if (!strcmp(argv[arg], "-v")) {

			cout << "Verbose" << "\n" << flush;
			verbose = true;

		} else if (!strcmp(argv[arg], "-super")) {

			cout << "Superimpose" << "\n" << flush;
			superimpose = true;

		} else if (!strcmp(argv[arg], "-apply")) {

			apply = true;
			cout << "Apply to third molecule" << flush;
			strcpy(to_apply, argv[arg + 1]);
			cout << to_apply << "\n" << flush;

		} else if (!strcmp(argv[arg], "-eps")) {

			epsilon = atof(argv[arg + 1]);
			if (verbose)
				cout << "Epsilon: " << epsilon << "\n" << flush;

		} else if (!strcmp(argv[arg], "-blosum")) {

			strcpy(blosum, argv[arg + 1]);
			if (verbose)
				cout << "Blosum matrix: " << blosum << "\n" << flush;

		} else if (!strcmp(argv[arg], "-dW")) {

			dW = atof(argv[arg + 1]);
			if (verbose)
				cout << "Sequence weight: " << dW << "\n" << flush;

		} else if (!strcmp(argv[arg], "-seed")) {

			seed = atoi(argv[arg + 1]);

			if (verbose)
				cout << "New seed: " << atoi(argv[arg + 1]) << "\n" << flush;

		} else if (!strcmp(argv[arg], "-ch")) {

			chainMult = atof(argv[arg + 1]);
			if (verbose)
				cout << "New chainMult: " << atof(argv[arg + 1]) << "\n"
				<< flush;

		}else if (!strcmp(argv[arg], "-T0")) {

			T0 = atof(argv[arg + 1]);
			if (verbose)
				cout << "New initial Temperature: " << atof(argv[arg + 1]) << "\n"
				<< flush;

		}else if (!strcmp(argv[arg], "-map")) {

			map = true;
			if (verbose)
				cout << "Output mapping of residues\n"
				<< flush;

		}

	}


	srand(seed);
	fast_srand(seed);

	Molecule* mol1 = new Molecule(argv[2], "PDB");
	Molecule* mol2 = new Molecule(argv[3], "PDB");

	int length1 = mol1->getLength();
	int length2 = mol2->getLength();

	nulls = round(double(min(length1, length2)) / 100.0 * percent);
	if (verbose)
		cout << "Max nulls: " << nulls << "\n";

	if (length1 > length2) {

		molA = mol2;
		molB = mol1;

	} else {

		molA = mol1;
		molB = mol2;

	}

	molB->shuffle();

	//load PDBs
	if (verbose) {

		cout << "Mol A length:" << " " << molA->getLength() << "\n" << flush;
		cout << "Mol B length:" << " " << molB->getLength() << "\n" << flush;

	}

	rs.Clear();

	if (verbose)
		cout << "Creating InterFaceComparison object..." << "\n";

	compare = new InterfaceComparison(molA, molB);
	compare->setD0(d0);

	if (verbose)
		cout << "Done." << "\n";

	for (int r = 0; r < rounds; r++) {

		anneal();

	}

	int pvalue_xindex = int(molA->getLength() / 5);

	pvalue_xindex = min(pvalue_xindex, 29);

	int pvalue_zindex = round(best * 1000);

	double best_pvalue = pvalues[pvalue_xindex][pvalue_zindex];

	if(map){
		cout << "Mapping of residue numbers: (Input 1 -> Input 2)\n";
		for(int i = 0; i < molA->getLength(); i++){
			cout << mol1->getResN(i) << " -> " << mol2->getResN(i) << "\n";
		}
	}

	cout << "Best" << ": " << best << " sequence score: "
			<< compare->getsScore() << " length A: "
			<< molA->getLength() - nulls << " length B: " << molB->getLength()
			<< flush;

	if (d0 == 0.5 && dW == 0.5 && epsilon == 0.5) {
		cout << " p-value: " << best_pvalue << "\n" << flush;
	} else {
		cout << "\n" << flush;
	}

	best = 0;

	if (superimpose) {
		//molA->printSeq(molA->getLength());
		//molB->printSeq(molA->getLength());

		FILE *fp1, *fp2, *fp3;
		fp1 = fopen("molA.pdb", "w+");
		fp2 = fopen("molB.pdb", "w+");
		//compare->printDiff();
		compare->superimpose(molA, molB);

		if (apply) {

			Molecule* mol3 = new Molecule(to_apply, "fullPDB");
			fp3 = fopen("molC.pdb", "w+");
			//mol3->center_molecule(molA->getLength());
			double* firstTrans = compare->getTranslation2();

			for (int i = 0; i < 3; i++) {
				firstTrans[i] = -firstTrans[i];
			}

			mol3->rototranslate(compare->getI(), firstTrans);
			mol3->rototranslate(compare->getRotation(),
					compare->getTranslation1());
			mol3->printPDB(fp3, mol3->getLength(), "fullPDB");
			molA->printPDB(fp1, molA->getLength(), "PDB",
					compare->getDiffCumulative());
			molB->printPDB(fp2, molA->getLength(), "PDB",
					compare->getDiffCumulative());
		} else {
			molA->printPDB(fp1, molA->getLength(), "PDB");
			molB->printPDB(fp2, molA->getLength(), "PDB");
		}

		fclose(fp1);
		fclose(fp2);

	}


	mol1->destroy();
	mol2->destroy();

}

std::vector<Molecule*> read_list_of_molecules(char* list){

	std::vector<Molecule*> molecules;

	std::string line;
	std::ifstream infile(list);

	while (std::getline(infile, line)){
		molecules.push_back(new Molecule(const_cast<char*>(line.c_str()), "PDB"));
	}

	return molecules;

}

void run_pdb_list(int argc, char* argv[]) {

	p_values(); //load the pvalues from file into array

	char to_apply[200];

	double percent = 0;

	for (int arg = 4; arg < argc; arg++) {
		if (!strcmp(argv[arg], "-nullP")) {

			percent = atoi(argv[arg + 1]);


		} else if (!strcmp(argv[arg], "-anneal")) {

			rounds = atoi(argv[arg + 1]);
			if (verbose)
				cout << "Annealing rounds: " << rounds << "\n";

		} else if (!strcmp(argv[arg], "-d0")) {

			if (!strcmp(argv[arg + 1], "auto")) {

				d0 = 0.5;

			} else {

				d0 = atof(argv[arg + 1]);

			}
			if (verbose)
				cout << "d0: " << d0 << "\n";

		} else if (!strcmp(argv[arg], "-v")) {

			cout << "Verbose" << "\n" << flush;
			verbose = true;

		} else if (!strcmp(argv[arg], "-super")) {

			cout << "Superimpose" << "\n" << flush;
			superimpose = true;

		} else if (!strcmp(argv[arg], "-apply")) {

			apply = true;
			cout << "Apply to third molecule" << flush;
			strcpy(to_apply, argv[arg + 1]);
			cout << to_apply << "\n" << flush;

		} else if (!strcmp(argv[arg], "-eps")) {

			epsilon = atof(argv[arg + 1]);
			if (verbose)
				cout << "Epsilon: " << epsilon << "\n" << flush;

		} else if (!strcmp(argv[arg], "-blosum")) {

			strcpy(blosum, argv[arg + 1]);
			if (verbose)
				cout << "Blosum matrix: " << blosum << "\n" << flush;

		} else if (!strcmp(argv[arg], "-dW")) {

			dW = atof(argv[arg + 1]);
			if (verbose)
				cout << "Sequence weight: " << dW << "\n" << flush;

		} else if (!strcmp(argv[arg], "-seed")) {

			seed = atoi(argv[arg + 1]);

			if (verbose)
				cout << "New seed: " << atoi(argv[arg + 1]) << "\n" << flush;

		} else if (!strcmp(argv[arg], "-ch")) {

			chainMult = atof(argv[arg + 1]);
			if (verbose)
				cout << "New chainMult: " << atof(argv[arg + 1]) << "\n"
				<< flush;

		}else if (!strcmp(argv[arg], "-T0")) {

			T0 = atof(argv[arg + 1]);
			if (verbose)
				cout << "New initial Temperature: " << atof(argv[arg + 1]) << "\n"
				<< flush;

		}

	}

	Molecule* mol1 = new Molecule(argv[2], "PDB");
	int length1 = mol1->getLength();

	std::vector<Molecule*> mols2 = read_list_of_molecules(argv[3]);
	for(std::vector<Molecule*>::iterator mol2 = mols2.begin(); mol2 != mols2.end(); ++mol2){

		srand(seed);
		fast_srand(seed);
		//Molecule* mol2 = new Molecule(argv[3], "PDB");
		int length2 = (*mol2)->getLength();

		nulls = round(double(min(length1, length2)) / 100.0 * percent);

		if (length1 > length2) {

			molA = (*mol2);
			molB = mol1;

		} else {

			molA = mol1;
			molB = (*mol2);
			molB->shuffle();

		}


		//load PDBs
		if (verbose) {

			cout << "Mol A length:" << " " << molA->getLength() << "\n" << flush;
			cout << "Mol B length:" << " " << molB->getLength() << "\n" << flush;

		}

		rs.Clear();

		if (verbose)
			cout << "Creating InterFaceComparison object..." << "\n";

		compare = new InterfaceComparison(molA, molB);
		compare->setD0(d0);

		if (verbose)
			cout << "Done." << "\n";

		for (int r = 0; r < rounds; r++) {

			anneal();

		}

		int pvalue_xindex = int(molA->getLength() / 5);

		pvalue_xindex = min(pvalue_xindex, 29);

		int pvalue_zindex = round(best * 1000);

		double best_pvalue = pvalues[pvalue_xindex][pvalue_zindex];

		cout << "Best" << ": " << best << " sequence score: "
				<< compare->getsScore() << " length A: "
				<< molA->getLength() - nulls << " length B: " << molB->getLength()
				<< flush;

		if (d0 == 0.5 && dW == 0.5 && epsilon == 0.5) {
			cout << " p-value: " << best_pvalue << "\n" << flush;
		} else {
			cout << "\n" << flush;
		}

		best = 0;

		if (superimpose) {
			//molA->printSeq(molA->getLength());
			//molB->printSeq(molA->getLength());

			FILE *fp1, *fp2, *fp3;
			fp1 = fopen("molA.pdb", "w+");
			fp2 = fopen("molB.pdb", "w+");
			//compare->printDiff();
			compare->superimpose(molA, molB);

			if (apply) {

				Molecule* mol3 = new Molecule(to_apply, "fullPDB");
				fp3 = fopen("molC.pdb", "w+");
				//mol3->center_molecule(molA->getLength());
				double* firstTrans = compare->getTranslation2();

				for (int i = 0; i < 3; i++) {
					firstTrans[i] = -firstTrans[i];
				}

				mol3->rototranslate(compare->getI(), firstTrans);
				mol3->rototranslate(compare->getRotation(),
						compare->getTranslation1());
				mol3->printPDB(fp3, mol3->getLength(), "fullPDB");
				molA->printPDB(fp1, molA->getLength(), "PDB",
						compare->getDiffCumulative());
				molB->printPDB(fp2, molA->getLength(), "PDB",
						compare->getDiffCumulative());
			} else {
				molA->printPDB(fp1, molA->getLength(), "PDB");
				molB->printPDB(fp2, molA->getLength(), "PDB");
			}

			fclose(fp1);
			fclose(fp2);

		}


		//mol1->destroy();
		(*mol2)->destroy();
	}
}

void run_mol(int argc, char* argv[]) {

	//according to algorithm with null correspondences, B molecule is always the largest of the two
	Molecule* mol1 = new Molecule(argv[2], "mol");

	int length1, length2;

	length1 = mol1->getLength();

	ifstream fin(argv[3]); //to count stuff
	ifstream fin2(argv[3]);	//to read stuff

	string line;
	string molecule("@<TRIPOS>MOLECULE");
	int nMolecules = 0;

	getline(fin, line);

	while (fin.good()) {

		if (line.compare(0, molecule.length(), molecule) == 0)
			nMolecules++;

		getline(fin, line);
	}

	if (verbose)
		cout << "reading " << nMolecules << " molecules...\n";

	fin.clear();
	fin.seekg(0, std::ios::beg);

	double percent = 0;

	for (int arg = 4; arg < argc; arg++) {

		if (!strcmp(argv[arg], "-nullP")) {

			percent = atoi(argv[arg + 1]);

		} else if (!strcmp(argv[arg], "-anneal")) {

			rounds = atoi(argv[arg + 1]);
			if (verbose)
				cout << "Annealing rounds: " << rounds << "\n";

		} else if (!strcmp(argv[arg], "-d0")) {

			if (!strcmp(argv[arg + 1], "auto")) {

				d0 = 0.0002745 * pow(min(length1, length2), 1.736);

			} else {

				d0 = atof(argv[arg + 1]);

			}
			if (verbose)
				cout << "d0: " << d0 << "\n";

		} else if (!strcmp(argv[arg], "-v")) {

			cout << "Verbose" << "\n" << flush;
			verbose = true;

		} else if (!strcmp(argv[arg], "-eps")) {

			epsilon = atof(argv[arg + 1]);
			cout << "Epsilon: " << epsilon << "\n" << flush;

		} else if (!strcmp(argv[arg], "-seed")) {

			seed = atoi(argv[arg + 1]);
			cout << "New seed: " << atoi(argv[arg + 1]) << "\n" << flush;

		} else if (!strcmp(argv[arg], "-ch")) {

			chainMult = atof(argv[arg + 1]);
			cout << "New chainMult: " << atoi(argv[arg + 1]) << "\n" << flush;

		}

	}

	//Molecule** mol2 = new Molecule*[nMolecules];

	Molecule* mol2;
	for (int i = 0; i < nMolecules; i++) {

		srand(seed);
		fast_srand(seed);
		mol2 = new Molecule(&fin, &fin2);

		length2 = mol2->getLength();

		nulls = round(double(min(length1, length2)) / 100.0 * percent);

		if (verbose)
			cout << "Nulls: " << nulls << "\n";

		if (length1 > length2) {

			molA = mol2;
			molB = mol1;

		} else {

			molA = mol1;
			molB = mol2;

		}

		molB->shuffle();

		//load PDBs
		if (verbose) {

			cout << "Mol A length:" << " " << molA->getLength() << "\n"
					<< flush;
			cout << "Mol B length:" << " " << molB->getLength() << "\n"
					<< flush;

		}

		rs.Clear();

		if (verbose)
			cout << "Creating InterFaceComparison object..." << "\n";

		compare = new InterfaceComparison(molA, molB);
		compare->setD0(d0);

		if (verbose)
			cout << "Done." << "\n";

		for (int r = 0; r < rounds; r++) {

			anneal();

		}

		cout << "Best for molecule number " << i << ": " << best << "\n"
				<< flush;
		best = 0;
		mol2->destroy();

	}

}

int main(int argc, char *argv[]) {

	fast_srand(time(NULL));
	srand(time(NULL));

	char type[4];

	if (!strcmp(argv[1], "-PDB")) {

		strcpy(type, "PDB");
		run_pdb(argc, argv);

	}else if (!strcmp(argv[1], "-list")) {

		strcpy(type, "lis");
		run_pdb_list(argc, argv);

	} else if (!strcmp(argv[1], "-mol2")) {

		strcpy(type, "mol");
		run_mol(argc, argv);

	} else {

		cout
		<< "Usage: InterfaceComparison [-PDB,-mol2] target1 target2 <options>\n";
		return EXIT_FAILURE;
	}

}

