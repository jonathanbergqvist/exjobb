/*
 * Molecule.h
 *
 *  Created on: Feb 10, 2016
 *      Author: claudio
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string>

#include <malloc.h>
#include <sstream>
#include <iostream>
#include <fstream>


#include <algorithm>
#include <cmath>        // std::abs

using namespace std;

class Molecule{

	private:

		char AAmap(string AA);
		char Emap(string AA);
		string AAmap2(char AA);
		void get_xyz(string line, double *x, double *y, double *z, char *resname, char *type);
		void get_xyz(string line, double *x, double *y, double *z, char *resname, char *chain, int *resnumber, char *atomname, char *type);
		void RemoveSpaces(char* source);
		int readPDB(char *filename, char *type);
		int readMol(ifstream *fin2);

		//private sets
		void setLength(char *filename, char const *type);
		void setLength(ifstream* fin);

		void setDist();
		void setDistSqrt();
		void setDiff();

		int length;

		double **coord;		//coordinates[length][3]
		int    *resnumber;
		char  **atomname;
		char   *seq;		//protein sequences
		char   *chain;		//protein sequences

		bool *null;

		double **dist;  //distance matrix

		//temporary stuff
		double* tmp;// = new double[size];
		double* tmpCoord;// = new double[3];

	public:

		Molecule(char* filename, char* type);
		Molecule(ifstream* fin, ifstream *fin2);

		void destroy();

		int getLength();

		char* getSeq();
		double** getDist();
		double* getCoord(int i);
		int getResN(int i);
		void setCoord(int i, int j, double newCoord);
		bool* getNull();

		void printSeq(int size);
		//void Molecule::printSeq(int size);
		void printSeq2();
		void printDist();
		void printPDB(FILE *name, int length=0, char *type="PDB");
		void printPDB(FILE *name, int length, char *type, double* distances);

		void swapA(int i, int j);
		void swapB(int i, int j);

		void resetNull();
		void setNull(int i);
		void shuffle();
		void shuffle(int rounds);
        void rototranslate(double** rot, double* trans);
        double* center_molecule(int length, double* distances);


};

Molecule::Molecule(char* filename, char* type){

	setLength(filename, type);

    coord = new double*[this->length];
    resnumber = new int[this->length];
    atomname = new char*[this->length];

	for(int i=0; i<this->length; i++){
		coord[i] = new double[3];
		resnumber[i] = 0;
		atomname[i] = new char[5];
	}

    this->seq   = new char[this->length+1];
    this->chain   = new char[this->length+1];
	null  = new bool[this->length];

	seq[this->length] = '\0';
	chain[this->length] = '\0';
	//load PDB coordinates
	this->readPDB(filename, type);


    dist = new double*[this->length];


	for(int x=0; x<this->length; x++){

		dist[x] = new double[this->length];
		null[x] = false;
	}

	if(!strcmp(type, "PDB") || !strcmp(type, "fullPDB")){
		this->setDistSqrt();
	}else{
		this->setDist();
	}

	tmpCoord = new double[3];
	tmp = new double[this->length];

}

Molecule::Molecule(ifstream* fin, ifstream* fin2){

	setLength(fin);

    coord = new double*[this->length];

	for(int i=0; i<this->length; i++){
		coord[i] = new double[3];
	}

    this->seq   = new char[this->length+1];
	null  = new bool[this->length];

	seq[this->length] = '\0';

	//load PDB coordinates
	this->readMol(fin2);


    dist = new double*[this->length];


	for(int x=0; x<this->length; x++){

		dist[x] = new double[this->length];
		null[x] = false;
	}

	this->setDist();

	tmpCoord = new double[3];
	tmp = new double[this->length];

}

void Molecule::rototranslate(double** rot, double* trans){
 
	for(int i=0; i<this->length; i++){

            double x = rot[0][0] * this->getCoord(i)[0] + rot[0][1] * this->getCoord(i)[1]  + rot[0][2] * this->getCoord(i)[2] + trans[0];
            double y = rot[1][0] * this->getCoord(i)[0] + rot[1][1] * this->getCoord(i)[1]	+ rot[1][2] * this->getCoord(i)[2] + trans[1];
            double z = rot[2][0] * this->getCoord(i)[0] + rot[2][1] * this->getCoord(i)[1]	+ rot[2][2] * this->getCoord(i)[2] + trans[2];

            this->setCoord(i, 0, x);
            this->setCoord(i, 1, y);
            this->setCoord(i, 2, z);
        
            
        }
}


double* Molecule::center_molecule(int length, double* distances){
    
  int	i;	/* Counter variables */
  int   natoms;  // Number of selected atoms
  double	xcen,ycen,zcen;		/* Temporary coordinates values */
  double* vec = new double[3];
  
  natoms=0;
  xcen=ycen=zcen=0;
  for (i=0; i<length; i++){
    //if (m->atm[i].selected){
	  if(distances[i] < 5){
      xcen+=this->getCoord(i)[0];
      ycen+=this->getCoord(i)[1];
      zcen+=this->getCoord(i)[2];
      natoms++;
	  }
    //}
  }
  /* Now center molecule */
  xcen/=(double)natoms;
  ycen/=(double)natoms;
  zcen/=(double)natoms;
  
  vec[0]=xcen;
  vec[1]=ycen;
  vec[2]=zcen;
  
  //  printf("TEST CEN natoms: %d %f %f %f \n",natoms,xcen,ycen,zcen);
  for (i=0;i<this->getLength();i++)
    {

      this->setCoord(i, 0, this->getCoord(i)[0]-xcen);
      this->setCoord(i, 1, this->getCoord(i)[1]-ycen);
      this->setCoord(i, 2, this->getCoord(i)[2]-zcen);

      // printf("TEST natoms: %d %f %f %f \n",i,m->atm[i].x,m->atm[i].y,m->atm[i].z);
    }
  //i--;
  //printf("TEST natoms: %d %f %f %f \n",i,m->atm[i].x,m->atm[i].y,m->atm[i].z);
  return vec;
}


char Molecule::AAmap(string AA){

    char A=' ';
    if(     AA.compare("BCK")==0)   A='X';
    else if(AA.compare("GLY")==0)   A='G';
    else if(AA.compare("ALA")==0)   A='A';
    else if(AA.compare("SER")==0)   A='S';
    else if(AA.compare("CYS")==0)   A='C';
    else if(AA.compare("VAL")==0)   A='V';
    else if(AA.compare("THR")==0)   A='T';
    else if(AA.compare("ILE")==0)   A='I';
    else if(AA.compare("PRO")==0)   A='P';
    else if(AA.compare("MET")==0)   A='M';
    else if(AA.compare("ASP")==0)   A='D';
    else if(AA.compare("ASN")==0)   A='N';
    else if(AA.compare("LEU")==0)   A='L';
    else if(AA.compare("LYS")==0)   A='K';
    else if(AA.compare("GLU")==0)   A='E';
    else if(AA.compare("GLN")==0)   A='Q';
    else if(AA.compare("ARG")==0)   A='R';
    else if(AA.compare("HIS")==0)   A='H';
    else if(AA.compare("PHE")==0)   A='F';
    else if(AA.compare("TYR")==0)   A='Y';
    else if(AA.compare("TRP")==0)   A='W';
    else if(AA.compare("CYX")==0)   A='C';
    else A='Z';

    return A;
}

char Molecule::Emap(string AA){

    char A=' ';
    if(     AA.compare("Br")==0)   A='A';
    else if(AA.compare("C.1")==0)   A='B';
    else if(AA.compare("C.2")==0)   A='C';
    else if(AA.compare("C.3")==0)   A='D';
    else if(AA.compare("C.ar")==0)   A='E';
    else if(AA.compare("Cl")==0)   A='F';
    else if(AA.compare("F")==0)   A='G';
    else if(AA.compare("F.3")==0)   A='H';
    else if(AA.compare("H")==0)   A='I';
    else if(AA.compare("N.1")==0)   A='L';
    else if(AA.compare("N.2")==0)   A='M';
    else if(AA.compare("N.3")==0)   A='N';
    else if(AA.compare("N.4")==0)   A='O';
    else if(AA.compare("N.am")==0)   A='P';
    else if(AA.compare("N.ar")==0)   A='Q';
    else if(AA.compare("N.pl3")==0)   A='R';
    else if(AA.compare("O.2")==0)   A='S';
    else if(AA.compare("O.3")==0)   A='T';
    else if(AA.compare("O.co2")==0)   A='U';
    else if(AA.compare("S.2")==0)   A='V';
    else if(AA.compare("S.3")==0)   A='W';
    else if(AA.compare("S.02")==0)   A='X';
    else A='Z';

    return A;
}

string Molecule::AAmap2(char AA){

    string A=" ";

    if(     AA=='X')                A="BCK";
    else if(AA=='G')                A="GLY";
    else if(AA=='A')                A="ALA";
    else if(AA=='S')                A="SER";
    else if(AA=='C')                A="CYS";
    else if(AA=='V')                A="VAL";
    else if(AA=='T')                A="THR";
    else if(AA=='I')                A="ILE";
    else if(AA=='P')                A="PRO";
    else if(AA=='M')                A="MET";
    else if(AA=='D')                A="ASP";
    else if(AA=='N')                A="ASN";
    else if(AA=='L')                A="LEU";
    else if(AA=='K')                A="LYS";
    else if(AA=='E')                A="GLU";
    else if(AA=='Q')                A="GLN";
    else if(AA=='R')                A="ARG";
    else if(AA=='H')                A="HIS";
    else if(AA=='F')                A="PHE";
    else if(AA=='Y')                A="TYR";
    else if(AA=='W')                A="TRP";
    else if(AA=='C')                A="CYX";
    else A="UNK";

    return A;
}

void Molecule::get_xyz(string line, double *x, double *y, double *z, char *resname, char *type){

    char cstr[50];

	if(!strcmp(type, "PDB") || !strcmp(type, "fullPDB")){
		strcpy(cstr, (line.substr(30, 8)).c_str());
		sscanf(cstr, "%lf", x);

		strcpy(cstr, (line.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", y);

		strcpy(cstr, (line.substr(46, 8)).c_str());
		sscanf(cstr, "%lf", z);

		strcpy(cstr, (line.substr(17, 3)).c_str());
		*resname=AAmap(cstr);

	}else if(!strcmp(type, "mol")){


		strcpy(cstr, (line.substr(18, 8)).c_str());
		sscanf(cstr, "%lf", x);

		strcpy(cstr, (line.substr(28, 8)).c_str());
		sscanf(cstr, "%lf", y);

		strcpy(cstr, (line.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", z);

		strcpy(cstr, (line.substr(47, 6)).c_str());
		RemoveSpaces(cstr);
		*resname=Emap(cstr);


	}
    //strcpy(cstr, (line.substr(22, 4)).c_str());
    //sscanf(cstr, "%d", no);
}

void Molecule::get_xyz(string line, double *x, double *y, double *z, char *resname, char *chain, int* resnumber,  char *atomname, char *type){

    char cstr[50];

	if(!strcmp(type, "PDB") || !strcmp(type, "fullPDB")){
		strcpy(cstr, (line.substr(30, 8)).c_str());
		sscanf(cstr, "%lf", x);

		strcpy(cstr, (line.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", y);

		strcpy(cstr, (line.substr(46, 8)).c_str());
		sscanf(cstr, "%lf", z);

		strcpy(cstr, (line.substr(17, 3)).c_str());
		*resname=AAmap(cstr);

		strcpy(cstr, (line.substr(22, 4)).c_str());
		sscanf(cstr, "%d", resnumber);

	}if(!strcmp(type, "fullPDB")){

		strcpy(cstr, (line.substr(30, 8)).c_str());
		sscanf(cstr, "%lf", x);

		strcpy(cstr, (line.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", y);

		strcpy(cstr, (line.substr(46, 8)).c_str());
		sscanf(cstr, "%lf", z);

		strcpy(cstr, (line.substr(17, 3)).c_str());
		*resname=AAmap(cstr);

		strcpy(cstr, (line.substr(22, 4)).c_str());
		sscanf(cstr, "%d", resnumber);

		strcpy(atomname, (line.substr(12, 4)).c_str());

		*chain = (line[21]);

	}else if(!strcmp(type, "mol")){


		strcpy(cstr, (line.substr(18, 8)).c_str());
		sscanf(cstr, "%lf", x);

		strcpy(cstr, (line.substr(28, 8)).c_str());
		sscanf(cstr, "%lf", y);

		strcpy(cstr, (line.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", z);

		strcpy(cstr, (line.substr(47, 6)).c_str());
		RemoveSpaces(cstr);
		*resname=Emap(cstr);


	}
    //strcpy(cstr, (line.substr(22, 4)).c_str());
    //sscanf(cstr, "%d", no);
}

void Molecule::RemoveSpaces(char* source){
  char* i = source;
  char* j = source;
  while(*j != 0)
  {
    *i = *j++;
    if(*i != ' ')
      i++;
  }
  *i = 0;
}

void Molecule::setLength(char *filename, char const *type){

    int i=0;
    string line="";



    ifstream fin(filename);

	if(!strcmp(type, "PDB")){

		string atom ("ATOM ");
		if (fin.is_open()){

		        while ( fin.good() ){

		            getline(fin, line);

		            if(line.compare(0, atom.length(), atom)==0){

		                if( line.compare(12, 4, "CA  ")==0 ||\
		                    line.compare(12, 4, " CA ")==0 ||\
		                    line.compare(12, 4, "  CA")==0 ||\
		                    line.compare(12, 4, "C1' ")==0 ||\
		                    line.compare(12, 4, " C1'")==0 ){

			                    if( line.compare(16, 1, " ")==0 ||\
	                        		line.compare(16, 1, "A")==0 ){

				                        i++;
	                    		    }
		                }
		            }
		        }
		        fin.close();
		    }else{
			char error[5000];
			sprintf(error, "Can not open file: %s", filename);
			cerr << error << endl;
			exit(1);

		    }
		}else if(!strcmp(type, "fullPDB")){

		string atom ("ATOM ");
		if (fin.is_open()){

		        while ( fin.good() ){

		            getline(fin, line);

		            if(line.compare(0, atom.length(), atom)==0){
		               
			                    if( line.compare(16, 1, " ")==0 ||\
	                        		line.compare(16, 1, "A")==0 ){

				                        i++;
	                    		    }
		            }
		        }
		        fin.close();
		    }else{
			char error[5000];
			sprintf(error, "Can not open file: %s", filename);
			cerr << error << endl;
			exit(1);

		    }
		}else if(!strcmp(type, "mol")){

			string atom("@<TRIPOS>ATOM");
			string bond("@<TRIPOS>BOND");

			if (fin.is_open()){

				while (fin.good()){

				        getline(fin, line);

					if(line.compare(0, atom.length(), atom)==0){

					        getline(fin, line);

						while(line.compare(0, bond.length(), bond)){

							if(line.at(47) != 'H'){
								i++;
							}
						        getline(fin, line);


						}
	           		 }
	        	}
	        fin.close();
	    }else{
		char error[5000];
		sprintf(error, "Can not open file: %s", filename);
		cerr << error << endl;
		exit(1);

	    }

	}

	if(i==0){
		char error[5000];
		sprintf(error, "Can not find CAs in file: %s", filename);
		cerr << error << endl;
		exit(1);

	}

    this->length = i;

}

void Molecule::setLength(ifstream* fin){

    int i=0;
    string line="";



    //ifstream fin(filename);

    //int place = 0;

    string molecule("@<TRIPOS>MOLECULE");
	string atom ("@<TRIPOS>ATOM");
	string bond ("@<TRIPOS>BOND");

	if ((*fin).is_open()){

		//place = (*fin).tellg();

		while (line.compare(0, bond.length(), bond) && (*fin).good()){

			getline(*fin, line);
			//cout << "in the while loop " << line << "\n";

			if(line.compare(0, atom.length(), atom)==0){

				getline(*fin, line);

				while(line.compare(0, bond.length(), bond)){

					if(line.at(47) != 'H'){
						i++;
					}

					getline(*fin, line);


				}
	        }
	   }
		//go back to where it started
		//(*fin).seekg(place, ios_base::cur);
		//fin.close();

	    }else{

	    	char error[5000];
	    	sprintf(error, "Read error");
	    	cerr << error << endl;

	    }



	if(i==0){

		char error[5000];
		sprintf(error, "Read error, size zero");
		cerr << error << endl;

	}
	cout << "current molecule size: " << i << "\n";
    this->length = i;

}

int Molecule::readPDB(char *filename, char *type){
    int i=0;
    string line, str;



    ifstream fin (filename);

	if(!strcmp(type, "PDB")){

	    //cout << "Reading PDB format file...\n";
	    string atom ("ATOM ");

	    if (fin.is_open()){

	        while ( fin.good() ){
	            getline(fin, line);

	            if(line.compare(0, atom.length(), atom)==0){

	                if( line.compare(12, 4, "CA  ")==0 ||\
	                    line.compare(12, 4, " CA ")==0 ||\
	                    line.compare(12, 4, "  CA")==0){

	                    if( line.compare(16, 1, " ")==0 ||\
	                        line.compare(16, 1, "A")==0 ){

		                        //get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i], type);
	                    		get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i], &chain[i], &resnumber[i], atomname[i], type);
	                        	i++;
	                    }
	                }
	            }
	        }
	        fin.close();
	    }
	    else{
		char error[5000];
		sprintf(error, "Can not open file: %s\n", filename);
		cerr << error << endl;
		exit(1);
    	}
    }else if(!strcmp(type, "fullPDB")){

	    //cout << "Reading PDB format file...\n";
	    string atom ("ATOM ");

	    if (fin.is_open()){

	        while ( fin.good() ){
	            getline(fin, line);

	            if(line.compare(0, atom.length(), atom)==0){           

	                    if( line.compare(16, 1, " ")==0 ||\
	                        line.compare(16, 1, "A")==0 ){

		                        get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i], &chain[i], &resnumber[i], atomname[i], type);

	                        	i++;
	                    }
	            }
	        }
	        fin.close();
	    }
	    else{
		char error[5000];
		sprintf(error, "Can not open file: %s\n", filename);
		cerr << error << endl;
		exit(1);
    	}
    }else if(!strcmp(type, "mol")){

			string atom ("@<TRIPOS>ATOM");
			string bond ("@<TRIPOS>BOND");
			cout << "Readingt mol2 format file...\n";
			if (fin.is_open()){

				while (fin.good()){

				        getline(fin, line);

					if(line.compare(0, atom.length(), atom)==0){
						printf("%s\n", line.c_str());
					        getline(fin, line);
						printf("%s\n", line.c_str());
						while(line.compare(0, bond.length(), bond)){

							if(line.at(47) != 'H'){
								//cout << line.at(47) << "\n";
				                        	get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i], type);
								i++;
							}
						        getline(fin, line);

						}
	           		 }
	        	}
	        fin.close();


	    }
	}

    return 0;
}


int Molecule::readMol(ifstream* fin){
    int i=0;
    string line, str;

    getline(*fin, line);

    char type[] = "mol";
   //%ifstream fin (filename);

    string molecule("@<TRIPOS>MOLECULE");
	string atom("@<TRIPOS>ATOM");
	string bond("@<TRIPOS>BOND");

	cout << "Readingt mol2 format file...\n";

	if ((*fin).is_open()){

		getline(*fin, line);
		//cout << line << "\n";

		while (line.compare(0, bond.length(), bond)!=0 && (*fin).good()){//(*fin).good()){

		        getline(*fin, line);
//cout << line << endl;
				if(line.compare(0, atom.length(), atom)==0){

						//printf("%s\n", line.c_str());

						getline(*fin, line);
						//printf("%s\n", line.c_str());

						while(line.compare(0, bond.length(), bond)){

							if(line.at(47) != 'H'){
								//cout << line.at(47) << "\n";
				                        	get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i], type);
								i++;
							}
						        getline(*fin, line);

						}
	           		 }
	        	}
	        //(*fin).close();


	    }


    return 0;
}




//fills the distance matrices the first time
void Molecule::setDist(){



	for(int i=0; i<this->length; i++){

		for(int j=i; j<this->length; j++){

				if(coord[i][0] != 0 && coord[j][0] != 0 && i!=j){
					this->dist[i][j] = this->dist[j][i] = (coord[i][0]-coord[j][0])*(coord[i][0]-coord[j][0])+(coord[i][1]-coord[j][1])*(coord[i][1]-coord[j][1])+(coord[i][2]-coord[j][2])*(coord[i][2]-coord[j][2]);
				}else{
					this->dist[i][j] = this->dist[j][i] = 0;
				}

		}


	}

}

//fills the distance matrices the first time
void Molecule::setDistSqrt(){



	for(int i=0; i<this->length; i++){

		for(int j=i; j<this->length; j++){

				if(coord[i][0] != 0 && coord[j][0] != 0 && i!=j){
					this->dist[i][j] = this->dist[j][i] = sqrt((coord[i][0]-coord[j][0])*(coord[i][0]-coord[j][0])+(coord[i][1]-coord[j][1])*(coord[i][1]-coord[j][1])+(coord[i][2]-coord[j][2])*(coord[i][2]-coord[j][2]));
				}else{
					this->dist[i][j] = this->dist[j][i] = 0;
				}

		}


	}

}

int Molecule::getLength(){

	return this->length;

}

char* Molecule::getSeq(){

	return this->seq;

}

double** Molecule::getDist(){

	return this->dist;

}

double* Molecule::getCoord(int i){

	return this->coord[i];

}

void Molecule::setCoord(int i, int j, double newCoord){

	this->coord[i][j] = newCoord;

}


int Molecule::getResN(int i){

	return this->resnumber[i];

}


bool* Molecule::getNull(){

	return this->null;

}


//swapping A is much simpler, it is only about setting a new null point
void Molecule::swapA(int c, int r){

	//int toTrue = toFalse = 0;
	//
	//while(nullA[toTrue] != false && nullA[toFalse] != true){
	//	toTrue  = rand() % lengthA;
	//	toFalse = rand() % lengthA;
	//}

	if(!null[r] and null[c]){
		null[r] = true;
		null[c]= false;
	}else if(null[r] and !null[c]){
		null[r] = false;
		null[c] = true;
	}else{
		cout << "Something's wrong in swapA\n";
	}
}

//swap two columns/rows on the distance matrix
//(use only when you have accepted the move)
void Molecule::swapB(int i, int j){


	char tmpC;
	int m;

	//swap coordinates

	//this->tmpCoord = this->coord[i];
	//this->tmpCoord[0] = this->coord[i][0];
	//this->tmpCoord[1] = this->coord[i][1];
	//this->tmpCoord[2] = this->coord[i][2];

	//this->coord[i] = this->coord[j];
	//this->coord[i][0] = this->coord[j][0];
	//this->coord[i][1] = this->coord[j][1];
	//this->coord[i][2] = this->coord[j][2];

	//this->coord[j] = this->tmpCoord;

	std::swap(this->coord[i], this->coord[j]);
	//this->coord[j][0] = this->tmpCoord[0];
	//this->coord[j][1] = this->tmpCoord[1];
	//this->coord[j][2] = this->tmpCoord[2];

	//normal D swapping, no need to recalculate distances
	//swap rows in the distance matrix
	std::swap(this->dist[i], this->dist[j]);

	std::swap(this->dist[i][j], this->dist[j][j]);
	std::swap(this->dist[j][i], this->dist[i][i]);

	for(m=0; m<this->length; m++){

		//this->tmp[m] = this->dist[i][m];

		if(m!=i && m!=j){

			//this->dist[i][m] = this->dist[m][i] = this->dist[j][m];
			//this->dist[j][m] = this->dist[m][j] = this->tmp[m];

			this->dist[m][i] = this->dist[i][m];
			this->dist[m][j] = this->dist[j][m];
		}
	}

	//swap letters in the sequence
	std::swap(this->seq[i], this->seq[j]);
	std::swap(this->resnumber[i], this->resnumber[j]);
	//tmpC   = this->seq[i];
	//this->seq[i] = this->seq[j];
	//this->seq[j] = tmpC;


}

void Molecule::shuffle(){

	//cout << "Shuffling... ";

	int i = 0;
	int j = 0;

	for(int r = 0; r<10000; r++){

		do{
			i = rand() % this->length;
			j = rand() % this->length;
		}while(i == j);

		//cout << "Swapping " << i << " " << j << "\n";
		this->swapB(i, j);

	}


	//cout << "Done.\n";

}

void Molecule::shuffle(int rounds){

	//cout << "Shuffling... ";

	int i = 0;
	int j = 0;

	for(int r = 0; r<rounds; r++){

		do{
			i = rand() % this->length;
			j = rand() % this->length;
		}while(i == j);

		//cout << "Swapping " << i << " " << j << "\n";
		this->swapB(i, j);

	}


	//cout << "Done.\n";

}

void Molecule::resetNull(){


	for(int i = 0; i< this->length; i++){

		this->null[i] = false;

	}

}

void Molecule::setNull(int i){

	this->null[i] = true;

}

/*void Molecule::printSeq(){

	cout << "> Molecule";

	cout << "\n";

	for(int a=0; a<this->length; a++){
		if(!this->null[a]){
			cout << seq[a] << "";
		}else{
			cout << "-" << "";
		}
	}

	cout << "\n";
}*/

void Molecule::printSeq(int size){

	cout << "> Molecule";

	cout << "\n";

	for(int a=0; a<size; a++){
		if(!this->null[a]){
			cout << seq[a] << "";
		}else{
			cout << "-" << "";
		}
	}

	cout << "\n";
}

void Molecule::printSeq2(){

	for(int a=0; a<this->length; a++){
		if(this->null[a]){
			cout << " ";
		}else{
			cout << "-";
		}
	}
	cout << "\n";

	for(int a=0; a<this->length; a++){

			cout << seq[a] << "";
	}

	cout << "\n";
}

void Molecule::printDist(){

	for(int a=0; a<this->length; a++)
		cout << seq[a] << "\t";

	cout << "\n";

	for(int i=0; i<this->length; i++){
		for(int j=0; j<this->length; j++){

			if(!this->null[i] && !this->null[j]){
				cout << dist[i][j] << "\t";
			}else{
				cout << "0\t";
			}
		}
		cout << "\n";
	}
}

void Molecule::printPDB(FILE *name, int length, char *type){

	for(int i=0; i<length; i++){

		if(!this->null[i] && !strcmp(type, "PDB"))
			fprintf(name, "ATOM  %5d  CA  %3s B%4d    %8.3f%8.3f%8.3f\n", i, AAmap2(this->seq[i]).c_str(), i, this->coord[i][0], this->coord[i][1], this->coord[i][2]);
		if(!this->null[i] && !strcmp(type, "fullPDB"))
			fprintf(name, "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n", i, this->atomname[i], AAmap2(this->seq[i]).c_str(), this->chain[i] , this->resnumber[i], this->coord[i][0], this->coord[i][1], this->coord[i][2]);

	}
}

void Molecule::printPDB(FILE *name, int length, char *type, double* distances){

	for(int i=0; i<length; i++){
		if(distances[i] < 5){
		if(!this->null[i] && !strcmp(type, "PDB"))
			fprintf(name, "ATOM  %5d  CA  %3s B%4d    %8.3f%8.3f%8.3f\n", i, AAmap2(this->seq[i]).c_str(), i, this->coord[i][0], this->coord[i][1], this->coord[i][2]);
		if(!this->null[i] && !strcmp(type, "fullPDB"))
			fprintf(name, "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n", i, this->atomname[i], AAmap2(this->seq[i]).c_str(), this->chain[i] , this->resnumber[i], this->coord[i][0], this->coord[i][1], this->coord[i][2]);
		}
	}
}

void Molecule::destroy(){

	//fix this
	delete[] seq;
	delete[] null;

	//delete[] tmp;
	//delete[] tmpCoord;

	for(int i=0; i<this->length; i++){
		//cout << "deleting row: " << i << " /" << length << "\n" << flush;
		delete[] dist[i];
		delete[] coord[i];
	}

	delete[] dist;
	delete[] coord;

}

#endif /*  */
