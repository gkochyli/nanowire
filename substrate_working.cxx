#include <iostream>
#include <vector>
#include "particle.h"
#include <cmath>
#include <unordered_set>
#include <iomanip>
#include <fstream>
#include <random>
#include <iomanip>
using namespace std;

void find_closest_neighbors(vector<int> &, Particle);
void check_neighborhood_boundaries(int&);
void find_set_coordinates(int&, int&, Particle);
void PBC(Particle&);
double calculate_energy(int, Particle, vector<int>, vector< unordered_set<int> > &);
bool check_Boltzmann(long double);

constexpr int L=100;					//linear dimension of lattice
constexpr int N=L/10;					//linear dimension of neighbors' vector.
constexpr double R = 0.5;				// radius of substrate atom
constexpr double Rc=3*R + 0.000001;		// cutoff distance when checking for neighbors
constexpr double PI = 3.14159265;
constexpr double ax=2*R, ay=2*sqrt(3)*R, az=2*sqrt(2.0/3)*2*R;	// distance between the centers of two atoms only for x!!! when on y, ay is the distance between two atoms with equal (i/L)%2, respectively for az.
constexpr double Lx = L*ax;
constexpr double Ly = L*ay/2;

ofstream print_1("1st_layer.txt"), print_2("2nd_layer.txt"), print_3("3rd_layer.txt");
ofstream print_energy("energy.txt");

vector<Particle> substrate(3*L*L);

mt19937_64 gen(13518);
uniform_real_distribution<> random_phi(0.0, 2*PI);
uniform_real_distribution<> random_theta(0.0, PI);
uniform_real_distribution<> rn(0.0, 1.0);

int main()
{
	vector<int> movable_substrate(L*L,0);					// periexei ta kinoumena swmatia tou 3ou layer.
	vector< unordered_set<int> > neighborhood(N*N);			// periexei set apo indeces pou vriskontai sthn idia "geitonia"
//----------------------------------------------------> Substrate Formation
	
	for(int i=0; i<L*L; i++)
	{
		if( (i/L)%2==0 ) 									
		{ 
			substrate[i].x = (i%L)*ax;
			substrate[i+L*L].x = (i%L)*ax + R;
			substrate[i+2*L*L].x = (i%L)*ax; 
		}	
		else 				
		{ 
			substrate[i].x = (i%L)*ax + R; 
			substrate[i+L*L].x = (i%L)*ax; 
			substrate[i+2*L*L].x = (i%L)*ax + R; 
		}	
			substrate[i].y = (i/L)*ay/2;
			substrate[i+L*L].y = (i/L)*ay/2 + sqrt(3)/3*R;
			substrate[i+2*L*L].y = (i/L)*ay/2;
		
			substrate[i].z = 0;
			substrate[i+L*L].z = az/2;
			substrate[i+2*L*L].z = az;
		
			substrate[i].type = 1;
		print_1 << '<' << substrate[i].x << ',' << substrate[i].z << ',' << substrate[i].y <<'>'<<endl;
		print_2 << '<' << substrate[i+L*L].x << ',' << substrate[i+L*L].z << ',' << substrate[i+L*L].y <<'>'<<endl;
		print_3 << '<' << substrate[i+2*L*L].x << ',' << substrate[i+2*L*L].z << ',' << substrate[i+2*L*L].y <<'>'<<endl;
//---------------------------------------------------->
		
		movable_substrate[i] = i+2*L*L;					// contains indeces of the 3rd layer. thats the moving substrate.
		
		int x_set, y_set;
		find_set_coordinates(x_set, y_set, substrate[i]);

		substrate[i].set = y_set*N+ x_set;
		substrate[i+L*L].set = y_set*N+ x_set;
		substrate[i+2*L*L].set = y_set*N+ x_set;
		
		neighborhood[y_set*N+ x_set].insert(i);			// vazei sto swsto set to swmatidio.
		neighborhood[y_set*N + x_set].insert(i+L*L);
		neighborhood[y_set*N + x_set].insert(i+2*L*L);
	}
//-----------------------------------------------------> Calcuting Distance of two atoms, for every close neighboring atom
	
	double sum_dE=0;										//gia sunolikh energeia
	for(int r=0; r<1000; r++)
	{
//		double sum_dE=0;									//gia diafora energeias ana time-step
		cout<< "Runs: " << r+1 << endl;
		for(int p=0; p<movable_substrate.size(); p++)		// Gia ka8e swmatio tou kinoumenou substrate
		{			
			//cout<< p << ')' << endl;
			int i = movable_substrate[p];					// i = label tou sugkekrimenou swmatidiou pou tsekarw

	//		cout << "Trexon Swmatidio: " << i << endl<<endl;
	//		cout << "Syntetagmnenes Trexontos Swmatidiou\t" << x_i << '\t' << y_i << '\t' << z_i << endl<<endl;
			vector<int> closest_neighbors(9,0);
			find_closest_neighbors(closest_neighbors, substrate[i]);	//	epistrefei mia lista mege8ous 9 me tous indeces twn geitonwn
			 	 
			long double E_old = calculate_energy(i, substrate[i], closest_neighbors, neighborhood);
			//cout << E_old << '\t' << substrate[i].x << '\t' << substrate[i].y << '\t' << substrate[i].z <<endl;
	//-----------------------------------------------------> Trial Move
			double phi = random_phi(gen);
			double theta = random_theta(gen);

			Particle trial_particle(substrate[i]);
			trial_particle.x += 0.1*ax*cos(phi)*sin(theta);
			trial_particle.y += 0.1*ax*sin(phi)*sin(theta);
			trial_particle.z += 0.1*ax*cos(theta);	

			PBC(trial_particle);

			find_closest_neighbors(closest_neighbors, trial_particle); 		

			long double E_new = calculate_energy(i, trial_particle, closest_neighbors, neighborhood);
			//cout<< E_new << '\t' << trial_particle.x << '\t' << trial_particle.y << '\t' << trial_particle.z << '\t' << endl;
			long double dE = E_new - E_old;

			double dr = sqrt(pow(substrate[i].x-trial_particle.x,2)+pow(substrate[i].y-trial_particle.y,2)+pow(substrate[i].z-trial_particle.z,2));
			//cout<< "dr = " << dr << endl;
			
			if(E_new<E_old || check_Boltzmann(dE))
			{
				//MOVE
				int x_set, y_set;
				find_set_coordinates(x_set, y_set, substrate[i]);

				int i_set_old = y_set*N + x_set; //cout<< i_set_old << '\t' << substrate[i].set<< endl<<endl;

				substrate[i] = trial_particle;	//cout<< "MOVED" <<endl;

				find_set_coordinates(x_set, y_set, trial_particle);
				trial_particle.set = y_set*N + x_set;
				int i_set_trial = y_set*N + x_set; //cout<< i_set_trial << '\t' << trial_particle.set<< endl<<endl;

				if(i_set_old!=i_set_trial)
				{
					auto result = neighborhood[i_set_old].erase(i);
					if(result == 0) cout << "DEN DIAGRAFHKE" <<endl;

					auto result1 = neighborhood[i_set_trial].insert(i);
					if(!result1.second) cout << "DEN MPHKE" <<endl;
				}
				
				sum_dE+=dE;
			}
			//cout<< "-------------------------------------------------------------------" << endl;		
		}
		print_energy<<sum_dE<<endl;
	}
}


//-----------------------------------------------------> Functions

void find_closest_neighbors(vector<int> &closest_neighbors, Particle s)  		// finds neighboring boxes in "sudoku box"
{
	int x,y;
	find_set_coordinates(x, y, s);
	
	int x_p= x+1, x_m= x-1, y_p= y+1, y_m= y-1;		// needed variables for sudoku box
	
	check_neighborhood_boundaries(x_p);
	check_neighborhood_boundaries(x_m);
	check_neighborhood_boundaries(y_p);
	check_neighborhood_boundaries(y_m);
	
	closest_neighbors[0] = y_m*N+x_m; 
	closest_neighbors[1] = y_m*N+x;
	closest_neighbors[2] = y_m*N+x_p;
	closest_neighbors[3] = y*N+x_m;
	closest_neighbors[4] = y*N+x;
	closest_neighbors[5] = y*N+x_p;
	closest_neighbors[6] = y_p*N+x_m;
	closest_neighbors[7] = y_p*N+x;
	closest_neighbors[8] = y_p*N+x_p;
}
//-----------------------------------------------------
void check_neighborhood_boundaries(int& x)					// PARAMETERS (x+ or x- or y+ or y-)
{
	if(x<0) 	x= N-1;
	if(x>N-1)	x= 0;
}
//----------------------------------------------------
void find_set_coordinates(int& x_set, int& y_set, Particle s)
{
	x_set = (int) s.x/(L/10);				 
	y_set = (int) s.y/(10*ay/2);
}
//----------------------------------------------------
void PBC(Particle& s)
{
	if(s.x<0) 	s.x += Lx;
	if(s.x>Lx) 	s.x -= Lx;
	if(s.y<0) 	s.y += Ly;
	if(s.y>Ly) 	s.y -= Ly;
}
//----------------------------------------------------
double calculate_energy(int i, Particle particle_cur, vector<int> closest_neighbors, vector< unordered_set<int> > &neighborhood)
{
	double sum_Er=0;
	double sum_Eb=0;
	
	double x_i, y_i, z_i;
	x_i = particle_cur.x;
	y_i = particle_cur.y;
	z_i = particle_cur.z;							// ta upologizw gia na ta xrhsimopoihsw parakatw.

	for(int k=0; k<closest_neighbors.size(); k++)	// For every set_index of the list
	{
		int i_set = closest_neighbors[k];			// i_set --> o index tou set pou prepei na tsekarw.

		for(int j:neighborhood[i_set])				// Gia ka8e swmatio (j) tou sugkekrimenou set.
		{
			double x_j, y_j, z_j;
			x_j = substrate[j].x;
			y_j = substrate[j].y;
			z_j = substrate[j].z;

			double dx = abs(x_i-x_j); if(dx>L/2) dx=L-dx;
			double dy = abs(y_i-y_j); if(dy>L*ay/4) dy=L*ay/2-dy;
			double dz = abs(z_i-z_j); 

			double r = sqrt(dx*dx + dy*dy + dz*dz);

			if(i!=j && r<Rc)
			{				
				if(r<=0) cout << "SHIT" <<endl;

				sum_Er += 0.0855*exp(-10.96*(r-1));
				sum_Eb += 1.224*1.224*exp(-2*2.278*(r-1));	
			}
		}
	}	
	double E_i = sum_Er - sqrt(sum_Eb);
	return E_i;
}
//-----------------------------------------------------------
bool check_Boltzmann(long double dE)
{
	long double T = 0.05;			// in kT units
	long double m = dE/T;
	long double w = exp(-m);
	
	long double R = rn(gen);
	if(R<w) 
	{
		//cout<< "Boltzmann me R = " << R << endl;
		return true;
	}
	else return false;
}
