#include <iostream>
#include <vector>
#include "particle.h"
#include <cmath>
#include <unordered_set>
#include <iomanip>
#include <fstream>
#include <random>
#include <sstream>
using namespace std;

void find_closest_neighbors(vector<int> &, Particle);
void check_neighborhood_boundaries(int&);
void find_set_coordinates(int&, int&, Particle);
void PBC(Particle&);
double calculate_energy(int, Particle, vector<int>, vector< unordered_set<int> > &);
bool check_Boltzmann(double);
void place_substrate();
void place_nanowire();
void print2files(int r);

constexpr int L=100;					//linear dimension of lattice
constexpr double R = 0.5;				// radius of substrate atom
constexpr double Rc=3*R + 0.000001;		// cutoff distance when checking for neighbors
constexpr double PI = 3.14159265;
constexpr double ax=2*R, ay=2*sqrt(3)*R, az=2*sqrt(2.0/3)*2*R;	// distance between the centers of two atoms only for x!!! when on y, ay is the distance between two atoms with equal (i/L)%2, respectively for az.
constexpr double Lx = L*ax;
constexpr double Ly = L*ay/2;

// Neighborhood Variables
constexpr int N=10;					//linear dimension of neighbors' vector.

// Sta8eres gia Pb-nanowire
constexpr double R_pb = 0.66;
constexpr double ax_pb = 2*R_pb;
constexpr int num_of_pb = 68; 		//number of Pb particles

int snap_at = 10;					// take snapshots every 100 steps.
int Runs = 1000;
	
stringstream filename;
ofstream print_energy("energy.txt");

vector<Particle> particles(3*L*L + num_of_pb);			// periexei ola ta particles tou susthmatos.
vector<int> movable_particles(L*L + num_of_pb,0);		// periexei ta labels gia ta kinoumena swmatia tou 3ou layer kai to nanowire
vector< unordered_set<int> > neighborhood(N*N);			// periexei set apo indeces pou vriskontai sthn idia "geitonia"

mt19937_64 gen(13518);
uniform_real_distribution<> random_phi(0.0, 2*PI);
uniform_real_distribution<> random_theta(0.0, PI);
uniform_real_distribution<> rn(0.0, 1.0);

int main()
{
//----------------------------------------------------> System Formation
	place_substrate();	
	place_nanowire();
	print2files(0);
//-----------------------------------------------------> Calcuting Distance of two atoms, for every close neighboring atom
	double sum_dE=0;
	for(int r=0; r<Runs; r++)
	{
		cout << "Run: " <<r+1 << endl;
		for(int p=0; p<movable_particles.size(); p++)		// Gia ka8e kinoumeno swmatio
		{			
			//sum_dE=0;
			//cout <<"Checking particle: " << p << ')' << endl;
			int i = movable_particles[p];					// i = label tou sugkekrimenou swmatidiou pou tsekarw

			vector<int> closest_neighbors(9,0);
			find_closest_neighbors(closest_neighbors, particles[i]);	//	epistrefei mia lista mege8ous 9 me tous indeces twn geitonwn		

			double E_old = calculate_energy(i, particles[i], closest_neighbors, neighborhood);
	//-----------------------------------------------------> Trial Move
			double phi = random_phi(gen);
			double theta = random_theta(gen);

			Particle trial_particle(particles[i]);
			trial_particle.x += 0.1*ax*cos(phi)*sin(theta);
			trial_particle.y += 0.1*ax*sin(phi)*sin(theta);
			trial_particle.z += 0.1*ax*cos(theta);	

			PBC(trial_particle);

			find_closest_neighbors(closest_neighbors, trial_particle);

			double E_new = calculate_energy(i, trial_particle, closest_neighbors, neighborhood);

			double dE = E_new - E_old;

			//double dr = sqrt(pow(particles[i].x-trial_particle.x,2)+pow(particles[i].y-trial_particle.y,2)+pow(particles[i].z-trial_particle.z,2));

			if(E_new<E_old || check_Boltzmann(dE))
			{
				//MOVE
				int x_set, y_set;
				find_set_coordinates(x_set, y_set, particles[i]);

				int i_set_old = y_set*N + x_set; //cout << i_set_old << '\t' << particles[i].set<< endl<<endl;

				particles[i] = trial_particle;	//cout << "MOVED" <<endl;

				find_set_coordinates(x_set, y_set, trial_particle);
				//trial_particle.set = y_set*N + x_set;
				int i_set_trial = y_set*N + x_set; //cout << i_set_trial << '\t' << trial_particle.set<< endl<<endl;

				if(i_set_old!=i_set_trial)
				{
					auto result = neighborhood[i_set_old].erase(i);
					if(result == 0) 
					{
						cout << "DEN DIAGRAFHKE" <<endl;
						return 1053;
					}

					auto result1 = neighborhood[i_set_trial].insert(i);
					if(!result1.second)
					{
						cout << "DEN MPHKE" <<endl;
						return 4930;
					}
				}
				sum_dE+=dE;
			}
			//cout << "-------------------------------------------------------------------" << endl;		
		}	
		print_energy << sum_dE <<endl;
		
		if((r+1)%snap_at==0)
		{
			print2files(r);
		}	
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
	x_set = (int) s.x/(Lx/N);				 
	y_set = (int) s.y/(Ly/N);
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
	z_i = particle_cur.z;

	for(int k=0; k<closest_neighbors.size(); k++)	// For every set_index of the list
	{
		int i_set = closest_neighbors[k];			// i_set --> o index tou set pou prepei na tsekarw.

		for(int j:neighborhood[i_set])				// Gia ka8e swmatio (j) tou sugkekrimenou set.
		{
			double x_j, y_j, z_j;
			x_j = particles[j].x;
			y_j = particles[j].y;
			z_j = particles[j].z;

			double dx = abs(x_i-x_j); if(dx>L/2) dx=L-dx;
			double dy = abs(y_i-y_j); if(dy>L*ay/4) dy=L*ay/2-dy;
			double dz = abs(z_i-z_j); 

			double r = sqrt(dx*dx + dy*dy + dz*dz);

			double A_ab, j_ab, p_ab, q_ab;
			if(particle_cur.type == particles[j].type)		// kai ta duo Cu OR kai ta duo PB
			{
				if(particle_cur.type==1)					// einai Cu
				{
					A_ab = 0.0855;
					j_ab = 1.224;
					p_ab = 10.96;
					q_ab = 2.278;
				}
				else										// einai Pb
				{
					A_ab = 0.098;
					j_ab = 0.914;
					p_ab = 9.576;
					q_ab = 3.648;
				}
			}
			else											// PREPEI NA TSEKARW if type!= 1 h 2
			{
				A_ab = sqrt(pow(0.0855,2) + pow(0.098,2));
				j_ab = sqrt(pow(1.224,2) + pow(0.914,2));
				p_ab = sqrt(pow(10.96,2) + pow(9.576,2));
				q_ab = sqrt(pow(2.278,2) + pow(3.648,2));
			}
			
			
			if(i!=j && r<Rc)
			{				
				if(r<=0) 
				{
					cout << "SHIT" <<endl;
					return 482;
				}

				sum_Er += A_ab*exp(-p_ab*(r-1));
				sum_Eb += j_ab*exp(-2*q_ab*(r-1));	
			}
		}
	}	
	double E_i = sum_Er - sqrt(sum_Eb);
	return E_i;
}
//-----------------------------------------------------------
bool check_Boltzmann(double dE)
{
	double T = 0.05;			// in kT units
	double w = exp(-dE/T);
	
	double R = rn(gen);
	if(R<w) 
	{
		//cout << "Boltzmann me R = " << R << endl;
		return true;
	}
	else return false;
}
//-----------------------------------------------------------
void place_substrate()
{
	for(int i=0; i<L*L; i++)
	{
		if( (i/L)%2==0 ) 									
		{ 
			particles[i].x = (i%L)*ax;
			particles[i+L*L].x = (i%L)*ax + R;
			particles[i+2*L*L].x = (i%L)*ax; 
		}	
		else 				
		{ 
			particles[i].x = (i%L)*ax + R; 
			particles[i+L*L].x = (i%L)*ax; 
			particles[i+2*L*L].x = (i%L)*ax + R; 
		}	
			particles[i].y = (i/L)*ay/2;
			particles[i+L*L].y = (i/L)*ay/2 + sqrt(3)/3*R;
			particles[i+2*L*L].y = (i/L)*ay/2;
		
			particles[i].z = 0;
			particles[i+L*L].z = az/2;
			particles[i+2*L*L].z = az;
		
			particles[i].type = 1;
			particles[i+L*L].type = 1;
			particles[i+2*L*L].type = 1;
		
//---------------------------------------------------->
		
		movable_particles[i] = i+2*L*L;					// ta prwta L*L stoixeia tou vector periexoun ta labels twn kinoumenwn swmatidiwn tou 3ou layer.
		
		int x_set, y_set;
		find_set_coordinates(x_set, y_set, particles[i]);	//ISWS EINAI LA8OS AUTO!!!!
		neighborhood[y_set*N+ x_set].insert(i);			// vazei sto swsto set to swmatidio.
		
		find_set_coordinates(x_set, y_set, particles[i+L*L]);	//ISWS EINAI LA8OS AUTO!!!!
		neighborhood[y_set*N + x_set].insert(i+L*L);
		
		find_set_coordinates(x_set, y_set, particles[i+2*L*L]);	//ISWS EINAI LA8OS AUTO!!!!
		neighborhood[y_set*N + x_set].insert(i+2*L*L);		
	}	
}
//-----------------------------------------------------------
void place_nanowire()
{
	for(int i=0; i<num_of_pb; i++)
	{
		particles[3*L*L+i].x = i*ax_pb + 10.0;		// 3ekinaei apto x = 10
		particles[3*L*L+i].y = 5*ay/2 + 4*Ly/10;	//topo8etw to particles sthn pempth seira tou 4 set
		particles[3*L*L+i].z = az + R + R_pb;		//topo8etw to nanowire se uyos mia aktinas Cu kai mias aktinas Pb
		
		particles[3*L*L+i].type = 2;
		
		movable_particles[L*L+i] = 3*L*L+i;	//ta swmatidia tou nanowire 8a exoun labels 30000, 30001, 30002, ..., 30067 
				
		int x_set, y_set;
		find_set_coordinates(x_set, y_set, particles[3*L*L+i]);
		neighborhood[y_set*N+ x_set].insert(3*L*L+i);			// vazei sto swsto set to swmatidio.		
	}
}
//------------------------------------------------------------
void print2files(int r)
{
	if (r==0) 	filename << r << "steps";
	else		filename << r+1 << "steps";
	
	ofstream print_1("./coordinates/1st_layer_"+filename.str()+".txt"), print_2("./coordinates/2nd_layer_"+filename.str()+".txt"), print_3("./coordinates/3rd_layer_"+filename.str()+".txt");
	ofstream print_nanowire("./coordinates/nanowire_coordinates_"+filename.str()+".txt");

	for(int m=0; m<L*L; m++)
	{

		print_1 << '<' << particles[m].x << ',' << particles[m].z << ',' << particles[m].y <<'>'<<endl;
		print_2 << '<' << particles[m+L*L].x << ',' << particles[m+L*L].z << ',' << particles[m+L*L].y <<'>'<<endl;
		print_3 << '<' << particles[m+2*L*L].x << ',' << particles[m+2*L*L].z << ',' << particles[m+2*L*L].y <<'>'<<endl;
	}

	for(int m=0; m<num_of_pb; m++)
	{
		print_nanowire << '<' << particles[m+3*L*L].x << ',' << particles[m+3*L*L].z << ',' << particles[m+3*L*L].y <<">,"<<endl;
	}

	filename.str(string());			//empty filename stream
}
//made a comment at the end