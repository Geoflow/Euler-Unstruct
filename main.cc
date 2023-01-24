#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string data_file_name = argv[1];

  // ----------------------- Fichier de données --------------------------------
  DataFile* data_file = new DataFile(data_file_name);
  data_file->ReadDataFile();
  // ---------------------------------------------------------------------------

  // ------------------Définition du nombre d'itérations------------------------
  int nb_iterations = int(ceil((data_file->Get_tfinal()-data_file->Get_t0())/data_file->Get_dt()));
 // data_file->Adapt_dt(data_file->Get_tfinal() / nb_iterations);
  // ---------------------------------------------------------------------------

  // ---------------------------- Résolution  ----------------------------------
  Mesh2D* mesh = new Mesh2D();
  mesh->ReadMesh(data_file->Get_mesh_name());
  Function* source_fct = new Function(data_file);
  Advection* adv = new Advection(source_fct, data_file, mesh);
  TimeScheme* time_scheme = NULL;

  if (data_file->Get_scheme() == "RungeKutta")
    time_scheme = new RungeKuttaScheme();
  else
    time_scheme = new EulerScheme();

  cout << "-------------------------------------------------" << endl;
  cout << "Search W such that : " << endl;
  cout << "dt u + grad_x F(W)+grad_y G(W) = 0,  B" << endl;
  cout << "-------------------------------------------------" << endl;

  // Démarrage du chrono
  
  auto start = chrono::high_resolution_clock::now();
  time_scheme->Initialize(data_file, adv); // Initialisation
  time_scheme->SaveSolution(0); // Sauvegarde condition initiale

  double Tmax(data_file->Get_tfinal());
  int n(1);


  while (time_scheme->get_t() <Tmax)
  {
    time_scheme->Advance();
    time_scheme->SaveSolution(n);
    n++;

  }

//time_scheme->SaveSolution(n);
  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< t << " milliseconds" << endl;
  // ---------------------------------------------------------------------------




  delete time_scheme;
  delete adv;
  delete data_file;
  delete source_fct;

  return 0;
}
