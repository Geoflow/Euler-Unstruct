#ifndef _ADVECTION_CPP
#include "Advection.h"
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;
using namespace Eigen;

// Constructeur
Advection::Advection(Function* source_fct, DataFile* data_file, Mesh2D* mesh) :
_function(source_fct), _vertices(mesh->GetVertices()),_triangles(mesh->GetTriangles()),
_edges(mesh->GetEdges()), _tri_center(mesh->GetTrianglesCenter()), _tri_area(mesh->GetTrianglesArea()),
_edg_normal(mesh->GetEdgesNormal()), _edg_center(mesh->GetEdgesCenter()), _edg_length(mesh->GetEdgesLength()),
_results(data_file->Get_results()), _dx(mesh->Getdx())
{
	system(("mkdir -p ./" + _results).c_str());
	system(("rm -f ./" + data_file->Get_results() + "/*.vtk").c_str());
	// Copier le fichier de données dans le dossier résultats
	system(("cp -r ./" + data_file->Get_file_name() + " ./" + _results + "/" + data_file->Get_file_name()).c_str());

	if (data_file->Get_numerical_flux_choice() == "centered")
    _numerical_flux = Centered;
  else
    _numerical_flux = Upwind;

    W_prev.resize(5,_triangles.size());
    W_next.resize(5,_triangles.size());
    InitialEulerCondition();

}



double Advection::Set_dt()
{
    std::vector<double> Wave;
    double sl,sr,sd,su,step;
    double cfl(0.5);


    for (int i = 0; i < _triangles.size(); i++)
    {
        sl=abs(W_prev.col(i)[1]/W_prev.col(i)[0]-sqrt(1.4*W_prev.col(i)[4]/W_prev.col(i)[0]));
        sr=abs(W_prev.col(i)[1]/W_prev.col(i)[0])+sqrt(1.4*W_prev.col(i)[4]/W_prev.col(i)[0]);

        sd=abs(W_prev.col(i)[2]/W_prev.col(i)[0]-sqrt(1.4*W_prev.col(i)[4]/W_prev.col(i)[0]));
        su=abs(W_prev.col(i)[2]/W_prev.col(i)[0])+sqrt(1.4*W_prev.col(i)[4]/W_prev.col(i)[0]);

        const auto tmp={sr,su,sl,sd};
        const auto max = std::max_element(tmp.begin(),tmp.end());
        Wave.push_back(*max);
    }

        const auto h=_edg_length.minCoeff();
        const auto max = std::max_element(begin(Wave), end(Wave));
        step=cfl*h/ *max;
        Wave.clear();
        if (step==0) std::cout<<" PB  dt=0. "<<std::endl;
        return step;

}

void Advection::SetW_prev(Eigen::MatrixXd X){W_prev=X;};
void Advection::SetW_next(Eigen::MatrixXd X){W_next=X;};




// Construit la condition initiale au centre des triangles
void Advection::InitialEulerCondition()
{
   W_prev.resize(5,_triangles.size());
    //Construire sol0 au centre des triangles
    for (int i = 0; i < _triangles.size(); i++)
  {
    W_prev.col(i)=_function->InitialEuler(0,_tri_center(i,0),_tri_center(i,1));
   }

}

void Advection::BuildFlux(int axe )

{
    Flux=MatrixXd::Zero(5,_triangles.size());
    VectorXd tmp(5),Wx(5),Wy(5),Wr(5),Wl(5);
    double theta;




    for (int i = 0; i < _edges.size(); i++)
    {
       if(_edges[i].GetT2()!=-1)
         // internal cells of the mesh
        {
           Eigen::MatrixXd R_theta(5,5);
           R_theta=Eigen::MatrixXd::Identity(5,5);
           R_theta(1,1)=_edg_normal(i,0);
           R_theta(1,2)=-_edg_normal(i,1);
           R_theta(2,1)=_edg_normal(i,1);
           R_theta(2,2)=_edg_normal(i,0);
        
          if (axe==0){
           Wx=_function->FHLL(R_theta*W_prev.col(_edges[i].GetT1()),W_prev.col(_edges[i].GetT2()));
            tmp=_edg_length[i]*Wx*_edg_normal(i,0);
          }
          else{
           Wy=_function->GHLL(R_theta.transpose()*W_prev.col(_edges[i].GetT1()),R_theta.transpose()*W_prev.col(_edges[i].GetT2()));
          tmp=_edg_length[i]*Wy*_edg_normal(i,1);
          }
          


            //tmp=_edg_length[i]*(Wx*_edg_normal(i,0)+Wy*_edg_normal(i,1));
            Flux.col(_edges[i].GetT1())-=tmp;
            Flux.col(_edges[i].GetT2())+=tmp;

         }
         else
         { // boundary condition
                 //Eigen::VectorXd Wx(5);
                 Wx=Eigen::VectorXd::Zero(5);
                 Flux.col(_edges[i].GetT1())=Wx;
          }
    }

}

void Advection::Update()
{

    W_next=W_prev-0.5*Set_dt()*Flux;
    Pressure();
    W_prev=W_next;

}




void Advection::Pressure()
{
    for (int i = 0; i < _triangles.size(); i++)
    {
   W_next.col(i)[4]=std::max(0.,0.4*(W_next.col(i)[3]/W_next.col(i)[0]-0.5*W_next.col(i)[0]*(pow(W_next.col(i)[1]/W_next.col(i)[0],2)+pow(W_next.col(i)[2]/W_next.col(i)[0],2))));
    }

}



// Sauvegarde la solution
void Advection::SaveSol(const Eigen::VectorXd& sol, int n)
{
	string name_file = _results + "/solution_" + std::to_string(n) + ".vtk";

  int nb_vert = _vertices.size();

  assert((W_prev.cols() == _triangles.size()) && "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  solution << "# vtk DataFile Version 3.0 " << endl;
  solution << "2D Unstructured Grid" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET UNSTRUCTURED_GRID" << endl;

  solution << "POINTS " << nb_vert << " float " << endl;
  for (int i = 0 ; i < nb_vert ; ++i)
  {
    solution << ((_vertices[i]).GetCoor())[0] << " " << ((_vertices[i]).GetCoor())[1] << " 0." << endl;
  }
  solution << endl;

  solution << "CELLS " << _triangles.size() << " " << _triangles.size()*4 << endl;
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << 3 << " " << ((_triangles[i]).GetVertices())[0] << " " << ((_triangles[i]).GetVertices())[1]
    << " " << ((_triangles[i]).GetVertices())[2] << endl;
  }
  solution << endl;

  solution << "CELL_TYPES " << _triangles.size() << endl;
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << 5 << endl;
  }
  solution << endl;

  solution << "CELL_DATA " << _triangles.size() << endl;
  solution << "SCALARS W_prev float 5" << endl;
  solution << "LOOKUP_TABLE default" << endl;

  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << float(GetW_prev().col(i)[0]) <<' '<< float(GetW_prev().col(i)[1]/GetW_prev().col(i)[0]) <<' '<< float(GetW_prev().col(i)[2]/GetW_prev().col(i)[0]) <<' '<< float(GetW_prev().col(i)[3]/GetW_prev().col(i)[0])<<' '<< float(GetW_prev().col(i)[4])  << endl;
  }
  solution << endl;

	solution.close();
}

#define _ADVECTION_CPP
#endif
