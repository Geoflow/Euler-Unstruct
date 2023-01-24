#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme() : _adv(0)
{}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{}

// Initialisation de vos différentes variables
void TimeScheme::Initialize(DataFile* data_file, Advection* adv)
{
  _t = data_file->Get_t0();
  _adv = adv;
}




void TimeScheme::SaveSolution(int n)
{
  _adv->SaveSol(_sol, n);
}



// Euler Explicite
void EulerScheme::Advance()
{
    //construire les flux via Hll 1d
    _adv->BuildFlux(0);
    //update
    _adv->Update();
//construire les flux via Hll 1d
    _adv->BuildFlux(1);
_adv->Update();
    _t += _dt;
}

// RungeKutta
void RungeKuttaScheme::Advance()
{
  VectorXd k1, k2, k3, k4;
 /* _adv->BuildF(_t, _sol);
  k1 = _adv->GetF();
  _adv->BuildF(_t+_dt/2., _sol+_dt/2.*k1);
  k2 = _adv->GetF();
  _adv->BuildF(_t+_dt/2., _sol+_dt/2.*k2);
  k3 = _adv->GetF();
  _adv->BuildF(_t+_dt, _sol+_dt*k3);
  k4 = _adv->GetF();
  _sol += _dt/6.*(k1 + 2.*k2 + 2.*k3 + k4);*/
  _t += _dt;
}

#define _TIME_SCHEME_CPP
#endif
