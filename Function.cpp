#ifndef _FUNCTION_CPP

#include "Function.h"
#define _USE_MATH_DEFINES
#include <cmath>

Function::Function(DataFile* data_file) :
_x0(data_file->Get_param_x0()),_y0(data_file->Get_param_y0()),_a(data_file->Get_param_a()),
_b(data_file->Get_param_b()),_c(data_file->Get_param_c()),_d(data_file->Get_param_d()),
_e(data_file->Get_param_e()),_f(data_file->Get_param_f())
{
  if (data_file->Get_initial_condition_choice() == "gaussian")
  {
    _initial_condition_choice = Gaussian;
  }
  else if (data_file->Get_initial_condition_choice() == "rectangular")
  {
    _initial_condition_choice = Rectangular;
  }
  else
  {
    std::cout << "Invalid initial condition" << std::endl;
    abort();
  }

  if (data_file->Get_velocity_choice() == "uniform")
  {
    _velocity_choice = Uniform;
  }
  else if (data_file->Get_velocity_choice() == "rotational")
  {
    _velocity_choice = Rotational;
  }
  else if (data_file->Get_velocity_choice() == "sinusoidal")
  {
    _velocity_choice = Sinusoidal;
  }
  else
  {
    std::cout << "Invalid velocity" << std::endl;
    abort();
  }
}



double Function::InitialCondition(const double x, const double y) const
{
  switch (_initial_condition_choice)
  {
    case Gaussian:
      return exp(-_a*(pow(x-_x0,2)+pow(y-_y0,2))); //DONE
      break;
    case Rectangular:
       if (pow(x-_x0,2)+pow(y-_y0,2)<_b) return 1.;
       else return 0.; // DONE
      break;
    default:
      abort();
  }
}


Eigen::VectorXd Function::Explosion_test(const double x, const double y) const
{

 Eigen::VectorXd W0(5);
if(pow(x,2)+pow(y,2)  <0.25)
    {
        W0[0]=1.;
        W0[1]=0.0;
        W0[2]=0.;
        W0[4]=1.;
        W0[3]=(W0[4]+0.5*(pow(W0[1],2)+pow(W0[2],2)));// We need to compute E for the initial data
    }
    else
    {
        W0[0]=0.125;
        W0[1]=0.;
        W0[2]=0.;
        W0[4]=0.1;
        W0[3]=(W0[4]+0.5*(pow(W0[1],2)+pow(W0[2],2)));// We need to compute E for the initial data
    }


return W0;
}


Eigen::VectorXd Function::case_test_3(const double x, const double y) const
{
 Eigen::VectorXd W0(5);
if(x <=0. && y<= 0.) // Bottom Left
    {
        W0[0]=0.138;
        W0[1]=1.206*W0[0];
        W0[2]=1.206*W0[0];
        W0[4]=0.029;
        W0[3]=(W0[4]/0.4+0.5*(pow(W0[1],2)+pow(W0[2],2)));// We need to compute E for the initial data          
    }
if(x <=0. && 0.<y) // Top Left
    {
        W0[0]=0.5323;
        W0[1]=1.206*W0[0];
        W0[2]=0.;
        W0[4]=0.3;
        W0[3]=(W0[4]/0.4+0.5*(pow(W0[1],2)+pow(W0[2],2)));// We need to compute E for the initial data
}

if(0.<x && y< 0.) // Bottom Right
    {
        W0[0]=0.5323;
        W0[1]=0.0;
        W0[2]=1.206*W0[0];
        W0[4]=0.0;
        W0[3]=(W0[4]/0.4+0.5*(pow(W0[1],2)+pow(W0[2],2)));// We need to compute E for the initial data       
    }
if( 0<x &&  0.<y) // Top Right
    {
        W0[0]=1.5;
        W0[1]=0.;
        W0[2]=0.;
        W0[4]=1.5;
        W0[3]=(W0[4]/0.4+0.5*(pow(W0[1],2)+pow(W0[2],2)));// We need to compute E for the initial data
    }

return W0;
}


Eigen::VectorXd Function::InitialEuler(int cas_test, const double x, const double y) const
{
  
    Eigen::VectorXd W0(5);
    //switch (cas_test) {
    //case 0:
    W0=Explosion_test(x,y);
    //case 3:
    //W0=case_test_3(x,y);

//}
return W0;
}


Eigen::VectorXd Function::F(Eigen::VectorXd& W)
{
    assert(W.size()==5);
 Eigen::VectorXd F(5);
 F[0]=W[1];
 F[1]=W[0]*pow(W[1],2)/W[0]+W[4];
 F[2]=W[1]*W[2]/W[0];
 F[3]=W[1]/W[0]*(W[3]+W[4]);
 F[4]=W[4];
 return F;
}

Eigen::VectorXd Function::G(Eigen::VectorXd& W)
{
    assert(W.size()==5);
 Eigen::VectorXd G(5);
 G[0]=W[2];
 G[1]=W[1]*W[2]/W[0];
 G[2]=W[0]*pow(W[2],2)/W[0]+W[4];
 G[3]=W[2]/W[0]*(W[3]+W[4]);
 G[4]=W[4];
 return G;
}

Eigen::VectorXd Function::FHLL(Eigen::VectorXd WL,Eigen::VectorXd WR)
{


 double SL,SR,u,a;
 Eigen::VectorXd W(5);


  u=((WR[1]/sqrt(WR[0]))+(WL[1]/sqrt(WL[0])))/(sqrt(WL[0])+sqrt(WR[0]));
  a=0.4*(((WL[3]+WL[4])/sqrt(WL[0]))+(WR[3]+WR[4])/sqrt(WR[0]))/(sqrt(WL[0])+sqrt(WR[0]))-0.2*(pow(WR[1],2)+pow(WL[1],2));
  SL=std::min(u-a, WL[1]/WL[0]-sqrt(1.4*WL[4]/WL[0]));
  SR=std::max(u+a, WR[1]/WR[0]+sqrt(1.4*WR[4]/WR[0]));
  //SL=WL[1]/WL[0]-sqrt(1.4*WL[4]/WL[0]);
  //SR=WR[1]/WR[0]+sqrt(1.4*WR[4]/WR[0]);
 if (SL<0. && 0.< SR)
 {
     W= (SR*F(WL)-SL*F(WR)+SR*SL*(WR-WL))/(SR-SL);

 }
 else{
 if (0.<= SL)  W=F(WL);
 if (SR<= 0.) W=F(WR);
 }


 return W;
}


 Eigen::VectorXd Function::GHLL(Eigen::VectorXd WL,Eigen::VectorXd WR)
 {

  double SL,SR,v,a;
  Eigen::VectorXd W(5);


  v=((WR[2]/sqrt(WR[0]))+(WL[2]/sqrt(WL[0])))/(sqrt(WL[0])+sqrt(WR[0]));
  a=0.4*(((WL[3]+WL[4])/sqrt(WL[0]))+((WR[3]+WR[4])/sqrt(WR[0])))/(sqrt(WL[0])+sqrt(WR[0]))-0.2*(pow(WL[2],2)+pow(WR[2],2));
  SL=std::min(v-a, WL[2]/WL[0]-sqrt(1.4*WL[4]/WL[0]));
  SR=std::max(v+a, WR[2]/WR[0]+sqrt(1.4*WR[4]/WR[0]));
 // SL=WL[2]/WL[0]-sqrt(1.4*WL[4]/WL[0]);
 // SR=WR[2]/WR[0]+sqrt(1.4*WR[4]/WR[0]);

  
  if (SL<0. && 0.< SR)
  {
      W=(SR*G(WL)-SL*G(WR)+SR*SL*(WR-WL))/(SR-SL);
  }
  else{
  if (0. <= SL) W=G(WL);
  if (SR <= 0.) W=G(WR);
 }
return W;
}




#define _FUNCTION_CPP
#endif
