#ifndef _FUNCTION_H

#include "DataFile.h"
#include "Eigen/Eigen/Dense"
class Function {
  private:
    int _initial_condition_choice; // Gaussian or Rectangular
    enum {Gaussian, Rectangular};
    int _velocity_choice; // Uniform, Rotational or Sinusoidal
    enum {Uniform, Rotational, Sinusoidal};

    // parameters x0, y0, a, b, c, d, e and f
    double _x0, _y0, _a, _b, _c, _d, _e, _f;
    bool _in_flow;
	public: // Méthodes et opérateurs de la classe
    Function(DataFile* data_file);
    double InitialCondition(const double x, const double y) const;



 //-------------------Euler system functions------------------------------
    Eigen::VectorXd InitialEuler(int cas_test,const double x, const double y) const;
    Eigen::VectorXd Explosion_test(const double x, const double y) const;
    Eigen::VectorXd case_test_3(const double x, const double y) const;
    Eigen::VectorXd F(Eigen::VectorXd& W) ;
    Eigen::VectorXd G(Eigen::VectorXd& W) ;
    Eigen::VectorXd FHLL(Eigen::VectorXd WR,Eigen::VectorXd WL) ;
    Eigen::VectorXd GHLL(Eigen::VectorXd WR,Eigen::VectorXd WL) ;


};

#define _FUNCTION_H
#endif
