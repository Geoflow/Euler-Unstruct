#ifndef _ADVECTION_H
#include <string>

#include "Eigen/Eigen/Sparse"
#include "Function.h"
#include "Mesh2D.h"

class Advection {
	private:
		// Pointeur de la classe Function (accès à la condition initiale, la solution
		// exacte et à la vitesse)
        Function* _function;
		// Variables géométriques
		const std::vector<Triangle>& _triangles;
		const std::vector<Vertex>& _vertices;
		const std::vector<Edge>& _edges;
		const std::vector<double>& _dx;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _tri_center;
		const Eigen::VectorXd & _tri_area;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _edg_normal;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _edg_center;
		const Eigen::VectorXd & _edg_length;
		// Dossier qui contient les résultats
    const std::string _results;
		// Le choix du flux
		int _numerical_flux; // Centered ou Upwind
  	enum{Centered, Upwind};

		// Écoulement V aux centres des triangles
		Eigen::Matrix<double, Eigen::Dynamic, 2> _V;

		// Vecteur f = g(u,t) (EDO : du/dt = g(u))
		Eigen::VectorXd _f;
        Eigen::MatrixXd Flux;

	public:

        Eigen::MatrixXd W_prev,W_next;
		// Constructeur
        Advection(Function* source_fct, DataFile* data_file, Mesh2D* mesh);


       void InitialEulerCondition();
        double Set_dt();
        // Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
        void BuildFlux(int axe);
        void Update();
        Eigen::MatrixXd GetFlux() {return Flux;};
        void Pressure();
        const Eigen::MatrixXd  & GetW_prev() const{return W_prev;};
        const Eigen::MatrixXd  & GetW_next() const{return W_next;}
        void SetW_prev(Eigen::MatrixXd X);
        void SetW_next(Eigen::MatrixXd X);


        // Sauvegarde la solution
        void SaveSol(const Eigen::VectorXd& sol, int n);

};

#define _ADVECTION_H
#endif
