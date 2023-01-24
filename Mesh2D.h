#ifndef _MESH_2D_H

#include<vector>
#include<string>
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
class Vertex
{
  private:
    Eigen::Vector2d _v_coor;
    int _ref;
  public:
    Vertex();
    Vertex(double x, double y, int ref);
    void Print() const;
    const Eigen::Vector2d GetCoor() const {return _v_coor;};
};

class Edge
{
  private:
    Eigen::Vector2i _v_edge;
    int _ref;
    int _t1, _t2;
  public:
    Edge();
    Edge(int vertex1, int vertex2, int ref);
    void Print() const;
    void AddTriangle(int t)
    {
      if (_t1 == -1)
        _t1 = t;
      else
        _t2 = t;
    }
    const Eigen::Vector2i& GetVertices() const { return _v_edge;}
    int GetT1() const { return _t1; };
    int GetT2() const { return _t2; };
    int GetReference() const { return _ref;};
};

class Triangle
{
  private:
    Eigen::Vector3i _v_triangle;
    int _ref;
  public:
    Triangle();
    Triangle(int vertex1, int vertex2, int vertex3, int ref);
    void Print() const;
    const Eigen::Vector3i& GetVertices() const { return _v_triangle; }
};

class Mesh2D
{
private:
  // liste de tous les sommets
  std::vector<Vertex> _vertices;
  // liste de tous les triangles
  std::vector<Triangle> _triangles;
  // centre de tous les triangles
  Eigen::Matrix<double, Eigen::Dynamic, 2> _tri_center;
  // aire de tous les triangles
  Eigen::VectorXd _tri_area;
  // liste de toutes les arêtes
  std::vector<Edge> _edges;
  
    // liste de toutes les arêtes
  std::vector<double> _dx;
  // liste de toutes les normales unitaires !!!
  Eigen::Matrix<double, Eigen::Dynamic, 2> _edg_normal;
  // liste de toutes les longueurs d'arêtes
  Eigen::VectorXd _edg_length;
  // centre des aretes
  Eigen::Matrix<double, Eigen::Dynamic, 2> _edg_center;

public:
  Mesh2D();
  void ReadMesh(std::string name_mesh);
  void BuildTrianglesCenterAndArea();
  void BuildEdgesNormalLengthAndCenter();

  const std::vector<Vertex> & GetVertices() const {return _vertices;};

  const std::vector<Triangle> & GetTriangles() const {return _triangles;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2> & GetTrianglesCenter() const {return _tri_center;};
  const Eigen::VectorXd & GetTrianglesArea() const  {return _tri_area;};

  const std::vector<Edge> & GetEdges() const {return _edges;};
  const Eigen::VectorXd & GetEdgesLength() const {return _edg_length;};
  const std::vector<double> & Getdx() const {return _dx;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2> & GetEdgesNormal() const {return _edg_normal;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2> & GetEdgesCenter() const {return _edg_center;};

protected:
  void AddSingleEdge(const Edge& edge, int ne, std::vector<int>& head_minv,
		     std::vector<int>& next_edge, int& nb_edges);
};

#define _MESH_2D_H
#endif
