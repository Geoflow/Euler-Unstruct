#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name),  _if_mesh_name(false), _if_t0(false), _if_tfinal(false), _if_dt(false),
_if_scheme(false), _if_velocity_choice(false), _if_initial_condition_choice(false),
_if_numerical_flux_choice(false), _if_results(false)
{}

  void DataFile::ReadDataFile()
  {
    ifstream data_file(_file_name.data());
    if (!data_file.is_open())
    {
      cout << "Unable to open file " << _file_name << endl;
      abort();
    }
    else
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Reading data file " << _file_name << endl;
    }

    string file_line;

    while (!data_file.eof())
    {
      getline(data_file, file_line);
      if (file_line.find("mesh") != std::string::npos)
      {
        data_file >> _mesh_name; _if_mesh_name = true;
      }

      if (file_line.find("velocity") != std::string::npos)
      {
        data_file >> _velocity_choice; _if_velocity_choice = true;
        if (_velocity_choice == "uniform")
        {
          data_file >> _c >> _d;
        }
        else if (_velocity_choice == "rotational")
        {
          data_file >> _c >> _d;
        }
        else if (_velocity_choice == "sinusoidal")
        {
          data_file >> _e >> _f;
        }
        else
        {
          cout << "Only uniform, sinusoidal and rotational velocities are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("initial_condition") != std::string::npos)
      {
        data_file >> _initial_condition_choice; _if_initial_condition_choice = true;
        if (_initial_condition_choice == "gaussian")
        {
          data_file >> _x0 >> _y0 >> _a;
        }
        else if (_initial_condition_choice == "rectangular")
        {
          data_file >> _x0 >> _y0 >> _b;
        }
        else
        {
          cout << "Only gaussian or rectangular initial conditions are implemented." << endl;
          abort();
        }
      }
      if (file_line.find("numerical_flux") != std::string::npos)
      {
        data_file >> _numerical_flux_choice; _if_numerical_flux_choice = true;
        if ((_numerical_flux_choice != "centered") && (_numerical_flux_choice != "upwind"))
        {
          cout << "Only centered and upwind numerical flows are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("t0") != std::string::npos)
      {
        data_file >> _t0; _if_t0 = true;
      }

      if (file_line.find("tfinal") != std::string::npos)
      {
        data_file >> _tfinal; _if_tfinal = true;
      }

      if (file_line.find("dt") != std::string::npos)
      {
        data_file >> _dt; _if_dt = true;
      }

      if (file_line.find("scheme") != std::string::npos)
      {
        data_file >> _scheme; _if_scheme = true;
        if ((_scheme != "ExplicitEuler") && (_scheme != "RungeKutta"))
        {
          cout << "Only Explicit Euler and Runge Kutta 4 are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("results") != std::string::npos)
      {
        data_file >> _results; _if_results = true;
      }
    }

    if (!_if_t0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (0.) is used for t0." << endl;
      _t0 = 0.;
    }
    if (!_if_tfinal)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (0.1) is used for tfinal." << endl;
      _tfinal = 0.1;
    }
    if (!_if_dt)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (0.001) is used for dt." << endl;
      _dt = 0.001;
    }
    if (!_if_scheme)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default scheme (Runge Kutta 4 scheme) is used." << endl;
      _scheme = "RungeKutta";
    }
    if (!_if_results)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default results folder name (results) is used." << endl;
      _results = "results";
    }
    cout << "-------------------------------------------------" << endl;
    if (!_if_mesh_name)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Do not forget to give the mesh name in the data file." << endl;
      abort();
    }
    if (!_if_velocity_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (uniforme) is used for the velocity." << endl;
      _velocity_choice = "uniform";
    }
    if (!_if_initial_condition_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (gaussian) is used for the initial condition." << endl;
      _initial_condition_choice = "gaussian";
    }
    if (!_if_numerical_flux_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (centered) is used for the numerical flow." << endl;
      _numerical_flux_choice = "centered";
    }
    // Conditions
    if (_scheme == "ExplicitEuler")
    {
      cout << "Beware to the CFL condition." << endl;
    }
  }

  #define _DATA_FILE_CPP
  #endif
