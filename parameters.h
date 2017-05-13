// Define parameters in SBE
#ifndef PARAMETERS_H
#define PARAMETERS_H
#include<fstream>
#include<string>
#include<Eigen/Dense>

class Parameters
{
public:
	Parameters(const std::string &);
protected:
	double kspacing;
	unsigned int NKPoints;
	double gamma;
	double sigma_wave; // temporal width of E-field
	double height_wave;
	double mv;
	double mc;
	double Eg;
	double Omega0;
	double Vq;
	Eigen::VectorXd evgrid;
  Eigen::VectorXd ecgrid;
	Eigen::VectorXd ksquared;
};
#endif
