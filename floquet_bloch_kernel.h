#ifndef FLOQUET_BLOCH_KERNEL_H
#define FLOQUET_BLOCH_KERNEL_H
#include"parameters.h"
#include"e_grid.h"
#include"dynamics.h"

class Floquet_Bloch_Kernel : public Parameters {
public:
  Floquet_Bloch_Kernel(const std::string &fname, const int ndressed, const int t_steps, const int gd):Parameters(fname){
    TOTAL_STEPS=t_steps;
    Population.resize(TOTAL_STEPS+1);
    Num_Dressed=ndressed;
    Grid_Density=gd;
    e_grid=new E_Grid(evgrid.minCoeff(),ecgrid.maxCoeff(),Omega0,height_wave,Num_Dressed,Grid_Density);
    for(int i=0; i<TOTAL_STEPS+1; i++){
      Population[i]=Eigen::VectorXd::Zero(e_grid->GetSize());
    }
  };
  void calc(void);
  Eigen::VectorXd GetPopulation(const int i){return Population[i];};
  double GetEnergy(const int i){return e_grid->GetGrid(i);};
  double GetDT(void){return 12*M_PI/Omega0/TOTAL_STEPS;};
  ~Floquet_Bloch_Kernel(){delete e_grid;};
protected:
  std::vector<Eigen::VectorXd> Population;
  int TOTAL_STEPS;
  E_Grid *e_grid;
  int Num_Dressed;
  int Grid_Density;
};

#endif
