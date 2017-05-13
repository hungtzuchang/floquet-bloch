#ifndef E_GRID_H
#define E_GRID_H
#include"floquet.h"

class E_Grid : public Floquet_Hamiltonian {
public:
  E_Grid():Floquet_Hamiltonian(){};
  E_Grid(const double var11, const double var12, const double var2, const double var3, const int var4, const int var5):Floquet_Hamiltonian(var11, var12, var2, var3, var4){
    grid_density=var5;
    grid.resize((5+Num_Dressed)*grid_density);
    const Eigen::VectorXd evaltmp=diagonalize().eval;
    const double energy_difference=evaltmp(2*Num_Dressed-(Num_Dressed/2)*2)-evaltmp((Num_Dressed/2)*2-1);
    energy_range=energy_difference+Omega*(3+Num_Dressed);
    dE=energy_range/grid.size();
    for(int i=0; i<grid.size(); i++){
      grid[i]=dE*i-energy_range/2;
    }
  };
  double GetE_Range(void){return energy_range;}
  double GetGrid(const int i){return grid[i];};
  int GetSize(void){return grid.size();};
  Eigen::VectorXd MapGrid(const Eigen::VectorXcd &WFN, const std::vector<double> &E_diff, const double omega_in);
protected:
  std::vector<double> grid;
  int grid_density;
  double energy_range;
  double dE;
};
#endif
