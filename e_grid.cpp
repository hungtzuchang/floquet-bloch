#include"e_grid.h"
#include<cmath>
#include<cstdlib>

Eigen::VectorXd E_Grid::MapGrid(const Eigen::VectorXcd &WFN, const std::vector<double> &E_diff, const double omega_in){
  Eigen::VectorXd map=Eigen::VectorXd::Zero(grid.size());
  Eigen::VectorXd population=WFN.array().abs2();
  const int wfn_num_dressed=WFN.size()/2;
  const int wfn_max_dressed=wfn_num_dressed/2;
  for (int i = 0; i < wfn_num_dressed; i++) {
    const double E_dressing=(i-wfn_max_dressed)*omega_in;
    const double E_valence=E_diff[0]+E_dressing;
    const double E_conduction=E_diff[1]+E_dressing;
    const int gn_valence=std::round((E_valence+energy_range/2)/dE);
    const int gn_conduction=std::round((E_conduction+energy_range/2)/dE);
    if((gn_valence<0) || (gn_conduction>=grid.size())){
      std::cerr<<"Grid energy range too small !!!!!"<<std::endl;
      std::abort();
    }
    map(gn_valence)+=population(i);
    map(gn_conduction)+=population(i+wfn_num_dressed);
  }
  return map;
}

#ifdef TEST
#include"dynamics.h"
int main(int argc, char const *argv[]) {
  const double DeltaE=std::atof(argv[1]);
  const double Omega=std::atof(argv[2]);
  const double FieldStrength=std::atof(argv[3]);
  const unsigned int Num_Dressed=std::atoi(argv[4]);
  const double pulse_width=std::atof(argv[5]);
  Dynamics dcalc(-DeltaE/2,DeltaE/2,Omega,FieldStrength,Num_Dressed,pulse_width,1,M_PI/Omega/20);
  dcalc.propagate();
  E_Grid energy_grid(-DeltaE/2,DeltaE/2,Omega,FieldStrength,Num_Dressed,20);
  for (int i = 0; i < dcalc.GetSteps()+1; i++) {
    Eigen::VectorXd population_e=energy_grid.MapGrid(dcalc.GetWFN(i),dcalc.GetE_diff(i),Omega);
    for(int j=0; j<population_e.size(); j++){
      std::cout<<M_PI/Omega/20*i<<'\t'<<energy_grid.GetGrid(j)<<'\t'<<population_e(j)<<std::endl;
    }
    std::cout<<'\n';
  }
  return 0;
  return 0;
}
#endif
