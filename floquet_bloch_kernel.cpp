#include"floquet_bloch_kernel.h"

void Floquet_Bloch_Kernel::calc(void){
  for (int i = 0; i < NKPoints; i++) {
    Dynamics dcalc(evgrid(i),ecgrid(i),Omega0,height_wave,Num_Dressed,sigma_wave,1.,TOTAL_STEPS);
    dcalc.propagate();
    for (int j = 0; j < TOTAL_STEPS+1; j++) {
      Population[j]+=e_grid->MapGrid(dcalc.GetWFN(j),dcalc.GetE_diff(j),Omega0);
    }
  }
}

#ifdef TEST
int main(int argc, char const *argv[]) {
  const std::string fname(argv[1]);
  const int Num_Dressed=std::atoi(argv[2]);
  const int TOTAL_STEPS=std::atoi(argv[3]);
  const int Grid_Size=std::atoi(argv[4]);
  Floquet_Bloch_Kernel fb_calc(fname,Num_Dressed,TOTAL_STEPS,Grid_Size);
  fb_calc.calc();
  const double DT=fb_calc.GetDT();
  for(int i=0; i<TOTAL_STEPS; i++){
    Eigen::VectorXd population=fb_calc.GetPopulation(i);
    for(int j=0; j<population.size(); j++){
      std::cout<<DT*(i-TOTAL_STEPS/2)<<'\t'<<fb_calc.GetEnergy(j)<<'\t'<<population(j)<<std::endl;
    }
    std::cout<<std::endl;
  }
  return 0;
}
#endif
