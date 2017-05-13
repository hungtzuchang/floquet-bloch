#include"floquet.h"
#include<cstdlib>

#ifdef TEST
int main(int argc, char const *argv[]) {
  const double DeltaE=std::atof(argv[1]);
  const double Omega=std::atof(argv[2]);
  const double FieldStrength=std::atof(argv[3]);
  const unsigned int Num_Dressed=std::atoi(argv[4]);
  Floquet_Hamiltonian H_F(-DeltaE/2,DeltaE/2,Omega,FieldStrength,Num_Dressed);
  std::cout<<H_F.Get()<<"\n\n";
  std::cout<<H_F.diagonalize().eval<<'\n'<<std::endl;

  for(int i=0; i<30; i++){
    const double f=i*FieldStrength/30;
    H_F.Set(-DeltaE/2,DeltaE/2,Omega,f,Num_Dressed);
    std::cerr<<f<<'\t'<<H_F.diagonalize().eval.transpose()<<std::endl;
    if(i==0)
      std::cout<<H_F.diagonalize().eval<<std::endl;
  }
  std::cout<<H_F.diagonalize().evec<<std::endl;
  return 0;
}
#endif
