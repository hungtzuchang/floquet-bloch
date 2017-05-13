#include"parameters.h"
#include<cmath>
#include<iostream>
#include<cstdlib>

Eigen::Vector3i Get_k_vector(const unsigned int NKpoints, const int index)
{
  Eigen::Vector3i a;
  int ktmp1=index/NKpoints;
  a(2)=index%NKpoints;
  a(0)=ktmp1/NKpoints;
  a(1)=ktmp1%NKpoints;
  return a;
}

Parameters::Parameters(const std::string &filename)
{
  std::string command;
  std::ifstream file_in(filename.c_str());
  if(!file_in.is_open()){
    std::cerr<<"ERROR OPENING FILE"<<std::endl;
    std::abort();
  }

  bool chk_NKPoints=0;
  bool chk_mv=0;
  bool chk_mc=0;
  bool chk_Eg=0;
  bool chk_Omega0=0;
  bool chk_sigma_wave=0;
  bool chk_height_wave=0;
  bool chk_gamma=0;
	while(!file_in.eof()){
		std::getline(file_in, command);
		if(command=="$NKPOINTS"){
			file_in>>NKPoints;
      chk_NKPoints=1;
			continue;
		}
    if(command=="$mv"){
			file_in>>mv;
      chk_mv=1;
			continue;
		}
    if(command=="$mc"){
			file_in>>mc;
      chk_mc=1;
			continue;
		}
    if(command=="$Eg"){
      file_in>>Eg;
      chk_Eg=1;
      continue;
    }
    if(command=="$Omega0"){
			file_in>>Omega0;
      chk_Omega0=1;
			continue;
		}
    if(command=="$PULSE_WIDTH"){
			file_in>>sigma_wave;
      chk_sigma_wave=1;
			continue;
		}
    if(command=="$FIELD_STRENGTH"){
			file_in>>height_wave;
      chk_height_wave=1;
			continue;
		}
    if(command=="$gamma"){
			file_in>>gamma;
      chk_gamma=1;
			continue;
		}
  }

  file_in.close();

  if (!chk_NKPoints) {
      NKPoints=10;
      std::cerr<<"#WARNING: NKPoints NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: NKPoint=10"<<std::endl;
  }
  if (!chk_mv) {
      mv=1;
      std::cerr<<"#WARNING: mv NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: mv=1"<<std::endl;
  }
  if (!chk_mc) {
      mc=1;
      std::cerr<<"#WARNING: mc NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: mc=1"<<std::endl;
  }
  if (!chk_Eg) {
      Eg=1;
      std::cerr<<"#WARNING: Eg NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: Eg=1"<<std::endl;
  }
  if (!chk_Omega0) {
      Omega0=Eg;
      std::cerr<<"#WARNING: Omega0 NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: Omega0=Eg"<<std::endl;
  }
  if (!chk_sigma_wave) {
      sigma_wave=5/Omega0;
      std::cerr<<"#WARNING: PULSE_WIDTH NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: PULSE_WIDTH=5/Omega0"<<std::endl;
  }
  if (!chk_height_wave) {
      height_wave=0.1;
      std::cerr<<"#WARNING: FIELD_STRENGTH NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: FIELD_STRENGTH=0.1"<<std::endl;
  }
  if (!chk_gamma) {
      gamma=0;
      std::cerr<<"#WARNING: gamma NOT FOUND"<<std::endl;
      std::cerr<<"#USING DEFAULT: gamma=0"<<std::endl;
  }
  const double UppermostE=Omega0-Eg+0.1;
  const double msmall= (mv<mc) ? mv : mc;
  kspacing=std::sqrt(2*msmall*UppermostE)/NKPoints;
  Vq=M_PI*4*std::pow(kspacing,3);

  ksquared=Eigen::VectorXd::Zero(NKPoints);
  for(int i=0; i<NKPoints; i++){
    ksquared(i)=Get_k_vector(NKPoints,i).array().square().sum()*kspacing*kspacing;
  }
  ecgrid=ksquared/(2*mc)+Eigen::VectorXd::Constant(NKPoints,Eg);
  evgrid=-ksquared/(2*mv);
}

#ifdef TEST
int main(int argc, char **argv)
{
  std::string a(argv[1]);
  Parameters p(a);
  return 0;
}
#endif
