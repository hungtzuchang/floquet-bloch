#ifndef DYNAMICS_H
#define DYNAMICS_H
#include"floquet.h"
#include<cmath>
#include<complex>
#include<vector>

class Dynamics : public Floquet_Hamiltonian {
public:
  Dynamics():Floquet_Hamiltonian(){};
  Dynamics(const double var11, const double var12, const double var2, const double var3, const int var4, const double var5,  const double var6, const double var7)
  :Floquet_Hamiltonian(var11,var12,var2,var3,var4){
    peak_field=var3;
    dipole_strength=var6;
    pulse_duration=var5;
    DT=var7;
    TOTAL_STEPS=0;
    Init();
  };
  Dynamics(const double var11, const double var12, const double var2, const double var3, const int var4, const double var5,  const double var6, const int var7)
  :Floquet_Hamiltonian(var11,var12,var2,var3,var4){
    peak_field=var3;
    dipole_strength=var6;
    pulse_duration=var5;
    TOTAL_STEPS=var7;
    Init();
  };
  void propagate(void){
    for (int i = 0; i < TOTAL_STEPS; i++) {
      const bool chk_converge=cnupdate(i);
      if(chk_converge)
        std::cerr<<"Not Converged Step: "<<i<<std::endl;
    }
  };
  Eigen::VectorXcd GetWFN(const int i){return WFN[i];};
  std::vector<double> GetE_diff(const int i){return E_diff[i];};
  int GetSteps(void){return TOTAL_STEPS;};
//protected:
  double dipole_strength;
  double pulse_duration;
  double peak_field;
  std::vector<Eigen::VectorXcd> WFN;
  double DT;
  double t_window;
  int TOTAL_STEPS;
  std::vector<Eigen::MatrixXcd> TMatrix;
  std::vector<std::vector<double> > E_diff;
  void Init(void){
    t_window=6*M_PI/Omega;
    if(TOTAL_STEPS==0){
      TOTAL_STEPS=2*std::ceil(t_window/DT);
    } else {
      DT=2*t_window/TOTAL_STEPS;
    }

    WFN.resize(TOTAL_STEPS+1);
    TMatrix.resize(TOTAL_STEPS+1);
    E_diff.resize(TOTAL_STEPS+1);
    for (int i = 0; i < TOTAL_STEPS+1; i++) {
      WFN[i]=Eigen::VectorXcd::Zero(2*Num_Dressed);
      TMatrix[i]=Eigen::MatrixXcd::Zero(2*Num_Dressed,2*Num_Dressed);
      E_diff[i].resize(2);
    }
    WFN[0]((Num_Dressed/2)*2-1)=std::complex<double>(1.,0.);
    TMatrix_calc();
  };
  double Field(const double t){
    const double sigma2=(pulse_duration*pulse_duration);
    return peak_field*std::exp(-t*t/(2*sigma2));
  };
  double dField(const double t){
    const double sigma2=(pulse_duration*pulse_duration);
    return -peak_field*t/sigma2*std::exp(-t*t/(2*sigma2));
  };
  void TMatrix_calc(void);
  void Reset(const double var3){
    FieldStrength=var3;
    SetOffDiag();
  };
  Eigen::VectorXcd dtdWFN(const int step, const Eigen::VectorXcd &wfnt){
    return TMatrix[step]*wfnt;
  };
  bool cnupdate(const int step){
    Eigen::VectorXcd wfnt1=Eigen::VectorXcd::Zero(2*Num_Dressed);
    Eigen::VectorXcd wfnt=WFN[step]+DT*dtdWFN(step,WFN[step]);
    for (int i = 0; i < 1000; i++) {
      wfnt1=WFN[step]+DT/2*(dtdWFN(step,WFN[step])+dtdWFN(step,wfnt));
      if((wfnt-wfnt1).array().abs().sum()<1e-8){
        WFN[step+1]=wfnt1;
        return 0;
      }
      wfnt=wfnt1;
    }
    return 1;
  }
};
#endif
