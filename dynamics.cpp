#include"dynamics.h"

void Dynamics::TMatrix_calc(void){
  double t0=-t_window;
  double fm1=Field(t0-DT);
  Reset(fm1);
  EigenContainer r=diagonalize();
  //TMatrix[i]=r.evec.transpose()*DipoleMatrix*r.evec;
  Eigen::MatrixXd evaltmp=r.eval.replicate(1,2*Num_Dressed);
  Eigen::MatrixXd tmp=(evaltmp-evaltmp.transpose());
  Eigen::MatrixXd evec_pre=r.evec;

  double f0=Field(t0);
  Reset(f0);
  r=diagonalize();
  //EigenContainer eig=
  Eigen::MatrixXd evec=r.evec;
  evaltmp=r.eval.replicate(1,2*Num_Dressed);
  Eigen::MatrixXd tmp0=(evaltmp-evaltmp.transpose());
  Eigen::MatrixXd prediff=(evec-evec_pre)/(f0-fm1);
  for (int i=0;  i < TOTAL_STEPS+1; i++) {
    const double t = DT*(i)-t_window;
    const double fp1=Field(t+DT);
    Reset(fp1);
    r=diagonalize();
    const Eigen::MatrixXd evec_post=r.evec;
    evaltmp=r.eval.replicate(1,2*Num_Dressed);
    tmp=(evaltmp-evaltmp.transpose());
    const Eigen::MatrixXd postdiff=(evec_post-evec)/(fp1-f0);//*dField(t-DT);
    TMatrix[i]=evec.transpose()*(prediff+postdiff)/2;
    TMatrix[i].array()*=-((t*tmp0).array().cos()+std::complex<double>(0,1)*(t*tmp0).array().sin())*dField(t);
    //((t*tmp).array().cos()+std::complex<double>(0,1)*(t*tmp).array().sin())*dField(t)
    std::cerr<<TMatrix[i]<<std::endl;
    prediff=postdiff;
    evec=evec_post;
    f0=fp1;
    tmp0=tmp;
  }
  //std::cout<<2*Num_Dressed-(Num_Dressed/2)*2<<'\t'<<(Num_Dressed/2)*2-1<<std::endl;
}


#ifdef TEST
int main(int argc, char const *argv[]) {
  const double DeltaE=std::atof(argv[1]);
  const double Omega=std::atof(argv[2]);
  const double FieldStrength=std::atof(argv[3]);
  const unsigned int Num_Dressed=std::atoi(argv[4]);
  const double pulse_width=std::atof(argv[5]);
  Dynamics dcalc(-DeltaE/2,DeltaE/2,Omega,FieldStrength,Num_Dressed,pulse_width,1,M_PI/Omega/20);
  dcalc.propagate();
  for (int i = 0; i < dcalc.GetSteps()+1; i++) {
    const double t=M_PI/Omega/20*(i-dcalc.GetSteps()/2);
    std::cout<<t<<'\t'<<dcalc.Field(t)<<'\t'<<dcalc.GetWFN(i).transpose().array().abs2()<<std::endl;
  }
  return 0;
}
#endif
