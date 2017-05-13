#ifndef Floquet_H
#define Floquet_H
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include<iostream>

struct EigenContainer {
  //EigenContainer(const EigenContainer &a){evec=a.evec; eval=a.eval;};
  Eigen::MatrixXd evec;
  Eigen::VectorXd eval;
};

class Floquet_Hamiltonian {
public:
  Floquet_Hamiltonian(){};
  Floquet_Hamiltonian(const double var11, const double var12, const double var2, const double var3, const int var4){Set(var11,var12,var2,var3,var4);};
  void Set(const double var11, const double var12, const double var2, const double var3, const int var4){
    Ev=var11;
    Ec=var12;
    Omega=var2;
    FieldStrength=var3;
    Num_Dressed=var4;
    H_F=Eigen::MatrixXd::Zero(2*Num_Dressed,2*Num_Dressed);
    Eigen::VectorXd H_diag=Eigen::VectorXd::Zero(Num_Dressed);
    for (int i = 0; i < Num_Dressed; i++) {
      H_diag(i)=Omega*(i-((double) Num_Dressed-1)/2);
    }
    H_F.topLeftCorner(Num_Dressed,Num_Dressed).diagonal()=H_diag+Ev*Eigen::VectorXd::Ones(Num_Dressed);
    H_F.bottomRightCorner(Num_Dressed,Num_Dressed).diagonal()=H_diag+Ec*Eigen::VectorXd::Ones(Num_Dressed);
    //H_F.diagonal()<<H_diag+Ev*Eigen::VectorXd::Ones(Num_Dressed),H_diag+Ec*Eigen::VectorXd::Ones(Num_Dressed);
    SetOffDiag();
  };
  Eigen::MatrixXd Get(void){return H_F;};
  EigenContainer diagonalize(void){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_F);
    EigenContainer r;
    r.eval=es.eigenvalues();
    r.evec=es.eigenvectors();
    return r;
  };
protected:
  double Ec;
  double Ev;
  double Omega;
  double FieldStrength;
  int Num_Dressed;
  Eigen::MatrixXd H_F;
  void SetOffDiag(void){
    Eigen::MatrixXd H_offdiag=Eigen::MatrixXd::Zero(Num_Dressed,Num_Dressed);
    H_offdiag.diagonal(-1)=FieldStrength*Eigen::VectorXd::Ones(Num_Dressed-1);
    H_offdiag.diagonal(+1)=H_offdiag.diagonal(-1);
    H_F.topRightCorner(Num_Dressed,Num_Dressed)=H_offdiag;
    H_F.bottomLeftCorner(Num_Dressed,Num_Dressed)=H_offdiag;
  };
};
#endif
