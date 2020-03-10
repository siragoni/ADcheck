// #include "Riostream.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TVectorD.h"

void EMDcorrectionMatrix()
{

  /* -
   * - Pileup and ZDC efficiencies
   * -
   */
  Double_t pA = 0.034;
  Double_t pC = 0.036;
  Double_t eA = 0.93;
  Double_t eC = 0.93;


  Double_t det1;
  TMatrixD CorrectionMatrix(4, 4);
  TArrayD  NullValues(16);
  for (Int_t i = 0; i < 16; i++) {
     NullValues[i] = 0;
  }
  CorrectionMatrix.SetMatrixArray(NullValues.GetArray());
  CorrectionMatrix.Print();

  /* -
   * -
   * - N_0N0N           M_0N0N
   * - N_0NXN = c_{ij}  M_0NXN
   * - N_XN0N           M_XN0N
   * - N_XNXN           M_XNXN
   * -
   */

  TArrayD  c_ij(16);
  /* -
   * - 0N0N
   */
  c_ij[0]  = (1. - pA)*(1. - pC);
  c_ij[1]  = (1. - eA)*(1. - pA)*(1. - pC);
  c_ij[2]  = (1. - eC)*(1. - pA)*(1. - pC);
  c_ij[3]  = (1. - eA)*(1. - eC)*(1. - pA)*(1. - pC);
  /* -
   * - 0NXN
   */
  // c_ij[4]  = (pA)*(1. - pC);
  // c_ij[5]  = (pC + eA * pA - pA * pC - eA * pC);
  // c_ij[6]  = (1. - eC)*pA*eA;
  // c_ij[7]  = (1. - eC)*(1. - pC)*( eA + (1. - eA) * pA );
  c_ij[4]  = (pA)*(1. - pC);
  c_ij[5]  = 1. - (   (1. - eA)*(1. - pA)*(1. - pC) + eA*pC +  (1. - eA)*pA*pC  );
  c_ij[6]  = (1. - eC)*pA*eA;
  c_ij[7]  = (1. - eC)*(1. - pC)*( eA + (1. - eA) * pA );
  /* -
   * - XN0N
   */
  // c_ij[8]  = (pC)*(1. - pA);
  // c_ij[9]  = (1. - eA)*pC*eC;
  // c_ij[10] = (pA + eC * pC - pC * pA - eC * pA);
  // c_ij[11] = (1. - eA)*(1. - pA)*( eC + (1. - eC) * pC );
  c_ij[8]  = (pC)*(1. - pA);
  c_ij[9]  = (1. - eA)*pC*eC;
  c_ij[10] = 1. - (   (1. - eC)*(1. - pC)*(1. - pA) + eC*pA +  (1. - eC)*pC*pA  );
  c_ij[11] = (1. - eA)*(1. - pA)*( eC + (1. - eC) * pC );
  /* -
   * - XNXN
   */
  c_ij[12] = (pC)*(pA);
  c_ij[13] = ( eA * pC + (1. - eA) * pC * pA );
  c_ij[14] = ( eC * pA + (1. - eC) * pA * pC );
  c_ij[15] = 1 - (    (1. - eA)*(1. - pA)*( eC + (1. - eC) * pC ) + (1. - eC)*(1. - pC)*( eA + (1. - eA) * pA ) + (1. - eA)*(1. - eC)*(1. - pA)*(1. - pC)   );

  CorrectionMatrix.SetMatrixArray( c_ij.GetArray() );
  CorrectionMatrix.Print();


  CorrectionMatrix.InvertFast(&det1);
  CorrectionMatrix.Print();

  // Get the maximum off-diagonal matrix value . One way to do this is to set the
  // diagonal to zero .

  // TMatrixD U1(CorrectionMatrix,TMatrixD::kMult,H_square);
  // TMatrixDDiag diag1(U1); diag1 = 0.0;
  // const Double_t U1_max_offdiag = (U1.Abs()).Max();
  // cout << "  Maximum off-diagonal = " << U1_max_offdiag << endl;
  // cout << "  Determinant          = " << det1 <<endl;

  // cout << "correction Matrix after inversion" << endl << CorrectionMatrix << endl;
}
