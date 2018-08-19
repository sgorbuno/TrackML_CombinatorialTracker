#ifndef PolynomFit_H
#define PolynomFit_H


/**
 * The class PolynomFit allows one to fit polynomial coefficients c0..cN
 *
 * with given measurements m_i, i=0, ... :
 *
 * m_i = c0*f0_i + c1*f1_i + .. + cN*fN_i + measurement_error
 *
 * The matrix inversion mathematics is home-made and can be replaced by ROOT TMatrixDSym utility 
 *
 */

class PolynomFit
{
public:

  PolynomFit(): fN(0), fA(0), fB(0)
  { }
  
  PolynomFit( int nCoefficients ): fN(0), fA(0), fB(0)
  {
    Reset(nCoefficients);
  }
  
  ~PolynomFit(){
    delete[] fA;
    delete[] fB;
  }
  
  static int invS( long double M[], int N );

  void Reset( int nCoefficients = -1 );

  void AddMeasurement( double f[], double m);
  
  int Fit( double Coefficients[] );    

private:

  int fN;
  long double *fA;
  long double *fB;
};

#endif
