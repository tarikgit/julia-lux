### Brock-Mirman Model
### Method: Projection method
### by Tarik Ocaktan (October 2008 - Amsterdam) 	### 

# Define some constants
param pi  := 3.14159265358979;
param eps := 1e-6;

# Parameterization #
param dAlpha 	:= 0.3;
param dBeta 	:= 0.99;
param dDelta 	:= 1;
param dGamma 	:= 1;
param dRho 	:= 0.9;
param dSigma 	:= 0.01;

param dKss := ( 1/(dAlpha*dBeta) )^(1/(dAlpha-1));
param dCss := dKss^dAlpha - dKss;
param dThetass := 1;

# run time
param runTime default 0;

# Remember that m > n
param m := 11;  # number of Chebyshev nodes
param n := 2;  # degree of Chebyshev polynomial

# Bounds of capital 
param dKmin := 0.8*dKss;  # minimum state value
param dKmax := 1.2*dKss;  # maximum state value
# Bounds of productivity  
param dThetamin := 0.8*dThetass;  # minimum state value
param dThetamax := 1.2*dThetass;  # maximum state value

# Define parameters values
set   SM := 1..m;       # index set of states
set   SN := 0..(n-1);

#Calculating the Chebyshev nodes and polynomials at these nodes
 # Chebyshev nodes on [0,1]
 param z {k in SM} 	  := -cos( ((2*k-1)/(2*m))*pi ); 	
 # Chebyshev nodes on capital and productivity
 param dK {k in SM} 	  := ( z[k]+1 )*( (dKmax-dKmin)/2 ) + dKmin;		
 param dTheta {k in SM}	  := ( z[k]+1 )*( (dThetamax-dThetamin)/2 ) + dThetamin;
 # Ordinary Polynomial (Tensor Product Bases - see Judd 237)
 param mKZ {p in SN, k in SM} := ( log(dK[k]) )^p * ( dZ[k] )^p;

# This is the function that we want to approximate
 # current value function's values at KNodes
  param f {k in SM, l in SM} := dAlpha*dBeta*;
#initialize variables
# var coefs {i in SN};
# let {i in SN} coefs[i] := ( sum {h in SM} (f[h]*T[i,h]) ) 
#     	       		   / ( sum {u in SM} (T[i,u]*T[i,u]) );

#display 2*coefs[0];
#display {i in 1..(n-1)} coefs[i];

#### END ###

#for {index in SN}{
#printf "%f\n", coefs[index] >> coefficients.txt;
#}

display dTheta;

