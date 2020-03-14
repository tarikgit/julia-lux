
# Define Chebyshev state values
param 	pi   := 3.141592654;
param  	kmin := 0.5;         # minimum state value
param  	kmax := 1.5;         # maximum state value
param  	nk   :=  20;         # number of states
# grid of Chebyshev nodes
set     N    := 1..nk;       # index set of states

# set of Chebyshev state variables
param  K {i in N} := kmin + (kmax-kmin)/(2)*(1 - cos((2*i-1)/(2*nk)*pi));  

param  uk   :=  101;         	# number of states
set    UN   := 1..uk;       	# index set of states
param  U {i in UN} := kmin + (kmax-kmin)/(uk-1)*(i-1);  # set of states 

# Define parameters values
param alpha := 0.33;
param  beta := 0.99;
param delta := 0.023;

param F {i in N} 	:= K[i] + ((1-beta)/(alpha*beta))*K[i]^alpha;
param UF {i in UN} 	:= U[i] + ((1-beta)/(alpha*beta))*U[i]^alpha;

# Define degree of polynomial
param na  := 6;
set A     := 1..(na+1);
var a {i in A};
var c {i in N} >= 0;

maximize TV : sum {i in N} ( log(c[i]) + beta*( sum {j in A} a[j]*( F[i]-c[i] )^(j-1) ) );

 subject to  cL {i in N} : F[i] - c[i] >= kmin;

 subject to  cU {i in N} : F[i] - c[i] <= kmax;


param Vi {i in N} default 0;

param Va_plus_1 {k in UN}; 

param Va {k in UN}; 



# least squares of a (coefficients of a)

minimize LeastSQ_a: sum {i in N} ( sum {j in A} a[j]*K[i]^(j-1) - Vi[i] )^2;

 subject to Vgrad {i in N} : sum {j in A: j >= 2} a[j]*(j-1)*(K[i]^(j-2)) >= 0;

 subject to Vhess {i in N} : sum {j in A: j >= 3} a[j]*(j-1)*(j-2)*(K[i]^(j-3)) <= 0;


######################################

problem Bellman: TV, c, cL, cU; 				# a is not a variable 

problem Fitting: a, LeastSQ_a, Vgrad, Vhess; 	# c is not a variable

######################################

param _mu default 1e-6;
param MAX_count default 5;
param counter, default 0;

let {i in A} a[i] := 0; 
let {i in N} c[i] := 2;

repeat { option solver_msg  0;

  solve Bellman;

   let {i in N} Vi[i] := log(c[i]) + beta*(sum {j in A} a[j]*(F[i]-c[i])^(j-1)); 
   let {k in UN} Va[k] := sum {j in A} a[j]*U[k]^(j-1); # uniform grid to test when to stop
  
  option solver_msg  0;

  solve Fitting; # Vi are parameters in the Fitting

   let {k in UN} Va_plus_1[k] := sum {j in A} a[j]*U[k]^(j-1);
   let counter := counter + 1;

} while ( max {k in UN} abs(Va_plus_1[k]-Va[k]) > _mu) && (counter < MAX_count);


########## END





