# # Minimize a differentiable multivariate function using conjugate gradients.
# #
# # Usage: (X, fX, i) = mincg(X, f, length )
# # 
# # X       initial guess; may be of any type, including struct and cell array
# # f       the name or pointer to the function to be minimized. The function
# #         f must return two arguments, the value of the function, and it's
# #         partial derivatives wrt the elements of X. The partial derivative  
# #         must have the same type as X.
# # length  length of the run; if it is positive, it gives the maximum number of
# #         line searches, if negative its absolute gives the maximum allowed
# #         number of function evaluations. Optionally, length can have a second
# #         component, which will indicate the reduction in function value to be
# #         expected in the first line-search (defaults to 1.0).
# #
# # X       the returned solution
# # fX      vector of function values indicating progress made
# # i       number of iterations (line searches or function evaluations, 
# #         depending on the sign of "length") used at termination.
# %
# # The function returns when either its length is up, or if no further progress
# # can be made (ie, we are at a (local) minimum, or so close that due to
# # numerical problems, we cannot get any closer). NOTE: If the function
# # terminates within a few iterations, it could be an indication that the
# # function values and derivatives are not consistent (ie, there may be a bug in
# # the implementation of your "f" function).
# #
# # The Polack-Ribiere flavour of conjugate gradients is used to compute search
# # directions, and a line search using quadratic and cubic polynomial
# # approximations and the Wolfe-Powell stopping criteria is used together with
# # the slope ratio method for guessing initial step sizes. Additionally a bunch
# # of checks are made to make sure that exploration is taking place and that
# # extrapolation will not be unboundedly large.
#
# 
#  Adapted from original MATLAB version: 
#  Copyright (C) 2001 - 2010 by Carl Edward Rasmussen, 2010-01-03
#  http://www.gaussianprocess.org/gpml/code/matlab/Copyright
#  
#  Julia Version 2012
#     MIT License
#  

# Read options

function mincg (f::Function, X::Matrix, length)


  RHO = 0.01;                            # a bunch of constants for line searches
  SIG = 0.5;       # RHO and SIG are the constants in the Wolfe-Powell conditions
  INT = 0.1;    # don't reevaluate within 0.1 of the limit of the current bracket
  EXT = 3.0;                    # extrapolate maximum 3 times the current bracket
  MAX = 20;                         # max 20 function evaluations per line search
  RATIO = 100;                                      # maximum allowed slope ratio

  red=1


  i = 0;                                            # zero the run length counter
  ls_failed = 0;                             # no previous line search has failed
  fX = [];
  (f1,df1) = apply(f, (X, ));                   # get function value and gradient
  i = i + ((length<0)?1:0);                                      # count epochs?!
  s = -df1;                                        # search direction is steepest
  d1::Float64 = (-s'*s)[1];                                   # this is the slope
  z1::Float64 = (red/(1-d1));                       # initial step is red/(|s|+1)

  ls_failed = 0;
  z2::Float64 = 0;
  while i < abs(length)                                      # while not finished
    i = i + ((length>0)?1:0);                                # count iterations?!

    X0 = X; f0 = f1; df0 = df1;                   # make a copy of current values
    # println("Z2 $z2 $(s[1]) $(s[2])")
    X = X + z1*s;                                             # begin line search
    (f2, df2) = apply(f, (X, ));
    # println ("Result of Apply: $f2")
    i = i + (length<0);                                          # count epochs?!
    d2::Float64 = (df2'*s)[1];
    f3 = f1; d3 = d1; z3 = -z1;             # initialize point 3 equal to point 1
    if length>0 M = MAX; else M = min(MAX, -length-i); end
    success = 0; limit = -1;                     # initialize quanteties
    while true
      while (((f2 > (f1+z1*RHO*d1))|| (d2 > (-SIG*d1)) ) && (M>0))
        limit = z1;                                         # tighten the bracket  
        # println("Line 100: $f1 $f2")      
        
        if (f2 > f1)
          z2 = (z3 - (0.5*d3*z3*z3)/(d3*z3+f2-f3));               # quadratic fit
        else
          A = 6*(f2-f3)/z3+3*(d2+d3);                                 # cubic fit
          B = 3*(f3-f2)-z3*(d3+2*d2);
          z2 = ((sqrt(B*B-A*d2*z3*z3)-B)/A);      # numerical error possible - ok!
        end
        
        if isnan(z2) || isinf(z2)
          z2 = z3/2;                  # if we had a numerical problem then bisect
        end
        z2 = max(min(z2, INT*z3),(1-INT)*z3);  # don't accept too close to limits
        z1 = z1 + z2;                                           # update the step
        X = X + z2*s;
        # println("Inner z2 $z2 $(s[1]) $(s[2])")
        (f2, df2) = apply(f, (X, ));
        # println ("Result of Inner Apply: $f2")
        M = M - 1; i = i + (length<0);                           # count epochs?!
        d2::Float64 = (df2'*s)[1];
        z3 = z3-z2;                    # z3 is now relative to the location of z2
      end
      # println("Line 118: $f2 $d2  $((f1+z1*RHO*d1)[1])  $((-SIG*d1)[1])  rhs: $( d2 > ((-SIG*d1)[1]) )  lhs: $( f2 > ((f1+z1*RHO*d1)[1]) )"); flush(stdout_stream);
      if ( ( f2 > ((f1+z1*RHO*d1)) ) || ( d2 > ((-SIG*d1)) ) )
        break;                                                # this is a failure
      elseif (d2 > (SIG*d1))
        success = 1; break;                                             # success
      elseif M == 0
        break;                                                          # failure
      end
      A = 6*(f2-f3)/z3+3*(d2+d3);                      # make cubic extrapolation
      B = 3*(f3-f2)-z3*(d3+2*d2);
      z2::Float64 = (-d2*z3*z3/(B+sqrt(B*B-A*d2*z3*z3))); # num. error possible - ok!
      if !isreal(z2) || isnan(z2) || isinf(z2) || z2 < 0   # num prob or wrong sign?
        if limit < -0.5                               # if we have no upper limit
          z2 = z1 * (EXT-1);                 # the extrapolate the maximum amount
        else
          z2 = (limit-z1)/2;                                   # otherwise bisect
        end
      elseif (limit > -0.5) & (z2+z1 > limit)          # extraplation beyond max?
        z2 = (limit-z1)/2;                                               # bisect
      elseif (limit < -0.5) & (z2+z1 > z1*EXT)       # extrapolation beyond limit
        z2 = z1*(EXT-1.0);                           # set to extrapolation limit
      elseif z2 < -z3*INT
        z2 = -z3*INT;
      elseif (limit > -0.5) & (z2 < (limit-z1)*(1.0-INT))   # too close to limit?
        z2 = (limit-z1)*(1.0-INT);
      end
      f3 = f2; d3 = d2; z3 = -z2;                  # set point 3 equal to point 2
      z1 = z1 + z2; X = X + z2*s;                      # update current estimates
      (f2, df2) = apply(f, (X,));
      M = M - 1; i = i + (length<0);                             # count epochs?!
      d2::Float64 = (df2'*s)[1];
    end                                                      # end of line search
    
    if (success == 1)                                         # if line search succeeded
      f1 = f2; fX = [fX' f1]';

      # println("Iteration $i | Cost: $f1")
      s = (((df2'*df2-df1'*df2)[1])/((df1'*df1)[1]))*s - df2;      # Polack-Ribiere direction
      tmp = df1; df1 = df2; df2 = tmp;                         # swap derivatives
      d2 = (df1'*s)[1];
      if d2 > 0                                      # new slope must be negative
        s = -df1;                              # otherwise use steepest direction
        d2 = (-s'*s)[1];    
      end
      z1 = z1 * min(RATIO, d1/(d2-realmin()));          # slope ratio but max RATIO
      d1 = d2;
      ls_failed = 0;                              # this line search did not fail
    else
      # println("Line Search Failed, $ls_failed")
      X = X0; f1 = f0; df1 = df0;  # restore point from before failed line search
      if (ls_failed ==1) || (i > abs(length))          # line search failed twice in a row
        break;                             # or we ran out of time, so we give up
      end
      tmp = df1; df1 = df2; df2 = tmp;                         # swap derivatives
      s = -df1;                                                    # try steepest
      d1 = (-s'*s)[1];
      z1 = 1/(1-d1);                     
      ls_failed = 1;                                    # this line search failed
    end
    flush(stdout_stream)
  end
  print('\n');
  return (X, fX, i)
end

# (/){T<:Union(Float64,Float32,Complex128,Complex64)}(A::VecOrMat{T}, B::StridedMatrix{T}) = (B'\A')';

# (/){T1<:Real, T2<:Real}(A::StridedVecOrMat{T1}, B::StridedMatrix{T2}) = (/)(float64(A), float64(B))



