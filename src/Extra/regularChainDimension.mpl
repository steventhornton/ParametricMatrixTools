# ======================================================================= #
# ======================================================================= #
#                                                                         #
# regularChainDimension.mpl                                               #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 1/2017                                                 #
#                                                                         #
# Compute the dimension of a regular chain. The dimension is              #
# nVars - (# of equations containing one or more of the largest nVars     #
# variables of R). The input regular chain is assumed to be linear in     #
# each of the nVars largest variables of R.                               #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   regularChainDimension(rc, n, R)                                       #
#                                                                         #
# INPUT                                                                   #
#   rc ... Regular chain                                                  #
#   n .... Number of variables                                            #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a positive integer corresponding to the dimension of the      #
#   input regular chain.                                                  #
# ======================================================================= #
# ======================================================================= #
regularChainDimension := proc(rc::TRDrc, n::nonnegint, R::TRDring, $)

    local x :: set(name),
          eqns :: list(polynom),
          eqn :: polynom,
          i :: nonnegint;

    # Get the largest 'nVars' variables of R
    x := {op((R['variables'])[1..n])};

    # Get the equations
    eqns := RC:-Equations(rc, R);

    # Initialize the dimension at nVars
    i := n;

    for eqn in eqns do
        if not evalb(nops(indets(eqn) intersect x) = 0) and not evalb(eqn = 0) then
            i := i - 1;
        end if;
    end do;

    return i;

end proc;