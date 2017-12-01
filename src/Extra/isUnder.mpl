# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isUnder.mpl                                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if all polynomials in a regular system or constructible set   #
# (including inequations) only contain variables that are strictly less   #
# than v in R.                                                            #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   isUnder(rs, v, R)                                                     #
#   isUnder(cs, v, R)                                                     #
#                                                                         #
# INPUT                                                                   #
#   rs ... Regular system                                                 #
#   cs ... Constructible set                                              #
#   v .... Variable                                                       #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   True if the equations and inequations of rs or cs only contain the    #
#   variables of R  strictly less that v, false otherwise.                #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a, b]):                                     #
#   >                                                                     #
#   > cs := GeneralConstruct([a*x+b^2], [a, b], R):                       #
#   > isUnder(cs, x, R);                                                  #
#         false                                                           #
#   > cs := GeneralConstruct([a-1], [a-1], R):                            #
#   > isUnder(cs, x, R);                                                  #
#         true                                                            #
# ======================================================================= #
# ======================================================================= #
isUnder := overload(
    [
        proc(rs::TRDrs, v::name, R::TRDring, $) :: truefalse;
            
            option overload;
            
            local rc :: TRDrc,
                  vars :: set;
            
            # Ensure v is the largest variable in R
            if not isGreatestVariable(v, R) then
                error "%1 is not the greatest variable of R", v;
            end if;
            
            rc := RC_CST:-RepresentingChain(rs, R);
            vars := indets(RC_CST:-RepresentingInequations(rs, R));
            vars := `union`(vars, indets(RC:-Inequations(rc, R)));
            vars := `union`(vars, indets(RC:-Equations(rc, R)));
            
           return not evalb(v in vars);
            
        end,
        
        
        proc(cs::TRDcs, v::name, R::TRDring, $) :: truefalse;
            
            option overload;
            
            local lrs :: TRDlrs,
                  rs :: TRDrs;
            
            lrs := RC_CST:-RepresentingRegularSystems(cs, R);
            
            for rs in lrs do
                if not isUnder(rs, v, R) then
                    return false;
                end if;
            end do;
            
            return true;
            
        end
    ]
);