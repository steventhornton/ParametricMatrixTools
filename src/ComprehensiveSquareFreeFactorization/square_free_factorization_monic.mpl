# ======================================================================= #
# ======================================================================= #
#                                                                         #
# square_free_factorization_monic.mpl                                     #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 29/2017                                                #
#                                                                         #
# Compute the square-free factorization of a parametric, univariate       #
# polynomial that is monic in its main variable. The computation is can   #
# return a complete case disussion such that each branch forms a          #
# partition of the input regular system. Alternatively, the computation   #
# can be done in the sense of kalkbrener.                                 #
#                                                                         #
# INPUT                                                                   #
#   p ............ Polynomial                                             #
#   v ............ Variable                                               #
#   rs ........... Regular system                                         #
#   R ............ Polynomial ring                                        #
#   opt_lazard ... Output option, if true the square-free factorization   #
#                  is computed in the sense of lazard, otherwise it is in #
#                  the sense of kalkbrener.                               #
#                                                                         #
# OUTPUT                                                                  #
#   A list of with elements of the form                                   #
#       [lp_i, rs_i]                                                      #
#   where lp_i is a list with elements of the form                        #
#       [p_j, n_j]                                                        #
#   and rs_i is a regular system.                                         #
#   p = m_i*product(p_j^n_j) and p_j are the square-free factors in the   #
#   zero set of rs_i and m_i is some rational function of the parameters. #
#                                                                         #
# LICENSE                                                                 #
#   This program is free software: you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation, either version 3 of the License, or     #
#   any later version.                                                    #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program.  If not, see http://www.gnu.org/licenses/.   #
# ======================================================================= #
# ======================================================================= #
square_free_factorization_monic := module()

    export ModuleApply;

    local implementation,
          monic_sqf_mod_rs_lazard,
          monic_sqf_mod_rs_kalkbrener,
          monic_sqf_mod_rc_h;

    ModuleApply := proc(p::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, opt_lazard::truefalse, $)
        return implementation(p, v, rs, R, opt_lazard);
    end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Check the input and process the input such that it matches the          #
# specifications of the SquarefreeFactorization_monic method.             #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as square_free_factorization_monic                               #
# ----------------------------------------------------------------------- #
implementation := proc(p_in::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, opt_lazard, $)
    
    local p :: polynom;
    
    p := expand(p_in);
    
    if not R['variables'][1] = v then
        error "v must be the greatest variable of R";
    end if;
    
    if not RC:-TRDuniv_lcoeff(p, v) = 1 then
        error "p must be a monic polynomial in v";
    end if;
    
    # If p is constant w.r.t. v
    if degree(p, v) = 0 then
        return [[sqrfree(p)[2], rs]];
    end if;
    
    # Case where p does not contain any parameters
    if v in indets(p) and nops(indets(p)) = 1 then
        return [[sqrfree(p, v)[2], rs]];
    end if;
    
    if opt_lazard then
        monic_sqf_mod_rs_lazard(p, v, rs, R);
    else
        return monic_sqf_mod_rs_kalkbrener(p, v, rs, R);
    end if;
    
end proc;


# ----------------------------------------------------------------------- #
# monic_sqf_mod_rs_lazard                                                 #
#                                                                         #
# Computes the square-free factorization of a polynomial. The result is   #
# in the sense of lazard.                                                 #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as square_free_factorization_monic                               #
#                                                                         #
# ASSUMPTIONS                                                             #
#   - p is expanded                                                       #
#   - v is the greatest variable of R                                     #
#   - deg(p, v) > 0                                                       #
#   - p is a monic polynomial in v                                        #
# ----------------------------------------------------------------------- #
monic_sqf_mod_rs_lazard := proc(p::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, $)
    
    local numParam :: posint,
          result :: list([list([polynom, posint]), TRDrs]),
          newResult :: list([list([polynom, posint]), TRDrs]),
          item :: [list([polynom, posint]), TRDrs],
          sqf_list :: list([polynom, posint]),
          rs_i :: TRDrs,
          p_i :: polynom;
    
    numParam := nops(R['variables']) - 1;
    
    result := monic_sqf_mod_rs_kalkbrener(p, v, rs, R);
    
    to numParam-1 do
        newResult := [];
        for item in result do
            sqf_list, rs_i := op(item);
            p_i := expand(mul(map(x -> x[1]^x[2], sqf_list)));
            newResult := [op(newResult), op(monic_sqf_mod_rs_kalkbrener(p_i, v, rs_i, R))];
        end do;
        result := newResult;
    end do;
    
    return result;
    
end proc;


# ----------------------------------------------------------------------- #
# monic_sqf_mod_rs_kalkbrener                                             #
#                                                                         #
# Computes the square-free factorization of a polynomial. The result is   #
# in the sense of kalkbrener.                                             #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as square_free_factorization_monic                               #
#                                                                         #
# ASSUMPTIONS                                                             #
#   - p is expanded                                                       #
#   - v is the greatest variable of R                                     #
#   - deg(p, v) > 0                                                       #
#   - p is a monic polynomial in v                                        #
# ----------------------------------------------------------------------- #
monic_sqf_mod_rs_kalkbrener := proc(p::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, $)
    
    local result :: list([list([polynom, posint]), TRDrs]),
          rc :: TRDrc,
          h :: polynom,
          p_sqf :: polynom,
          p_sqf_disc :: polynom,
          h_disc :: polynom,
          lrc :: TRDlrc,
          rc_i :: TRDrc;
    
    rc := RC:-TRDregular_chain_rs(rs, R);
    h := RC:-TRDlist_mul_polys(RC:-TRDinequations_rs(rs, R), R);
    
    # Get the square-free part of p when viewed as a univariate polynomial 
    # in v
    p_sqf := RC:-TRDuniv_square_free_part(p, v, R);

    # Compute the discriminant of the square-free part of p w.r.t. v
    p_sqf_disc := expand(discrim(p_sqf, v));
    
    h_disc := p_sqf_disc*h;
    
    # Factorization when p_sqf_disc <> 0
    result := monic_sqf_mod_rc_h(p, v, rc, h_disc, R);
    
    # Factorization when p_sqf_disc = 0
    lrc := RC:-Intersect(p_sqf_disc, rc, R);
    for rc_i in lrc do
        result := [op(result), op(monic_sqf_mod_rc_h(p, v, rc_i, h, R))];
    end do;
    
    return result;

end proc;


# ----------------------------------------------------------------------- #
# monic_sqf_mod_rc_h                                                      #
#                                                                         #
# Computes the square-free factorization over a regular chain in with an  #
# inequation constraint.                                                  #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   v .... Variable                                                       #
#   rc ... Regular chain                                                  #
#   h .... Inequation                                                     #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Same as square_free_factorization_monic                               #
#                                                                         #
# ASSUMPTIONS                                                             #
#   - p is expanded                                                       #
#   - v is the greatest variable of R                                     #
#   - deg(p, v) > 0                                                       #
#   - p is a monic polynomial in v                                        #
# ----------------------------------------------------------------------- #
monic_sqf_mod_rc_h := proc(p, v, rc, h, R)
    
    local result :: list([list([polynom, posint]), TRDrs]),
          sqf :: list([list([polynom, posint]), TRDrc]),
          item :: [list([polynom, posint]), TRDrc],
          lp :: list([polynom, posint]),
          rcnew :: TRDrc,
          reg :: TRDlrc,
          rc_reg :: TRDrc,
          rs :: TRDrs;
    
    result := [];
    
    sqf := RC:-TRDsqf_mod_rc(p, v, rc, R);
    
    for item in sqf do
        lp, rcnew := op(item);
        
        # h must be regular w.r.t. rcnew
        if not RC:-TRDis_regular(h, rcnew, R) then
            reg := RC_CT:-Regularize(h, rcnew, R)[1];
            for rc_reg in reg do
                rs := RC_CST:-RegularSystem(rc_reg, [h], R);
                result := [op(result), [lp, rs]];
            end do;
        else
            rs := RC_CST:-RegularSystem(rcnew, [h], R);
            result := [op(result), [lp, rs]];
        end if;
    end do;
    
    return result;
    
end proc;

end module;