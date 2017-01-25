# ======================================================================= #
# ======================================================================= #
#                                                                         #
# comprehensive_square_free_factorization_yun.mpl                         #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Computes the square free factorization of a parametric univariate       #
# polynomial in the sense of Lazard. Computation is done over a           #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   v .... Variable                                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list of with elements of the form                                   #
#       [m_i, lp_i, rs_i]                                                 #
#   where                                                                 #
#       - rs_i is a regular system                                        #
#       - m_i is a rational function of the paramters                     #
#       - lp_i is a list with elements of the form:                       #
#             [p_j, n_j]                                                  #
#   such that p = m_i*product(p_j^n_j) and p_j are the square-free        #
#   factors in the zero set of rs_i.                                      #
#                                                                         #
# REFERENCES                                                              #
#   - Yun, D. Y. (1976, August). On square-free decomposition algorithms. #
#     In Proceedings of the third ACM symposium on Symbolic and algebraic #
#     computation (pp. 26-35). ACM.                                       #
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
comprehensive_square_free_factorization_yun := module()

    export ModuleApply;

    local implementation,
          computeFactors,
          oneYunIteration,
          computeMultipliers,
          csff_gcd,
          cleanResult1,
          cleanResult2;

    ModuleApply := proc(p::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, $)
        return implementation(p, v, cs, R);
    end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Parametric implementation of Yun's algorithm.                           #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as comprehensive_square_free_factorization_yun                  #
# ----------------------------------------------------------------------- #
implementation := proc(p::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, $)

    local dp :: polynom,
          gcdResult,
          result, a, b, c, rs, task;

    dp := diff(p, v);
    
    gcdResult := csff_gcd(p, dp, v, cs, R);

    result := [];

    for task in gcdResult do
        a, b, c, rs := op(task);
        result := [op(result), op(computeFactors(a, b, c, v, rs, R))];
    end do;
    
    result := cleanResult1(result, R);
    result := computeMultipliers(p, result, v);
    result := cleanResult2(result, v, R);
    return result;

end proc;


# ----------------------------------------------------------------------- #
# csff_gcd                                                                #
#                                                                         #
# Compute the gcd of g of p1 and p2 and compute the polynomials c1 and c2 #
# such that c1*g = p1, and c2*g = p2.                                     #
#                                                                         #
# INPUT                                                                   #
#   p1 ... Polynomial                                                     #
#   p2 ... Polynomial                                                     #
#   v .... Variable                                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, c1, c2, rs]                                                   #
#   where g = gcd(p1, p2) in the zero set of rs, and c1*g = p1, and       #
#   c2*g = p2.                                                            #
# ----------------------------------------------------------------------- #
csff_gcd := proc(p1, p2, v, cs, R)
    
    local result,
          cs_zero :: TRDcs,
          item :: [polynom, TRDrs],
          g :: polynom,
          rs :: TRDrs,
          rc :: TRDrc,
          out :: list([polynom, ratpoly, ratpoly, TRDrs]),
          c1 :: ratpoly,
          c2 :: ratpoly,
          l :: polynom,
          s :: integer;
    
    result, cs_zero := ComprehensiveGcd:-comprehensive_gcd_src(p1, p2, v, cs, R);
    
    ASSERT(RC:-TRDis_empty_constructible_set(cs_zero, R));
    
    result := ComprehensiveGcd:-convertToRS(result, R);
    
    # Compute the cofactors
    out := [];
    for item in result do
        
        g, rs := op(item);
        
        rc := RC_CST:-RepresentingChain(rs, R);
        
        c1 := ComprehensiveGcd:-pseudo_cofactor(p1, g, v, rc, R);
        c2 := ComprehensiveGcd:-pseudo_cofactor(p2, g, v, rc, R);
        
        # Clean the cofactors
        s := sign(lcoeff(g,v));
        g := s*g;
        
        c1 := s*c1;
        c2 := s*c2;
        
        l := lcm(denom(c1), denom(c2));
        c1 := normal(c1*l);
        c2 := normal(c2*l);
        
        out := [op(out), [g, c1, c2, rs]];
        
    end do;
    
end proc;


# ----------------------------------------------------------------------- #
# computeMultipliers                                                      #
#                                                                         #
# For each element of the result list, compute the multiplier m such that #
# the product of the square-free terms is equal to the input polynomial.  #
#                                                                         #
# INPUT                                                                   #
#   p ........ Polynomial                                                 #
#   result ... List with elements of the form                             #
#                  [lp, rs]                                               #
#              where lp is a list with elements of the form               #
#                  [p_i, n_i]                                             #
#   v ........ Variable                                                   #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [m_i, lp_i, rs_i]                                                 #
#   where lp is a list with elements of the form                          #
#       [p_j, n_j]    (same as the input lp)                              #
#   such that p = m_i*product(p_j^n_j) where m_i is a rational function   #
#   in the paramters (i.e. doesn't contain v).                            #
# ----------------------------------------------------------------------- #
computeMultipliers := proc(p::polynom, result, v::name, $)
    
    local out, task, a, rs, q, m;
    
    out := [];
    
    for task in result do
        a, rs := op(task);
        
        q := mul(map(x -> x[1]^x[2], a));
        m := normal(lcoeff(p, v)/lcoeff(q, v));
        
        out := [op(out), [m, a, rs]];
    end do;
    
    return out;
    
end proc;


# ----------------------------------------------------------------------- #
# cleanResult1                                                            #
#                                                                         #
# Remove the first element of each square-free factors list from the      #
# result and remove any factors that are equal to 1.                      #
# ----------------------------------------------------------------------- #
cleanResult1 := proc(result, R, $)

    local out, i, polyList, rs;

    out := [];

    # Clean up factors that = 1
    # Remove first factor
    for i to nops(result) do

        polyList, rs := op(result[i]);
        
        polyList := polyList[2..-1];
        polyList := remove(x -> isZeroOverRS(x[1]-1, rs, R), polyList);
        
        out := [op(out), [polyList, rs]]
    end do;
    
    return out;
    
end proc;


# ----------------------------------------------------------------------- #
# cleanResult2                                                            #
#                                                                         #
# - Remove the content from each factor (move into the multiplier)        #
# - Call SparsePseudoRemainder on each factor                             #
# - Remove the content again (move into multiplier)                       #
# - Call SparsePseudoRemainder on multiplier                              #
# ----------------------------------------------------------------------- #
cleanResult2 := proc(result, v, R, $)
    
    local out, item, lp, rs, rc, lp_tmp, task, p, n, m, m_n, m_d, co, h, h_d, h_n;
    
    out := [];
    
    for item in result do
        
        m, lp, rs := op(item);
        
        rc := RC_CST:-RepresentingChain(rs, R);
        
        # Remove content from each factor
        lp_tmp := [];
        for task in lp do
            p, n := op(task);
            p := primpart(p, v, 'co');
            m := m*co^n;
            lp_tmp := [op(lp_tmp), [p, n]];
        end do;
        lp := lp_tmp;
        
        # Call SparsePseudoRemainder on each factor
        lp_tmp := [];
        for task in lp do
            p, n := op(task);
            p := RC:-SparsePseudoRemainder(p, rc, R, 'h');
            m := m/(h^n);
            lp_tmp := [op(lp_tmp), [p, n]];
        end do;
        lp := lp_tmp;
        
        # Remove content from each factor (again)
        lp_tmp := [];
        for task in lp do
            p, n := op(task);
            p := primpart(p, v, 'co');
            m := m*co^n;
            lp_tmp := [op(lp_tmp), [p, n]];
        end do;
        lp := lp_tmp;
        
        # Call SparsePseudoRemainder on multiplier
        m_n := RC:-SparsePseudoRemainder(numer(m), rc, R, 'h_n');
        m_d := RC:-SparsePseudoRemainder(denom(m), rc, R, 'h_d');
        m := normal((m_n*h_d)/(m_d*h_n));
        
        out := [op(out), [m, lp, rs]];
        
    end do;
    
    return out;

end proc;


# ----------------------------------------------------------------------- #
# computeFactors                                                          #
#                                                                         #
# Yun's algorithm                                                         #
#                                                                         #
# INPUT/OUTPUT                                                            #
# ----------------------------------------------------------------------- #
computeFactors := proc(a_in::depends(polyInRing(R)), b_in::polynom, c_in::polynom, v::name, rs_in::TRDrs, R::TRDring, $)

    local a :: list([polynom, posint]),
          b :: polynom,
          d :: polynom,
          rs :: TRDrs,
          newTasks :: list([list([polynom, posint]), polynom, polynom, TRDrs, posint]),
          out :: list([list([polynom, posint]), TRDrs]),
          tasks :: list([list([polynom, posint]), polynom, polynom, TRDrs, posint]),
          result,
          i :: posint;

    # a is a list, b, d are polynomials
    d := expand(c_in - diff(b_in, v));
    tasks := [[[[a_in, 1]], b_in, d, rs_in, 1]];
    result := [];

    while nops(tasks) > 0 do

        a, b, d, rs, i := op(tasks[1]);

        newTasks, out := oneYunIteration(a, b, d, rs, v, R, i);

        tasks := [op(tasks[2..-1]), op(newTasks)];

        result := [op(result), op(out)];

    end do;

    return result;

end proc;


# ----------------------------------------------------------------------- #
# oneYunIteration                                                         #
#                                                                         #
# Run one iteration of Yun's algorithm.                                   #
#                                                                         #
# INPUT                                                                   #
#   a_in .... List of polynomials                                         #
#   b_in .... Polynomial                                                  #
#   d_in .... Polynomial                                                  #
#   rs_in ... Regular system                                              #
#   v ....... Variable                                                    #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   Two lists:                                                            #
#       out, task                                                         #
#   where task is a list with elements of the form                        #
#       [a, b, c, d, rs]                                                  #
#   where a is a list of polynomials, b is a polynomial, c is a           #
#   polynomial, d is a polynomial and rs is a regular system. out is a    #
#   list with elements of the form                                        #
#       [a, rs]                                                           #
#   where a is a list of polynomials and rs is a regular system.          #
# ----------------------------------------------------------------------- #
oneYunIteration := proc(a_in::depends(list([polyInRing(R), posint])), b_in::depends(polyInRing(R)), d_in::depends(polyInRing(R)), rs_in::TRDrs, v::name, R::TRDring, i::posint, $)

    local a :: polynom,
          b :: ratpoly,
          c :: ratpoly,
          d :: polynom,
          out :: list([list([polynom, posint]), TRDrs]),
          task :: list([list([polynom, posint]), polynom, polynom, TRDrs, posint]),
          g :: [polynom, ratpoly, ratpoly, rs],
          gcdResult :: list([polynom, ratpoly, ratpoly, TRDrs]),
          rs :: TRDrs;

    out, task := [], [];

    if isConstant(b_in, R) then
        return [], [[a_in, rs_in]];
    end if;

    # gcdResult is a list with elements of the form
    #   [g, c1, c2, rs]
    gcdResult := csff_gcd(b_in, d_in, v, RC_CST:-ConstructibleSet([rs_in], R), R);

    for g in gcdResult do
        a, b, c, rs := op(g);
        
        if not v in indets(b) then
            out := [op(out), [[op(a_in), [a, i]], rs]];
        else
            d := expand(c - (diff(b, v)));
            task := [op(task), [[op(a_in), [a, i]], b, d, rs, i+1]];
        end if;
    end do;

    return task, out;

end proc;

end module;