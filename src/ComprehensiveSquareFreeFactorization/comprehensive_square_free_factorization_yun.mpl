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
#   A list of lists of the form                                           #
#       [lp_i, cs_i]                                                      #
#   Where                                                                 #
#       - cs_i is a constructible set                                     #
#       - lp_i is a list of the form:                                     #
#             [m_i, llp_i]                                                #
#         where llp_i is a list of the form                               #
#             [p_j, n_j]                                                  #
#         such that p = m_i*product(p_j^n_j) and p_j are the square-free  #
#         factors in the zero set of cs_i. Each m_i is a rational         #
#         function in the parameters.                                     #
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
          cleanCofactors,
          computeMultipliers,
          cleanResult;

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
          cs_zero :: TRDcs,
          result, a, b, c, rs, task;

    dp := diff(p, v);

    gcdResult, cs_zero := ComprehensiveGcd(p, dp, v, cs, R, 'outputType' = 'RS', 'cofactors'=true);

    ASSERT(RC:-TRDis_empty_constructible_set(cs_zero, R));

    result := [];

    for task in gcdResult do
        a, b, c, rs := op(task);
        b, c := cleanCofactors(b, c);
        result := [op(result), op(computeFactors(a, b, c, v, rs, R))];
    end do;
    
    result := cleanResult(result, v, R);
    result := computeMultipliers(p, result, v);
    return result;

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
#                  [lp, cs]                                               #
#              where lp is a list with elements of the form               #
#                  [p_i, n_i]                                             #
#   v ........ Variable                                                   #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [lp, cs]                                                          #
#   where lp is a list with elements of the form                          #
#       [m, llp]                                                          #
#   where llp is a list with elements of the form                         #
#       [p_i, n_i]    (same as the input lp)                              #
#   such that p = m*product(p_i^n_i) where m is a rational function.      #
# ----------------------------------------------------------------------- #
computeMultipliers := proc(p::polynom, result, v::name, $)
    
    local out, task, a, cs, q, m;
    
    out := [];
    
    for task in result do
        a, cs := op(task);
        q := mul(map(x -> x[1]^x[2], a));
        m := normal(lcoeff(p, v)/lcoeff(q, v));
        out := [op(out), [[m, a], cs]];
    end do;
    
    return out;
    
end proc;


# ----------------------------------------------------------------------- #
# cleanCofactors                                                          #
#                                                                         #
# Convert the input rational functions to polynomials by multiplying both #
# by the lcm of their denominators.                                       #
#                                                                         #
# INPUT                                                                   #
#   b ... Rational function                                               #
#   c ... Rational function                                               #
#                                                                         #
# OUTPUT                                                                  #
#   Two polynomials:                                                      #
#       b_out = b*lcm(denom(b),denom(c))                                  #
#       c_out = c*lcm(denom(b),denom(c))                                  #
# ----------------------------------------------------------------------- #
cleanCofactors := proc(b::ratpoly, c::ratpoly, $)
    
    local d_b :: polynom,
          d_c :: polynom,
          l :: polynom;
    
    # Compute the cofactors (b = p1/g)
    d_b := denom(b);
    
    # Compute the cofactors (c = p2/g)
    d_c := denom(c);
    
    l := lcm(d_c, d_b);
    
    return normal(b*l), normal(c*l);
    
end proc;


# ----------------------------------------------------------------------- #
# cleanResult                                                             #
#                                                                         #
# Remove the first element of each square-free factors list from the      #
# result and remove any factors that are equal to 1.                      #
# ----------------------------------------------------------------------- #
cleanResult := proc(result, v, R, $)

    local out, i, polyList, rs;

    out := [];

    for i to nops(result) do

        polyList, rs := op(result[i]);
        polyList := polyList[2..-1];
        polyList := remove(x -> isZeroOverRS(x[1]-1, rs, R), polyList);

        out := [op(out), [polyList, rs]]
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
          cs_zero :: TRDcs,
          rs :: TRDrs;

    out, task := [], [];

    if isConstant(b_in, R) then
        return [], [[a_in, rs_in]];
    end if;

    # pGcd is a list with elements of the form
    #   [g, c1, c2, rs]
    # where g = gcd(b_in, d_in), c1*g = b_in, c2*g = d_in
    gcdResult, cs_zero := ComprehensiveGcd(b_in, d_in, v, rs_in, R, 'outputType' = 'RS', 'cofactors'=true);

    ASSERT(RC:-TRDis_empty_constructible_set(cs_zero, R));

    for g in gcdResult do
        a, b, c, rs := op(g);
        b, c := cleanCofactors(b, c);
        
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