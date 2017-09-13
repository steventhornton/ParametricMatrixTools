# ======================================================================= #
# ======================================================================= #
#                                                                         #
# square_free_factorization_monic.mpl                                     #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Mar. 9/2017                                                 #
#                                                                         #
# Compute the square-free factorization of a parametric, univariate       #
# polynomial that is monic in its main variable. A complete case          #
# discussion forming a partition of the input reular system is returned.  #
# The partition of the input regular system is such that over each branch #
# the square-free factorization of the input polynomial is continuous.    #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   v .... Variable                                                       #
#   rs ... Regular system                                                 #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list of with elements of the form                                   #
#       [lp_i, rs_i]                                                      #
#   where lp_i is a list with elements of the form                        #
#       [p_j, n_j]                                                        #
#   and rs_i is a regular system.                                         #
#   There exists a rational function m_i such that                        #
#       p = m_i*product(p_j^n_j)                                          #
#   where p_j are the square-free factors of p such that the              #
#   factorization is continuous over rs_i.                                #
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
          sqf_mod_rs,
          getNewTasks,
          getDiscrims,
          getResultants;

    ModuleApply := proc(p::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, $)
        return implementation(p, v, rs, R);
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
implementation := proc(p_in::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, $)
    
    local p :: polynom,
          rc :: TRDrc,
          h :: list(polynom);
    
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
    
    rc := RC_CST:-RepresentingChain(rs, R);
    h := RC_CST:-RepresentingInequations(rs, R);
    
    return sqf_mod_rs(p, v, rc, h, R);
    
end proc;


# ----------------------------------------------------------------------- #
# sqf_mod_rs                                                              #
#                                                                         #
# Compute the square-free factorization of a polynomial given a regular   #
# chain and a list of inequations that form a regular system.             #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   v .... Variable                                                       #
#   rc ... Regular chain                                                  #
#   H .... List of polynomials representing inequations                   #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Same as square_free_factorization_monic                               #
# ----------------------------------------------------------------------- #
sqf_mod_rs := proc(p::depends(polyInRing(R)), v::name, rc::TRDrc, H::list(polynom), R::TRDring, $)

    local out, tasks, new_tasks, task, tmp_tasks, new_out;
    
    # Get the first task
    tasks, out := getNewTasks(p, v, rc, H, R);
    
    while nops(tasks) > 0 do
        new_tasks := [];
        for task in tasks do
            tmp_tasks, new_out := getNewTasks(p, v, task, H, R); 
            new_tasks := [op(new_tasks), op(tmp_tasks)];
            out := `union`(out, new_out);
        end do;
        tasks := new_tasks;
        tasks := LT:-MakeUnique(tasks, 1, proc (x, y) options operator, arrow; RC_CT:-EqualSaturatedIdeals(x, y, R) end proc);
    end do;
    return convert(out, list);
end proc;


# ----------------------------------------------------------------------- #
# getNewTasks                                                             #
#                                                                         #
# Compute one square-free factorization.                                  #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   v .... Variable                                                       #
#   rc ... Regular chain                                                  #
#   H .... List of polynomials representing inequations                   #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Two values:                                                           #
#       tasks, out                                                        #
#   - tasks is a list of regular chains that the square-free              #
#     factorization of p is not continuous over.                          #
#   - out is a set with elements of the form                              #
#          [lp_i, rs_i]                                                   #
#     where lp_i is a list with elements of the form                      #
#       [p_j, m_j]                                                        #
#    that gives the square-free factorization of p. The square-free       #
#    factorization of p is continous over rs_i.                           #
# ----------------------------------------------------------------------- #
getNewTasks := proc(p::depends(polyInRing(R)), v::name, rc::TRDrc, H::list(polynom), R::TRDring, $)
    
    local out, tasks, sqf, item, lp_s, rc_s, Hs, h, rs, rc_h, rc_h_i, inRad, h_i;
    
    out := {};
    tasks := {};

    sqf := RC:-TRDsqf_mod_rc(p, v, rc, R);
    
    for item in sqf do
        lp_s, rc_s := op(item);
        Hs := `union`(getDiscrims(lp_s, v), getResultants(lp_s, v));
        Hs := `minus`(Hs, {1});
        
        rs := RC_CST:-RegularSystem(rc_s, `union`(Hs, convert(H,set)), R);
        
        out := `union`(out, {[lp_s, rs]});
        
        for h in Hs do
            rc_h := RC:-Intersect(h, rc, R);
            for rc_h_i in rc_h do
                inRad := false;
                for h_i in H do
                    inRad := RC_CT:-IsInRadical(h_i, rc_h_i, R);
                    if inRad then break end if;
                end do;
                if not inRad then
                    tasks := `union`(tasks, {rc_h_i});
                end if;
            end do;
        end do;
    end do;
    tasks := convert(tasks, list);
    return tasks, out;
end proc;


# ----------------------------------------------------------------------- #
# getNewTasks                                                             #
#                                                                         #
# Computes the discriminant of all polynomials in a list                  #
#                                                                         #
# INPUT                                                                   #
#   lp ... List of pairs of the form                                      #
#            [p, i]                                                       #
#          where p is a polynomial and i is an integer.                   #
#   v .... Variable                                                       #
#                                                                         #
# OUTPUT                                                                  #
#   A set containing the discriminants of all the polynomials p.          #
# ----------------------------------------------------------------------- #
getDiscrims := proc(lp::list([polynom, posint]), v::name, $)
    return convert(map(x -> discrim(x[1], v), lp), set);
end proc;


# ----------------------------------------------------------------------- #
# getNewTasks                                                             #
#                                                                         #
# Computes all pairwise resultants of the polynomials in a list           #
#                                                                         #
# INPUT                                                                   #
#   lp ... List of pairs of the form                                      #
#            [p, i]                                                       #
#          where p is a polynomial and i is an integer.                   #
#   v .... Variable                                                       #
#                                                                         #
# OUTPUT                                                                  #
#   The set of all pairwise resultants the polynomials in the input list. #
# ----------------------------------------------------------------------- #
getResultants := proc(lp_in::list([polynom, posint]), v::name, $)
    local lp, i, j, out;
    out := {};
    lp := map(x -> x[1], lp_in);
    for i to nops(lp_in)-1 do
        for j from i+1 to nops(lp_in) do
            out := `union`(out, {resultant(lp[i], lp[j], v)});
        end do;
    end do;
    return out;
end proc;

end module;