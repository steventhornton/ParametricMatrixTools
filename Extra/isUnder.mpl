# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isUnder.mpl                                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 24/2016                                                #
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