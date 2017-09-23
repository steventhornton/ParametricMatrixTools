# This method is for printing the results of a call to the 
# ComprehensiveSquareFreeFactorization method with rs as the output type.
printResult_CSFF := proc(result, p, cs, R, k)
    
    uses RegularChains, RegularChains:-ConstructibleSetTools;
    
    local i, rc, eqn, rs;
    
    printf("============================================================\n");
    printf("ComprehensiveSquareFreeFactorization - Example %d\n", k);
    printf("============================================================\n\n");
    
    printf("Computing the Square-free Decomposition of\n");
    printf("\tp = %a\n\n", p);
    
    if not IsEmpty(Complement(cs, R),R) then
        printf("\tParameter constraints:\n");
        for rs in RepresentingRegularSystems(cs, R) do
            print_rs(rs, R,1);
        end do;
    end if;
    
    for i to nops(result) do
        printf("Case %d:\n", i);
        printf("\tm = %a\n", result[i][1]);
        printf("\tSquare-free Decomposition: %a\n", result[i][2]);
        printf("\tParameter constraints:\n");
        print_rs(result[i][3],R,2);
    end do;
    
end proc:


printResult_CGCD := proc(result, cs_zero, p1, p2, cs, R, k)
    
    uses RegularChains, RegularChains:-ConstructibleSetTools;
    
    local i, rc, eqn, rs;
    
    printf("============================================================\n");
    printf("ComprehensiveGcd - Example %d\n", k);
    printf("============================================================\n\n");
    
    printf("Computing the gcd of\n");
    printf("\tp1 = %a\n", p1);
    printf("\tp2 = %a\n\n", p2);
    
    if not IsEmpty(Complement(cs, R),R) then
        printf("Parameter constraints:\n");
        for rs in RepresentingRegularSystems(cs, R) do
            print_rs(rs, R,1);
        end do;
        printf("\n");
    end if;
    
    if not IsEmpty(cs_zero, R) then
        printf("Gcd does not exist when:\n");
        for rs in RepresentingRegularSystems(cs_zero, R) do
            print_rs(rs, R,1);
        end do;
        printf("\n");
    end if;
    
    for i to nops(result) do
        printf("Case %d:\n", i);
        printf("\tGcd = %a\n", result[i][1]);
        printf("\tParameter constraints:\n");
        print_rs(result[i][2],R,2);
        printf("\n");
    end do;
    
end proc:


printResult_SFFm := proc(result, p, rs, R, opt_lazard, k)

    uses RegularChains, RegularChains:-ConstructibleSetTools;

    local i, rc, eqn;

    printf("============================================================\n");
    printf("SquarefreeFactorization_monic - Example %d\n", k);
    printf("============================================================\n\n");
    
    printf("Computing the Square-free factorization in the sense of ");
    if opt_lazard then
        printf("Lazard");
    else
        printf("Kalkbrener");
    end if;
    printf(" of \n");
    printf("\tp = %a\n\n", p);
    
    if not IsEmpty(Complement(ConstructibleSet([rs],R),R),R) then
        printf("\tParameter constraints:\n");
        print_rs(rs, R,1);
    end if;

    for i to nops(result) do
        printf("Case %d:\n", i);
        printf("\tSquare-free Decomposition: %a\n", map(factor, result[i][1]));
        printf("\tParameter constraints:\n");
        print_rs(result[i][2],R,2);
    end do;

end proc:


printResult_CR := proc(result, A, cs, R, k)

    uses RegularChains, RegularChains:-ConstructibleSetTools;

    local i, rc, eqn;

    printf("============================================================\n");
    printf("ComprehensiveRank - Example %d\n", k);
    printf("============================================================\n\n");
    
    printf("\tA = %a\n\n", A);

    if not IsEmpty(Complement(cs,R),R) then
        printf("\tParameter constraints:\n");
        print_rs(cs, R,1);
    end if;

    for i to nops(result) do
        printf("Case %d:\n", i);
        print(result[i][1]);
        printf("\tRank: %d\n", result[i][1]);
        printf("\tParameter constraints:\n");
        print_rs(result[i][2], R, 2);
    end do;

end proc:


print_rs := proc(rs, R, n)
    
    local rc, eqn;
    
    rc := RepresentingChain(rs, R);
    for eqn in Equations(rc, R) do
        printf(cat(seq("\t", i=1..n), "{ %a = 0\n"), eqn);
    end do;
    for eqn in Inequations(rc, R) do
        printf(cat(seq("\t", i=1..n), "{ %a ≠ 0\n"), eqn);
    end do;
    for eqn in RepresentingInequations(rs, R) do
        printf(cat(seq("\t", i=1..n), "{ %a ≠ 0\n"), eqn);
    end do;
end proc: