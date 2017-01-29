test := module()

    export ModuleApply;

    local verifyResult_rs,
          test1,
          test2,
          test3,
          test4,
          test5;

    uses ParametricMatrixTools, RegularChains, RegularChains:-ConstructibleSetTools, RegularChains:-ChainTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList, t; 

        testList := ['test1', 'test2', 'test3', 'test4', 'test5'];

        printf("Testing SquarefreeFactorization_monic\n");

        passCount, failCount := 0, 0;

        for test in testList do
            printf("\t%a: ...", test);
            try
                t := test();
            catch:
                printf("\b\b\b");
                printf("Error in test!\n");
                t := false;
            end try;
            if t then
                passCount := passCount + 1
            else
                failCount := failCount + 1;
            end if;
        end do;

        printf("\n");

        return passCount, failCount;

    end proc;
    
    verifyResult_rs := proc(p, v, rs, result, R)
        
        local cs, cs_result, item, d1, d2, m, sqr, es, g, cs_zero, G, ds, i, j, g_item;
        
        # Check that the union of the regular systems in result equal rs
        cs := ConstructibleSet([rs], R);
        cs_result := RegularChains:-TRDempty_constructible_set(R);
        for item in result do
            cs_result := Union(cs_result, ConstructibleSet([item[2]], R),R);
        end do;
        if not ParametricMatrixTools:-TRDequal_cs(cs, cs_result, R) then
            printf("Cases missing");
            return false;
        end if;
        
        # Check that the total degree of the result matches the input 
        # polynomial
        d1 := degree(p, v);
        for item in result do
            d2 := add(x, x in map(x -> degree(x[1], v)*x[2], item[1]));
            if d1 <> d2 then
                printf("Degrees don't match\n");
                return false;
            end if;
        end do;
        
        # Check that each pairwise factor has gcd=1 over the given
        # regular system
        for item in result do
            sqr, es := op(item);
            sqr := ListTools:-Flatten(map(x -> x[1], sqr));
            if nops(sqr) = 1 then
                next;
            end if;
            for i to nops(sqr)-1 do
                for j from i+1 to nops(sqr) do
                    g, cs_zero := ComprehensiveGcd(sqr[i], sqr[j], v, es, R);
                    if not RegularChains:-TRDis_empty_constructible_set(cs_zero, R) then
                        printf("Two of the factors can be zero.\n");
                        return false;
                    end if;
                    # Verify the gcd is 1
                    for g_item in g do
                        G, ds := op(g_item);
                        if not ParametricMatrixTools:-isZeroOverCS(G-1, ds, R) then
                            printf("Gcd is not equal to 1.\n");
                            return false;
                        end if;
                    end do;
                end do;
            end do;
        end do;
        
        return true;
        
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 1                                                              #
    #                                                                     #
    # p = (x+1)^2*(x+a)                                                   #
    # rs = {}                                                             #
    # opt_lazard = true                                                   #
    # ------------------------------------------------------------------- #
    test1 := proc($)

        local p, R, rs, result, correct;
        
        p := (x+1)^2*(x+a);
        
        R := PolynomialRing([x, a]);
        rs := RegularSystem(R);
        
        try
            result := SquarefreeFactorization_monic(p, 'x', rs, R, true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;

        correct := verifyResult_rs(p, 'x', rs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 2                                                              #
    #                                                                     #
    # p = (x+a)*(x+b)^2*(x+1)^3                                           #
    # rs = {}                                                             #
    # opt_lazard = true                                                   #
    # ------------------------------------------------------------------- #
    test2 := proc($)
    
        local p, R, rs, result, correct;
        
        p := (x+a)*(x+b)^2*(x+1)^3;
        
        R := PolynomialRing([x, a, b]);
        rs := RegularSystem(R);
        
        try
            result := SquarefreeFactorization_monic(p, 'x', rs, R, true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_rs(p, 'x', rs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 3                                                              #
    #                                                                     #
    # p = (x+a)*(a*b+x)^2                                                 #
    # rs = {}                                                             #
    # opt_lazard = true                                                   #
    # ------------------------------------------------------------------- #
    test3 := proc($)
    
        local p, R, rs, result, correct;
        
        p := (x+a)*(a*b+x)^2;
        
        R := PolynomialRing([x, a, b]);
        rs := RegularSystem(R);
        
        try
            result := SquarefreeFactorization_monic(p, 'x', rs, R, true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_rs(p, 'x', rs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 4                                                              #
    #                                                                     #
    # p = (x+1)^2*(x^2+x+1)*(x+a)                                         #
    # rs = {}                                                             #
    # opt_lazard = true                                                   #
    # ------------------------------------------------------------------- #
    test4 := proc($)
    
        local p, R, rs, result, correct;
        
        p := (x+1)^2*(x^2+x+1)*(x+a);
        
        R := PolynomialRing([x, a]);
        rs := RegularSystem(R);
        
        try
            result := SquarefreeFactorization_monic(p, 'x', rs, R, true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_rs(p, 'x', rs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 5                                                              #
    #                                                                     #
    # p = (x+a)*(x+b)^2*(x-a-c)^3                                         #
    # rs = {}                                                             #
    # opt_lazard = true                                                   #
    # ------------------------------------------------------------------- #
    test5 := proc($)
    
        local p, R, rs, result, correct;
        
        p := (x+a)*(x+b)^2*(x-a-c)^3;
        
        R := PolynomialRing([x, a, b, c]);
        rs := RegularSystem(R);
        
        try
            result := SquarefreeFactorization_monic(p, 'x', rs, R, true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_rs(p, 'x', rs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
end module: