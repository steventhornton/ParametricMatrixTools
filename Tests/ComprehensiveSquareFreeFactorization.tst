test := module()

    export ModuleApply;

    local verifyResult_cs,
          test1,
          test2,
          test3,
          test4,
          test5,
          test6,
          test7,
          test8,
          test9,
          test10,
          test11,
          test12,
          test13,
          test14,
          test15,
          test16;

    uses ParametricMatrixTools, RegularChains, RegularChains:-ConstructibleSetTools, RegularChains:-ChainTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList, t; 

        testList := ['test1',  'test2',  'test3',  'test4',  'test5', 
                     'test6',  'test7',  'test8',  'test9',  'test10',
                     'test11', 'test12', 'test13', 'test14', 'test15',
                     'test16'];

        printf("Testing ComprehensiveSquareFreeFactorization\n");

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
    
    
    verifyResult_cs := proc(p::depends(polyInRing(R)), v::name, cs::TRDcs, result, R::TRDring, $)
        
        local csList :: TRDlcs,
              d1, d2, item, sqr, es, ds, i, j, g, G, g_item, cs_zero, q, m;
        
        # Check that the constructible sets in the result list form a 
        # partition of cs
        csList := map(x -> x[-1], result);
        if not ParametricMatrixTools:-TRDis_partition_cs(csList, cs, R) then
            printf("Not a partition!");
            return false;
        end if;
        
        # Check that the total degree of the result matches the input 
        # polynomial
        d1 := degree(p, v);
        for item in result do
            d2 := add(map(x -> degree(x[1], v)*x[2], item[2]));
            if d1 <> d2 then
                printf("Degrees don't match\n");
                return false;
            end if;
        end do;
        
        # Check that each pairwise factor has gcd=1 over the given
        # constructible set
        for item in result do
            m, sqr, es := op(item);
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
        
        # Check that the product of the square-free decomposition equals
        # the input polynomial over the given constructible set
        for item in result do
            m, sqr, es := op(item);
            q := m*mul(map(x -> x[1]^x[2], sqr));
            if not ParametricMatrixTools:-isZeroOverCS(numer(normal(p - q)), es, R) then
                printf("Square-free factorization not equal to input polynomial\n");
                return false;
            end if;
        end do;
        
        return true;
        
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 1-4                                                            #
    #                                                                     #
    # p = (x+1)^2 * (x+a)                                                 #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test1 := proc($)

        local p, R, cs, result, correct;
        
        p := (x+1)^2 * (x+a);
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', cs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;

        correct := verifyResult_cs(p, 'x', cs, result, R);

        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;
    
    
    test2 := proc($)
    
        local p, R, rs, cs, result, correct;
        
        p := (x+1)^2 * (x+a);
        
        R := PolynomialRing([x,a]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', rs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test3 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x+1)^2 * (x+a);
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test4 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x+1)^2 * (x+a);
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
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
    # TEST 5-8                                                            #
    #                                                                     #
    # p = (x-b)^4 + a*x^2                                                 #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test5 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x-b)^4+a*x^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', cs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test6 := proc($)
    
        local p, R, rs, cs, result, correct;
        
        p := (x-b)^4+a*x^2;
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', rs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test7 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x-b)^4+a*x^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test8 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x-b)^4+a*x^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
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
    # TEST 9-12                                                           #
    #                                                                     #
    # p = (x+1)^2*(x^2+x+1)*(x+a)                                         #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test9 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x+1)^2*(x^2+x+1)*(x+a);
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', cs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test10 := proc($)
    
        local p, R, rs, cs, result, correct;
        
        p := (x+1)^2*(x^2+x+1)*(x+a);
        
        R := PolynomialRing([x,a]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', rs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test11 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x+1)^2*(x^2+x+1)*(x+a);
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test12 := proc($)
    
        local p, R, cs, result, correct;
        
        p := (x+1)^2*(x^2+x+1)*(x+a);
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
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
    # TEST 13-16                                                          #
    #                                                                     #
    # p = x^3 + a*x^2 - x - a*b - a                                       #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test13 := proc($)
    
        local p, R, cs, result, correct;
        
        p := x^3 + a*x^2 - x - a*b - a;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', cs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test14 := proc($)
    
        local p, R, rs, cs, result, correct;
        
        p := x^3 + a*x^2 - x - a*b - a;
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', rs, R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test15 := proc($)
    
        local p, R, cs, result, correct;
        
        p := x^3 + a*x^2 - x - a*b - a;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test16 := proc($)
    
        local p, R, cs, result, correct;
        
        p := x^3 + a*x^2 - x - a*b - a;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            result := ComprehensiveSquareFreeFactorization(p, 'x', [], [], R, 'outputType'='CS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := verifyResult_cs(p, 'x', cs, result, R);
    
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