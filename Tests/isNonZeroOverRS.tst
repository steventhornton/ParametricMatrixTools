test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4,
          test5,
          test6,
          test7,
          test8,
          test9;

    uses ParametricMatrixTools, 
         RegularChains, 
         RegularChains:-ConstructibleSetTools,
         RegularChains:-ChainTools;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1','test2', 'test3', 'test4', 'test5', 'test6', 'test7', 'test8', 'test9'];
        
        printf("Testing isNonZeroOverRS\n");
        
        passCount, failCount := 0, 0;
        
        for test in testList do
            printf("\t%a: ...", test);
            if test() then
                passCount := passCount + 1
            else
                failCount := failCount + 1;
            end if;
        end do;
        
        printf("\n");
        
        return passCount, failCount;
        
    end proc;


    test1 := proc($)

        local R, p, cs, lrs, rs, result, correct;

        R := PolynomialRing([x, y]);
        cs := GeneralConstruct([x+y], [], R);
        lrs := RepresentingRegularSystems(cs, R);
        rs := lrs[1];
        p := x+y+1; 
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;
        
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test2 := proc($)

        local R, p, cs, lrs, rs, result, correct;

        R := PolynomialRing([a, b, c]);
        cs := GeneralConstruct([a+b], [a], R);
        lrs := RepresentingRegularSystems(cs, R);
        rs := lrs[1];
        p := (a+b+1)*a^3; 
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test3 := proc($)

        local R, p, rc, rs, result, correct;

        R := PolynomialRing([x, y]);
        rc := Chain([(x+y)^2], Empty(R), R);
        rs := RegularSystem(rc, [y], R);
        p := y^2+x+y;
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test4 := proc($)

        local R, p, rc, rs, result, correct;

        R := PolynomialRing([x, y]);
        rc := Chain([(x+y)^2], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := x+y+1;
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test5 := proc($)

        local R, p, rs, result, correct;

        R := PolynomialRing([x, a]);
        rs := RegularSystem(Empty(R), [x+1], R); 
        p := (x+1)*(x^2-a+2*x);
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := false;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test6 := proc($)

        local R, p, rs, result, correct;

        R := PolynomialRing([x, a]);
        rs := RegularSystem(Empty(R), [(x+1)*(x^2-a+2*x)*(x^2+4)], R);
        p := (x+1)^2*(x^2-a+2*x);
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test7 := proc($)

        local R, p, rc, rs, result, correct;

        R := PolynomialRing([x, y, z]); 
        rc := Chain([(z+1)*(z+2), y^2+z, (x-z)*(x-y)], Empty(R), R); 
        rs := RegularSystem(rc, R); 
        p := z-x; 
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := false;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test8 := proc($)

        local R, p, rc, rs, result, correct;

        R := PolynomialRing([x, y, z]);
        rc := Chain([(z+1)*(z+2), y^2+z, (x-z)*(x-y)], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := (z+1)*(x^3+5);
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := false;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test9 := proc($)

        local R, p, rc, rs, result, correct;

        R := PolynomialRing([x, y, z]);
        rc := Chain([(z+1)*(z+2), y^2+z, (x-z)*(x-y)], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := x+y+z;
        
        try
            result := ParametricMatrixTools:-isNonZeroOverRS(p, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    end proc;


end module: