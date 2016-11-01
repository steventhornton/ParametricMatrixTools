test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4,
          test5;

    uses ParametricMatrixTools,
         RegularChains,
         RegularChains:-ConstructibleSetTools;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1', 'test2', 'test3', 'test4', 'test5'];
        
        printf("Testing ListIntersection\n");
        
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

        local R, p, q, cs1, cs2, cs, cs3, result;
        
        # Example from the Intersection help page
        R := PolynomialRing([x, y, t]);
        p := (5*t+5)*x-y-10*t-7;
        q := (5*t-5)*x-(t+2)*y-7*t+11;
        cs1 := GeneralConstruct([p, q], [x-t], R);
        cs2 := GeneralConstruct([p, q], [x+t], R);
        
        cs3 := GeneralConstruct([p, q], [x+t, x-t], R);
        
        try
            cs := ParametricMatrixTools:-ListIntersection([cs1, cs2], R);
            result := IsContained(cs3, cs, R) and IsContained(cs, cs3, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test2 := proc($)

        local R, p, q, cs1, cs2, cs, cs3, result;
        
        # Example from the Intersection help page
        R := PolynomialRing([x, y, t]);
        p := (5*t+5)*x-y-10*t-7;
        q := (5*t-5)*x-(t+2)*y-7*t+11;
        cs1 := GeneralConstruct([p, q], [x-t], R);
        cs2 := GeneralConstruct([p, q], [x+t], R);
        cs3 := GeneralConstruct([p, q], [x+t, x-t], R);
        
        try
            cs := ParametricMatrixTools:-ListIntersection({cs1, cs2}, R);
            result := IsContained(cs3, cs, R) and IsContained(cs, cs3, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test3 := proc($)

        local R, p, q, cs1, cs2, cs, cs3, result;
        
        # Example from the Intersection help page
        R := PolynomialRing([x, y, t]);
        p := (5*t+5)*x-y-10*t-7;
        q := (5*t-5)*x-(t+2)*y-7*t+11;
        cs1 := GeneralConstruct([p, q], [x-t], R);
        cs2 := GeneralConstruct([p, q], [x+t], R);
        cs3 := GeneralConstruct([p, q], [x+t, x-t], R);
        
        try
            cs := ParametricMatrixTools:-ListIntersection(Array([cs1, cs2]), R);
            result := IsContained(cs3, cs, R) and IsContained(cs, cs3, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;


    test4 := proc($)
    
        local R, cs1, cs2, cs3, cs4, cs, result;
        
        R := PolynomialRing([a,b,c]);
        
        cs1 := GeneralConstruct([], [], R);
        cs2 := GeneralConstruct([a-1], [a-1], R);
        cs3 := GeneralConstruct([a+1], [b-2], R);
        cs4 := GeneralConstruct([a*c-b*a+2], [], R);
        
        try
            cs := ParametricMatrixTools:-ListIntersection([cs1, cs2, cs3, cs4], R);
            result := IsEmpty(cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test5 := proc($)
    
        local R, cs1, cs2, cs3, cs4, cs5, cs, result;
        
        R := PolynomialRing([a,b,c]);
        
        cs1 := GeneralConstruct([], [], R);
        cs2 := GeneralConstruct([a-1], [], R);
        cs3 := GeneralConstruct([a-2], [b-12], R);
        cs4 := GeneralConstruct([c+3, b], [b-2, a], R);
        
        cs5 := GeneralConstruct([a-1, a-2, c+3, b], [b-12, b-2, a], R);
        
        try
            cs := ParametricMatrixTools:-ListIntersection([cs1, cs2, cs3, cs4], R);
            result := IsContained(cs5, cs, R) and IsContained(cs, cs5, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;

end module: