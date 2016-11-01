test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4,
          test5;

    uses ParametricMatrixTools, 
         RegularChains;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1','test2', 'test3', 'test4', 'test5'];
        
        printf("Testing InitialByVar\n");
        
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

        local R, p, result, correct;

        R := PolynomialRing([x, a, b]);

        p := x^2 + 2*x + 1;

        try
            result := ParametricMatrixTools:-InitialByVar(p, x, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := 1;
    
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

        local R, p, result, correct;

        R := PolynomialRing([x, a, b]);

        p := 10;

        try
            result := ParametricMatrixTools:-InitialByVar(p, x, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := 10;

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

        local R, p, result, correct;

        R := PolynomialRing([x, a, b]);

        p := 10*a^2 + 13*b^3 + a*b - 4;

        try
            result := ParametricMatrixTools:-InitialByVar(p, x, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;

        correct := 10*a^2 + 13*b^3 + a*b - 4;

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

        local R, p, result, correct;

        R := PolynomialRing([x, a, b]);

        p := (a^2 + b)*x^3 + (a*b + a^2 - b^2 -1)*x - 20;

        try
            result := ParametricMatrixTools:-InitialByVar(p, x, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;

        correct := a^2 + b;

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
    
        local R, p, result, correct;
    
        R := PolynomialRing([x, a, b]);
    
        p := 0;
    
        try
            result := ParametricMatrixTools:-InitialByVar(p, x, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := 0;
    
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