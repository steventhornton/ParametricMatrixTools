test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4,
          test5,
          test6,
          test7;

    uses ParametricMatrixTools, 
         RegularChains, 
         RegularChains:-ConstructibleSetTools;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1','test2', 'test3', 'test4', 'test5', 'test6', 'test7'];
        
        printf("Testing isUnder\n");
        
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
        
        local R, cs, result, correct;
        
        R := PolynomialRing([x, a, b]);
        
        cs := GeneralConstruct([a^2 + a*b - 2], [b-1], R);
        
        try
            result := ParametricMatrixTools:-isUnder(cs, x, R);
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
        
        local R, cs, result, correct;
        
        R := PolynomialRing([x, a, b]);
        
        cs := GeneralConstruct([], [], R);
        
        try
            result := ParametricMatrixTools:-isUnder(cs, x, R);
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
        
        local R, cs, result, correct;
        
        R := PolynomialRing([x, a, b]);
        
        cs := GeneralConstruct([a-1], [a-1], R);
        
        try
            result := ParametricMatrixTools:-isUnder(cs, x, R);
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
        
        local R, cs, result, correct;
        
        R := PolynomialRing([x, a, b]);
        
        cs := GeneralConstruct([a*x+b^2], [a, b], R);
        
        try
            result := ParametricMatrixTools:-isUnder(cs, x, R);
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
    
    
    test5 := proc($)
        
        local R, cs, result, correct;
        
        R := PolynomialRing([x, a, b]);
        
        cs := GeneralConstruct([(x+1)*(x+2)*(x+a)*(x+b)], [], R);
        
        try
            result := ParametricMatrixTools:-isUnder(cs, x, R);
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
        
        local R, rc, rs, result, correct;
        
        R := PolynomialRing([x, a, b, c]);
        
        rc := Triangularize([a*x^2 + b*x + c], R)[1];
        
        rs := RegularSystem(rc, R);
        
        try
            result := ParametricMatrixTools:-isUnder(rs, x, R);
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
    
    
    test7 := proc($)
        
        local R, rc, rs, result, correct;
        
        R := PolynomialRing([x, a, b, c]);
        
        rc := Triangularize([a^2 + b^2 - 1], R)[1];
        
        rs := RegularSystem(rc, R);
        
        try
            result := ParametricMatrixTools:-isUnder(rs, x, R);
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