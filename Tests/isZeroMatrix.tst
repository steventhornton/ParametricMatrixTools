test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4,
          test5,
          test6;

    uses ParametricMatrixTools,
         RegularChains,
         RegularChains:-ChainTools,
         RegularChains:-ConstructibleSetTools;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1','test2', 'test3', 'test4', 'test5', 'test6'];
        
        printf("Testing isZeroMatrix\n");
        
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
        
        local A, result, correct;
        
        A := Matrix(2, 3);
        
        try
            result := ParametricMatrixTools:-isZeroMatrix(A);
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
        
        local A, R, cs, result, correct;
        
        R := PolynomialRing([a, b]);
        
        A := Matrix([[a+1, a-2], [b+2, b-2]]);
        
        cs := GeneralConstruct([a+1, a-2, b+2, b-2], [], R);
        
        try
            result := ParametricMatrixTools:-isZeroMatrix(A, cs, R);
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
        
        local A, R, cs, result, correct;
        
        R := PolynomialRing([a, b]);
        
        A := Matrix([[a+1, a-2], [b+2, b-2]]);
        
        cs := GeneralConstruct([a+1, b+2], [], R);
        
        try
            result := ParametricMatrixTools:-isZeroMatrix(A, cs, R);
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
    
    
    test4 := proc($)
        
        local A, R, rc, rs, result, correct;
        
        R := PolynomialRing([x, y]);
        
        A := Matrix(2, 4);
        A[2, 3] := x + y;
        
        rc := Chain([(x+y)^2], Empty(R), R); 
        rs := RegularSystem(rc, R);
        
        try
            result := ParametricMatrixTools:-isZeroMatrix(A, rs, R);
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
        
        local A, R, rc, rs, result, correct;
        
        R := PolynomialRing([x, y]);
        
        A := Matrix(2, 4);
        A[2, 3] := x + y;
        A[2, 4] := 1;
        
        rc := Chain([(x+y)^2], Empty(R), R); 
        rs := RegularSystem(rc, R);
        
        try
            result := ParametricMatrixTools:-isZeroMatrix(A, rs, R);
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
        
        local A, R, rc, rs, result, correct;
        
        R := PolynomialRing([x, y]);
        
        A := Matrix(2, 4);
        A[2, 3] := x + y;
        A[1, 1] := 1;
        
        rc := Chain([(x+y)^2], Empty(R), R); 
        rs := RegularSystem(rc, R);
        
        try
            result := ParametricMatrixTools:-isZeroMatrix(A, rs, R);
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
    
    
end module: