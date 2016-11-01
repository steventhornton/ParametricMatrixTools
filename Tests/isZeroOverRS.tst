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
    
    uses RegularChains, 
         RegularChains:-ConstructibleSetTools,
         RegularChains:-ChainTools,
         ParametricMatrixTools;
    
    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1', 'test2', 'test3', 'test4', 'test5', 'test6', 'test7', 'test8', 'test9'];
        
        printf("Testing isZeroOverRS\n");
        
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
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y]);
        rc := Chain([x+y], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := x+y+1;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
    
    
    test2 := proc($)
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([a, b, c]);
        rc := Chain([a+b], Empty(R), R);
        rs := RegularSystem(rc, [a], R);
        p := (a+b+1)*a^3;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
    
    
    test3 := proc($)
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y]);
        rc := Chain([(x+y)^2], Empty(R), R);
        rs := RegularSystem(rc, [y], R); 
        p := y^2+x+y;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y]);
        rc := Chain([(x+y)^2], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := x+y+1;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y]);
        rc := Chain([(x+y)^2], Empty(R), R);
        rs := RegularSystem(rc, [], R);
        p := x+y;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
    
    
    test6 := proc($)
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y]);
        rc := Chain([(x+y)^2-1], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := x+y;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y, z]);
        rc := Chain([(z+1)*(z+2), y^2+z, (x-z)*(x-y)], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := z-x;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y, z]);
        rc := Chain([(z+1)*(z+2), y^2+z, (x-z)*(x-y)], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := (z+1)*(x^3+5);
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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
        
        local R, rc, rs, p, result, correct;
        
        R := PolynomialRing([x, y, z]);
        rc := Chain([(z+1)*(z+2), y^2+z, (x-z)*(x-y)], Empty(R), R);
        rs := RegularSystem(rc, R);
        p := x+y+z;
        
        try
            result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
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