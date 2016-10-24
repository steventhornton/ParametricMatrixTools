test := module()
    
    export ModuleApply;
    
    local test1,
          test2,
          test3;
    
    uses RegularChains, 
         RegularChains:-ConstructibleSetTools,
         RegularChains:-ChainTools,
         ParametricMatrixTools;
    
    ModuleApply := proc($)
        printf("Testing isZeroOverRS\n");
        printf("\tTest 1: ...");
        test1();
        printf("\tTest 2: ...");
        test2();
        printf("\tTest 3: ...");
        test3();
        printf("\n");
    end proc;
    
    
    test1 := proc($)
        
        local R, rc, rs, p, result;
        
        R := PolynomialRing([x]);
        
        rc := Empty(R);
        rc := Chain([(x+1)^2], rc, R);
        
        rs := RegularSystem(rc, R);
        
        p := x+1;
        
        result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
        
    end proc;
    
    
    test2 := proc($)
        
        local R, cs, rs, p, result;
        
        R := PolynomialRing([x, y]);
        
        cs := GeneralConstruct([x*y+1, x^2-1], [x], R);
        
        rs := RepresentingRegularSystems(cs, R)[1];
        
        p := y+1;
        
        result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
        
    end proc;
    
    
    test3 := proc($)
        
        local R, cs, rs, p, result;
        
        R := PolynomialRing([x, y]);
        
        cs := GeneralConstruct([x*y+1, x^2-1], [x], R);
        
        rs := RepresentingRegularSystems(cs, R)[1];
        
        p := x;
        
        result := ParametricMatrixTools:-isZeroOverRS(p, rs, R);
        
        printf("\b\b\b");
        if not result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
        
    end proc;
    
end module: