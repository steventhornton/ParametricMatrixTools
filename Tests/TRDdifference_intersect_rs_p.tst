test := module()

    export ModuleApply;

    local test1,
          test2,
          test3;

    uses ParametricMatrixTools, RegularChains, RegularChains:-ConstructibleSetTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList; 

        testList := ['test1', 'test2', 'test3'];

        printf("Testing TRDdifference_intersect_rs_p\n");

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

        local p, R, cs, cs_nz, cs_z, correct, lrs, rs, lcs_nz, lcs_z;

        p := a - 1;

        R := PolynomialRing([a, b]);
        cs := GeneralConstruct([b-1], [], R);
        
        lrs := RepresentingRegularSystems(cs, R);

        try
            lcs_nz, lcs_z := [], [];
            for rs in lrs do
                cs_nz, cs_z := ParametricMatrixTools:-TRDdifference_intersect_rs_p(rs, p, R);
                lcs_nz := [op(lcs_nz), cs_nz];
                lcs_z := [op(lcs_z), cs_z];
            end do;
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;
        for cs_z in lcs_z do
            correct := correct and ParametricMatrixTools:-isZeroOverCS(p, cs_z, R);
        end do;
        for cs_nz in lcs_nz do
            correct := correct and ParametricMatrixTools:-isNonZeroOverCS(p, cs_nz, R);
        end do;

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

        local p, R, cs, cs_nz, cs_z, correct, lrs, rs, lcs_nz, lcs_z;

        p := a^2 - 1;

        R := PolynomialRing([a, b]);
        cs := GeneralConstruct([], [a*b - b^2 + a^2 - 10], R);

        lrs := RepresentingRegularSystems(cs, R);
        
        try
            lcs_nz, lcs_z := [], [];
            for rs in lrs do
                cs_nz, cs_z := ParametricMatrixTools:-TRDdifference_intersect_rs_p(rs, p, R);
                lcs_nz := [op(lcs_nz), cs_nz];
                lcs_z := [op(lcs_z), cs_z];
            end do;
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;
        for cs_z in lcs_z do
            correct := correct and ParametricMatrixTools:-isZeroOverCS(p, cs_z, R);
        end do;
        for cs_nz in lcs_nz do
            correct := correct and ParametricMatrixTools:-isNonZeroOverCS(p, cs_nz, R);
        end do;

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

        local p, R, cs, cs_nz, cs_z, correct, lrs, rs, lcs_nz, lcs_z;

        p := a*b + b^2 - a*b^2 + a^10;

        R := PolynomialRing([a, b]);
        cs := GeneralConstruct([], [a*b - b^2 + a^2 - 10], R);

        lrs := RepresentingRegularSystems(cs, R);
        
        try
            lcs_nz, lcs_z := [], [];
            for rs in lrs do
                cs_nz, cs_z := ParametricMatrixTools:-TRDdifference_intersect_rs_p(rs, p, R);
                lcs_nz := [op(lcs_nz), cs_nz];
                lcs_z := [op(lcs_z), cs_z];
            end do;
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := true;
        for cs_z in lcs_z do
            correct := correct and ParametricMatrixTools:-isZeroOverCS(p, cs_z, R);
        end do;
        for cs_nz in lcs_nz do
            correct := correct and ParametricMatrixTools:-isNonZeroOverCS(p, cs_nz, R);
        end do;
        
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