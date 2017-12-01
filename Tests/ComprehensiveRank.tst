test := module()

    export ModuleApply;

    local verifyResult_rs,
          verifyResult_cs,
          test1,
          test2,
          test3,
          test4;

    uses ParametricMatrixTools, RegularChains, RegularChains:-ConstructibleSetTools, RegularChains:-ChainTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList, t; 

        testList := ['test1',  'test2',  'test3',  'test4'];

        printf("Testing ComprehensiveRank\n");

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
    
    
    verifyResult_rs := proc(A, cs::TRDcs, result, R::TRDring, $)
        
        local rsList :: TRDlrs,
              csList :: TRDlcs, m, item;
        
        # Check that the regular systems in the result list form a 
        # partition of cs
        # rsList := map(x -> x[-1], result);
        # csList := map(x -> ConstructibleSet([x], R), rsList);
        # if not ParametricMatrixTools:-TRDis_partition_cs(csList, cs, R) then
        #     printf("Not a partition!");
        #     return false;
        # end if;
        
        # Check that the ranks are non-negative integers no greater than
        # the number of rows of A
        m := LinearAlgebra:-RowDimension(A);
        for item in result do
            if item[1] > m or item[1] < 0 then
                printf("Invalid rank!");
                return false;
            end if;
        end do;
        
        return true;
        
    end proc;
    
    
    verifyResult_cs := proc(A, cs::TRDcs, result, R::TRDring, $)
        
        local csList :: TRDlcs, m, item;
        
        # Check that the constructible sets in the result list form a 
        # partition of cs
        csList := map(x -> x[-1], result);
        if not ParametricMatrixTools:-TRDis_partition_cs(csList, cs, R) then
            printf("Not a partition!");
            return false;
        end if;
        
        # Check that the ranks are non-negative integers no greater than
        # the number of rows of A
        m := LinearAlgebra:-RowDimension(A);
        for item in result do
            if item[1] > m or item[1] < 0 then
                printf("Invalid rank!");
                return false;
            end if;
        end do;
        
        return true;
        
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 1                                                              #
    #                                                                     #
    # A = [1, 2, 3]                                                       #
    #     [4, 5, 6]                                                       #
    #     [7, 8, 9]                                                       #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test1 := proc($)
    
        local A, R, cs, rs, i, results, correct;
        
        A := Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
        
        R := PolynomialRing([a]);
        cs := GeneralConstruct([],[],R);
        rs := RepresentingRegularSystems(cs, R)[1];
        
        try
            results := [0$10];
            results[1]  := ComprehensiveRank(A, R,         'outputType'='CS');
            results[2]  := ComprehensiveRank(A, cs, R,     'outputType'='CS');
            results[3]  := ComprehensiveRank(A, rs, R,     'outputType'='CS');
            results[4]  := ComprehensiveRank(A, [], R,     'outputType'='CS');
            results[5]  := ComprehensiveRank(A, [], [], R, 'outputType'='CS');
            results[6]  := ComprehensiveRank(A, R,         'outputType'='RS');
            results[7]  := ComprehensiveRank(A, cs, R,     'outputType'='RS');
            results[8]  := ComprehensiveRank(A, rs, R,     'outputType'='RS');
            results[9]  := ComprehensiveRank(A, [], R,     'outputType'='RS');
            results[10] := ComprehensiveRank(A, [], [], R, 'outputType'='RS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        for i to 5 do
            correct := verifyResult_cs(A, cs, results[i], R);
            if not correct then
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        for i from 6 to 10 do
            correct := verifyResult_rs(A, cs, results[i], R);
            if not correct then
                printf("\b\b\b");
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        
        printf("\b\b\b");
        printf("Pass\n");
        return true;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 2                                                              #
    #                                                                     #
    # A = [1, 2, 3]                                                       #
    #     [4, 5, 6]                                                       #
    #     [7, 8, a]                                                       #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test2 := proc($)
    
        local A, R, cs, rs, i, results, correct;
        
        A := Matrix([[1, 2, 3], [4, 5, 6], [7, 8, a]]);
        
        R := PolynomialRing([a]);
        cs := GeneralConstruct([],[],R);
        rs := RepresentingRegularSystems(cs, R)[1];
        
        try
            results := [0$10];
            results[1]  := ComprehensiveRank(A, R,         'outputType'='CS');
            results[2]  := ComprehensiveRank(A, cs, R,     'outputType'='CS');
            results[3]  := ComprehensiveRank(A, rs, R,     'outputType'='CS');
            results[4]  := ComprehensiveRank(A, [], R,     'outputType'='CS');
            results[5]  := ComprehensiveRank(A, [], [], R, 'outputType'='CS');
            results[6]  := ComprehensiveRank(A, R,         'outputType'='RS');
            results[7]  := ComprehensiveRank(A, cs, R,     'outputType'='RS');
            results[8]  := ComprehensiveRank(A, rs, R,     'outputType'='RS');
            results[9]  := ComprehensiveRank(A, [], R,     'outputType'='RS');
            results[10] := ComprehensiveRank(A, [], [], R, 'outputType'='RS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        for i to 5 do
            correct := verifyResult_cs(A, cs, results[i], R);
            if not correct then
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        for i from 6 to 10 do
            correct := verifyResult_rs(A, cs, results[i], R);
            if not correct then
                printf("\b\b\b");
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        
        printf("\b\b\b");
        printf("Pass\n");
        return true;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 3                                                              #
    #                                                                     #
    # z3 = [0, 0, 0]                                                      #
    #      [0, 0, 0]                                                      #
    #      [0, 0, 0]                                                      #
    # z1 = [0]                                                            #
    #      [0]                                                            #
    #      [0]                                                            #
    # E = [1, 3, 1]                                                       #
    #     [3, 1, 1]                                                       #
    #     [0, 0, 0]                                                       #
    # A1 = [1, 1, 3]                                                      #
    #      [1, 3, 1]                                                      #
    #      [0, 0, 0]                                                      #
    # A2 = [     lambda,  3*lambda,      lambda]                          #
    #      [3*lambda+mu, lambda+mu, lambda+3*mu]                          #
    #      [          0,         0,           0]                          #
    # B = [0]                                                             #
    #     [0]                                                             #
    #     [1]                                                             #
    # A = [ -E, z3, z3, z3, B, z1, z1, z1, z1, z1]                        #
    #     [-A1, -E, z3, z3, z1, B, z1, z1, z1, z1]                        #
    #     [ A2, A1, -E, z3, z1, z1, B, z1, z1, z1]                        #
    #     [ z3, A2, -A1, -E, z1, z1, z1, B, z1, z1]                       #
    #     [ z3, z3, A2, -A1, z1, z1, z1, z1, B, z1]                       #
    #     [ z3, z3, z3, A2, z1, z1, z1, z1, z1, B]                        #
    # cs1 = {lambda = 0}                                                  #
    # cs2 = {lambda <> 0}                                                 #
    # ------------------------------------------------------------------- #
    test3 := proc($)
    
        local A, z3, z1, A1, A2, B, E, R, i, cs1, cs2, rs1, rs2, results, correct;
        
        z3 := LinearAlgebra:-ZeroMatrix(3);
        z1 := LinearAlgebra:-ZeroMatrix(3, 1);
        E := Matrix([[1, 3, 1],
                     [3, 1, 1],
                     [0, 0, 0]]);
        A1 := Matrix([[1, 1, 3],
                      [1, 3, 1],
                      [0, 0, 0]]);
        A2 := Matrix([[lambda, 3*lambda, lambda],
                      [3*lambda+mu, lambda+mu, lambda+3*mu],
                      [0, 0, 0]]);
        B := Matrix([[0],
                     [0],
                     [1]]);
        A := Matrix([[-E, z3, z3, z3, B, z1, z1, z1, z1, z1],
                     [-A1, -E, z3, z3, z1, B, z1, z1, z1, z1],
                     [A2, A1, -E, z3, z1, z1, B, z1, z1, z1],
                     [z3, A2, -A1, -E, z1, z1, z1, B, z1, z1],
                     [z3, z3, A2, -A1, z1, z1, z1, z1, B, z1],
                     [z3, z3, z3, A2, z1, z1, z1, z1, z1, B]]);
        
        R := PolynomialRing([mu, lambda]);
        
        cs1 := GeneralConstruct([lambda],[],R);
        cs2 := GeneralConstruct([],[lambda],R);
        rs1 := RepresentingRegularSystems(cs1, R)[1];
        rs2 := RepresentingRegularSystems(cs2, R)[1];
        
        try
            results := [0$15];
            
            results[1]  := ComprehensiveRank(A, cs1, R,          'outputType'='CS');
            results[2]  := ComprehensiveRank(A, rs1, R,          'outputType'='CS');
            results[3]  := ComprehensiveRank(A, [lambda], R,     'outputType'='CS');
            results[4]  := ComprehensiveRank(A, [lambda], [], R, 'outputType'='CS');
            
            results[5]  := ComprehensiveRank(A, cs2, R,          'outputType'='CS');
            results[6]  := ComprehensiveRank(A, rs2, R,          'outputType'='CS');
            results[7]  := ComprehensiveRank(A, [], [lambda], R, 'outputType'='CS');
            
            results[8]  := ComprehensiveRank(A, cs1, R,           'outputType'='RS');
            results[9]  := ComprehensiveRank(A, rs1, R,           'outputType'='RS');
            results[10]  := ComprehensiveRank(A, [lambda], R,     'outputType'='RS');
            results[11]  := ComprehensiveRank(A, [lambda], [], R, 'outputType'='RS');
            
            results[12]  := ComprehensiveRank(A, cs2, R,          'outputType'='RS');
            results[13]  := ComprehensiveRank(A, rs2, R,          'outputType'='RS');
            results[14]  := ComprehensiveRank(A, [], [lambda], R, 'outputType'='RS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        for i from 1 to 4 do
            print(i);
            correct := verifyResult_cs(A, cs1, results[i], R);
            if not correct then
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        for i from 5 to 7 do
            print(i);
            correct := verifyResult_cs(A, cs2, results[i], R);
            if not correct then
                printf("\b\b\b");
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        for i from 8 to 11 do
            print(i);
            correct := verifyResult_rs(A, cs1, results[i], R);
            if not correct then
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        for i from 12 to 14 do
            print(i);
            correct := verifyResult_rs(A, cs2, results[i], R);
            if not correct then
                printf("\b\b\b");
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        
        
        printf("\b\b\b");
        printf("Pass\n");
        return true;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 4                                                              #
    #                                                                     #
    # A = [-4*z[1, 1]-4*z[1, 2], -4*z[1, 2]-4*z[1, 3], 20*z[1, 3]+...     #
    #        24*z[1, 1]+44*z[1, 2]],                                      #
    #     [-7*z[1, 1]-6*z[1, 2]+z[1, 3], -18*z[1, 2]-12*z[1, 3]-...       #
    #       6*z[1, 1], 54*z[1, 3]+72*z[1, 1]+126*z[1, 2]                  #
    #     [-z[2, 1]+z[2, 3], -12*z[2, 2]-6*z[2, 1]-6*z[2, 3],             #
    #       24*z[2, 3]+60*z[2, 2]+36*z[2, 1]]                             #
    # cs := {}                                                            #
    # ------------------------------------------------------------------- #
    test4 := proc($)
    
        local A, R, cs, rs, i, results, correct;
        
        A := Matrix([[-4*z[1, 1]-4*z[1, 2], -4*z[1, 2]-4*z[1, 3], 20*z[1, 3]+24*z[1, 1]+44*z[1, 2]], [-7*z[1, 1]-6*z[1, 2]+z[1, 3], -18*z[1, 2]-12*z[1, 3]-6*z[1, 1], 54*z[1, 3]+72*z[1, 1]+126*z[1, 2]], [-z[2, 1]+z[2, 3], -12*z[2, 2]-6*z[2, 1]-6*z[2, 3], 24*z[2, 3]+60*z[2, 2]+36*z[2, 1]]]);
        
        R := PolynomialRing([z[1, 1], z[1, 2], z[1, 3], z[2, 1], z[2, 2], z[2, 3]]);
        
        cs := GeneralConstruct([],[],R);
        rs := RepresentingRegularSystems(cs, R)[1];
        
        try
            results := [0$10];
            results[1]  := ComprehensiveRank(A, R,         'outputType'='CS');
            results[2]  := ComprehensiveRank(A, cs, R,     'outputType'='CS');
            results[3]  := ComprehensiveRank(A, rs, R,     'outputType'='CS');
            results[4]  := ComprehensiveRank(A, [], R,     'outputType'='CS');
            results[5]  := ComprehensiveRank(A, [], [], R, 'outputType'='CS');
            results[6]  := ComprehensiveRank(A, R,         'outputType'='RS');
            results[7]  := ComprehensiveRank(A, cs, R,     'outputType'='RS');
            results[8]  := ComprehensiveRank(A, rs, R,     'outputType'='RS');
            results[9]  := ComprehensiveRank(A, [], R,     'outputType'='RS');
            results[10] := ComprehensiveRank(A, [], [], R, 'outputType'='RS');
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        for i to 5 do
            correct := verifyResult_cs(A, cs, results[i], R);
            if not correct then
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        for i from 6 to 10 do
            correct := verifyResult_rs(A, cs, results[i], R);
            if not correct then
                printf("\b\b\b");
                printf("FAIL: Incorrect result\n");
                return false;
            end if;
        end do;
        
        printf("\b\b\b");
        printf("Pass\n");
        return true;
    
    end proc;
    
end module: