SetClassGroupBounds("GRH");
SetNthreads(4);
AttachSpec("proj");

procedure TestMaximality()

    print "Testing non-maximal cubics...";
    assert CheckMaximalAtP(4, 2, 1, 1, 2) eq false;
    assert CheckMaximalAtP(25, 5, 1, 2, 5) eq false;
    assert CheckMaximalAtP(1, 1, 3, 9, 3) eq false;
    assert CheckMaximalAtP(1, 0, 0, -17, 3) eq false;
    assert CheckMaximalAtP(1, 0, -4, -16, 2) eq false;
    print "Passed non-maximal.";

    print "Testing maximal cubics...";
    assert CheckMaximalAtP(1, 0, -1, -1, 23) eq true;

    //factorisation 3 7 463
    assert CheckMaximalAtP(15, 11, 5, 2, 2) eq true;
    
    assert CheckMaximalAtP(19, 2, 3, 1, 9743) eq true;

    
    assert CheckMaximalAtP(-4, -9, 5, 1, 9749) eq true;
    // factorisation[ <7, 1>, <11, 1>, <127, 1> ]
    assert CheckMaximalAtP(22, 5, 3, 1, 7) eq true;
    assert CheckMaximalAtP(22, 5, 3, 1, 11) eq true;
    assert CheckMaximalAtP(22, 5, 3, 1, 127) eq true;
    print "Passed maximal checks.";
    
end procedure;

procedure TestMapForClassGroup(X, is_real)
    
    Q<x> := PolynomialRing(Rationals());
    
    coeff := GenerateCubics(X, is_real);
    
    g := Q!coeff;
    print "Cubic Polynomial: ", g;
    L := SplittingField(g);
    f := DefiningPolynomial(L);
    print "Sextic Polynomial: ", f;
    L_coeff := Coefficients(f);
    f := f / LeadingCoefficient(f);
    L := NumberField(f);
    L := OptimisedRepresentation(L);
    ZL := MaximalOrder(L);
    Cl, m := ClassGroup(ZL);

    print "Class group ", AbelianInvariants(Cl);

    G, Aut, mG := AutomorphismGroup(L);
    sigma := mG([ g : g in G | Order(g) eq 3 ][1]);

    
    subs := Subfields(L);
    
    K1_data := [ s : s in subs | Degree(s[1]) eq 3 ][1];
    F_data  := [ s : s in subs | Degree(s[1]) eq 2 ][1];

    K1, mK1 := Explode(K1_data);
    F, mF := Explode(F_data);
    
    
    Cl1, m1 := ClassGroup(K1);
    Cl4, m4 := ClassGroup(F);
    
    Remove(~subs, 1);
    LoverK1 := RelativeField(K1, L);
    LoverF := RelativeField(F, L);
    
    Zk1 := MaximalOrder(LoverK1);
    Zf := MaximalOrder(LoverF);
    
    
    autsF := Automorphisms(F);
    tau := [ a : a in autsF | a(F.1) ne F.1 ][1];

    Phi := function(I)
        // We need to compute the relative norms.
        // We just view the generators of I in the field over F and over K.

        // In the case its quadratic, I need tau(Nm(I)) to do that, first view it as an ideal in LoverF
        I_F := ideal< Zf | [ LoverF | g : g in Generators(I) ] >;
        N_F := Norm(I_F);
        N_F_tau := ideal< MaximalOrder(F) | [ tau(g) : g in Generators(N_F) ] >;
        
        // This is the 2nd cubic field, i.e. sigma(L)
        I_sig1 := ideal< ZL | [ sigma(L!g) : g in Generators(I) ] >;       // sigma^1
        // This is the 3nd cubic field, i.e. sigma(sigma(L))
        I_sig2 := ideal< ZL | [ sigma(sigma(L!g)) : g in Generators(I) ] >; // sigma^2

        // View them as ideals
        I_K1_1 := ideal< Zk1 | [ LoverK1 | g : g in Generators(I) ] >;
        I_K1_2 := ideal< Zk1 | [ LoverK1 | g : g in Generators(I_sig2) ] >;
        I_K1_3 := ideal< Zk1 | [ LoverK1 | g : g in Generators(I_sig1) ] >;
        // Take the norm
        return <Norm(I_K1_1), Norm(I_K1_2), Norm(I_K1_3), N_F_tau>;
    end function;
    
    Cl_gens := [ Cl.i : i in [1..Ngens(Cl)] ];
    ideal_gens := [ m(g) : g in Cl_gens ];
    
    D, inc, proj := DirectSum([Cl1, Cl1 , Cl1, Cl4]);
    images_in_D := [];
    for I in ideal_gens do
        // For every generator of Cl, we need to view the images in the map Phi.
        // K x K x K x F
        J := Phi(I);
        
        // The 3 cubic fields all have the same class group.
        c1 := J[1] @@ m1;
        c2 := J[2] @@ m1;
        c3 := J[3] @@ m1;

        // Quadratic field
        c4 := J[4] @@ m4;

        d_elem := inc[1](c1) + inc[2](c2) + inc[3](c3) + inc[4](c4);
        Append(~images_in_D, d_elem);
    end for;
    
    // Now create the map phi
    phi := hom< Cl -> D | images_in_D>;
    
    im_phi := Image(phi);
    ker_phi := Kernel(phi);
    coker_phi := quo < D | im_phi >;
    
    
    Psi := function(x)
        // Given an element in D, we want to view it in the field L. 
        c_1 := proj[1](x);
        c_2 := proj[2](x);
        c_3 := proj[3](x);
        c_4 := proj[4](x);
        
        
        I_1 := m1(c_1);
        I_2 := m1(c_2);
        I_3 := m1(c_3);
        I_4 := m4(c_4);
        
        J_1 := ideal< ZL | [ mK1(g) : g in Generators(I_1) ] >;
        J_2 := ideal< ZL | [ sigma(mK1(g)) : g in Generators(I_2) ] >;
        J_3 := ideal< ZL | [ sigma(sigma(mK1(g))) : g in Generators(I_3) ] >;
        J_4 := ideal< ZL | [ mF(g) : g in Generators(I_4) ] >;
    
        // Need -J_4 as a_4 = -1
        
        return (J_1 @@ m) + (J_2 @@ m) + (J_3 @@ m) - (J_4 @@ m);
    end function;
    
    psi := hom< D -> Cl | [ Psi(D.i) : i in [1..Ngens(D)] ] >;
    
    im_psi := Image(psi);
    ker_psi := Kernel(psi);
    coker_psi := quo < Cl | im_psi>;

    // If it were successful, P(x) = 3 * x for any x;
    P := phi * psi;
    
    for cl in Generators(Cl) do
        assert 3*cl eq P(cl);
    end for;
    // If there is no 3 primary part, the map psi must be surjective.
    if #Cl mod 3 ne 0 then
        assert #coker_psi eq 1;
    end if;
    print "Passed";

end procedure;

procedure RunAllTests()
    print "Starting tests...";

    TestMaximality();
    
    SEED := 42;
    SetSeed(SEED);
    print "Running Cases for Class group...";

    // We just do some basic sanity tests, composing phi then psi should be eq to multiplying by 3.
    // In the case the map has no 3 sylow part, the map Psi should ALWAYS be surjective.
    for _ in [1..100] do
        TestMapForClassGroup(100, true);
    end for;


end procedure;


RunAllTests();


quit;