SetClassGroupBounds("GRH");

SetNthreads(4);
// Same as the manager process.
host := "192.168.0.3";
port := 10000;
SetNthreads(4);
SetClassGroupBounds("GRH");


function foo(arg)
    // Input is a cubic
    Q<x> := PolynomialRing(Rationals());

    g := Q!arg;
    k := NumberField(g);
    K := SplittingField(k);

    ZK := MaximalOrder(K);

    Cl, m := ClassGroup(ZK);

     P := Subfields(K);
    // quadratic sub;
    P1 := P[2];
    f1 := DefiningPolynomial(P1[1]);
    // lmfdb labelling
    D_K := Discriminant(ZK);
    r_K, c_K := Signature(K);
    label := Sprintf("%o.%o.%o.1", 6, r_K, Abs(D_K));   
    subs := Subfields(K);
    // class num must be divisible by 3
    if #subs ne 5 or #Cl mod 3 ne 0 then
        result := <
            label,
            <Coefficients(DefiningPolynomial(ZK)), Signature(K), D_K, Factorisation(D_K), Coefficients(f1), arg, AbelianInvariants(Cl), #Cl>,
            < [-1], -1 >, // 1. Sylow3
            < [-1], -1 >, // 2. D
            < [-1], -1 >, // 3. im_phi
            < [-1], -1 >, // 4. ker_phi
            < [-1], -1 >, // 5. coker_phi
            < [-1], -1 >, // 6. im_psi
            < [-1], -1 >, // 7. ker_psi
            < [-1], -1 >  // 8. coker_psi
        >;
        return result;
    end if;


    // Cubics
    F1, mF1 := subs[3][1], subs[3][2];
    F2, mF2 := subs[4][1], subs[4][2];
    F3, mF3 := subs[5][1], subs[5][2];

    // Quadratic
    F, mF := subs[2][1], subs[2][2];

    Remove(~subs, 1);
    KoverF1 := RelativeField(F1, K);
    KoverF2 := RelativeField(F2, K);
    KoverF3 := RelativeField(F3, K);

    KoverF := RelativeField(F, K);

    Zf1 := MaximalOrder(KoverF1);
    Zf2 := MaximalOrder(KoverF2);
    Zf3 := MaximalOrder(KoverF3);

    Zf := MaximalOrder(KoverF);

    tau := Automorphisms(F)[2];

    Phi := function(I)
        gens := Generators(I);
    
        I_F := ideal< Zf |  gens>;
    
    
        I_F1 := ideal< Zf1 | gens>;
        I_F2 := ideal< Zf2 | gens >;
        I_F3 := ideal< Zf3 | gens >;

        return <Norm(I_F1), Norm(I_F2), Norm(I_F3), tau(Norm(I_F))^-1>;
    
    end function;

    Sylow3 := SylowSubgroup(Cl, 3);





    gens_Sylow3 := [ Sylow3.i : i in [1..Ngens(Sylow3)] ];

    ideal_gens := [ m(g) : g in gens_Sylow3 ];

    phis := [Phi(I) : I in ideal_gens];

    Cl1, m1 := ClassGroup(F1);
    Cl2, m2 := ClassGroup(F2);
    Cl3, m3 := ClassGroup(F3);
    Cl4, m4 := ClassGroup(F);
    D, inc, proj := DirectSum([SylowSubgroup(Cl1, 3),SylowSubgroup(Cl2,3 ),SylowSubgroup(Cl3, 3),SylowSubgroup(Cl4, 3)]);




    images_in_D := [];
    for I in ideal_gens do
        J := Phi(I);
    
        c1 := Inverse(m1)(J[1]); 
        c2 := Inverse(m2)(J[2]); 
        c3 := Inverse(m3)(J[3]); 
        c4 := Inverse(m4)(J[4]); 

    
        d_elem := inc[1](c1) + inc[2](c2) + inc[3](c3) + inc[4](c4);
    
        Append(~images_in_D, d_elem);
    
    end for;


    phi := hom< Sylow3 -> D | images_in_D>;

    im_phi := Image(phi);
    ker_phi := Kernel(phi);
    coker_phi := quo < D | im_phi >;

    ms := [m1, m2, m3, m4];

    Psi := function(x)

        P := ideal< ZK | 1 >; 
    

        for i in [1..#proj] do

            c_i := proj[i](x);
        

            a_i := ms[i](c_i); 
        
            P *:= ideal< ZK | a_i >;
        end for;
    
        return Inverse(m)(P);
    end function;

    psi := hom< D -> Sylow3 | [ Psi(D.i) : i in [1..Ngens(D)] ] >;

    im_psi := Image(psi);
    ker_psi := Kernel(psi);
    coker_psi := quo < Sylow3 | im_psi>;
     
    result := <
        label,
        <Coefficients(DefiningPolynomial(ZK)), Signature(K), D_K, Factorisation(D_K), Coefficients(f1), arg, AbelianInvariants(Cl), #Cl>,
        < Invariants(Sylow3),    #Sylow3 >,   
        < Invariants(D),         #D >,       
        < Invariants(im_phi),    #im_phi >, 
        < Invariants(ker_phi),   #ker_phi >, 
        < Invariants(coker_phi), #coker_phi >,
        < Invariants(im_psi),    #im_psi >,  
        < Invariants(ker_psi),   #ker_psi >,
        < Invariants(coker_psi), #coker_psi >
    >;
    return result;
end function;

DistributedWorker(host, port, foo);

quit;
