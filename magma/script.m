SetClassGroupBounds("GRH");
//SetNthreads(4);
//Alarm(300);

X := StringToInteger(X);
is_real := StringToInteger(is_real);
Q<x> := PolynomialRing(Rationals());


check_maximal := function(a,b,c,d, Disc)

    if Disc eq 0 then
        return false;
    end if;

    // 3. Find candidate primes (p^2 | Disc)
    facts := Factorization(Abs(Disc));
    cand_primes := [ f[1] : f in facts | f[2] ge 2 ];

    // 4. Test each candidate prime
    for p in cand_primes do


        if (a mod p eq 0) and (b mod p eq 0) then
            if a mod p^2 eq 0 then
                return false;
            end if;

        else
            Fp := GF(p);
            P<x> := PolynomialRing(Fp);
            g := Fp!a*x^3 + Fp!b*x^2 + Fp!c*x + Fp!d;
            dg := Derivative(g);

            // The multiple root modulo p is a root of gcd(g, g')
            common_poly := GCD(g, dg);
            mult_roots := Roots(common_poly);

            for r in mult_roots do
                u := Integers()!r[1]; // Lift root back to Z

                // Evaluate f(u, 1) mod p^2
                val := a*u^3 + b*u^2 + c*u + d;
                if val mod p^2 eq 0 then
                    return false; // Fails maximality at p
                end if;
            end for;
        end if;
    end for;

    return true; // Passed tests for all candidate primes
end function;

generate_cubics := function(X, N, is_real)
    cubics := [];
    num_gen := 0;
    while num_gen lt N do
        a := Random(X);
        b := Random(X);
        c := -X + Random(2*X);
        d := -X + Random(2*X);
        if a eq 0 then continue; end if;
    
        is_reduced := false;
    
        P := b^2-3*a*c;
        Q := b*c-9*a*d;
        R := c^2 - 3*b*d;
        // Check if reduced
        if (b gt 0 or d lt 0) and
           (Q ne 0 or d lt 0) and
           (P ne Q or b lt Abs(3*a - b)) and
           (P ne R or (a lt Abs(d) and
               (a ne Abs(d) or b lt Abs(c)))) then
            is_reduced := true;
        end if;
        // If not reduced, return
        if not is_reduced then continue; end if;
        // 1. Primitivity check
        if GCD([a, b, c, d]) gt 1 then
            continue;
        end if;
        Disc := b^2*c^2 - 4*a*c^3 - 4*b^3*d - 27*a^2*d^2 + 18*a*b*c*d;
        if is_real eq 1 and Disc lt 0 then continue; end if;
        if is_real eq 0 and Disc gt 0 then continue; end if;

        // Check square free
        if IsSquare(Disc) then continue; end if;
    
        f := a*x^3 + b*x^2 +c*x + d;
    
        if not IsIrreducible(f) then continue; end if;
    
        is_maximal := check_maximal(a,b,c,d,Disc);
        if not is_maximal then continue; end if;
    
        num_gen +:= 1;
        Append(~cubics, [a,b,c,d]);

    end while;
    return cubics;
end function;



coeff := generate_cubics(X, 1, is_real)[1];

g := Q!coeff;

//k := NumberField(g);
L := SplittingField(g);
L := OptimisedRepresentation(L);
ZL := MaximalOrder(L);
auts := Automorphisms(L);
Cl, m := ClassGroup(ZL);



subs := Subfields(L);

// Cubics
K1, mK1 := subs[3][1], subs[3][2];

// Quadratic
F, mF := subs[2][1], subs[2][2];

f_quad := DefiningPolynomial(F);

D_k := Discriminant(MaximalOrder(K1));
D_f := Discriminant(MaximalOrder(F));

f := Isqrt(Integers()!(D_k / FundamentalDiscriminant(D_k)));


D_L := Discriminant(ZL);
r_L, c_L := Signature(L);
label := Sprintf("%o.%o.%o.1", 6, r_L, Abs(D_L));

L_coeff := Coefficients(DefiningPolynomial(L));

sextic_info := <label, L_coeff, c_L, D_L, Factorisation(D_L), Coefficients(f_quad), coeff, AbelianInvariants(Cl), #Cl, f>;
print sextic_info;

Cl1, m1 := ClassGroup(K1);
Cl4, m4 := ClassGroup(F);

Sylow_Cl1 := SylowSubgroup(Cl1, 3);
Sylow_Cl4 := SylowSubgroup(Cl4, 3);

// cubic info

r_k, c_k := Signature(K1);
label_cubic := Sprintf("%o.%o.%o.1", 3, r_k, Abs(D_k));
cubic_info := <label_cubic, coeff, c_k, D_k, Factorisation(D_k), AbelianInvariants(Cl1), #Cl1, AbelianInvariants(Sylow_Cl1), #Sylow_Cl1>;
print cubic_info;

// quad info

r_f, c_f := Signature(F);
label_quad := Sprintf("%o.%o.%o.1", 2, r_f, Abs(D_f));
quad_info := <label_quad, Coefficients(f_quad), c_f, D_f, Factorisation(D_f), AbelianInvariants(Cl4), #Cl4, AbelianInvariants(Sylow_Cl4), #Sylow_Cl4>;
print quad_info;

if #Cl mod 3 ne 0 then
    print <[], 1, [], 1>;
    print <[], 1, [], 1, [], 1>;
    print <[], 1, [], 1, [], 1>;
    quit;
end if;

Remove(~subs, 1);
LoverK1 := RelativeField(K1, L);
LoverF := RelativeField(F, L);

Zk1 := MaximalOrder(LoverK1);
Zf := MaximalOrder(LoverF);


sigma := auts[1];

for aut in auts do
    if aut(aut(aut(L.1))) eq L.1 and aut(L.1) ne L.1 then
        sigma := aut;
        break;
    end if;
end for;

autsF := Automorphisms(F);
tau := autsF[1];
for aut in autsF do
    if aut(F.1) ne F.1 then
        tau := aut;
    end if;
end for;

Phi := function(I)
    gens := Generators(I);
    
    gens_F := [ LoverF | g : g in gens ];
    I_F := ideal< Zf | gens_F >;
    N_F := Norm(I_F);
    N_F_tau := ideal< MaximalOrder(F) | [ tau(g) : g in Generators(N_F) ] >;
    
    gens_K1 := [ LoverK1 | g : g in gens ];
    I_K1 := ideal< Zk1 | gens_K1 >;
    
    gens_K2 := [ LoverK1 | sigma(sigma(g)) : g in gens ];
    I_K2 := ideal< Zk1 | gens_K2 >;
    
    gens_K3 := [ LoverK1 | sigma(g) : g in gens ];
    I_K3 := ideal< Zk1 | gens_K3 >;

    return <Norm(I_K1), Norm(I_K2), Norm(I_K3), N_F_tau^-1>;
end function;

Sylow3 := SylowSubgroup(Cl, 3);
gens_Sylow3 := [ Sylow3.i : i in [1..Ngens(Sylow3)] ];
ideal_gens := [ m(g) : g in gens_Sylow3 ];
phis := [Phi(I) : I in ideal_gens];



D, inc, proj := DirectSum([Sylow_Cl1, Sylow_Cl1 , Sylow_Cl1, Sylow_Cl4]);
print <Invariants(Sylow3), #Sylow3, Invariants(D), #D>;

images_in_D := [];
for I in ideal_gens do
    J := Phi(I);
    
    c1 := Inverse(m1)(J[1]); 
    c2 := Inverse(m1)(J[2]); 
    c3 := Inverse(m1)(J[3]); 
    c4 := Inverse(m4)(J[4]); 
    
    d_elem := inc[1](c1) + inc[2](c2) + inc[3](c3) + inc[4](c4);
    Append(~images_in_D, d_elem);
end for;

phi := hom< Sylow3 -> D | images_in_D>;

im_phi := Image(phi);
ker_phi := Kernel(phi);
coker_phi := quo < D | im_phi >;

print <Invariants(im_phi), #im_phi, Invariants(ker_phi), #ker_phi, Invariants(coker_phi), #coker_phi>;

Psi := function(x)
    P := ideal< ZL | 1 >; 
    
    a_1 := m1(proj[1](x)); 
    P *:= ideal< ZL | [ mK1(g) : g in Generators(a_1) ] >;
    
    a_2 := m1(proj[2](x)); 
    P *:= ideal< ZL | [ sigma(mK1(g)) : g in Generators(a_2) ] >;
    
    a_3 := m1(proj[3](x)); 
    P *:= ideal< ZL | [ sigma(sigma(mK1(g))) : g in Generators(a_3) ] >;
    
    a_4 := m4(proj[4](x)); 
    P *:= ideal< ZL | [ mF(g) : g in Generators(a_4) ] >;
    
    return Inverse(m)(P);
end function;

psi := hom< D -> Sylow3 | [ Psi(D.i) : i in [1..Ngens(D)] ] >;

im_psi := Image(psi);
ker_psi := Kernel(psi);
coker_psi := quo < Sylow3 | im_psi>;

print <Invariants(im_psi), #im_psi, Invariants(ker_psi), #ker_psi, Invariants(coker_psi), #coker_psi>;

exit;
