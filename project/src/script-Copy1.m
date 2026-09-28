SetClassGroupBounds("GRH");
//SetNthreads(4);
AttachSpec("proj");

X := StringToInteger(X);
is_real := (is_real eq "1");
Q<x> := PolynomialRing(Rationals());

//coeff := GenerateCubics(X, is_real);


//g := Q!coeff;

//g := -66*x^3 - 27*x^2 + 73*x + 10;
g := -87*x^3 - 65*x^2 + 98*x + 6; 
L := SplittingField(g);


f_def := DefiningPolynomial(L);
L_coeff := Coefficients(f_def);
f_def := f_def / LeadingCoefficient(f_def);

L := OptimisedRepresentation(NumberField(f_def));
ZL := MaximalOrder(L);
Cl, m := ClassGroup(ZL);


G, Aut, mG := AutomorphismGroup(L);
sigma := mG([ g : g in G | Order(g) eq 3 ][1]);


subs := Subfields(L);
K1_data := [ s : s in subs | Degree(s[1]) eq 3 ][1];
F_data  := [ s : s in subs | Degree(s[1]) eq 2 ][1];

K1, mK1 := Explode(K1_data);
F, mF := Explode(F_data);

f_quad := DefiningPolynomial(F);
autsF := Automorphisms(F);
tau := [ a : a in autsF | a(F.1) ne F.1 ][1];

D_k := Discriminant(MaximalOrder(K1));
D_f := Discriminant(MaximalOrder(F));
f := Isqrt(Integers()!(D_k / FundamentalDiscriminant(D_k)));

D_L := Discriminant(ZL);
r_L, c_L := Signature(L);
label := Sprintf("6.%o.%o.1", r_L, Abs(D_L));

coeff := Coefficients(g); 
sextic_info := <label, L_coeff, c_L, D_L, Factorisation(D_L), Coefficients(f_quad), coeff, AbelianInvariants(Cl), #Cl, f>;
print sextic_info;

Cl1, m1 := ClassGroup(K1);
Cl4, m4 := ClassGroup(F);

Sylow_Cl1 := SylowSubgroup(Cl1, 3);
Sylow_Cl4 := SylowSubgroup(Cl4, 3);

r_k, c_k := Signature(K1);
label_cubic := Sprintf("3.%o.%o.1", r_k, Abs(D_k));
cubic_info := <label_cubic, coeff, c_k, D_k, Factorisation(D_k), AbelianInvariants(Cl1), #Cl1, AbelianInvariants(Sylow_Cl1), #Sylow_Cl1>;
print cubic_info;

LoverK1 := RelativeField(K1, L);
LoverF := RelativeField(F, L);
Zk1 := MaximalOrder(LoverK1);
Zf := MaximalOrder(LoverF);

D_rel := Discriminant(Zf);
rel_factors := [ < [ Eltseq(b) : b in Basis(x[1]) ], x[2] > : x in Factorization(D_rel) ];

r_f, c_f := Signature(F);
label_quad := Sprintf("2.%o.%o.1", r_f, Abs(D_f));
quad_info := <label_quad, Coefficients(f_quad), c_f, D_f, Factorisation(D_f), rel_factors, #rel_factors, AbelianInvariants(Cl4), #Cl4, AbelianInvariants(Sylow_Cl4), #Sylow_Cl4>;
print quad_info;



Phi := function(I)
    gens := Generators(I);
    gens_F := [ LoverF | g : g in gens ];
    I_F := ideal< Zf | gens_F >;
    N_F := Norm(I_F);
    N_F_tau := ideal< MaximalOrder(F) | [ tau(x) : x in Generators(N_F) ] >;
    
    gens_K1 := [ LoverK1 | g : g in gens ];
    norm_K1 := Norm(ideal< Zk1 | gens_K1 >);

    return <norm_K1, norm_K1, norm_K1, N_F_tau^-1>;
end function;


ideal_gens := [ m(Cl.i) : i in [1..Ngens(Cl)] ];

D, inc, proj := DirectSum([Cl1, Cl1, Cl1, Cl4]);
print <Invariants(Cl), #Cl, Invariants(D), #D>;

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

phi := hom< Cl -> D | images_in_D>;
im_phi := Image(phi);
ker_phi := Kernel(phi);
coker_phi := quo < D | im_phi >;
print <Invariants(im_phi), #im_phi, Invariants(ker_phi), #ker_phi, Invariants(coker_phi), #coker_phi>;

Psi := function(x)
    P := ideal< ZL | 1 >; 
    P *:= ideal< ZL | [ mK1(g) : g in Generators(m1(proj[1](x))) ] >;
    P *:= ideal< ZL | [ sigma(mK1(g)) : g in Generators(m1(proj[2](x))) ] >;
    P *:= ideal< ZL | [ sigma(sigma(mK1(g))) : g in Generators(m1(proj[3](x))) ] >;
    P *:= ideal< ZL | [ mF(g) : g in Generators(m4(proj[4](x))) ] >;
    return P @@ m;
end function;

psi := hom< D -> Cl | [ Psi(D.i) : i in [1..Ngens(D)] ] >;
im_psi := Image(psi);
ker_psi := Kernel(psi);
coker_psi := quo < Cl | im_psi>;

P := phi * psi;


// What we expect is, we have L(Cl.i) == 3 * Cl.i;
for cl in Generators(Cl) do
    print 3*cl, P(cl);
    assert 3*cl eq P(cl);
end for;



