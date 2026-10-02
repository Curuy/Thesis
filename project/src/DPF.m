SetClassGroupBounds("GRH");

AttachSpec("proj");
//X := StringToInteger(X);
//is_real := (is_real eq "1");

Q<x> := PolynomialRing(Rationals());

coeff := GenerateCubics(X, is_real);

//coeff :=  [6, 22, -71, -33]; // sgn + w^*^2 which is type alpha
//coeff := [16, 13, -85, -12]; //W* + P* which is type delta
//coeff := [7, 96, -24, -58]; // sgn + W^2 type gamma

//coeff := [41, 97, -61, -47]; // sgn + w + w* which is type beta

//coeff := [23, 97, -83, -75]; // W + P* which is type epsilon


g := Q!coeff;



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
K_data := [ s : s in subs | Degree(s[1]) eq 3 ][1];
F_data  := [ s : s in subs | Degree(s[1]) eq 2 ][1];

K, iK := Explode(K_data);
F, iF := Explode(F_data);

ZF := MaximalOrder(F);
ZK := MaximalOrder(K);

ClK, mK := ClassGroup(ZK);
ClF, mF := ClassGroup(ZF);



f_quad := DefiningPolynomial(F);
autsF := Automorphisms(F);
tau := [ a : a in autsF | a(F.1) ne F.1 ][1];

D_k := Discriminant(MaximalOrder(K));
D_f := Discriminant(MaximalOrder(F));
f := Isqrt(Integers()!(D_k / FundamentalDiscriminant(D_k)));

D_L := Discriminant(ZL);
r_L, c_L := Signature(L);
label := Sprintf("6.%o.%o.1", r_L, Abs(D_L));

coeff := Coefficients(g); 
sextic_info := <label, L_coeff, c_L, D_L, Factorisation(D_L), Coefficients(f_quad), coeff, AbelianInvariants(Cl), #Cl, f>;
print sextic_info;

Sylow_ClK := SylowSubgroup(ClK, 3);
Sylow_ClF := SylowSubgroup(ClF, 3);

r_k, c_k := Signature(K);
label_cubic := Sprintf("3.%o.%o.1", r_k, Abs(D_k));
cubic_info := <label_cubic, coeff, c_k, D_k, Factorisation(D_k), AbelianInvariants(ClK), #ClK, AbelianInvariants(Sylow_ClK), #Sylow_ClK>;
print cubic_info;

LoverK := RelativeField(K, L);
LoverF := RelativeField(F, L);

Zk := MaximalOrder(LoverK);
Zf := MaximalOrder(LoverF);

D_rel := Discriminant(Zf);
rel_factors := [ < [ Eltseq(b) : b in Basis(x[1]) ], x[2] > : x in Factorization(D_rel) ];

r_f, c_f := Signature(F);
label_quad := Sprintf("2.%o.%o.1", r_f, Abs(D_f));
quad_info := <label_quad, Coefficients(f_quad), c_f, D_f, Factorisation(D_f), rel_factors, #rel_factors, AbelianInvariants(ClF), #ClF, AbelianInvariants(Sylow_ClF), #Sylow_ClF>;
print quad_info;


function DPFInvariants(L)

    classImages := [Cl | ];

    for j in [1..Ngens(ClF)] do
        a := mF(ClF.j);
        aL := ideal< ZL | [iF(b) : b in Basis(a, F)] >;
        Append(~classImages, aL @@ m);
    end for;

    T := hom< ClF -> Cl | classImages >;
    capitulation := Kernel(T);

    C := Valuation(#capitulation, 3);
    assert #capitulation eq 3^C;

    ramPrimes := [
        pe[1] : pe in Factorization(Abs(Discriminant(ZK)))
    ];

    absoluteImages := [ClK | ];

    // Want to compute A which is dimension of P_L / P

    for q in ramPrimes do
        for de in Decomposition(ZK, q) do
            // This is the condition to be in the ambigious ideal, I^3 in Q
            if de[2] eq 3 then
                Append(~absoluteImages, de[1] @@ mK);
            end if;
        end for;
    end for;

    tA := #absoluteImages;

    // Now we want the principal ones.
    if tA eq 0 then
        A := 0;
    else
        // We know that the ambigious factors is isomorphic to (Z/Z3)^{tA}
        VA := AbelianGroup([3 : j in [1..tA]]);
        fA := hom< VA -> ClK | absoluteImages >;
        // It is principal if its in the kernal so,
        A := Valuation(#Kernel(fA), 3);
    end if;

    // We now want to compute the relative ambigious ideals, I_L| F / I_F

    // Cl(L) / I_K
    Qcl, pi := quo< Cl | classImages >;

    relativeImages := [Qcl | ];

    for q in ramPrimes do
        for de in Decomposition(ZL, q) do
            if de[2] mod 3 eq 0 then
                Append(~relativeImages, pi(de[1] @@ m));
            end if;
        end for;
    end for;

    tN := #relativeImages;

    if tN eq 0 then
        D := 0;
    else
        VN := AbelianGroup([3 : j in [1..tN]]);
        // This basically gives the map from relative ambigious class quo frac ideals into CL(L) / I_F
        // I_L|F / I_F -> I_L / P_L I_F.
        fN := hom< VN -> Qcl | relativeImages >;
        // The kernel of this map is basically all relative quotients that map to principal ideals and so,
        // I_L|F meets (P_L I_F) / I_F
        D := Valuation(#Kernel(fN), 3);
    end if;

    // From eq 3.8 in https://arxiv.org/pdf/2102.12187 we have D := A + R and so,
    R := D - A;
    U := D + C - 1;

    assert R ge 0;

    // U must be 0 or 1
    assert U in {0, 1};

    return A, R, C, U;
end function;



A_dpf, R_dpf, C_dpf, U_dpf := DPFInvariants(L);
printf "DPF: %o %o %o %o\n", A_dpf, R_dpf, C_dpf, U_dpf;



if #Cl mod 3 ne 0 then
    print <[], 1, [], 1>;
    print <[], 1, [], 1, [], 1>;
    print <[], 1, [], 1, [], 1>;
    quit;
end if;

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
    I_K_1 := ideal< Zk | [ LoverK | g : g in Generators(I) ] >;
    I_K_2 := ideal< Zk | [ LoverK | g : g in Generators(I_sig2) ] >;
    I_K_3 := ideal< Zk | [ LoverK | g : g in Generators(I_sig1) ] >;
    // Take the norm
    return <Norm(I_K_1), Norm(I_K_2), Norm(I_K_3), N_F_tau>;
end function;

Sylow3 := SylowSubgroup(Cl, 3);
ideal_gens := [ m(Sylow3.i) : i in [1..Ngens(Sylow3)] ];

D, inc, proj := DirectSum([Sylow_ClK, Sylow_ClK, Sylow_ClK, Sylow_ClF]);
print <Invariants(Sylow3), #Sylow3, Invariants(D), #D>;

images_in_D := [];
for I in ideal_gens do
    // For every generator of Cl, we need to view the images in the map Phi.
    // K x K x K x F
    J := Phi(I);
    
    // The 3 cubic fields all have the same class group.
    c1 := J[1] @@ mK;
    c2 := J[2] @@ mK;
    c3 := J[3] @@ mK;

    // Quadratic field
    c4 := J[4] @@ mF;

    d_elem := inc[1](c1) + inc[2](c2) + inc[3](c3) + inc[4](c4);
    Append(~images_in_D, d_elem);
end for;

phi := hom< Sylow3 -> D | images_in_D>;
im_phi := Image(phi);
ker_phi := Kernel(phi);
coker_phi := quo < D | im_phi >;
print <Invariants(im_phi), #im_phi, Invariants(ker_phi), #ker_phi, Invariants(coker_phi), #coker_phi>;

Psi := function(x)
    // Given an element in D, we want to view it in the field L. 
    c_1 := proj[1](x);
    c_2 := proj[2](x);
    c_3 := proj[3](x);
    c_4 := proj[4](x);
    
    
    I_1 := mK(c_1);
    I_2 := mK(c_2);
    I_3 := mK(c_3);
    I_4 := mF(c_4);
    
    J_1 := ideal< ZL | [ iK(g) : g in Generators(I_1) ] >;
    J_2 := ideal< ZL | [ sigma(iK(g)) : g in Generators(I_2) ] >;
    J_3 := ideal< ZL | [ sigma(sigma(iK(g))) : g in Generators(I_3) ] >;
    J_4 := ideal< ZL | [ iF(g) : g in Generators(I_4) ] >;

    // Need -J_4 as a_4 = -1
    
    return (J_1 @@ m) + (J_2 @@ m) + (J_3 @@ m) - (J_4 @@ m);
end function;
psi := hom< D -> Sylow3 | [ Psi(D.i) : i in [1..Ngens(D)] ] >;
im_psi := Image(psi);
ker_psi := Kernel(psi);
coker_psi := quo < Sylow3 | im_psi>;
print <Invariants(im_psi), #im_psi, Invariants(ker_psi), #ker_psi, Invariants(coker_psi), #coker_psi>;


function ComputeActionMatrices(QuotientGrp, QuotMap, G, mG, L, ZL, Cl_map)
    d := Ngens(QuotientGrp);
    F3 := GaloisField(3);
    matrices := [];
    DomainGrp := Domain(QuotMap); // Can be M_L or Sylow3
    
    for g in Generators(G) do
        aut := mG(g);
        M_g := ZeroMatrix(F3, d, d);
        for i in [1..d] do
            s := QuotientGrp.i @@ QuotMap;
            I := Cl_map(Cl ! s);
            J := ideal< ZL | [ aut(L ! x) : x in Generators(I) ] >;
            
            // Pullback to Cl, coerce to the subgroup domain, then map to quotient
            I_cl := DomainGrp ! (J @@ Cl_map); 
            coords := Eltseq(QuotMap(I_cl));
            for j in [1..d] do
                M_g[i, j] := F3 ! coords[j];
            end for;
        end for;
        Append(~matrices, M_g);
    end for;
    return matrices;
end function;

M_L := im_psi;
Cl_3 := sub< Cl | [ 3*x : x in Generators(Cl) ] >;

B, m_B := quo< M_L | M_L meet Cl_3 >;
matrices_B := ComputeActionMatrices(B, m_B, G, mG, L, ZL, m);

V, q_V := quo< Sylow3 | Sylow3 meet Cl_3 >;
matrices_V := ComputeActionMatrices(V, q_V, G, mG, L, ZL, m);

// 0 is trivial, 1 is sign, 2 is W, 3 is W dual, 4 is F3, 5 is F3 dual
function IdentifyIndecomposable(M)
    dim := Dimension(M);
    G := Group(M);
    
    t := [g : g in G | Order(g) eq 2][1];
    
    S := Socle(M);
    is_socle_trivial := (S.1 * t) eq S.1;
    
    if dim eq 1 then
        return is_socle_trivial select 0 else 1;
    elif dim eq 2 then
        return is_socle_trivial select 2 else  3;
    elif dim eq 3 then
        return is_socle_trivial select 4 else 5;
    else
        return -1;
    end if;
end function;

dB := Ngens(B);
dV := Ngens(V);
T := ZeroMatrix(GaloisField(3), dB, dV);

for i in [1..dB] do
    v_elem := q_V( Cl ! (B.i @@ m_B) );
    coords := Eltseq(v_elem);
    for j in [1..dV] do
        T[i,j] := coords[j];
    end for;
end for;

ModM_L := GModule(G, matrices_B);
ModC_L := GModule(G, matrices_V);
iota := hom< ModM_L -> ModC_L | T >;
E, proj := quo< ModC_L | Image(iota)>;

// Extract properties
get_decomps := func< M | Sort([ IdentifyIndecomposable(m) : m in Decomposition(M) ]) >;
get_comps := func< M | [ Sort([ IdentifyIndecomposable(a) : a in Decomposition(m) ]) : m in CompositionSeries(M) ] >;

print get_decomps(ModM_L);
print get_decomps(ModC_L);
print get_decomps(E);

print get_comps(ModM_L);
print get_comps(ModC_L);
if #E gt 1 then 
    print get_comps(E);
else
    print [];
end if;


//exit;