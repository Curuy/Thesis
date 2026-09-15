SetClassGroupBounds("GRH");
//SetNthreads(4);
AttachSpec("proj");

X := StringToInteger(X);
is_real := (is_real eq "1");
Q<x> := PolynomialRing(Rationals());

coeff := GenerateCubics(X, is_real);

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
    I_K1_1 := ideal< Zk1 | [ LoverK1 | g : g in Generators(I) ] >;
    I_K1_2 := ideal< Zk1 | [ LoverK1 | g : g in Generators(I_sig2) ] >;
    I_K1_3 := ideal< Zk1 | [ LoverK1 | g : g in Generators(I_sig1) ] >;
    // Take the norm
    return <Norm(I_K1_1), Norm(I_K1_2), Norm(I_K1_3), N_F_tau>;
end function;

Sylow3 := SylowSubgroup(Cl, 3);
ideal_gens := [ m(Sylow3.i) : i in [1..Ngens(Sylow3)] ];

D, inc, proj := DirectSum([Sylow_Cl1, Sylow_Cl1, Sylow_Cl1, Sylow_Cl4]);
print <Invariants(Sylow3), #Sylow3, Invariants(D), #D>;

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
psi := hom< D -> Sylow3 | [ Psi(D.i) : i in [1..Ngens(D)] ] >;
im_psi := Image(psi);
ker_psi := Kernel(psi);
coker_psi := quo < Sylow3 | im_psi>;
print <Invariants(im_psi), #im_psi, Invariants(ker_psi), #ker_psi, Invariants(coker_psi), #coker_psi>;


function ComputeActionMatrices(Q, mQ, G, mG, L, ZL, Cl_map)
    // Given Q with d generators, see how the generatos of G act on Q.
    // These objects are class groups however since Q is a quotient, X / 3*Cl, to get their view in the proper domain must take their preimage in mQ so that it is compatible with Cl map
    d := Ngens(Q);
    F3 := GaloisField(3);
    matrices := [];
    DomainQ := Domain(mQ); // Domain of the map.
    
    for g in Generators(G) do
        aut := mG(g);
        // Want to construct a Mat_F3(d) matrix
        M_g := ZeroMatrix(F3, d, d);
        for i in [1..d] do
            // Get the pre image of Q.i in mQ.
            s := Q.i @@ mQ;
            I := Cl_map(Cl ! s);
            // apply the action. First view the generators of x in L then we want an ideal so create the ideal over Zl.
            J := ideal< ZL | [ aut(L ! x) : x in Generators(I) ] >;
            
            // pullback to cl and construct the map.
            I_cl := DomainQ ! (J @@ Cl_map); 
            coords := Eltseq(mQ(I_cl));
            // we want a right action so,
            // sigma(q) = q sigma for a q in Q
            // This corresponds to,
            // (q.1, q.2, q.3) sigma = (sigma(q.1) \\ sigma(q.2) \\ sigma(q.3))
            // I.e. each row corresponds to how sigma acts on the ith generator.
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
print B, m_B;
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
// We now want to create an injective map from B to V. 
// We just need to know how the generators of B are "viewed" in V
// To do that, we just view each generator of B in V

for i in [1..dB] do
    v_elem := q_V( Cl ! (B.i @@ m_B) );
    coords := Eltseq(v_elem);
    for j in [1..dV] do
        T[i,j] := coords[j];
    end for;
end for;

print matrices_B;
print matrices_V;

print "";

ModM_L := GModule(G, matrices_B);
ModC_L := GModule(G, matrices_V);
iota := hom< ModM_L -> ModC_L | T >;
E, proj := quo< ModC_L | Image(iota)>;

print matrices_B[1] * T;
print "";
print T * matrices_V[1];

//print ActionGenerators(Image(iota));
exit;

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

exit;