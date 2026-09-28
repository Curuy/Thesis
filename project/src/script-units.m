SetClassGroupBounds("GRH");
//SetNthreads(4);
// Input is a string of coeff

if not assigned coeff then
    error "Supply coeff:=[a0,a1,a2,a3] on the command line.";
end if;
if Type(coeff) eq MonStgElt then
    coeff := [ StringToInteger(t) : t in Split(coeff, "[], \t\r\n") ];
end if;
assert #coeff eq 4;
assert coeff[4] ne 0;
input_coeff := [ Integers() | c : c in coeff ];

Q<x> := PolynomialRing(Rationals());

g := Q!input_coeff;

L := SplittingField(g);
f_def := DefiningPolynomial(L);
L_coeff := Coefficients(f_def);
f_def := f_def / LeadingCoefficient(f_def);

L := OptimisedRepresentation(NumberField(f_def));

ZL := MaximalOrder(L);

UL, umap_L := UnitGroup(ZL : GRH := true);


subs := Subfields(L);
D_L := Discriminant(ZL);
r_L, c_L := Signature(L);
label := Sprintf("6.%o.%o.1", r_L, Abs(D_L));

coeff := Coefficients(g); 
sextic_info := <label, L_coeff, c_L, D_L, Factorisation(D_L), coeff, AbelianInvariants(UL)>;
print sextic_info;


K1_data := [ s : s in subs | Degree(s[1]) eq 3 ][1];
F_data  := [ s : s in subs | Degree(s[1]) eq 2 ][1];




K1, mK1 := Explode(K1_data);
F, mF := Explode(F_data);

D_k := Discriminant(MaximalOrder(K1));
D_f := Discriminant(MaximalOrder(F));

autsF := Automorphisms(F);
tau := [ a : a in autsF | a(F.1) ne F.1 ][1];

G, Aut, mG := AutomorphismGroup(L);
sigma := mG([ g : g in G | Order(g) eq 3 ][1]);


UK, umap_K := UnitGroup(K1 : GRH := true);
UF, umap_F := UnitGroup(F : GRH := true);


r_k, c_k := Signature(K1);
label_cubic := Sprintf("3.%o.%o.1", r_k, Abs(D_k));
cubic_info := <label_cubic, coeff, c_k, AbelianInvariants(UK)>;
print cubic_info;

f_quad := DefiningPolynomial(F);
r_f, c_f := Signature(F);
label_quad := Sprintf("2.%o.%o.1", r_f, Abs(D_f));
quad_info := <label_quad, Coefficients(f_quad), c_f, AbelianInvariants(UF)>;
print quad_info;


// Just need to define Psi.


D, inc, proj := DirectSum([UK, UK, UK, UF]);

// For every generator of D, just want to define the map that maps it to UL.

Psi := function(x)
    // Given an element in D, we want to view it in the field L. 
    c_1 := proj[1](x);
    c_2 := proj[2](x);
    c_3 := proj[3](x);
    c_4 := proj[4](x);
    
    
    u_1 := K1 ! umap_K(c_1);
    u_2 := K1 ! umap_K(c_2);
    u_3 := K1 ! umap_K(c_3);
    u_4 := F  ! umap_F(c_4);

    v_1 := mK1(u_1);
    v_2 := sigma(mK1(u_2));
    v_3 := sigma(sigma(mK1(u_3)));
    v_4 := mF(u_4);

    v := v_1 * v_2 * v_3 / v_4;

    return (ZL ! v) @@ umap_L;

end function;

psi := hom< D -> UL | [ Psi(D.i) : i in [1..Ngens(D)] ] >;

U0 := Image(psi);




function ComputeActionMatrices(Q, mQ, G, mG, L, ZL, UL, umap_L)
    // Given Q with d generators, see how the generatos of G act on Q.
    // These objects are class groups however since Q is a quotient, X / 3*UL, to get their view in the proper domain must take their preimage in mQ so that it is compatible with Cl map
    d := Ngens(Q);
    F3 := GaloisField(3);
    matrices := [];
    DomainQ := Domain(mQ); // Domain of the map.
    
    for g in Generators(G) do
        aut := mG(g);
        // Want to construct a Mat_F3(d) matrix
        M_g := ZeroMatrix(F3, d, d);
        for i in [1..d] do
            v := UL ! (Q.i @@ mQ);
            
            u := L ! umap_L(v);

            u_image := ZL ! aut(u);

            v_in_UL := u_image @@ umap_L;

            // Pullback to B
            v_B := mQ(DomainQ ! v_in_UL);

            
            coords := Eltseq(v_B);
        
            for j in [1..d] do
                M_g[i, j] := F3 ! coords[j];
            end for;
        end for;
        Append(~matrices, M_g);
    end for;
    return matrices;
end function;


UL_3 := sub< UL | [ 3*x : x in Generators(UL) ] >;

B, m_B := quo< U0 | U0 meet UL_3 >;
matrices_B := ComputeActionMatrices(B, m_B, G, mG, L, ZL, UL, umap_L);


V, m_V := quo< UL | UL meet UL_3 >;
matrices_V := ComputeActionMatrices(V, m_V, G, mG, L, ZL, UL, umap_L);

ModU0 := GModule(G, matrices_B);
ModUL := GModule(G, matrices_V);

dB := Ngens(B);
dV := Ngens(V);
T := ZeroMatrix(GaloisField(3), dB, dV);

for i in [1..dB] do
    v_elem := m_V( UL ! (B.i @@ m_B) );
    coords := Eltseq(v_elem);
    for j in [1..dV] do
        T[i,j] := coords[j];
    end for;
end for;

iota := hom< ModU0 -> ModUL | T >;


UE, proj := quo< ModUL | Image(iota)>;

// Must be a SES so,

// Image(iota) == ker(proj) and iota must be injective and proj must be surjective.

// Check injective
assert 1 eq #Kernel(iota);

// Check they are equal
assert Image(iota) eq Kernel(proj);

// Check proj is surjective.

assert 1 eq #quo<UE | Image(proj)>;


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


get_decomps := func< M | Sort([ IdentifyIndecomposable(m) : m in Decomposition(M) ]) >;
get_comps := func< M | [ Sort([ IdentifyIndecomposable(a) : a in Decomposition(m) ]) : m in CompositionSeries(M) ] >;


print get_decomps(ModU0);
print get_decomps(ModUL);
print get_decomps(UE);

print get_comps(ModU0);
print get_comps(ModUL);
if #UE gt 1 then 
    print get_comps(UE);
else
    print [];
end if;




