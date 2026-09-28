// Usage: magma -b 'coeff:=[21,54,-30,-73]' unit_ses.m
// Ascending coefficients: [a0,a1,a2,a3] means a0+a1*x+a2*x^2+a3*x^3.
// Module codes: 0=T, 1=S, 2=W (socle T), 3=W* (socle S),
//               4=P_T, 5=P_S. Empty list means the zero module.
SetClassGroupBounds("GRH");
SetColumns(0);

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
g := Q ! input_coeff;
assert IsIrreducible(g);
assert not IsSquare(Discriminant(g));

L := SplittingField(g);
assert Degree(L) eq 6;
f_def := DefiningPolynomial(L);
f_def := f_def / LeadingCoefficient(f_def);
L := OptimisedRepresentation(NumberField(f_def));
ZL := MaximalOrder(L);
UL, umap_L := UnitGroup(ZL : GRH := true);

subs := Subfields(L);
K1_data := [ s : s in subs | Degree(s[1]) eq 3 ][1];
F_data := [ s : s in subs | Degree(s[1]) eq 2 ][1];
K1, mK1 := Explode(K1_data);
F, mF := Explode(F_data);
G, Aut, mG := AutomorphismGroup(L);
assert #G eq 6;
sigma := mG([ h : h in G | Order(h) eq 3 ][1]);
UK, umap_K := UnitGroup(K1 : GRH := true);
UF, umap_F := UnitGroup(F : GRH := true);

D, inc, projections := DirectSum([UK, UK, UK, UF]);
Psi := function(z)
    u1 := K1 ! umap_K(projections[1](z));
    u2 := K1 ! umap_K(projections[2](z));
    u3 := K1 ! umap_K(projections[3](z));
    u4 := F ! umap_F(projections[4](z));
    v := mK1(u1) * sigma(mK1(u2)) * sigma(sigma(mK1(u3))) / mF(u4);
    return (ZL ! v) @@ umap_L;
end function;
psi := hom< D -> UL | [Psi(D.i) : i in [1..Ngens(D)]] >;
U0 := Image(psi);
UL_3 := sub< UL | [3*u : u in Generators(UL)] >;
assert UL_3 subset U0;
B, m_B := quo< U0 | U0 meet UL_3 >;
V, m_V := quo< UL | UL_3 >;

function ComputeActionMatrices(A, q, G, mG, L, ZL, UL, umap_L)
    d := Ngens(A);
    assert #A eq 3^d;
    matrices := [];
    for k in [1..Ngens(G)] do
        aut := mG(G.k);
        R := ZeroMatrix(GF(3), d, d);
        for i in [1..d] do
            v := UL ! (A.i @@ q);
            u := L ! umap_L(v);
            v_image := (ZL ! aut(u)) @@ umap_L;
            coords := Eltseq(q(Domain(q) ! v_image));
            for j in [1..d] do
                R[i,j] := GF(3) ! coords[j];
            end for;
        end for;
        assert Rank(R) eq d;
        Append(~matrices, R);
    end for;
    return matrices;
end function;

matrices_B := ComputeActionMatrices(B, m_B, G, mG, L, ZL, UL, umap_L);
matrices_V := ComputeActionMatrices(V, m_V, G, mG, L, ZL, UL, umap_L);
ModU0 := GModule(G, matrices_B);
ModUL := GModule(G, matrices_V);

dB := Ngens(B);
dV := Ngens(V);
Tmat := ZeroMatrix(GF(3), dB, dV);
for i in [1..dB] do
    coords := Eltseq(m_V(UL ! (B.i @@ m_B)));
    for j in [1..dV] do
        Tmat[i,j] := GF(3) ! coords[j];
    end for;
end for;
assert Rank(Tmat) eq dB;
for k in [1..Ngens(G)] do
    assert matrices_B[k] * Tmat eq Tmat * matrices_V[k];
end for;
iota := hom< ModU0 -> ModUL | Tmat >;
UE, quotient_map := quo< ModUL | Image(iota) >;
assert Dimension(Kernel(iota)) eq 0;
assert Image(iota) eq Kernel(quotient_map);
assert Dimension(Image(quotient_map)) eq Dimension(UE);
Q_units := quo< UL | U0 >;
assert #Q_units eq 3^Dimension(UE);

function IdentifyIndecomposable(M)
    dim := Dimension(M);
    assert dim in {1,2,3};
    H := Group(M);
    t := [ h : h in H | Order(h) eq 2 ][1];
    soc := Socle(M);
    assert Dimension(soc) eq 1;
    trivial_socle := (soc.1 * t) eq soc.1;
    if dim eq 1 then
        return trivial_socle select 0 else 1;
    elif dim eq 2 then
        return trivial_socle select 2 else 3;
    else
        return trivial_socle select 4 else 5;
    end if;
end function;

function GetDecomps(M)
    if Dimension(M) eq 0 then
        return [ Integers() | ];
    end if;
    return Sort([IdentifyIndecomposable(N) : N in Decomposition(M)]);
end function;

// Preserve the original script's meaning: decompositions of the modules
// in a composition series, NOT merely a list of simple composition factors.
function GetComps(M)
    if Dimension(M) eq 0 then
        return [];
    end if;
    return [GetDecomps(N) : N in CompositionSeries(M)];
end function;

D_L := Discriminant(ZL);
D_k := Discriminant(MaximalOrder(K1));
D_f := Discriminant(MaximalOrder(F));
r_L, c_L := Signature(L);
r_k, c_k := Signature(K1);
r_f, c_f := Signature(F);
if r_L eq 6 then
    assert dV eq 5;
    assert #Q_units in {1,3,9};
end if;
assert forall{c : c in GetDecomps(UE) | c eq 1};

// Labelled records let Python ignore banners and diagnostic output.
print "@@UNIT_SES_BEGIN";
printf "schema_version=%o\n", 1;
printf "input_coeff=%o\n", input_coeff;
printf "sextic_label=6.%o.%o.1\n", r_L, Abs(D_L);
printf "sextic_coeff=%o\n", Coefficients(DefiningPolynomial(L));
printf "sextic_sig=%o\n", [r_L,c_L];
printf "sextic_disc=%o\n", D_L;
printf "sextic_disc_fact=%o\n", Factorisation(D_L);
printf "sextic_unit_inv=%o\n", AbelianInvariants(UL);
printf "cubic_label=3.%o.%o.1\n", r_k, Abs(D_k);
printf "cubic_coeff=%o\n", input_coeff;
printf "cubic_sig=%o\n", [r_k,c_k];
printf "cubic_disc=%o\n", D_k;
printf "cubic_unit_inv=%o\n", AbelianInvariants(UK);
printf "quad_label=2.%o.%o.1\n", r_f, Abs(D_f);
printf "quad_coeff=%o\n", Coefficients(DefiningPolynomial(F));
printf "quad_sig=%o\n", [r_f,c_f];
printf "quad_disc=%o\n", D_f;
printf "quad_unit_inv=%o\n", AbelianInvariants(UF);
printf "unit_index=%o\n", #Q_units;
printf "m_u_dim=%o\n", Dimension(ModU0);
printf "u_l_dim=%o\n", Dimension(ModUL);
printf "e_u_dim=%o\n", Dimension(UE);
printf "m_u_decomp=%o\n", GetDecomps(ModU0);
printf "u_l_decomp=%o\n", GetDecomps(ModUL);
printf "e_u_decomp=%o\n", GetDecomps(UE);
printf "m_u_comp=%o\n", GetComps(ModU0);
printf "u_l_comp=%o\n", GetComps(ModUL);
printf "e_u_comp=%o\n", GetComps(UE);
print "@@UNIT_SES_END";
quit;
