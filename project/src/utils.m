intrinsic CheckMaximalAtP(a::RngIntElt, b::RngIntElt,
    c::RngIntElt, d::RngIntElt, p::RngIntElt) -> BoolElt
    {Test p-maximality of the cubic ring associated to (a,b,c,d).}

    require IsPrime(p) : "p must be prime.";

    if (a mod p eq 0) and (b mod p eq 0) and
       (c mod p eq 0) and (d mod p eq 0) then
        return false;
    end if;

    // Check the projective point at infinity.
    if (a mod (p^2) eq 0) and (b mod p eq 0) then
        return false;
    end if;

    // Find repeated finite roots modulo p.
    Fp<t> := PolynomialRing(GF(p));
    fp := a*t^3 + b*t^2 + c*t + d;
    h := GCD(fp, Derivative(fp));

    if Degree(h) gt 0 then
        for rt in Roots(h) do
            r := Integers()!rt[1];

            if ((a*r + b)*r^2 + c*r + d) mod (p^2) eq 0 then
                return false;
            end if;
        end for;
    end if;

    return true;
end intrinsic;


intrinsic GenerateCubics(X::RngIntElt, is_real::BoolElt) -> SeqEnum
    {Uniformly sample reduced, irreducible, maximal, non-cyclic
     cubic forms with 1 <= a <= X, 0 <= b <= X,
     and -X <= c,d <= X.}

    require X ge 1 : "X must be at least 1.";
    require (not is_real) or (X ge 3) :
        "The real non-cyclic reduced box is empty for X < 3.";

    Zx<x> := PolynomialRing(Integers());

    while true do
        // Uniform proposal in the coefficient box.
        a := 1 + Random(X - 1);
        b := Random(X);
        c := -X + Random(2*X);
        d := -X + Random(2*X);

        Disc := b^2*c^2 - 4*a*c^3 - 4*b^3*d
                - 27*a^2*d^2 + 18*a*b*c*d;

        if is_real then
            if Disc le 0 then continue; end if;

            // Your real reduction criterion.
            P := b^2 - 3*a*c;
            Q := b*c - 9*a*d;
            R := c^2 - 3*b*d;

            if not (Abs(Q) le P and P le R) then
                continue;
            end if;

            if not (
                (b gt 0 or d lt 0) and
                (Q ne 0 or d lt 0) and
                (P ne Q or b lt Abs(3*a - b)) and
                (P ne R or (a le Abs(d) and
                    (a ne Abs(d) or b lt Abs(c))))
            ) then
                continue;
            end if;

            // Exclude cyclic cubic fields.
            if IsSquare(Disc) then continue; end if;

        else
            if Disc ge 0 then continue; end if;

            // Exact complex reduction: Belabas, Lemma 4.2.
            if b eq 0 and d le 0 then continue; end if;

            if d^2 - a^2 + a*c - b*d le 0 then
                continue;
            end if;

            if a*d - b*c le -(a-b)^2 - a*c then
                continue;
            end if;

            if a*d - b*c ge (a+b)^2 + a*c then
                continue;
            end if;

            // Negative discriminant is automatically non-square.
        end if;

        if GCD([a, b, c, d]) ne 1 then continue; end if;

        f := a*x^3 + b*x^2 + c*x + d;
        if not IsIrreducible(f) then continue; end if;

        is_maximal := true;

        for fact in Factorisation(Abs(Disc)) do
            if fact[2] ge 2 then
                if not CheckMaximalAtP(a, b, c, d, fact[1]) then
                    is_maximal := false;
                    break;
                end if;
            end if;
        end for;

        if is_maximal then
            return [a, b, c, d];
        end if;
    end while;
end intrinsic;