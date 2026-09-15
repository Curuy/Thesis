intrinsic CheckMaximalAtP(a::RngIntElt, b::RngIntElt, c::RngIntElt, d::RngIntElt, p::RngIntElt) -> BoolElt
    {Check if a given cubic ax^3 + bx^2 + cx + d is maximal at p}

    if a eq 0 or p lt 2 then return false; end if;

    P<x, y> := PolynomialRing(Integers(), 2);
    f := a*x^3 + b*x^2*y + c*x*y^2 + d*y^3;

    if (a mod p eq 0) and (b mod p eq 0) and (c mod p eq 0) and (d mod p eq 0) then
        return false;
    end if;
    for r in [0..p-1] do
        
        fp := Evaluate(f, [r*x - y, x]);
        if (a mod p^2 eq 0) and (b mod p eq 0) then
            return false;
        end if;
        
        ap := MonomialCoefficient(fp, x^3);
        bp := MonomialCoefficient(fp, x^2 * y);
        
        if (ap mod p^2 eq 0) and (bp mod p eq 0) then
            return false;
        end if;
    end for;
    
    return true;

end intrinsic;

intrinsic GenerateCubics(X::RngIntElt, is_real:: BoolElt) -> SeqEnum
    {Generate a random irreducible, reduced and maximal non-cyclic cubic in the range a,b in [0, X] and c, d in [-X,X]}
    Zx<x> := PolynomialRing(Integers());
    // Runs until it generates a valid cubic.
    while true do
        a := Random(X);
        b := Random(X);
        c := -X + Random(2*X);
        d := -X + Random(2*X);
        if a eq 0 then continue; end if;
    
        // Hessian
        P := b^2-3*a*c;
        Q := b*c-9*a*d;
        R := c^2 - 3*b*d;
    
        is_reduced := false;
    
        // Check if the hessian is reduced
    
        if not(Abs(Q) le P and P le R) then continue; end if;
    
        if  (b gt 0 or d lt 0) and
            (Q ne 0 or d lt 0) and
            (P ne Q or b lt Abs(3 * a - b)) and
            (P ne R or (a le Abs(d) and
                (a ne Abs(d) or b lt Abs(c)))) then
                is_reduced := true;
        end if;
    
        if not is_reduced then continue; end if;

        Disc := b^2*c^2 - 4*a*c^3 - 4*b^3*d - 27*a^2*d^2 + 18*a*b*c*d;

        // If is_real is true, Disc must be positive
        if is_real and Disc lt 0 then continue; end if;
        // If is_real is false (complex), Disc must be negative
        if not is_real and Disc gt 0 then continue; end if;

        // Check if the polynomial is primitive as its cheaper then irreduciblity check
        if GCD([a, b, c, d]) gt 1 then
            continue;
        end if;

        // We want the cubic to be non cyclic so its discriminant not be square.
        if IsSquare(Disc) then continue; end if;

        f := a*x^3 + b*x^2 +c*x + d;
        if not IsIrreducible(f) then continue; end if;
        factors := Factorisation(Disc);
        is_maximal := true;
        
        for fact in factors do
            p := fact[1];
            pow := fact[2];
            if pow ge 2 then
                if not CheckMaximalAtP(a, b, c, d, p) then
                    is_maximal := false;
                    break;
                end if;
            end if;
        end for;
        
        if not is_maximal then continue; end if;
        

        // Valid cubic
        return [a,b,c,d];

    end while;
    
    
end intrinsic;


intrinsic GenerateCubics(X::RngIntElt, is_real:: BoolElt) -> SeqEnum
    {Generate a random irreducible, reduced and maximal non-cyclic cubic in the range a,b in [0, X] and c, d in [-X,X]}
    Zx<x> := PolynomialRing(Integers());
    
    // Runs until it generates a valid cubic.
    while true do
        a := Random(X);
        b := Random(X);
        c := -X + Random(2*X);
        d := -X + Random(2*X);
        if a eq 0 then continue; end if;
    
        // Hessian
        P := b^2 - 3*a*c;
        Q := b*c - 9*a*d;
        R := c^2 - 3*b*d;
    
        is_reduced := false;
    
        // Check if the hessian is reduced
        if not(Abs(Q) le P and P le R) then continue; end if;
    
        if  (b gt 0 or d lt 0) and
            (Q ne 0 or d lt 0) and
            (P ne Q or b lt Abs(3 * a - b)) and
            (P ne R or (a le Abs(d) and
                (a ne Abs(d) or b lt Abs(c)))) then
                is_reduced := true;
        end if;
    
        if not is_reduced then continue; end if;

        Disc := b^2*c^2 - 4*a*c^3 - 4*b^3*d - 27*a^2*d^2 + 18*a*b*c*d;

        // If is_real is true, Disc must be positive
        if is_real and Disc lt 0 then continue; end if;
        // If is_real is false (complex), Disc must be negative
        if not is_real and Disc gt 0 then continue; end if;

        // Check if the polynomial is primitive as it's cheaper than irreducibility check
        if GCD([a, b, c, d]) gt 1 then
            continue;
        end if;

        // We want the cubic to be non-cyclic so its discriminant cannot be square.
        if IsSquare(Disc) then continue; end if;

        f := a*x^3 + b*x^2 + c*x + d;
        if not IsIrreducible(f) then continue; end if;
        
        factors := Factorisation(Disc);
        is_maximal := true;
        
        for fact in factors do
            p := fact[1];
            pow := fact[2];
            if pow ge 2 then
                // If it fails maximality at p, mark false and break out of the FOR loop
                if not CheckMaximalAtP(a, b, c, d, p) then
                    is_maximal := false;
                    break;
                end if;
            end if;
        end for;
        
        // If it was marked non-maximal, restart the WHILE loop
        if not is_maximal then continue; end if;

        // Valid cubic
        return [a, b, c, d];

    end while;
    
end intrinsic;