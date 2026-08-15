
AttachSpec("spec");
Z := "0";
X := "10^14";
SEED := 42;
//N := 10;
SetSeed(SEED);
cmd := Sprintf("./sextic %o %o 24", Z, X);
output := Pipe(cmd, "");

ParseTuple := function(str)
    parts := Split(str, ",");

    seq := [ StringToInteger(Join(Split(s, " "), "")) : s in parts ];

    return < seq[1], seq[2], seq[3] >;
    end function;

ParseList := function(str)
    parts := Split(str, ",");
    return [ StringToInteger(Join(Split(s, " "), "")) : s in parts ];
    end function;

lines := Split(output, "\n");
parsed_data := [];

for line in lines do

    if Index(line, ":") gt 0 and Index(line, "(") gt 0 and Index(line, "[") gt 0 then
        try

            colon_idx := Index(line, ":");
            int_str := Substring(line, 1, colon_idx - 1);
            int_str := Join(Split(int_str, " "), ""); // Strip any spaces
            first_elem := StringToInteger(int_str);


            open_paren := Index(line, "(");
            close_paren := Index(line, ")");
            tuple_str := Substring(line, open_paren + 1, close_paren - open_paren - 1);
            second_elem := ParseTuple(tuple_str);


            open_bracket := Index(line, "[");
            close_bracket := Index(line, "]");
            list_str := Substring(line, open_bracket + 1, close_bracket - open_bracket - 1);
            third_elem := ParseList(list_str);
            third_elem := Reverse(third_elem);

            Append(~parsed_data, <first_elem, second_elem, third_elem>);

        catch e

            continue;
        end try;
    end if;
end for;

PORT := 10000;
NUM_WORKERS := 5;
coeffs := [elem[3] : elem in parsed_data];
// randomly select 25% of the coeffs
N := Floor(0.25 * #coeffs);
random_subset := RandomSubset(Set(coeffs), N);

random_coeffs := (random_subset eq Set(coeffs)) select coeffs else [x : x in random_subset];


socket := Socket(:LocalHost := "0.0.0.0", LocalPort := PORT);

StartDistributedWorkers("worker.m", NUM_WORKERS);

results := DistributedManager(socket, random_coeffs);



file_sextics := Sprintf("sextics-%o-%o-%o-%o.csv", Z, X, SEED, N);
file_maps := Sprintf("sextics-map-%o-%o-%o-%o.csv", Z, X, SEED, N);

F_s := Open(file_sextics, "w");
F_m := Open(file_maps, "w");

fprintf F_s, "label,coeff,deg,disc,prime_div,quad_sub,cubic_sub,class_group_inv,class_number\n";
fprintf F_m, "label,coeff,disc,inv_sylow,size_sylow,inv_direct_sum,size_direct_sum,inv_im_phi,size_im_phi,inv_ker_phi,size_ker_phi,inv_coker_phi,size_coker_phi,inv_im_psi,size_im_psi,inv_ker_psi,size_ker_psi,inv_coker_psi,size_coker_psi\n";


for t in results do
    if #t ne 10 then continue; end if;
    label := t[1];
    sextic := t[2];
    disc := sextic[1];

    syl3 := t[3];
    direct := t[4];
    im_phi := t[5];
    ker_phi := t[6];
    cok_phi := t[7];
    im_psi := t[8];
    ker_psi := t[9];
    cok_psi := t[10];

    fprintf F_s, "\"%o\",\"%o\",\"%o\",\"%o\",\"%o\",\"%o\",\"%o\",\"%o\",%o\n",
        label,
        Sprint(sextic[1]), Sprint(sextic[2]), Sprint(sextic[3]),
        Sprint(sextic[4]), Sprint(sextic[5]), Sprint(sextic[6]), Sprint(sextic[7]), sextic[8];
        


    fprintf F_m, "\"%o\",\"%o\",\"%o\",\"%o\",%o,\"%o\",%o,\"%o\",%o,\"%o\",%o,\"%o\",%o,\"%o\",%o,\"%o\",%o,\"%o\",%o\n",
        label, Sprint(sextic[6]), Sprint(sextic[3]), 
        Sprint(syl3[1]),    syl3[2],
        Sprint(direct[1]),  direct[2],
        Sprint(im_phi[1]),  im_phi[2],
        Sprint(ker_phi[1]), ker_phi[2],
        Sprint(cok_phi[1]), cok_phi[2],
        Sprint(im_psi[1]),  im_psi[2],
        Sprint(ker_psi[1]), ker_psi[2],
        Sprint(cok_psi[1]), cok_psi[2];
end for;

delete F_s;
delete F_m;

delete socket;



quit;
