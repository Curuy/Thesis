import os.path
import re
import pandas as pd
import concurrent.futures
import ast

from utils.utils import magma_run

list_to_str = lambda x: ','.join(map(str, x))

def parse_magma_output(output: str):
    """
    Parses the Magma output format using regular expressions.
    Extracts the sextic info, cubic/quad info, group flags, and the new decomposition lists.
    """
    text = re.sub(r'\\\s+', '', output)
    text = text.replace('\n', ' ').replace('\r', ' ')
    
    INT = r"-?\d+"
    RAT = r"-?\d+(?:\/\d+)?"  
    LABEL = r"\"?\d+\.\d+\.\d+\.\d+\"?"
    
    LIST_INT = rf"(?<=<|,)\s*\[\s*(?:{INT}(?:\s*,\s*{INT})*)?\s*\]"
    LIST_RAT = rf"\[\s*(?:{RAT}(?:\s*,\s*{RAT})*)?\s*\]"
    
    MAGMA_TUP = rf"<\s*{INT}\s*,\s*{INT}\s*>"
    LIST_FACT = rf"\[\s*(?:{MAGMA_TUP}(?:\s*,\s*{MAGMA_TUP})*)?\s*\]"
    
    SEXTIC_PATTERN = rf"<\s*({LABEL})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({INT})\s*,\s*({LIST_FACT})\s*,\s*({LIST_INT})\s*,\s*({LIST_INT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({INT})\s*>"
    CUBIC_PATTERN = rf"<\s*({LABEL})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({INT})\s*,\s*({LIST_FACT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*>"
    
    LIST_THEN_INT = rf"({LIST_INT})\s*,\s*({INT})"
    FLAGS_PATTERN = rf"<\s*{LIST_THEN_INT}(?:\s*,\s*{LIST_THEN_INT})*\s*>"

    LIST_OF_LISTS = rf"\[\s*(?:{LIST_RAT}(?:\s*,\s*{LIST_RAT})*)?\s*\]"
    IDEAL_TUP = rf"<\s*{LIST_OF_LISTS}\s*,\s*{INT}\s*>"
    LIST_IDEAL_FACT = rf"\[\s*(?:{IDEAL_TUP}(?:\s*,\s*{IDEAL_TUP})*)?\s*\]"
    
    QUAD_PATTERN = rf"<\s*({LABEL})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({INT})\s*,\s*({LIST_FACT})\s*,\s*({LIST_IDEAL_FACT})\s*,\s*({INT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*>"

    def parse_fact(s):
        if not s or s == "[]": return []
        return ast.literal_eval(s.replace('<', '(').replace('>', ')'))

    def parse_ideal_fact(s):
        if not s or s == "[]": return []
        s = s.replace('<', '(').replace('>', ')')
        return eval(s)

    sextic = None
    s_m = re.search(SEXTIC_PATTERN, text)
    if s_m:
        sextic = (
            s_m.group(1).replace('"', ''),            
            ast.literal_eval(s_m.group(2)),          
            int(s_m.group(3)),                        
            int(s_m.group(4)),                        
            parse_fact(s_m.group(5)),                
            ast.literal_eval(s_m.group(6)),          
            ast.literal_eval(s_m.group(7)),          
            ast.literal_eval(s_m.group(8)),          
            int(s_m.group(9)),
            int(s_m.group(10))
        )
        
    c_matches = list(re.finditer(CUBIC_PATTERN, text))
    cubic = None
    if c_matches:
        c_m = c_matches[0]
        cubic = (
            c_m.group(1).replace('"', ''),
            ast.literal_eval(c_m.group(2)),
            int(c_m.group(3)),
            int(c_m.group(4)),
            parse_fact(c_m.group(5)),
            ast.literal_eval(c_m.group(6)),
            int(c_m.group(7)),
            ast.literal_eval(c_m.group(8)),
            int(c_m.group(9))
        )

    q_matches = list(re.finditer(QUAD_PATTERN, text))
    quad = None
    if q_matches:
        q_m = q_matches[0] 
        quad = (
            q_m.group(1).replace('"', ''),
            ast.literal_eval(q_m.group(2)),
            int(q_m.group(3)),
            int(q_m.group(4)),
            parse_fact(q_m.group(5)),
            parse_ideal_fact(q_m.group(6)),
            int(q_m.group(7)),
            ast.literal_eval(q_m.group(8)),
            int(q_m.group(9)),
            ast.literal_eval(q_m.group(10)),
            int(q_m.group(11))
        )

    flags_tuples = []
    for m in re.finditer(FLAGS_PATTERN, text):
        tup = []
        for k in re.finditer(LIST_THEN_INT, m.group()):
            list_str, int_str = k.groups()
            tup.append(ast.literal_eval(list_str))
            tup.append(int(int_str))
        flags_tuples.append(tuple(tup))

    
    # Base pattern for [1, 2, 3] or []
    SIMPLE_LST_INT = r"\[\s*(?:-?\d+(?:\s*,\s*-?\d+)*)?\s*\]"
    # Base pattern for [[1], [2, 3]] or []
    NESTED_LST_INT = rf"\[\s*(?:{SIMPLE_LST_INT}(?:\s*,\s*{SIMPLE_LST_INT})*)?\s*\]"
    
    # The Magma script prints exactly 3 simple lists followed by 3 nested lists
    DECOMP_COMP_PATTERN = rf"({SIMPLE_LST_INT})\s*({SIMPLE_LST_INT})\s*({SIMPLE_LST_INT})\s*({NESTED_LST_INT})\s*({NESTED_LST_INT})\s*({NESTED_LST_INT})"
    
    M_L_decomp = Cl_L_decomp = E_decomp = None
    M_L_comp = Cl_L_comp = E_comp = None
    ml_decomp = []
    cl_decomp = []
    e_decomp = []
    ml_comp = []
    cl_comp = []
    e_comp = []
    # Check the very end of the output first
    dc_match = re.search(DECOMP_COMP_PATTERN + r"\s*$", text)
    if not dc_match:
        # Fallback to the last matched block in the text just in case there's trailing garbage
        dc_matches = list(re.finditer(DECOMP_COMP_PATTERN, text))
        if dc_matches:
            dc_match = dc_matches[-1]
    
    DPF_matches = re.findall(
        r"\bDPF:[ \t]*([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)(?=\s|$)",
        output,
    )
    if not DPF_matches:
        DPF_matches = re.findall(
            r"^[ \t]*([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]*\r?$",
            output,
            flags=re.MULTILINE,
        )
    if not DPF_matches:
        raise ValueError(
            "No DPF invariants found. Use the updated src/DPF.m, "
            "which prints 'DPF: A R C U'."
        )
    A, R, C, U = map(int, DPF_matches[-1])
    if dc_match:
        ml_decomp = ast.literal_eval(dc_match.group(1))
        cl_decomp = ast.literal_eval(dc_match.group(2))
        e_decomp = ast.literal_eval(dc_match.group(3))
        ml_comp = ast.literal_eval(dc_match.group(4))
        cl_comp = ast.literal_eval(dc_match.group(5))
        e_comp = ast.literal_eval(dc_match.group(6))


    

    if not sextic or not cubic or not quad or len(flags_tuples) < 3:
        raise ValueError("Missing field information or class-group map output.")

    invs_1 = flags_tuples[0]
    invs_2 = flags_tuples[1]
    invs_3 = flags_tuples[2]

    if [len(invs_1), len(invs_2), len(invs_3)] != [4, 6, 6]:
        raise ValueError("Unexpected class-group map tuple lengths.")
    if not dc_match and invs_1[1] != 1:
        raise ValueError("Missing decomposition output for a nontrivial 3-class group.")
    
    return [sextic, cubic, quad, invs_1, invs_2, invs_3, ml_decomp, cl_decomp, e_decomp, ml_comp, cl_comp, e_comp, A, R, C, U]


def process_field(i):
    """
    worker
    """
    config = "X:=100 is_real:=1"
    out, err = magma_run(config, program_name="src/DPF.m")
    try:
        return i, parse_magma_output(out)
    except ValueError as exc:
        raise ValueError(
            f"{exc}\nMagma stderr: {err}\nMagma output (last 4000 characters):\n{out[-4000:]}"
        ) from exc


def write_to_csv(processed_data, filename):
    if not processed_data:
        return
        
    results_df = pd.DataFrame(processed_data)
    
    file_exists = os.path.isfile(filename)
    

    results_df.to_csv(
        filename, 
        mode='a', 
        header=not file_exists, 
        index=False
    )

def read_cubics(filename):
    extracted = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or '[' not in line:
                continue
            list_str = line.split('[')[1].replace(']', '')
            extracted.append(f"[{list_str}]")
    return extracted

def main():
    filename = "subfields-uniform-100-real-1-dpf.csv"

    n = 100
    start_idx = 0
    max_index = n

    MAX_WORKERS = 20

    processed_data = []
    count = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {}
        for i in range(start_idx, max_index):
            future = executor.submit(process_field, i)
            futures[future] = i

        for future in concurrent.futures.as_completed(futures):
            idx = futures[future]
            try:
                returned_idx, result = future.result()
                
                # Check for 16 elements now
                if len(result) != 16:
                    print(f"Index {returned_idx} returned {len(result)} entries instead of 16. Skipping.")
                    continue
                
                sextic, cubic, quad, invs_1, invs_2, invs_3, ml_decomp, cl_decomp, e_decomp, ml_comp, cl_comp, e_comp, A, R, C, U = result
                
                processed_data.append({
                    # Sextic mapping
                    'sextic_label': str(sextic[0]),
                    'sextic_coeff': list(sextic[1]),
                    'sextic_sig': sextic[2],
                    'sextic_disc': sextic[3],
                    'sextic_disc_fact': str(list(sextic[4])), 
                    'sextic_cl_inv': list(sextic[7]),
                    'sextic_cl_size': sextic[8],
                    'conductor': sextic[9],
                    
                    # Cubic mapping
                    'cubic_label': str(cubic[0]),
                    'cubic_coeff': list(cubic[1]),
                    'cubic_sig': cubic[2], 
                    'cubic_disc': cubic[3],
                    'cubic_disc_fact': str(list(cubic[4])),
                    'cubic_cl_inv': list(cubic[5]),
                    'cubic_cl_size': cubic[6],
                    'cubic_sylow_inv': list(cubic[7]),
                    'cubic_sylow_size': cubic[8],
                    
                    # Quadratic mapping
                    'quad_label': str(quad[0]),
                    'quad_coeff': list(quad[1]),
                    'quad_sig': quad[2], 
                    'quad_disc': quad[3],
                    'quad_disc_fact': str(list(quad[4])),
                    'quad_rel_disc_fact': str(list(quad[5])),
                    'quad_num_rel_primes': quad[6],
                    'quad_cl_inv': list(quad[7]),
                    'quad_cl_size': quad[8],
                    'quad_sylow_inv': list(quad[9]),
                    'quad_sylow_size': quad[10],
                    
                    # Core flags / maps
                    'inv_sylow': list(invs_1[0]),
                    'size_sylow': invs_1[1],
                    'inv_direct_sum': list(invs_1[2]),
                    'size_direct_sum': invs_1[3],
                    
                    'inv_im_phi': list(invs_2[0]),
                    'size_im_phi': invs_2[1],
                    'inv_ker_phi': list(invs_2[2]),
                    'size_ker_phi': invs_2[3],
                    'inv_coker_phi': list(invs_2[4]),
                    'size_coker_phi': invs_2[5],
                    
                    'inv_im_psi': list(invs_3[0]),
                    'size_im_psi': invs_3[1],
                    'inv_ker_psi': list(invs_3[2]),
                    'size_ker_psi': invs_3[3],
                    'inv_coker_psi': list(invs_3[4]),
                    'size_coker_psi': invs_3[5],

                    # Decompositions
                    'm_l_decomp': ml_decomp,
                    'cl_l_decomp': cl_decomp,
                    'e_decomp': e_decomp,
                    # Compositions
                    'm_l_comp': ml_comp,
                    'cl_l_comp': cl_comp,
                    'e_comp': e_comp,
                    'A': A,
                    'R': R,
                    'C': C,
                    'U': U,

                    
                })
                
                count += 1
                if count % 100 == 0:
                    print(f"{count}/{max_index}")
                    write_to_csv(processed_data, filename)
                
                    processed_data = []
                    
            except Exception as exc:
                print(f"Error processing index {idx}: {exc}")
                pass


        write_to_csv(processed_data, filename)

if __name__ == "__main__":
    main()
