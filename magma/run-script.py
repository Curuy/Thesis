import os.path
import re
import pandas as pd
import concurrent.futures
import ast


import time
import os
from subprocess import Popen, PIPE, TimeoutExpired
from functools import partial


PROGRAM = "magma"
PROJECT_DIR = "tests"
PROGRAM_RUN_CMD = "magma -b {config} {PROGRAM_NAME}"


Popen = partial(Popen, stdin=PIPE, stdout=PIPE)

timeout_sec = None

def magma_run(config, program_name="script.m"):
    cmd = ["magma", "-b"] + config.split() + [program_name]
    magma_instance = Popen(cmd)
    try:
        out, err = magma_instance.communicate(timeout=timeout_sec)
    except TimeoutExpired:
        magma_instance.kill()
        
        out, err = magma_instance.communicate()
        
        decoded_out = out.decode('utf-8') if out else ""
        return decoded_out, f"Error: Process timed out after {timeout_sec} seconds."

    decoded_out = out.decode('utf-8') if out else ""
    return decoded_out, err




list_to_str = lambda x: ','.join(map(str, x))

def parse_magma_output(output: str):
    """
    Parses the Magma output format using regular expressions.
    Extracts the sextic info, cubic/quad info, and the group flag invariant tuples separately.
    Handles multiple factorisation tuples and multiline outputs natively.
    """
    text = re.sub(r'\\\s+', '', output)
    text = text.replace('\n', ' ').replace('\r', ' ')
    
    INT = r"-?\d+"
    LABEL = r"\"?\d+\.\d+\.\d+\.\d+\"?"
    LIST_INT = rf"\[\s*(?:{INT}(?:\s*,\s*{INT})*)?\s*\]"
    MAGMA_TUP = rf"<\s*{INT}\s*,\s*{INT}\s*>"
    LIST_FACT = rf"\[\s*(?:{MAGMA_TUP}(?:\s*,\s*{MAGMA_TUP})*)?\s*\]"
    
    SEXTIC_PATTERN = rf"<\s*({LABEL})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({INT})\s*,\s*({LIST_FACT})\s*,\s*({LIST_INT})\s*,\s*({LIST_INT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({INT})\s*>"
    CUBIC_QUAD_PATTERN = rf"<\s*({LABEL})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({INT})\s*,\s*({LIST_FACT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*,\s*({LIST_INT})\s*,\s*({INT})\s*>"
    
    LIST_THEN_INT = rf"({LIST_INT})\s*,\s*({INT})"
    FLAGS_PATTERN = rf"<\s*{LIST_THEN_INT}(?:\s*,\s*{LIST_THEN_INT})*\s*>"

    def parse_fact(s):
        return ast.literal_eval(s.replace('<', '(').replace('>', ')'))
        
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
        
    cq_matches = list(re.finditer(CUBIC_QUAD_PATTERN, text))
    
    cubic = None
    if len(cq_matches) > 0:
        c_m = cq_matches[0]
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
        
    quad = None
    if len(cq_matches) > 1:
        q_m = cq_matches[1]
        quad = (
            q_m.group(1).replace('"', ''),
            ast.literal_eval(q_m.group(2)),
            int(q_m.group(3)),
            int(q_m.group(4)),
            parse_fact(q_m.group(5)),
            ast.literal_eval(q_m.group(6)),
            int(q_m.group(7)),
            ast.literal_eval(q_m.group(8)),
            int(q_m.group(9))
        )
        
    flags_tuples = []
    for m in re.finditer(FLAGS_PATTERN, text):
        tup = []
        for k in re.finditer(LIST_THEN_INT, m.group()):
            list_str, int_str = k.groups()
            tup.append(ast.literal_eval(list_str))
            tup.append(int(int_str))
        flags_tuples.append(tuple(tup))

    if not sextic or not cubic or not quad or len(flags_tuples) < 3:
        print(output)
        return []

    invs_1 = flags_tuples[0]
    invs_2 = flags_tuples[1]
    invs_3 = flags_tuples[2]
    
    return [sextic, cubic, quad, invs_1, invs_2, invs_3]


def process_field(i):
    """
    worker
    """
    config = "X:=100 is_real:=0"
    
    out, err = magma_run(config)
    return i, parse_magma_output(out)


def write_to_csv(processed_data, filename="subfields-uniform-100-complex-1.csv"):
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
            # Keeping the brackets so it matches the expected `coeff` string format for the Magma worker
            extracted.append(f"[{list_str}]")
    return extracted

def main():

    output_filename = "subfields-uniform-100-complex-1.csv"

    n = 7700
    start_idx = 0
    max_index = n

    MAX_WORKERS = 24

    processed_data = []
    count = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        
        futures = {}
        for i in range(start_idx, max_index):
            #coeff = coeffs[i]
            
            future = executor.submit(process_field, i)
            futures[future] = i

        for future in concurrent.futures.as_completed(futures):
            idx = futures[future]
            try:
                returned_idx, result = future.result()
                
                if len(result) != 6:
                    print(f"Index {returned_idx} returned {len(result)} tuples instead of 6. Skipping.")
                    continue
                
                sextic, cubic, quad, invs_1, invs_2, invs_3 = result
                
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
                    'quad_cl_inv': list(quad[5]),
                    'quad_cl_size': quad[6],
                    'quad_sylow_inv': list(quad[7]),
                    'quad_sylow_size': quad[8],
                    
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
                    'size_coker_psi': invs_3[5]
                })
                
                count += 1
                if count % 100 == 0:
                    print(f"{count}/{max_index}")
                    write_to_csv(processed_data, output_filename)
                    processed_data = []
            except Exception as exc:
                print(f"Error processing index {idx}: {exc}")
                pass
                
    write_to_csv(processed_data, output_filename)


if __name__ == "__main__":
    main()
