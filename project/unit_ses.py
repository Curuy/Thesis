"""Run: python unit_ses.py coeffs.csv --output units.csv --program src/unit_ses.m

Uses utils.utils.magma_run and the labelled output of the supplied unit_ses.m.
Input coefficients are in ascending order: [a0,a1,a2,a3].
"""

import argparse
import csv
import json
import os.path
import re
import pandas as pd
import concurrent.futures
import ast

from utils.utils import magma_run

list_to_str = lambda x: ','.join(map(str, x))


def parse_magma_output(output: str):
    """Parse the labelled output printed by unit_ses.m."""
    text = re.sub(r'\\\s*\n\s*', '', output)

    INT = r'-?\d+'
    RAT = r'-?\d+(?:/\d+)?'
    LABEL = r'\d+\.\d+\.\d+\.\d+'
    LIST_INT = rf'\[\s*(?:{INT}(?:\s*,\s*{INT})*)?\s*\]'
    LIST_RAT = rf'\[\s*(?:{RAT}(?:\s*,\s*{RAT})*)?\s*\]'
    LIST_FACT = rf'\[\s*(?:<\s*{INT}\s*,\s*{INT}\s*>(?:\s*,\s*<\s*{INT}\s*,\s*{INT}\s*>)*)?\s*\]'
    LIST_OF_LISTS = rf'\[\s*(?:{LIST_INT}(?:\s*,\s*{LIST_INT})*)?\s*\]'

    def match(key, pattern):
        found = re.search(rf'(?m)^\s*{key}=\s*({pattern})', text)
        if found is None:
            raise ValueError(f'Missing Magma output field: {key}')
        return found.group(1)

    def parse_list(key, pattern=LIST_INT):
        value = match(key, pattern).replace('<', '(').replace('>', ')')
        # Preserve nonintegral polynomial coefficients exactly as strings.
        value = re.sub(r'(-?\d+/\d+)', r'"\1"', value)
        return ast.literal_eval(value)

    result = {'input_coeff': parse_list('input_coeff')}

    for field in ('sextic', 'cubic', 'quad'):
        result[f'{field}_label'] = match(f'{field}_label', LABEL)
        result[f'{field}_coeff'] = parse_list(f'{field}_coeff', LIST_RAT)
        result[f'{field}_sig'] = parse_list(f'{field}_sig')
        result[f'{field}_disc'] = int(match(f'{field}_disc', INT))
        result[f'{field}_unit_inv'] = parse_list(f'{field}_unit_inv')

    result['sextic_disc_fact'] = parse_list('sextic_disc_fact', LIST_FACT)
    result['unit_index'] = int(match('unit_index', INT))

    for module in ('m_u', 'u_l', 'e_u'):
        result[f'{module}_dim'] = int(match(f'{module}_dim', INT))
        result[f'{module}_decomp'] = parse_list(f'{module}_decomp')
        result[f'{module}_comp'] = parse_list(f'{module}_comp', LIST_OF_LISTS)

    return result


def process_field(i, coeff, program_name):
    """Run one supplied cubic through Magma."""
    config = f'coeff:=[{list_to_str(coeff)}]'
    out, err = magma_run(config, program_name=program_name)
    if err:
        raise RuntimeError(err)
    try:
        result = parse_magma_output(out)
    except (ValueError, SyntaxError) as exc:
        raise RuntimeError(f'{exc}\n\nMagma output:\n{out}') from exc
    return i, result


def write_to_csv(processed_data, filename):
    if not processed_data:
        return

    results_df = pd.DataFrame(processed_data)
    file_exists = os.path.isfile(filename) and os.path.getsize(filename) > 0
    mode = 'a'
    if file_exists:
        with open(filename, newline='') as f:
            reader = csv.reader(f)
            header = next(reader, [])
            has_rows = next(reader, None) is not None
        if header != list(results_df.columns):
            if has_rows:
                raise ValueError('Existing CSV uses different columns; choose a new --output filename')
            # A previous failed run may have left only an obsolete header.
            mode, file_exists = 'w', False
    results_df.to_csv(filename, mode=mode, header=not file_exists, index=False)


def read_cubics(filename):
    """Read one [a0,a1,a2,a3] coefficient list per line."""
    extracted = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.split('#', 1)[0].strip()
            if not line or '[' not in line:
                continue
            list_str = re.search(r'\[[^\[\]]*\]', line).group()
            extracted.append(ast.literal_eval(list_str))
    return extracted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', default='cubics.txt')
    parser.add_argument('--output', default='unit-ses.csv')
    parser.add_argument('--program', default='src/unit_ses.m')
    parser.add_argument('--workers', type=int, default=48)
    args = parser.parse_args()

    filename = args.output
    cubics = read_cubics(args.input)

    start_idx = 0
    max_index = len(cubics)
    MAX_WORKERS = args.workers

    processed_data = []
    count = 0
    failed = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {}
        for i in range(start_idx, max_index):
            future = executor.submit(process_field, i, cubics[i], args.program)
            futures[future] = i

        for future in concurrent.futures.as_completed(futures):
            idx = futures[future]
            try:
                returned_idx, result = future.result()
            except Exception as exc:
                failed += 1
                with open(filename + '.errors.jsonl', 'a') as f:
                    f.write(json.dumps({'input_index': idx, 'coeff': cubics[idx], 'error': str(exc)}) + '\n')
                print(f'Error processing index {idx}: {exc}', flush=True)
                continue

            processed_data.append({'input_index': returned_idx, **result})
            count += 1
            if count % 100 == 0:
                print(f'{count}/{max_index - start_idx}', flush=True)
                write_to_csv(processed_data, filename)
                processed_data = []

        write_to_csv(processed_data, filename)

    print(f'Saved {count} rows to {filename}; {failed} failed.')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
