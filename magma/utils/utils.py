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

