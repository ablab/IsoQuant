
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import subprocess
import sys

def run_until_output(command, stop_string):
    """
    Run a command and terminate it as soon as its output contains stop_string.

    Args:
        command (list[str]): Command to run, e.g. ["python", "script.py"].
        stop_string (str): Substring to look for in the output.
    """
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,  # line-buffered
        universal_newlines=True
    )

    try:
        for line in process.stdout:
            print(line, end='')  # echo output live
            if stop_string in line:
                print(f"\n[INFO] Found '{stop_string}', terminating process...")
                process.terminate()
                process.wait(timeout=5)
                break
    except KeyboardInterrupt:
        print("\n[INFO] Interrupted by user.")
        process.terminate()
    except Exception as e:
        print(f"[ERROR] {e}")
        process.terminate()
    finally:
        if process.poll() is None:
            process.kill()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: python {sys.argv[0]} <command> <stop_string>")
        sys.exit(1)

    stop_string = sys.argv[2]
    command = sys.argv[1].split()
    run_until_output(command, stop_string)
