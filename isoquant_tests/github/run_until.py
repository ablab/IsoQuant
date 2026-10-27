############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Run a command and kill it, with every process it started, once its output says so.

Used to produce partial IsoQuant runs for the --resume tests. The kill has to be exact:

- the command runs in its own process group, and the whole group is frozen (SIGSTOP)
  the moment the stop line is read, then killed (SIGKILL) -- IsoQuant's worker
  processes would otherwise keep writing into the folder that becomes test data;
- --after lines (in order) arm the stop, and --occurrence picks the N-th match after
  that, so a line printed by several stages can still pin one of them;
- a run that finishes, or never prints the stop line, is an error: the caller must not
  store a complete run as a partial one.

Usage:
    run_until.py --stop "Finished processing chromosome" --after "Processing assigned reads" \\
                 --occurrence 2 [--log run.log] -- python3 isoquant.py -o out ...

Exit codes: 0 stopped as requested, 2 bad usage, 3 the command ended before the stop line.
"""

import argparse
import os
import signal
import subprocess
import sys
import time
from typing import List, Optional, TextIO


EXIT_STOPPED = 0
EXIT_NOT_STOPPED = 3


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stop", required=True, help="kill the command when a line contains this text")
    parser.add_argument("--after", action="append", default=[],
                        help="arm the stop only after lines containing these texts were seen, in order "
                             "(may be repeated)")
    parser.add_argument("--occurrence", type=int, default=1,
                        help="stop at the N-th matching line after arming [default: 1]")
    parser.add_argument("--log", help="also write the command output to this file")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="command to run, after --")
    args = parser.parse_args(argv)
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("no command given (put it after --)")
    if args.occurrence < 1:
        parser.error("--occurrence must be positive")
    return args


class StopCondition:
    def __init__(self, stop: str, after: List[str], occurrence: int) -> None:
        self.stop = stop
        self.pending_after = list(after)
        self.remaining = occurrence

    def armed(self) -> bool:
        return not self.pending_after

    def feed(self, line: str) -> bool:
        """True when this line is the one to stop at."""
        if self.pending_after:
            if self.pending_after[0] in line:
                self.pending_after.pop(0)
            return False
        if self.stop in line:
            self.remaining -= 1
            return self.remaining == 0
        return False


def group_alive(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def kill_group(process: subprocess.Popen, timeout: float = 30.0) -> None:
    """Freeze, then kill the whole process group and wait until it is gone."""
    pgid = process.pid  # start_new_session makes the child its own group leader
    for sig in (signal.SIGSTOP, signal.SIGKILL):
        try:
            os.killpg(pgid, sig)
        except ProcessLookupError:
            break
    process.wait()
    deadline = time.time() + timeout
    while group_alive(pgid) and time.time() < deadline:
        # orphaned workers are reparented and reaped by init; just wait for them
        time.sleep(0.1)
    if group_alive(pgid):
        print("[run_until] WARNING: processes of group %d are still alive" % pgid, file=sys.stderr)


def run_until(command: List[str], condition: StopCondition, log: Optional[TextIO] = None) -> int:
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                               bufsize=1, start_new_session=True)
    stopped = False
    try:
        for line in process.stdout:
            sys.stdout.write(line)
            if log:
                log.write(line)
            if condition.feed(line):
                kill_group(process)
                stopped = True
                break
    except BaseException:
        kill_group(process)
        raise
    finally:
        process.stdout.close()

    if stopped:
        print("\n[run_until] Stopped at: %s" % line.rstrip("\n"))
        if log:
            log.write("[run_until] Stopped at: %s" % line)
        return EXIT_STOPPED

    returncode = process.wait()
    reason = ("--after lines were not all seen" if not condition.armed()
              else "stop line seen only %d time(s) too few" % condition.remaining)
    print("\n[run_until] ERROR: the command ended (exit code %d) before the stop point: %s"
          % (returncode, reason), file=sys.stderr)
    return EXIT_NOT_STOPPED


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    condition = StopCondition(args.stop, args.after, args.occurrence)
    if args.log:
        with open(args.log, "w") as log:
            return run_until(args.command, condition, log)
    return run_until(args.command, condition)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
