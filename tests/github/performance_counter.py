#!/usr/bin/env python3
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
# Author: Andrey Prjibelski
############################################################################

import argparse
import os
import sys
import psutil
import subprocess
from traceback import print_exc
from time import time
from time import sleep
from collections import defaultdict
from threading import Thread
from threading import Event
from datetime import datetime


task_over = Event()


def to_gb(val):
    return float(val) / float(1024 * 1024 * 1024)


def human_readable_time(time_secs):
    time_secs = int(time_secs)
    h = time_secs // 3600
    remain_secs = time_secs - h * 3600
    m = remain_secs // 60
    s = remain_secs % 60
    m_prefix = "0" if m < 10 else ""
    s_prefix = "0" if s < 10 else ""
    return "%d:%s%d:%s%d" % (h, m_prefix, m, s_prefix, s)


def get_current_stats(pids):
    total_vms = 0
    total_rss = 0
    total_uss = 0
    max_shared = 0
    total_cpu = 0.0
    for pid in pids.keys():
        proc = pids[pid]
        try:
            mem_info = proc.memory_info()
            total_vms += mem_info.vms
            total_rss += mem_info.rss
            total_cpu += proc.cpu_percent(interval=None)
            # USS (Unique Set Size) = memory unique to this process, excluding shared
            # This avoids overcounting shared memory across processes
            try:
                full_info = proc.memory_full_info()
                total_uss += full_info.uss
                # Track max shared memory from any process (shared mem is the same across processes)
                max_shared = max(max_shared, full_info.shared)
            except (psutil.AccessDenied, AttributeError):
                # memory_full_info may require elevated privileges or not be available
                total_uss += mem_info.rss
        except psutil.NoSuchProcess:
            pass

    # Effective RAM = USS (unique to each process) + shared memory (counted once)
    effective_ram = total_uss + max_shared
    return total_rss, total_vms, total_uss, max_shared, effective_ram, total_cpu


def get_children(pid):
    pids = {pid}
    for p in psutil.Process(pid).children(recursive=True):
        pids.add(p.pid)
    return pids


def track_ram_and_cpu(task_subprocess, interval, outfile):
    outf = open(outfile, "w")
    outf.write("Time\tRSS\tUSS\tShared\tEffective\tVMS\tCPU %\n")

    task_pid = task_subprocess.pid
    task_process_obj = psutil.Process(task_pid)
    all_pids = {task_pid : task_process_obj}
    start_time = time()
    max_rss = 0
    max_uss = 0
    max_effective = 0
    cpu_times = defaultdict(float)
    flush_interval = 0

    while task_process_obj.status() != 'zombie':
        children_pids = get_children(task_pid)
        for child in children_pids:
            if child not in all_pids:
                all_pids[child] = psutil.Process(child)
        total_rss, total_vms, total_uss, max_shared, effective_ram, total_cpu = get_current_stats(all_pids)
        max_rss = max(max_rss, total_rss)
        max_uss = max(max_uss, total_uss)
        max_effective = max(max_effective, effective_ram)
        current_time = time() - start_time
        outf.write("%d\t%d\t%d\t%d\t%d\t%d\t%.1f\n" % (current_time, total_rss, total_uss, max_shared, effective_ram, total_vms, total_cpu))
        if flush_interval > 60:
            flush_interval = 0
            outf.flush()
        for p in all_pids.keys():
            try:
                cpu_times[p] = max(cpu_times[p], all_pids[p].cpu_times().user)
            except psutil.NoSuchProcess:
                pass
        sleep(interval)
        flush_interval += interval

    outf.close()
    cpu_time = sum(cpu_times.values())
    wall_clock_time = time() - start_time
    human_readable_stats = ("Max RSS: %.3f GB\n"
                            "Max USS: %.3f GB (excluding shared memory)\n"
                            "Max Effective: %.3f GB (USS + shared memory once)\n"
                            "CPU time:\t%s\n"
                            "Wall clock time:\t%s") % (to_gb(max_rss),
                                                       to_gb(max_uss),
                                                       to_gb(max_effective),
                                                       human_readable_time(cpu_time),
                                                       human_readable_time(wall_clock_time))
    tsv_stats = ("max_rss\t%d\t%.3fGB\n"
                 "max_uss\t%d\t%.3fGB\n"
                 "max_effective\t%d\t%.3fGB\n"
                 "cpu_time\t%.0f\t%s\n"
                 "clock_time\t%.0f\t%s") % (max_rss, to_gb(max_rss),
                                             max_uss, to_gb(max_uss),
                                             max_effective, to_gb(max_effective),
                                             cpu_time, human_readable_time(cpu_time),
                                             wall_clock_time, human_readable_time(wall_clock_time))

    return human_readable_stats, tsv_stats


def track_disk_usage(folder, interval, outfile, res_stats):
    outf = open(outfile, "w")
    outf.write("time\tdisk\n")
    start_time = time()
    max_usage = 0
    while not task_over.is_set():
        if not os.path.exists(folder):
            sleep(interval)
            continue

        disk_usage = int(subprocess.check_output(['du', '-s', folder]).split()[0].decode('utf-8'))
        current_time = time() - start_time
        max_usage =max(max_usage, disk_usage)
        outf.write("%d\t%d\n" % (current_time, disk_usage))
        sleep(interval)
    outf.close()
    res_stats = "Maximum disk space taken: %.3f GB" % to_gb(max_usage)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="A tiny util for measuring RAM, CPU and disk consumption. "
                                                 "All children processes will be taken into account and summed up.")
    input_args = parser.add_mutually_exclusive_group()
    input_args.add_argument("--cmd", "-c", help="command line to execute (from current dir); "
                                                "~ alias for $HOME is not supported", type=str)
    input_args.add_argument("--cmd_file", help="file with command line to execute (from current dir); "
                                               "~ alias for $HOME is not supported", type=str)

    parser.add_argument("--interval", "-i", help="time interval between measurements (seconds)", type=int, default=2)
    parser.add_argument("--out_dir", help="command line output folder to monitor disk usage; "
                                          "will not be monitored if not set; "
                                          "note, that disk monitoring may affect perfromance", type=str)
    parser.add_argument("--disk_interval", help="time interval between disk measurements (seconds)", type=int, default=10)
    parser.add_argument("--output", "-o", help="output folder, default performance_stats/YYYY_MM_DD_HH_MM_SS/", type=str)
    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    if args.cmd:
        cmd = args.cmd.split()
    elif args.cmd_file:
        cmd = open(args.cmd_file).readline().strip().split()
    else:
        sys.stderr.write("Provide command line to run with --cmd or --cmd_file")
        exit(-1)
    print("Running %s" % " ".join(cmd))

    if not args.output:
        args.output = "performance_stats/" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    with open(os.path.join(args.output, "command.txt"), "w") as cmd_file:
        cmd_file.write(" ".join(cmd) + "\n")

    thread = None
    disk_stats = ""
    if args.out_dir:
        disk_stat_file = os.path.join(args.output, "disk_usage.tsv")
        thread = Thread(target=track_disk_usage, args=(args.out_dir, args.disk_interval, disk_stat_file, disk_stats))
        thread.start()

    task_subprocess = subprocess.Popen(cmd, shell=False)
    ram_cpu_stat_file = os.path.join(args.output, "ram_cpu_usage.tsv")
    human_readable_stats, tsv_stats = track_ram_and_cpu(task_subprocess, args.interval, ram_cpu_stat_file)

    task_over.set()
    if thread:
        thread.join()

    print(human_readable_stats)
    print(disk_stats)
    with open(os.path.join(args.output, "stats.tsv"), "w") as stat_file:
        stat_file.write(tsv_stats + "\n")
        stat_file.write(disk_stats + "\n")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)


