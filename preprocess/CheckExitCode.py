import os
import sys
import argparse

WARNING = '\033[93m'
ERROR = '\033[91m'
ENDC = '\033[0m'

def log_error(log):
    return ERROR + log + ENDC

def log_warning(log):
    return WARNING + log + ENDC

def CheckExitCode(args):

    parallel_log_fn = args.parallel_log_fn
    if not os.path.exists(parallel_log_fn):
        print(log_error("[ERROR] Parallel log file {} not found!".format(parallel_log_fn)))
        sys.exit(1)
    with open(parallel_log_fn, 'r') as log_file:
        log_lines = log_file.readlines()
    if len(log_lines) == 0:
        print(log_warning("[WARNING] Parallel log file {} is empty!".format(parallel_log_fn)))
        sys.exit(0)
    header = log_lines[0].strip().split('\t')
    try:
        exitval_idx = header.index("Exitval")
        signal_idx = header.index("Signal")
        command_idx = header.index("Command")
    except ValueError as e:
        print(log_error(f"[ERROR] Required columns not found in parallel log header: {e}"))
        sys.exit(1)

    has_failures = False

    for line_num, line in enumerate(log_lines[1:], start=2):
        fields = line.strip().split('\t')
        if len(fields) <= max(exitval_idx, signal_idx, command_idx):
            print(log_warning(f"[WARNING] Line {line_num} has insufficient columns, skipping"))
            continue

        exitval = fields[exitval_idx]
        signal = fields[signal_idx]
        command = fields[command_idx]

        if exitval != '0' or signal != '0':
            has_failures = True
            print(log_error(
                f"\n[ERROR] Command failed: {command}\n"
                f"Exit status: {exitval}, Signal: {signal}\n"
                f"{'-' * 80}"
            ))

    if has_failures:
        print(log_error("\n[ERROR] Some parallel jobs failed (see above for details)"))
        sys.exit(1)
    else:
        sys.exit(0)


def main():
    parser = argparse.ArgumentParser(
        description="Check the exit code of parallel jobs.")

    parser.add_argument('--parallel_log_fn', type=str, default=None,
                        help="BAM file input, default: %(default)s")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    CheckExitCode(args)


if __name__ == "__main__":
    main()
