# (c) 2026 Giorgio Gonnella
"""Run a shell command and report its time and memory.

The command is run with bash -c (pipes allowed); reported: the wall
time, the CPU time (user + system, of all its processes) and the peak
resident memory of its largest process (the maximum of ru_maxrss of the
children: for a pipeline, that of its largest stage, not their sum).

Usage: python3 -I benchmark/measure.py <label> <command>
Output: a TSV line on the standard output: label, wall seconds, CPU
seconds, peak memory (MB); the exit status of the command is returned.
"""

import resource
import subprocess
import sys
import time


def main():
  if len(sys.argv) != 3:
    sys.exit(__doc__)
  label, command = sys.argv[1:]
  start = time.monotonic()
  status = subprocess.call(["bash", "-c", "set -euo pipefail; " + command])
  wall = time.monotonic() - start
  usage = resource.getrusage(resource.RUSAGE_CHILDREN)
  cpu = usage.ru_utime + usage.ru_stime
  print("%s\t%.2f\t%.2f\t%.0f" % (label, wall, cpu, usage.ru_maxrss / 1024))
  sys.exit(status)


if __name__ == "__main__":
  main()
