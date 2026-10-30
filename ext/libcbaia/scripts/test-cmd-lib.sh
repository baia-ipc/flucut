#
# Bash library for the program tests, sourced by the test targets of the
# Makefile of a project using libcbaia (the test goal of the libcbaia
# Makefile does not use it), see docs/manual/testing.rst
#

#
# Compare a results file with the expected results (if any, optionally
# after sorting both); on differences the results are copied to
# <expected>.found
#
function _run_diff { local RES=$1 EXP=$2 SOR=$3
  local RESC EXPC
  if [ "$EXP" != "" ]; then
    if [ "$SOR" -eq 1 ]; then
      echo "  Comparing (sorted) output ($RES) to expected results ($EXP)"
      RESC=$RES.s
      sort "$RES" > "$RESC"
      EXPC=$EXP.s
      sort "$EXP" > "$EXPC"
    else
      echo "  Comparing output ($RES) to expected results ($EXP)"
      RESC=$RES
      EXPC=$EXP
    fi
    diff "$RESC" "$EXPC"
    local CODE=$?
    if [ "$SOR" -eq 1 ]; then
      rm -f "$EXPC" "$RESC"
    fi
    if [ $CODE -ne 0 ]; then
      echo "  Comparison failed, differences found!"
      cp "$RES" "$EXP.found"
      echo "  Wrong file copied to $EXP.found"
      return 1
    fi
  fi
  return 0
}

#
# Run a command, check the exit code
# and compare up to three files with known results
#
# CPX and CMD are not quoted, as they are split into words
#
function run_cmp { local LBL=$1  # test label to display
                   local CPX=$2  # command prefix (e.g. valgrind)
                   local CMD=$3  # command to run and display
                   local RED=$4  # redirect output to RES1? 0: no; 1: yes
                   local RES1=$5 # first results file (opt; required if RED)
                   local EXP1=$6 # first expected results file (if RES1)
                   local SOR1=$7 # sort RES1/EXP1 before comparing? 0: no 1: y
                   local RES2=$8 # second results file (optional)
                   local EXP2=$9 # second expected results file (if RES2)
                   local RES3=${10} # third results file (optional)
                   local EXP3=${11} # third expected results file (if RES3)
  echo -e "\nTest $LBL ..."
  echo "  Command: $CMD"

  local FAILED=0

  if [ "$RED" -eq 0 ]
    then $CPX $CMD
    else $CPX $CMD > "$RES1"
  fi
  local CODE=$?
  if [ $CODE -ne 0 ]; then
    echo "  Command failed, exit code: $CODE"
    FAILED=1
  fi

  if [ $FAILED -eq 0 ]; then
    _run_diff "$RES1" "$EXP1" "${SOR1:-0}" || FAILED=1
  fi
  if [ $FAILED -eq 0 ]; then
    _run_diff "$RES2" "$EXP2" 0 || FAILED=1
  fi
  if [ $FAILED -eq 0 ]; then
    _run_diff "$RES3" "$EXP3" 0 || FAILED=1
  fi

  if [ $FAILED -eq 0 ]; then
    if [[ -t 1 ]]; then echo -n -e "\x1b[92m"; fi
    echo -e "\n[Test $LBL] SUCCESS"
    if [[ -t 1 ]]; then echo -n -e "\x1b[0m"; fi
  else
    if [[ -t 1 ]]; then echo -n -e "\x1b[91m"; fi
    echo -e "\n[Test $LBL] FAILED"
    if [[ -t 1 ]]; then echo -n -e "\x1b[0m"; fi
  fi

  return $FAILED
}
