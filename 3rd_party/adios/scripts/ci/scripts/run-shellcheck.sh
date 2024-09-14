#!/usr/bin/env bash

# We need the same sort order
export LC_ALL=C
EXCLUDE_SCRIPTS=.shellcheck_exclude_paths

echo "---------- Begin ENV ----------"
env | sort
echo "----------  End ENV  ----------"

if [ -n "${SOURCE_DIR}" ]
then
  cd "${SOURCE_DIR}" || exit
fi

# Give me a sorted list of the project scripts
found_scripts="$({
  find scripts source testing -regextype posix-extended -iregex '.*\.(sh|bash)' -print;
  grep -rnlE -e '#!/(/usr)?/bin/(bash|sh)' -e '#!(/usr)?/bin/env\s+(bash|sh)' scripts;
} | sort -u)"

echo "[I] Found the following files:"
echo "$found_scripts"

# Give me the list of scripts without $EXCLUDE_SCRIPTS
if [ -f "$EXCLUDE_SCRIPTS" ]
then
  if ! sort -uc "$EXCLUDE_SCRIPTS"
  then
    echo "[E] file: $EXCLUDE_SCRIPTS is not sorted."
    echo "    To fix this sort with: LC_ALL=C sort -u $EXCLUDE_SCRIPTS"
    exit 1
  fi

  check_scripts=$(comm -2 -3 <(echo "$found_scripts") "$EXCLUDE_SCRIPTS" )
fi

echo "[I] Checking the following files:"
echo "$check_scripts"

if ! xargs -n1 shellcheck <<<"$check_scripts"
then
  echo "[E] shellcheck:"
  echo "    Code format checks failed."
  echo "    Please run shellcheck on your changes before committing."
  exit 2
fi

exit 0
