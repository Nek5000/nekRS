#!/bin/bash
# author: Vicente Bolea <vicente.bolea@kitware.com>

function require_dependency()
{
  if ! command -v "$1" &> /dev/null
  then
    echo "[E] Missing dependencies: $1"
    exit 1
  fi
}

require_dependency "jq"
require_dependency "gh"
require_dependency "csvlook"

if [ "$#" != "2" ]
then
  echo "[E] Wrong arguments. Invoke as:"
  echo "scripts/developer/create-changelog.sh <new_release> <old_release>"
  exit 2
fi

new_release="$1"
old_release="$2"

prs="$(git log --oneline  --grep='Merge pull request ' "^${old_release}" "${new_release}" | grep -Po '\s#\K[0-9]+')"
num_prs="$(wc -w <<<"$prs")"

echo "[I] Found $num_prs PRs"

# Header first
output=("PR, Title")
i=0
for pr in ${prs}
do
  printf "\r[I] Processing: PR=%06d progress=%05d/%05d" "$pr" "$i" "$num_prs"
  output+=("$(gh api "/repos/ornladios/ADIOS2/pulls/$pr" | jq -r '["#\(.number)", .title] | @csv')")
  ((i++))
done
echo ""

printf '%s\n' "${output[@]}" | csvlook
