#/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CMD="uncrustify --no-backup --replace -c $DIR/uncrustify.cfg"
if [ "$#" -ne 0 ]; then
  $CMD $1
  exit 0
fi
find . -type f -regextype posix-extended -regex '.*.(cpp|hpp|c|h|okl)' -not -path "*/libP/*" -exec $CMD {} \;
