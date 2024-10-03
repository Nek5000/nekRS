#/bin/bash

if [ $# -ne 1 ]; then
  echo "usage: $0 <filename|directory>"
  exit 1
fi

if [ -d $1 ]; then
  find $1 -regex '.*\.\(cpp\)' -exec clang-format -style=file -i {} \;
  find $1 -regex '.*\.\(c\)' -exec clang-format -style=file -i {} \;
  find $1 -regex '.*\.\(udf\)' -exec clang-format -style=file -i {} \;
  find $1 -regex '.*\.\(hpp\)' -exec clang-format -style=file -i {} \;
  find $1 -regex '.*\.\(h\)' -exec clang-format -style=file -i {} \;
  find $1 -regex '.*\.\(okl\)' -exec clang-format -style=file -i {} \;
fi

if [ -f $1 ]; then
  clang-format -style=file -i $1;
fi
