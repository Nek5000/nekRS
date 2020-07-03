# file(ARCHIVE_CREATE ... ) doesn't support wildcards, so we have to exclude
# unwanted files in this temporary directory
file(MAKE_DIRECTORY ${TARNAME})
file(COPY 
  makenrs 
  CMakeLists.txt
  LICENSE
  README.md
  RELEASE.md
  3rd_party 
  examples 
  okl 
  scripts 
  src 
  DESTINATION ${TARNAME}
  REGEX ".git" EXCLUDE)
