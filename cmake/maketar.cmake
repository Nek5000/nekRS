# file(ARCHIVE_CREATE ... ) doesn't support wildcards, so we have to exclude
# unwanted files in this temporary directory
file(REMOVE_RECURSE ${DEST_DIR})
file(MAKE_DIRECTORY ${DEST_DIR})
file(COPY
  3rd_party
  build.sh 
  CMakeLists.txt
  LICENSE
  README.md
  RELEASE.md
  CONTRIBUTING.md
  cmake 
  doc
  examples 
  scripts 
  src 
  DESTINATION ${DEST_DIR}
  PATTERN ".gitignore" EXCLUDE
  PATTERN ".git" EXCLUDE
  PATTERN ".cache" EXCLUDE)
