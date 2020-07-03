# file(ARCHIVE_CREATE ... ) doesn't support wildcards, so we have to exclude
# unwanted files in this temporary directory
file(REMOVE_RECURSE ${TARNAME})
file(MAKE_DIRECTORY ${TARNAME})
file(COPY 
  makenrs 
  CMakeLists.txt
  LICENSE
  README.md
  RELEASE.md
  config
  examples 
  okl 
  scripts 
  src 
  DESTINATION ${TARNAME}
  PATTERN ".git" EXCLUDE)

file(COPY 
  ${LIBP_SOURCE_DIR}/
  DESTINATION ${TARNAME}/3rd_party/libparanumal
  PATTERN ".git" EXCLUDE)

file(COPY
  ${HYPRE_SOURCE_DIR}/
  DESTINATION ${TARNAME}/3rd_party/hypre
  PATTERN ".git" EXCLUDE)

