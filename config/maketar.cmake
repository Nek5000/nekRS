# file(ARCHIVE_CREATE ... ) doesn't support wildcards, so we have to exclude
# unwanted files in this temporary directory
file(REMOVE_RECURSE ${DEST_DIR})
file(MAKE_DIRECTORY ${DEST_DIR})
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
  DESTINATION ${DEST_DIR}
  PATTERN ".git" EXCLUDE
  PATTERN ".cache" EXCLUDE)

file(COPY 
  ${LIBP_SOURCE_DIR}/
  DESTINATION ${DEST_DIR}/3rd_party/libparanumal
  PATTERN ".git" EXCLUDE)

file(COPY
  ${HYPRE_SOURCE_DIR}/
  DESTINATION ${DEST_DIR}/3rd_party/hypre
  PATTERN ".git" EXCLUDE)

file(COPY
  ${NEK5000_SOURCE_DIR}/
  DESTINATION ${DEST_DIR}/3rd_party/nek5000
  PATTERN ".git" EXCLUDE
  REGEX ".*\\.o$" EXCLUDE
  REGEX ".*\\.a$" EXCLUDE)
