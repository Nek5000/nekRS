add_executable(axhelm-bin src/bench/axHelm/main.cpp)
set_target_properties(axhelm-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-axhelm)
target_link_libraries(axhelm-bin PRIVATE nekrs-lib)

add_executable(advsub-bin src/bench/advsub/main.cpp)
set_target_properties(advsub-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-advsub)
target_link_libraries(advsub-bin PRIVATE nekrs-lib)

add_executable(fdm-bin src/bench/fdm/main.cpp)
set_target_properties(fdm-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-fdm)
target_link_libraries(fdm-bin PRIVATE nekrs-lib)

set(BENCH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/bench)
set(BENCH_SOURCES
        ${BENCH_SOURCE_DIR}/fdm/benchmarkFDM.cpp
        ${BENCH_SOURCE_DIR}/axHelm/benchmarkAx.cpp
        ${BENCH_SOURCE_DIR}/advsub/benchmarkAdvsub.cpp
        ${BENCH_SOURCE_DIR}/core/kernelBenchmarker.cpp
)
