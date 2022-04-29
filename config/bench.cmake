set(BENCH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/bench)

set(BENCH_SOURCES
        ${BENCH_SOURCE_DIR}/fdm/benchmarkFDM.cpp
        ${BENCH_SOURCE_DIR}/axHelm/benchmarkAx.cpp
        ${BENCH_SOURCE_DIR}/advsub/benchmarkAdvsub.cpp
        ${BENCH_SOURCE_DIR}/core/kernelBenchmarker.cpp
)