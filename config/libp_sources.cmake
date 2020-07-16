set(LIBP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/libP)
set(OGS_SOURCE_DIR ${LIBP_SOURCE_DIR}/libs/gatherScatter)
set(PARALMOND_SOURCE_DIR ${LIBP_SOURCE_DIR}/libs/parAlmond)
set(ELLIPTIC_SOURCE_DIR ${LIBP_SOURCE_DIR}/solvers/elliptic)
set(BLASLAPACK_LIBP_DIR ${LIBP_SOURCE_DIR}/3rdParty/BlasLapack)

set(BLASLAPACK_LIBP_SOURCES 
        ${BLASLAPACK_LIBP_DIR}/dasum.f
        ${BLASLAPACK_LIBP_DIR}/daxpy.f
        ${BLASLAPACK_LIBP_DIR}/dbdsqr.f
        ${BLASLAPACK_LIBP_DIR}/dcabs1.f
        ${BLASLAPACK_LIBP_DIR}/dcopy.f
        ${BLASLAPACK_LIBP_DIR}/ddot.f
        ${BLASLAPACK_LIBP_DIR}/dgebak.f
        ${BLASLAPACK_LIBP_DIR}/dgebal.f
        ${BLASLAPACK_LIBP_DIR}/dgebd2.f
        ${BLASLAPACK_LIBP_DIR}/dgebrd.f
        ${BLASLAPACK_LIBP_DIR}/dgecon.f
        ${BLASLAPACK_LIBP_DIR}/dgeev.f
        ${BLASLAPACK_LIBP_DIR}/dgehd2.f
        ${BLASLAPACK_LIBP_DIR}/dgehrd.f
        ${BLASLAPACK_LIBP_DIR}/dgelq2.f
        ${BLASLAPACK_LIBP_DIR}/dgelqf.f
        ${BLASLAPACK_LIBP_DIR}/dgels.f
        ${BLASLAPACK_LIBP_DIR}/dgelss.f
        ${BLASLAPACK_LIBP_DIR}/dgemm.f
        ${BLASLAPACK_LIBP_DIR}/dgemv.f
        ${BLASLAPACK_LIBP_DIR}/dgeqr2.f
        ${BLASLAPACK_LIBP_DIR}/dgeqrf.f
        ${BLASLAPACK_LIBP_DIR}/dger.f
        ${BLASLAPACK_LIBP_DIR}/dgesvd.f
        ${BLASLAPACK_LIBP_DIR}/dgesv.f
        ${BLASLAPACK_LIBP_DIR}/dgetf2.f
        ${BLASLAPACK_LIBP_DIR}/dgetrf.f
        ${BLASLAPACK_LIBP_DIR}/dgetri.f
        ${BLASLAPACK_LIBP_DIR}/dgetrs.f
        ${BLASLAPACK_LIBP_DIR}/dhseqr.f
        ${BLASLAPACK_LIBP_DIR}/dlabad.f
        ${BLASLAPACK_LIBP_DIR}/dlabrd.f
        ${BLASLAPACK_LIBP_DIR}/dlacon.f
        ${BLASLAPACK_LIBP_DIR}/dlacpy.f
        ${BLASLAPACK_LIBP_DIR}/dladiv.f
        ${BLASLAPACK_LIBP_DIR}/dlae2.f
        ${BLASLAPACK_LIBP_DIR}/dlaev2.f
        ${BLASLAPACK_LIBP_DIR}/dlahqr.f
        ${BLASLAPACK_LIBP_DIR}/dlahrd.f
        ${BLASLAPACK_LIBP_DIR}/dlaln2.f
        ${BLASLAPACK_LIBP_DIR}/dlamch.f
        ${BLASLAPACK_LIBP_DIR}/dlange.f
        ${BLASLAPACK_LIBP_DIR}/dlanhs.f
        ${BLASLAPACK_LIBP_DIR}/dlanst.f
        ${BLASLAPACK_LIBP_DIR}/dlansy.f
        ${BLASLAPACK_LIBP_DIR}/dlanv2.f
        ${BLASLAPACK_LIBP_DIR}/dlapy2.f
        ${BLASLAPACK_LIBP_DIR}/dlarfb.f
        ${BLASLAPACK_LIBP_DIR}/dlarf.f
        ${BLASLAPACK_LIBP_DIR}/dlarfg.f
        ${BLASLAPACK_LIBP_DIR}/dlarft.f
        ${BLASLAPACK_LIBP_DIR}/dlarfx.f
        ${BLASLAPACK_LIBP_DIR}/dlartg.f
        ${BLASLAPACK_LIBP_DIR}/dlas2.f
        ${BLASLAPACK_LIBP_DIR}/dlascl.f
        ${BLASLAPACK_LIBP_DIR}/dlaset.f
        ${BLASLAPACK_LIBP_DIR}/dlasq1.f
        ${BLASLAPACK_LIBP_DIR}/dlasq2.f
        ${BLASLAPACK_LIBP_DIR}/dlasq3.f
        ${BLASLAPACK_LIBP_DIR}/dlasq4.f
        ${BLASLAPACK_LIBP_DIR}/dlasq5.f
        ${BLASLAPACK_LIBP_DIR}/dlasq6.f
        ${BLASLAPACK_LIBP_DIR}/dlasr.f
        ${BLASLAPACK_LIBP_DIR}/dlasrt.f
        ${BLASLAPACK_LIBP_DIR}/dlassq.f
        ${BLASLAPACK_LIBP_DIR}/dlasv2.f
        ${BLASLAPACK_LIBP_DIR}/dlaswp.f
        ${BLASLAPACK_LIBP_DIR}/dlatrd.f
        ${BLASLAPACK_LIBP_DIR}/dlatrs.f
        ${BLASLAPACK_LIBP_DIR}/dnrm2.f
        ${BLASLAPACK_LIBP_DIR}/dorg2l.f
        ${BLASLAPACK_LIBP_DIR}/dorg2r.f
        ${BLASLAPACK_LIBP_DIR}/dorgbr.f
        ${BLASLAPACK_LIBP_DIR}/dorghr.f
        ${BLASLAPACK_LIBP_DIR}/dorgl2.f
        ${BLASLAPACK_LIBP_DIR}/dorglq.f
        ${BLASLAPACK_LIBP_DIR}/dorgql.f
        ${BLASLAPACK_LIBP_DIR}/dorgqr.f
        ${BLASLAPACK_LIBP_DIR}/dorgtr.f
        ${BLASLAPACK_LIBP_DIR}/dorm2r.f
        ${BLASLAPACK_LIBP_DIR}/dormbr.f
        ${BLASLAPACK_LIBP_DIR}/dorml2.f
        ${BLASLAPACK_LIBP_DIR}/dormlq.f
        ${BLASLAPACK_LIBP_DIR}/dormqr.f
        ${BLASLAPACK_LIBP_DIR}/dpbtf2.f
        ${BLASLAPACK_LIBP_DIR}/dpbtrf.f
        ${BLASLAPACK_LIBP_DIR}/dpbtrs.f
        ${BLASLAPACK_LIBP_DIR}/dposv.f
        ${BLASLAPACK_LIBP_DIR}/dpotf2.f
        ${BLASLAPACK_LIBP_DIR}/dpotrf.f
        ${BLASLAPACK_LIBP_DIR}/dpotrs.f
        ${BLASLAPACK_LIBP_DIR}/drot.f
        ${BLASLAPACK_LIBP_DIR}/drscl.f
        ${BLASLAPACK_LIBP_DIR}/dscal.f
        ${BLASLAPACK_LIBP_DIR}/dsteqr.f
        ${BLASLAPACK_LIBP_DIR}/dsterf.f
        ${BLASLAPACK_LIBP_DIR}/dstev.f
        ${BLASLAPACK_LIBP_DIR}/dswap.f
        ${BLASLAPACK_LIBP_DIR}/dsyev.f
        ${BLASLAPACK_LIBP_DIR}/dsygs2.f
        ${BLASLAPACK_LIBP_DIR}/dsygst.f
        ${BLASLAPACK_LIBP_DIR}/dsygv.f
        ${BLASLAPACK_LIBP_DIR}/dsymm.f
        ${BLASLAPACK_LIBP_DIR}/dsymv.f
        ${BLASLAPACK_LIBP_DIR}/dsyr2.f
        ${BLASLAPACK_LIBP_DIR}/dsyr2k.f
        ${BLASLAPACK_LIBP_DIR}/dsyr.f
        ${BLASLAPACK_LIBP_DIR}/dsyrfs.f
        ${BLASLAPACK_LIBP_DIR}/dsyrk.f
        ${BLASLAPACK_LIBP_DIR}/dsytd2.f
        ${BLASLAPACK_LIBP_DIR}/dsytrd.f
        ${BLASLAPACK_LIBP_DIR}/dsytrs.f
        ${BLASLAPACK_LIBP_DIR}/dtbsv.f
        ${BLASLAPACK_LIBP_DIR}/dtrevc.f
        ${BLASLAPACK_LIBP_DIR}/dtrmm.f
        ${BLASLAPACK_LIBP_DIR}/dtrmv.f
        ${BLASLAPACK_LIBP_DIR}/dtrsm.f
        ${BLASLAPACK_LIBP_DIR}/dtrsv.f
        ${BLASLAPACK_LIBP_DIR}/dtrti2.f
        ${BLASLAPACK_LIBP_DIR}/dtrtri.f
        ${BLASLAPACK_LIBP_DIR}/dzasum.f
        ${BLASLAPACK_LIBP_DIR}/dzsum1.f
        ${BLASLAPACK_LIBP_DIR}/idamax.f
        ${BLASLAPACK_LIBP_DIR}/ieeeck.f
        ${BLASLAPACK_LIBP_DIR}/ilaenv.f
        ${BLASLAPACK_LIBP_DIR}/izamax.f
        ${BLASLAPACK_LIBP_DIR}/izmax1.f
        ${BLASLAPACK_LIBP_DIR}/lsame.f
        ${BLASLAPACK_LIBP_DIR}/xerbla.f
        ${BLASLAPACK_LIBP_DIR}/zaxpy.f
        ${BLASLAPACK_LIBP_DIR}/zcopy.f
        ${BLASLAPACK_LIBP_DIR}/zdotc.f
        ${BLASLAPACK_LIBP_DIR}/zdotu.f
        ${BLASLAPACK_LIBP_DIR}/zdrscl.f
        ${BLASLAPACK_LIBP_DIR}/zdscal.f
        ${BLASLAPACK_LIBP_DIR}/zgecon.f
        ${BLASLAPACK_LIBP_DIR}/zgemm.f
        ${BLASLAPACK_LIBP_DIR}/zgeru.f
        ${BLASLAPACK_LIBP_DIR}/zgetf2.f
        ${BLASLAPACK_LIBP_DIR}/zgetrf.f
        ${BLASLAPACK_LIBP_DIR}/zlacon.f
        ${BLASLAPACK_LIBP_DIR}/zladiv.f
        ${BLASLAPACK_LIBP_DIR}/zlaswp.f
        ${BLASLAPACK_LIBP_DIR}/zlatrs.f
        ${BLASLAPACK_LIBP_DIR}/zscal.f
        ${BLASLAPACK_LIBP_DIR}/zswap.f
        ${BLASLAPACK_LIBP_DIR}/ztrsm.f
        ${BLASLAPACK_LIBP_DIR}/ztrsv.f)

set(OGS_SOURCES
        ${OGS_SOURCE_DIR}/src/ogsGather.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherMany.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherScatter.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherScatterMany.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherScatterVec.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherVec.cpp
        ${OGS_SOURCE_DIR}/src/ogsHostGather.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherMany.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherScatter.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherScatterMany.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherScatterVec.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherVec.c
        ${OGS_SOURCE_DIR}/src/ogsHostScatter.c
        ${OGS_SOURCE_DIR}/src/ogsHostScatterMany.c
        ${OGS_SOURCE_DIR}/src/ogsHostScatterVec.c
        ${OGS_SOURCE_DIR}/src/ogsHostSetup.c
        ${OGS_SOURCE_DIR}/src/ogsKernels.cpp
        ${OGS_SOURCE_DIR}/src/ogsMappedAlloc.cpp
        ${OGS_SOURCE_DIR}/src/ogsScatter.cpp
        ${OGS_SOURCE_DIR}/src/ogsScatterMany.cpp
        ${OGS_SOURCE_DIR}/src/ogsScatterVec.cpp
        ${OGS_SOURCE_DIR}/src/ogsSetup.cpp
        ${OGS_SOURCE_DIR}/src/oogs.cpp)

set(LIBP_SOURCES
        ${LIBP_SOURCE_DIR}/src/hash.c
        ${LIBP_SOURCE_DIR}/src/matrixConditionNumber.c
        ${LIBP_SOURCE_DIR}/src/matrixInverse.c
        ${LIBP_SOURCE_DIR}/src/meshApplyElementMatrix.c
        ${LIBP_SOURCE_DIR}/src/meshConnect.c
        ${LIBP_SOURCE_DIR}/src/meshConnectBoundary.c
        ${LIBP_SOURCE_DIR}/src/meshConnectFaceNodes2D.c
        ${LIBP_SOURCE_DIR}/src/meshConnectFaceNodes3D.c
        ${LIBP_SOURCE_DIR}/src/meshConnectPeriodicFaceNodes2D.c
        ${LIBP_SOURCE_DIR}/src/meshConnectPeriodicFaceNodes3D.c
        ${LIBP_SOURCE_DIR}/src/meshFree.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricPartition2D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricPartition3D.c
        ${LIBP_SOURCE_DIR}/src/meshHaloExchange.c
        ${LIBP_SOURCE_DIR}/src/meshHaloExtract.c
        ${LIBP_SOURCE_DIR}/src/meshHaloSetup.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshOccaSetup2D.c
        ${LIBP_SOURCE_DIR}/src/meshOccaSetup3D.c
        ${LIBP_SOURCE_DIR}/src/meshOccaSetupQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelConnectOpt.c
        ${LIBP_SOURCE_DIR}/src/meshParallelConsecutiveGlobalNumbering.c
        ${LIBP_SOURCE_DIR}/src/meshParallelGatherScatterSetup.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshPartitionStatistics.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshPlotVTU2D.c
        ${LIBP_SOURCE_DIR}/src/meshPlotVTU3D.c
        ${LIBP_SOURCE_DIR}/src/meshPrint2D.c
        ${LIBP_SOURCE_DIR}/src/meshPrint3D.c
        ${LIBP_SOURCE_DIR}/src/meshSetup.c
        ${LIBP_SOURCE_DIR}/src/meshSetupBoxHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshSetupBoxQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshSetupHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshSetupQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshSetupQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshSetupTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshSetupTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshSetupTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshVTU2D.c
        ${LIBP_SOURCE_DIR}/src/meshVTU3D.c
        ${LIBP_SOURCE_DIR}/src/mysort.c
        ${LIBP_SOURCE_DIR}/src/occaHostMallocPinned.c
        ${LIBP_SOURCE_DIR}/src/parallelSort.c
        ${LIBP_SOURCE_DIR}/src/readArray.c
        ${LIBP_SOURCE_DIR}/src/setupAide.c
        ${LIBP_SOURCE_DIR}/src/timer.c)

set_source_files_properties(${LIBP_SOURCES} PROPERTIES LANGUAGE CXX)

set(PARALMOND_SOURCES
        ${PARALMOND_SOURCE_DIR}/hypre/hypre.c
        ${PARALMOND_SOURCE_DIR}/src/SpMV.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgLevel.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/agmgSetup.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/constructProlongation.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/formAggregates.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/galerkinProd.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/strongGraph.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/transpose.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSmoother.cpp
        ${PARALMOND_SOURCE_DIR}/src/coarseSolver.cpp
        ${PARALMOND_SOURCE_DIR}/src/kernels.cpp
        ${PARALMOND_SOURCE_DIR}/src/level.cpp
        ${PARALMOND_SOURCE_DIR}/src/matrix.cpp
        ${PARALMOND_SOURCE_DIR}/src/multigrid.cpp
        ${PARALMOND_SOURCE_DIR}/src/parAlmond.cpp
        ${PARALMOND_SOURCE_DIR}/src/pcg.cpp
        ${PARALMOND_SOURCE_DIR}/src/pgmres.cpp
        ${PARALMOND_SOURCE_DIR}/src/solver.cpp
        ${PARALMOND_SOURCE_DIR}/src/timer.cpp
        ${PARALMOND_SOURCE_DIR}/src/utils.cpp
        ${PARALMOND_SOURCE_DIR}/src/vector.cpp)

# ---------------------------------------------------------
# libelliptic
# ---------------------------------------------------------

set(ELLIPTIC_SOURCES
        ${ELLIPTIC_SOURCE_DIR}/src/NBFPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/NBPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/PCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildContinuous.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildContinuousGalerkin.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildIpdg.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildJacobi.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildLocalPatches.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildMultigridLevel.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticHaloExchange.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticKernelInfo.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMixedCopy.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridLevel.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticResidualProjection.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridLevelSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridSchwarz.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticOperator.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticPlotVTUHex3D.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticPreconditioner.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticPreconditionerSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticSEMFEMSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticScaledAdd.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticSetScalar.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticSolve.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticSolveSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticThinOas.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticThinOasSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticUpdateNBFPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticUpdateNBPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticUpdatePCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticVectors.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticWeightedInnerProduct.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticWeightedNorm2.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticZeroMean.c)

set_source_files_properties(${ELLIPTIC_SOURCES} PROPERTIES LANGUAGE CXX)

# ---------------------------------------------------------
# install
# ---------------------------------------------------------

set(file_pattern "\.okl$|\.c$|\.hpp$|\.tpp$|\.h$|hex.*\.dat$")

install(DIRECTORY 
  ${LIBP_SOURCE_DIR}/include 
  ${LIBP_SOURCE_DIR}/nodes
  ${LIBP_SOURCE_DIR}/okl 
  DESTINATION libparanumal
  FILES_MATCHING REGEX ${file_pattern})

install(DIRECTORY
  ${OGS_SOURCE_DIR}/include 
  ${OGS_SOURCE_DIR}/okl 
  DESTINATION gatherScatter
  FILES_MATCHING REGEX ${file_pattern})
install(FILES ${OGS_SOURCE_DIR}/ogs.hpp DESTINATION gatherScatter)

install(DIRECTORY
  ${PARALMOND_SOURCE_DIR}/include
  ${PARALMOND_SOURCE_DIR}/okl
  DESTINATION parAlmond
  FILES_MATCHING REGEX ${file_pattern})
install(FILES ${PARALMOND_SOURCE_DIR}/parAlmond.hpp DESTINATION gatherScatter)

install(DIRECTORY
  ${ELLIPTIC_SOURCE_DIR}/data
  ${ELLIPTIC_SOURCE_DIR}/okl
  DESTINATION elliptic
  FILES_MATCHING REGEX ${file_pattern})
install(FILES 
  ${ELLIPTIC_SOURCE_DIR}/elliptic.h
  ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGrid.h
  ${ELLIPTIC_SOURCE_DIR}/ellipticResidualProjection.h
  ${ELLIPTIC_SOURCE_DIR}/ellipticPrecon.h
  DESTINATION gatherScatter)
