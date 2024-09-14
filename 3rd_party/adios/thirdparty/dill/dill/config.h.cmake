/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define with the build version of dill */
#cmakedefine DILL_VERSION "@DILL_VERSION@"

/* Define if bfd functions use bfd_byte arguments. */
#undef BFD_BYTE

/* Define to 1 if the emulator should be built */
#cmakedefine BUILD_EMULATOR

/* Define if dill should attempt to DCG */
#cmakedefine EMULATION_ONLY

/* Define to 1 if you have the <dis-asm.h> header file. */
#cmakedefine HAVE_DIS_ASM_H

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H

/* Define to 1 if you have __clear_cache is defined  */
#cmakedefine CLEAR_CACHE_DEFINED

/* Define to 1 if you have the <malloc.h> header file. */
#cmakedefine HAVE_MALLOC_H

/* Define to 1 if you have the <sys/mman.h> header file. */
#cmakedefine HAVE_SYS_MMAN_H

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H

/* Define if you have the `print_insn_arm' function. */
#cmakedefine HAVE_PRINT_INSN_ARM

/* Define if you have the `print_insn_big_powerpc' function. */
#cmakedefine HAVE_PRINT_INSN_BIG_POWERPC

/* Define if you have the `print_insn_i386' function. */
#cmakedefine HAVE_PRINT_INSN_I386

/* Define if you have the `print_insn_ia64' function. */
#cmakedefine HAVE_PRINT_INSN_IA64

/* Define if you have the `print_insn_little_arm' function. */
#cmakedefine HAVE_PRINT_INSN_LITTLE_ARM

/* Define if you have the `print_insn_sparc' function. */
#cmakedefine HAVE_PRINT_INSN_SPARC

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H

/* Define if the host processor is ARM5 */
#cmakedefine HOST_ARM5

/* Define if the host processor is ARM6 */
#cmakedefine HOST_ARM6 @HOST_ARM6@

/* Define if the host processor is ARM7 */
#cmakedefine HOST_ARM7 @HOST_ARM7@

/* Define if the host processor is ARM8 */
#cmakedefine HOST_ARM8 @HOST_ARM8@

/* Define if the host processor is an ia64 */
#cmakedefine HOST_IA64

/* Define if the host processor is a powerpc */
#cmakedefine HOST_POWERPC

/* Define if the host processor is a powerpc64 */
#cmakedefine HOST_POWERPC64

/* Define if the host processor is a ppc64le */
#cmakedefine HOST_PPC64LE

/* Define if the host processor is a sparc */
#cmakedefine HOST_SPARC

/* Define if the host processor is a sparcv9 */
#cmakedefine HOST_SPARCV9

/* Define if the host processor is an x86 */
#cmakedefine HOST_X86

/* Define if the host processor is an x86_64 */
#cmakedefine HOST_X86_64

/* Define if INIT_DISASSEMBLE_INFO takes three arguments instead of two */
#cmakedefine INIT_DISASSEMBLE_INFO_THREE_ARG

/* Define if integrating with kernel plugins */
#cmakedefine KPLUGINS_INTEGRATION

/* Define if compiling for linux kernel */
#cmakedefine LINUX_KERNEL_MODULE

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#cmakedefine LT_OBJDIR

/* Define if you want more than just native target */
#cmakedefine MULTI_TARGET

/* Define for the host architecture type */
#cmakedefine NATIVE_ARCH "@NATIVE_ARCH@"

/* Define if there is no disassembler */
#cmakedefine NO_DISASSEMBLER

/* Define if we should not use inlined procedures from BFD */
#cmakedefine NO_INLINED_BFD_PROCS

/* The number of bytes in type long */
#define SIZEOF_LONG @SIZEOF_LONG@

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS

/* Define if the membar instruction should be used to sync Icache and Dcache
   */
#cmakedefine USE_MEMBAR

/* Define this if mmap should be used instead of malloc() for code memory */
#cmakedefine USE_MMAP_CODE_SEG

/* Define this if VirtualProtect should be used to change memory protections */
#cmakedefine USE_VIRTUAL_PROTECT

/* Define this if windows calling convention should be used */
#cmakedefine USE_WINDOWS_CALLS

/* Define if byteorder is bigendian */
#cmakedefine WORDS_BIGENDIAN

/* Define if using ARM hardware float ABIs for parameter passing */
#define ARM_HARD_FLOAT  @ARM_HARD_FLOAT@

/* Define if using ARMv7 hardware */
#cmakedefine ARMV7_AVAILABLE

/* Define if we shouldn't use native DCG */
#cmakedefine DILL_IGNORE_NATIVE

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#undef inline
#endif
