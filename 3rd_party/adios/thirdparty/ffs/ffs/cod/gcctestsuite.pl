#perl
use File::Basename;
#$TESTSUITE_PATH="$ENV{HOME}/prog/gcc-3.3.1-3/gcc/testsuite/gcc.c-torture/execute";

$TESTSUITE_PATH="$ENV{HOME}/prog/tinycc/tests/tests2";

$BUILD_PATH="$ENV{HOME}/";

@files = <$TESTSUITE_PATH/*>;

my $compile_failed_count = 0;
my $execute_failed_count = 0;
my $expected_failure_count = 0;

%tcc_exceptions =  ( "06_case" => "No CASE!",
		     "12_hashdefine" => "Uses #define",
		     "15_recursion" => "CoD can't do recursion",
		     "18_include" => "CoD has special #include",
		     "25_quicksort" => "CoD can't do recursion",
		     "30_hanoi" => "Uses #define, recurses",
		     "31_args" => "CoD doesn't have ARGV",
		     "32_led" => "Uses #define",
		     "40_stdio" => "No stdio in CoD",
		     "41_hashif" => "Uses #if",
		     "42_function_pointer" => "No function pointers in CoD",
		     "46_grep" => "Grep requires ARGV, opening files, etc.  Too ambitious",
		     "47_switch_return" => "No CASE!",
		     "59_function_array" => "No function pointers",
		     "55_lshift_type" => "Uses #define",
		     "64_macro_nesting" => "No Macros",
		     "65_macro_concat_start" => "No Macros",
		     "66_macro_concat_end" => "No Macros",
		     "67_macro_concat" => "No Macros",
		     "68_macro_param_list_err_1" => "No Macros",
		     "69_macro_param_list_err_2" => "No Macros",
		     "test" => "reason");

%gcc_exceptions =  ( "20000205-1" => "Recurses",
		     "20000113-1" => "Uses bitfields",
		     "20000223-1" => "Uses #define",
		     "20000402-1" => "Uses #if",
		     "20000403-1" => "forward declarations",
		     "20000412-1" => "Recurses",
		     "20000818-1" => "Uses #if",
		     "20000822-1" => "Uses #if",
		     "20001229-1" => "Uses #if",
		     "20010119-1" => "Uses #if",
		     "20010122-1" => "Uses #define",
		     "20010724-1" => "Uses #define",
		     "20011008-3" => "Uses #define",
		     "20011115-1" => "Uses #if",
		     "20020108-1" => "Uses #define",
		     "20020226-1" => "Uses #define",
		     "20020307-1" => "Uses #define",
		     "20020402-1" => "Uses #define",
		     "20020506-1" => "Uses #define",
		     "20020508-1" => "Uses #define",
		     "20020508-2" => "Uses #define",
		     "20020508-3" => "Uses #define",
		     "20020619-1" => "Uses #define",
		     "20020720-1" => "Uses #if",
		     "20021120-1" => "Uses #define",
		     "20021120-3" => "Uses #define",
		     "20030613-1" => "Uses #define",
		     "920302-1" => "Uses #if",
		     "920410-1" => "Uses #define",
		     "920415-1" => "Uses #if",
		     "920428-2" => "Uses #if",
		     "920501-3" => "Uses #if",
		     "920501-4" => "Uses #if",
		     "920501-5" => "Uses #if",
		     "920501-7" => "Uses #define",
		     "920612-2" => "Uses #if",
		     "920721-4" => "Uses #if",
		     "921007-1" => "Uses #define",
		     "921017-1" => "Uses #if",
		     "921113-1" => "Uses #define",
		     "921202-1" => "Uses #define",
		     "921208-2" => "Uses #define",
		     "921215-1" => "Uses #if",
		     "930106-1" => "Uses #define",
		     "930406-1" => "Uses #if",
		     "931002-1" => "Uses #if",
		     "950221-1" => "Uses #if",
		     "960311-1" => "Uses #define",
		     "960311-2" => "Uses #define",
		     "960311-3" => "Uses #define",
		     "960405-1" => "Uses #define",
		     "960416-1" => "Uses #define",
		     "960419-2" => "Uses #define",
		     "960521-1" => "Uses #define",
		     "960830-1" => "Uses #if",
		     "970214-1" => "Uses #define",
		     "970214-2" => "Uses #define",
		     "980526-1" => "Uses #if",
		     "980605-1" => "Uses #define",
		     "981001-1" => "Uses #define",
		     "990128-1" => "Uses #define",
		     "990208-1" => "Uses #if",
		     "990211-1" => "Uses #define",
		     "990923-1" => "Uses #define",
		     "991014-1" => "Uses #define",
		     "991216-1" => "Uses #define",
		     "991216-2" => "Uses #define",
		     "991216-3" => "Uses #define",
		     "991228-1" => "Uses #define",
		     "alloca-1" => "Uses #define",
		     "arith-rand-ll" => "Uses #define",
		     "arith-rand" => "Uses #define",
		     "ashldi-1" => "Uses #define",
		     "ashrdi-1" => "Uses #define",
		     "bcp-1" => "Uses #define",
		     "builtin-abs-1" => "Uses #define",
		     "builtin-abs-2" => "Uses #define",
		     "builtin-complex-1" => "Uses #if",
		     "builtin-constant" => "Uses #define",
		     "builtin-noret-1" => "Uses #if",
		     "builtin-prefetch-1" => "Uses #define",
		     "builtin-prefetch-4" => "Uses #define",
		     "builtin-prefetch-6" => "Uses #define",
		     "cbrt" => "Uses #if",
		     "cmpdi-1" => "Uses #define",
		     "comp-goto-2" => "Uses #define",
		     "compare-3" => "Uses #if",
		     "conversion" => "Uses #if",
		     "comp-goto-1" => "Uses #if",
		     "complex-6" => "Uses #define",
		     "eeprof-1" => "Uses #define",
		     "ffs-2" => "Uses #define",
		     "loop-13" => "Uses #define",
		     "loop-2f" => "Uses #define",
		     "loop-2g" => "Uses #define",
		     "lshrdi-1" => "Uses #define",
		     "memcpy-1" => "Uses #define",
		     "memcpy-2" => "Uses #define",
		     "memcpy-bi" => "Uses #define",
		     "memset-1" => "Uses #define",
		     "memset-2" => "Uses #define",
		     "memset-3" => "Uses #define",
		     "nestfunc-4" => "Uses #define",
		     "nestfunc-1" => "Uses #if",
		     "nestfunc-2" => "Uses #if",
		     "nestfunc-3" => "Uses #if",
		     "pure-1" => "Uses #if",
		     "stdio-opt-1" => "Uses #if",
		     "stdio-opt-2" => "Uses #if",
		     "stdio-opt-3" => "Uses #if",
		     "strcmp-1" => "Uses #define",
		     "strcpy-1" => "Uses #define",
		     "string-opt-1" => "Uses #if",
		     "string-opt-10" => "Uses #if",
		     "string-opt-11" => "Uses #if",
		     "string-opt-12" => "Uses #if",
		     "string-opt-13" => "Uses #if",
		     "string-opt-14" => "Uses #if",
		     "string-opt-15" => "Uses #if",
		     "string-opt-16" => "Uses #if",
		     "string-opt-17" => "Uses #if",
		     "string-opt-2" => "Uses #if",
		     "string-opt-3" => "Uses #if",
		     "string-opt-4" => "Uses #if",
		     "string-opt-6" => "Uses #if",
		     "string-opt-7" => "Uses #if",
		     "string-opt-8" => "Uses #if",
		     "string-opt-9" => "Uses #if",
		     "strlen-1" => "Uses #define",
		     "strncmp-1" => "Uses #define",
		     "tstdi-1" => "Uses #define",
		     "va-arg-8" => "Uses #if",
		     "va-arg-10" => "Uses #define",
		     "va-arg-21" => "Uses #define",
		     "va-arg-22" => "Uses #define",
		     "widechar-1" => "Uses #define",
		     "alpha" => "beta" 
                    );

foreach $file (@files) {
  my($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
  next if ("$suffix" ne ".c");
  if (defined $gcc_exceptions{$filename}) {
      print "$filename skipped : $gcc_exceptions{$filename}\n";
      $expected_failure_count++;
      next;
  }
  if (defined $tcc_exceptions{$filename}) {
      print "$filename skipped : $tcc_exceptions{$filename}\n";
      $expected_failure_count++;
      next;
  }
  system("perl ./gen_tests.pl $file /tmp > /dev/null");
  my $return = system("cd /tmp ; gcc -o $filename $filename.c -I$BUILD_PATH/include -L$BUILD_PATH/lib -lffs  >& /dev/null");
  if (($return >> 8) != 0) {
      print "$filename compilation failed\n";
      $compile_failed_count += 1;
      next;
  }
  $return = system("/tmp/$filename >& /dev/null");
  if (($return >> 8) != 0) {
      print "$filename test failed, return code $return\n";
      $execute_failed_count += 1;
  } else {
      print "$filename test SUCCEEDED\n";
      $execute_success_count += 1;
  }
}

print "\n\n$execute_success_count tests succeeded!\n$compile_failed_count tests failed to compile, $execute_failed_count test failed upon execution\n$expected_failure_count tests not run (expected failures)\n";
