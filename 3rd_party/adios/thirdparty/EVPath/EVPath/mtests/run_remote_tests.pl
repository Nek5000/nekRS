#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $target;
my $hostname = `hostname`;
my $build;
my $quiet = 0;
my (%build_for_host, %source_dir, %build_name, %external_hostname, %local_port_range, %external_ssh_arg);
sub load_db();
sub build_env_scripts($);

load_db();
chomp($hostname);

die "No build for this host" unless defined $build_for_host{$hostname};

$build = $build_for_host{$hostname};
die "Usage: run_remote_tests.pl <target>" unless defined $ARGV[0];
if ($ARGV[0] eq "-build_scripts") {
    build_env_scripts($ARGV[1]);
    exit(0);
}
if ($ARGV[0] eq "-q") {
    shift(@ARGV);
    $quiet = 1;
}
$target = $ARGV[0];
die "No build for target \"$target\"" unless defined $source_dir{$target};
print "Using build $build on host $hostname\n" unless $quiet;
opendir(DIR, "$source_dir{$build}");
my @files = grep(/\.tsts$/,readdir(DIR));
closedir(DIR);
my $tsts_dir = "$source_dir{$build}";

if (! @files ) {
    opendir(DIR, "$source_dir{$build}/..");
    @files = grep(/\.tsts$/,readdir(DIR));
    closedir(DIR);
    $tsts_dir = "$source_dir{$build}/..";
}

open CTEST_CONFIG, ">$source_dir{$build}/CTestConfig.cmake" or die $1;
my $dartsubmitinfo = <<'END';
SET (CTEST_NIGHTLY_START_TIME "10:00:00 GMT")
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

SET(CTEST_PROJECT_NAME "Evpath")
SET(CTEST_NIGHTLY_START_TIME "10:49:51 GMT")

IF(NOT DEFINED CTEST_DROP_METHOD)
  SET(CTEST_DROP_METHOD "http")
ENDIF(NOT DEFINED CTEST_DROP_METHOD)

IF(CTEST_DROP_METHOD STREQUAL "http")
  SET(CTEST_DROP_SITE "evpath.net")
  SET(CTEST_DROP_LOCATION "/CDash/submit.php?project=Evpath")
  SET(CTEST_DROP_SITE_CDASH TRUE)
  SET(CTEST_TRIGGER_SITE "cdash.cercs.gatech.edu")
ENDIF(CTEST_DROP_METHOD STREQUAL "http")
END
print CTEST_CONFIG "$dartsubmitinfo\n";
close CTEST_CONFIG;

open CTEST, ">$source_dir{$build}/CTestTestfile.cmake" or die $1;

print CTEST "# \n";
print CTEST "# generated from files: ";
print CTEST join( ', ', @files );
print CTEST "\n";
print CTEST "# \n";
print CTEST "message (\"SSH params are \$ENV{SSH_PARAMS}\")" unless $quiet;

foreach my $file (@files) {
    open(TSTS, "$tsts_dir/$file") or die("Could not open  $tsts_dir/$file.");
    print CTEST "# \n";
    print CTEST "# tests from $file\n";
    print CTEST "# \n";
    foreach my $line (<TSTS>)  {
	chomp $line;
	print CTEST "ADD_TEST($line \"$line\" \"-ssh\" \"\$ENV{SSH_PARAMS}\")\n";
    }
}

close CTEST;

open TEST, ">$source_dir{$build}/$target.cmake" or die $1;
my $TEST_BODY = <<'END';
## -- Set hostname
## --------------------------
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
set(CTEST_SITE                          "${HOSTNAME}")
set(MODEL                               "Experimental")
set(CMAKE_C_COMPILER "/usr/bin/gcc")
SET(CTEST_NIGHTLY_START_TIME "11:56:44 GMT")
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)
set(CTEST_BINARY_DIRECTORY              "${CTEST_SOURCE_DIRECTORY}")

set(CTEST_TIMEOUT           "7200")
set( $ENV{LC_MESSAGES}      "en_EN" )
ctest_start(${MODEL} TRACK ${MODEL})

END
$TEST_BODY .= 'message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
' unless $quiet;
$TEST_BODY .= 'ctest_test(     BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
';

$TEST_BODY .= 'message(" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
' unless $quiet;
$TEST_BODY .= "ctest_submit(  RETURN_VALUE res)\n";
$TEST_BODY .= 'message(" -- Finished ${MODEL}  - ${CTEST_BUILD_NAME} --")
' unless $quiet;

print TEST "set(CTEST_SOURCE_DIRECTORY 	\"$source_dir{$build}\")\n";
print TEST "set(CTEST_BUILD_NAME 	\"$build_name{$build}-to-$build_name{$target}\")\n";
print TEST "SET( ENV{SSH_PARAMS}    \"$external_ssh_arg{$target}:$source_dir{$target}\" )\n";
if (defined $external_hostname{$build}) {
    print TEST "SET( ENV{CM_HOSTNAME} \"$external_hostname{$build}\")\n";
}
if (defined $local_port_range{$build}) {
    print TEST "SET( ENV{CM_PORT_RANGE} \"$local_port_range{$build}\")\n";
}
print TEST $TEST_BODY;
close TEST;

chdir ("$source_dir{$build}");
if ($quiet) {
    system("ctest -S $target.cmake");
} else {
    system("ctest -V -S $target.cmake");
}

sub build_env_scripts($) {
    my $exe_dir = shift;
    my $file;
    print "Building scripts for executable directory $exe_dir\n";
    mkdir("$exe_dir/scripts");
    chdir ("$exe_dir");
    @files = <*>;
    foreach $file (@files) {
	if (-x $file && ! -d $file) {
	    print $file . "\n";
	    open SCRIPT, ">scripts/$file";
	    print SCRIPT "#!/bin/bash\n";
	    if (defined $external_hostname{$build}) {
		print SCRIPT "export CM_HOSTNAME=$external_hostname{$build}\n";
	    }
	    if (defined $local_port_range{$build}) {
		print SCRIPT "export CM_PORT_RANGE=$local_port_range{$build}\n";
	    }
	    print SCRIPT "$exe_dir/$file \"\$\@\"\n";
	    close SCRIPT;
	    chmod 0755, "scripts/$file"
	}
    }
}

sub load_db() {
    %build_for_host = (
        'jedi080', 'jedi',
        'chaos', 'chaos',
	'raspberrypi', 'raspberry',
	);
    $source_dir{"jedi"} = "/users/c/chaos/evpath_tests/rhe6-64";
    $source_dir{"chaos"} = "/home/greg/evpath_tests/";
    $source_dir{"raspberry"} = "/home/eisen/evpath_tests/scripts";
    $build_name{"jedi"} = "Jedi";
    $build_name{"chaos"} = "UB18";
    $build_name{"raspberry"} = "Pi";
    $external_hostname{"raspberry"} = "eisenhauer.dyndns.org";
    $local_port_range{"raspberry"} = "62000:62100";
    $external_ssh_arg{"jedi"} = "chaos\@jedi080.cc.gatech.edu";
    $external_ssh_arg{"chaos"} = "greg\@chaos.cc.gatech.edu";
    $external_ssh_arg{"raspberry"} = "eisen\@eisenhauer.dyndns.org:4022";
}
