#perl
use strict;
use warnings;
use Text::Balanced qw/extract_multiple extract_bracketed/;
use File::Basename;
use Getopt::Std;
my %options=();
getopts("v",\%options);

print "Verbose\n" if defined $options{v};
print "-v $options{v}\n" if defined $options{v};
sub parse_c_test($);
sub generate_cod_test($$$);

my ($filename, $outputdir )  =  @ARGV;

my $outputname = basename($filename);
my $sourcedir = dirname($filename);

my ($decls, $subroutines) = parse_c_test($filename) ;
generate_cod_test($outputdir . "/" . $outputname , $decls, $subroutines);

sub generate_cod_test($$$) 
{
    my ($filename, $decls, $subroutines) = @_;

    foreach my $sub (@$subroutines) {
        my ($name, $decl, $body) = @$sub;
	print "$name,   $decl,  $body\n" if ($options{v});
    }

    unless (open (INT, ">$outputdir/$outputname")) { die "Failed to open $outputdir/$outputname";}
print INT<<EOF;
#include "cod.h"
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <setjmp.h>

int exit_value = 0; /* success */
jmp_buf env;

void
my_abort()
{
    exit_value = 1; /* failure */
    longjmp(env, 1);
}

void
test_exit(int value)
{
    if (value != 0) {
	exit_value = value;
    }
    longjmp(env, 1);
}

FILE *test_output = NULL;
int verbose = 0;

int test_printf(const char *format, ...)
{
    int ret;
    va_list args;
    va_start(args, format);

    if (test_output == NULL) {
	test_output = fopen("$filename.output", "w");
    }
    if (verbose) vprintf(format, args);
    ret = vfprintf(test_output, format, args);

    va_end(args);
    return ret;
}

int
main(int argc, char**argv)
{
    while (argc > 1) {
	if (strcmp(argv[1], "-v") == 0) {
	    verbose++;
        }
	argc--; argv++;
    }
EOF
    my $externs_table= "    cod_extern_entry externs[] = \n    {\n";
    my $externs_string = "    char extern_string[] = \"\\n\\\n";
    my $function_bodies = "    char *func_bodies[] = {\n";
    my $function_decls = "    char *func_decls[] = {\n";
    my $global_decls = "    char *global_decls[] = {\n";
    foreach my $sub (@$subroutines) {
        my ($name, $decl, $body) = @$sub;
	chomp($decl);
	$decl =~ s/\n/ /s;
	$externs_table .= "	{\"$name\", (void*)(long)-1},\n";
	$externs_string .= "	$decl;\\n\\\n";
	$body =~ s/\\/\\\\/g;
	$body =~ s/\n/\\n\\\n/g;
	$body =~ s/"/\\"/g;
	$body =~ s/(\W)printf\(/$1test_printf\(/g;
	$body =~ s/ = NULL/ = \(void*\)0/g;
	$function_bodies .= "\n/* body for $name */\n\"$body\",\n";
	$function_decls .= "\t\"$decl;\",\n";
    }
    $externs_table .= "	{\"abort\", (void*)my_abort},\n	{\"exit\", (void*)test_exit},\n	{\"test_printf\", (void*)test_printf},\n	{(void*)0, (void*)0}\n    };\n";
    $externs_string .= "    	void exit(int value);\\n\\\n        void abort();\\n\\\n        int test_printf(const char *format, ...);\";";
    $externs_string =~ s/\(void\)/\(\)/g;
    $function_bodies .= "\"\"};\n";
    $function_decls .= "\t\"\"};\n";
    $function_decls =~ s/\(void\)/\(\)/g;
    foreach my $decl (@$decls) {
	chomp($decl);
	$decl =~ s/\n/\\n\\\n/g;
	$decl =~ s/\n/ /g;
	$decl =~ s/\\n\\ /\\n\\\n/g;
	$global_decls .=  "\t\"$decl\",\n";
    }
    $global_decls .= "\"\"};\n";
    print INT "$externs_table\n";
    print INT "$externs_string\n";
    print INT "$function_bodies\n";
    print INT "$function_decls\n";
    print INT "$global_decls\n";
    print INT "    int i;\n";
    print INT "    cod_code gen_code[". scalar @$subroutines. "];\n";
    print INT "    for (i=0; i < ". scalar @$subroutines."; i++) {\n";
    print INT "        int j;\n";
    print INT "        if (verbose) {\n";
    print INT "             printf(\"Working on subroutine %s\\n\", externs[i].extern_name);\n";
    print INT "        }\n";
    print INT "        cod_parse_context context = new_cod_parse_context();\n";
    print INT "        cod_assoc_externs(context, externs);\n";
    print INT "        for (j=0; j < " . scalar @$decls . "; j++) {\n";
    print INT "            cod_parse_for_globals(global_decls[j], context);\n";
    print INT "        }\n";
    print INT "        cod_parse_for_context(extern_string, context);\n";
    print INT "        cod_subroutine_declaration(func_decls[i], context);\n";
    print INT "        gen_code[i] = cod_code_gen(func_bodies[i], context);\n";
    print INT "        externs[i].extern_value = (void*) gen_code[i]->func;\n";
    print INT "        if (i == ". scalar @$subroutines - 1 . ") {\n";
    print INT "            int (*func)() = (int(*)()) externs[i].extern_value;\n";
    print INT "            if (setjmp(env) == 0) {\n";
    print INT "                func();\n";
    print INT "            }\n";
    print INT "            if (exit_value != 0) {\n";
    print INT "                printf(\"Test $filename failed\\n\");\n";
    print INT "                exit(exit_value);\n";
    print INT "            }\n";
    print INT "        }\n";
    print INT "    }\n";
    print INT "    if (test_output) {\n";
    print INT "        /* there was output, test expected */\n";
    my $expectedname = "$sourcedir/$outputname";
    $expectedname =~ s/\.c/\.expect/g;
    print INT "        fclose(test_output);\n";
    print INT "        int ret = system(\"cmp $filename.output $expectedname\");\n";
    print INT "        ret = ret >> 8;\n";
    print INT "        if (ret == 1) {\n";
    print INT "            printf(\"Test $filename failed, output differs\\n\");\n";
    print INT "            exit(1);\n";
    print INT "        }\n";
    print INT "        if (ret != 0) {\n";
    print INT "            printf(\"Test $filename failed, output missing\\n\");\n";
    print INT "            exit(1);\n";
    print INT "        }\n";
    print INT "    }\n";
    print INT "    if (verbose) printf(\"Test $filename Succeeded\\n\");\n";
    print INT "    return 0;\n";
    print INT "}\n";
}

sub parse_c_test($) {
    my ($filename) = @_;
    local $_;
    my $file;
    
    my @decls = ();
    my @subroutines = ();
    open my $fileHandle, '<', $filename;
    
    { 
	local $/ = undef; # or use File::Slurp
	$file = <$fileHandle>;
    }
    
    close $fileHandle;
    
    my @array = extract_multiple(
				 $file,
				 [ sub{extract_bracketed($_[0], '{}')},],
				 undef,
				 0
				);
    
    
    while (@array) {
	my $non_bracket_item = shift(@array);
	$non_bracket_item =~ s/^\s+//;
	next if ($non_bracket_item eq "");
	
	next if (substr($non_bracket_item, 0, 2) eq "//");
	if (substr($non_bracket_item, 0, 2) eq "/*") {
	  $non_bracket_item =~ s/.*\*\///g;
	  $non_bracket_item =~ s/^\s+//;
	  next if ($non_bracket_item eq "");
	}
	if (substr($non_bracket_item, 0, 1) eq "{") {
	    print "Failure in item {$non_bracket_item }\n";
	    next;
	}
	my @bits = split(/\n\n/, $non_bracket_item);
	my $last_segment = pop @bits;
	if (@bits) {
	    # everything prior to the last empty lines must be declarations
	    my $item = join ("\n\n", @bits);
	    $item =~ s/\s+$//;
	    next if ($item eq "");
	    print "This is a 1 declaration set [\n $item\n]\n\n\n" if ($options{v});
	    push @decls, $item;
	}
	#
	@bits = split(/ /, $last_segment);
	my $last = pop @bits;
	my $before_last = "";
	if (@bits) {
	    $before_last = pop @bits;
	    push @bits, $before_last;
	}
	push @bits, $last;
	
	if (($before_last eq "struct") || ($before_last eq "enum")) {
	    my $decl_set = $last_segment . shift(@array);
	    my $next = shift(@array);
	    if ($next =~ m/(.*;)(.*)/) { 
	        $decl_set .= $1;
		$next = $2;
	    }
	    unshift @array, $next;
	    print "This is a declaration set [\n $decl_set\n]\n\n\n" if ($options{v});
	    push @decls, $decl_set;
	    next;
	}
	$last_segment =~ s/\s+$//;
	if (substr($last_segment, -1) eq ")") {
	    my $subroutine_name;
	    if ($last_segment =~ /^\W*\w+\W*$/)   {
	      $last_segment = "void $last_segment";
	    }
	    if ($last_segment =~ /(\w+)\W*\(/)   {
		$subroutine_name = $1;
	    }
	    print "Ha, must be a subroutine 1, name $subroutine_name, profile $last_segment\n" if ($options{v});
	    push @subroutines, [($subroutine_name, $last_segment, shift(@array))];
	    next;
	}
	# Well, maybe K&R style decls.
	my @lines = split (/\n/, $last_segment);
	my $count = 0;
	my $last_line;
	while (@lines) {
	    $last_line = pop @lines;
	    $last_line =~ s/\s+$//;
	    if (substr($last_line, -1) eq ";") {
		$count++;
		next;
	    }
	    if (substr($last_line, -1) ne ")") {
		print "Didn't get it, last line $last_line\n";
	    }
	    while(@lines) {
		$last_line = pop(@lines) . "\n" . $last_line;
	    }
	}
	my ($subroutine_prefix, $params) = split /\(/, $last_line;
	if (!defined $params) {
	  print "This is a declaration set [\n $subroutine_prefix\n]\n\n\n" if ($options{v});
	  push @decls, $subroutine_prefix;
	  next;
	}
	my @params = split(/,/, $params);
	if (($count != 0) && (@params != $count)) {
	    print "Param count didn't match " . scalar @params . " vs $count\n";
	}
	my $subroutine_name;
	if ($subroutine_prefix =~ /^\W*\w+\W*$/)   {
	    $subroutine_prefix = "void $subroutine_prefix";
	}
	if ($subroutine_prefix =~ /(\w+)\W*$/)   {
	    $subroutine_name = $1;
	}
	$subroutine_prefix =~ s/\n/ /g;
	my $subroutine_header = $subroutine_prefix . "(";
	@lines = split (/\n/, $last_segment);
	while ($count) {
	    my $param = pop @lines;
	    $param =~ s/;\s*$//;
	    $param =~ s/^\s*//;
	    $subroutine_header .= "$param";
	    if ($count != 1) {
		$subroutine_header .= ", ";
	    }
	    $count--;
	}
	$subroutine_header .= ")";
	print "Ha, must be a subroutine 2, name $subroutine_name, profile $subroutine_header\n" if ($options{v});
	push @subroutines, [($subroutine_name, $subroutine_header, shift(@array))];
	next;
	
    }  
    return (\@decls, \@subroutines);
}

