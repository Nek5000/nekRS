#!perl
sub upperc {
    local($_) = pop(@_);
    tr/[a-z]/[A-Z]/;
    return $_;
}

%c_types = ('c', 'signed char', 'uc', 'unsigned char', 's', 'short', 'us', 'unsigned short', 'i', 'int', 'u', 'unsigned int', 'l', 'intptr_t', 'ul', 'uintptr_t', 'p', 'char*', 'f', 'float', 'd', 'double');
%cod_types = ('c', 'signed char', 'uc', 'unsigned char', 's', 'short', 'us', 'unsigned short', 'i', 'int', 'u', 'unsigned int', 'l', 'long', 'ul', 'unsigned long', 'p', 'char*', 'f', 'float', 'd', 'double');
%rand_types = ('c', 'l', 'uc', 'l', 's', 'l', 'us', 'l', 'i', 'l', 'u', 'l', 'l', 'l', 'ul', 'l', 'p', '(intptr_t)l', 'f', 'd', 'd', 'd');
%c_print_formats = ('c', '%d', 'uc', '%u', 's', '%d', 'us', '%u', 'i', '%d', 'u', '%u', 'l', '%zd', 'ul', '%zu', 'p', '%lx', 'f', '%g', 'd', '%g');
%drisc_arg_formats = ('c', '%i', 'uc', '%u', 's', '%i', 'us', '%u', 'i', '%i', 'u', '%u', 'l', '%l', 'ul', '%ul', 'p', '%p', 'f', '%f', 'd', '%d');

%type_masks = ('c', "f", 'uc', "f", 's', "f ff fff", 'us', "f ff fff", 'i', "f ff fff ffff fffff ffffff fffffff", 'u', "f ff fff ffff fffff ffffff fffffff", 'l', "f ff fff ffff fffff ffffff fffffff", 'ul', "f ff fff ffff fffff ffffff fffffff");
%type_max = ('c', "ff", 'uc', "ff", 's', "ffff", 'us', "ffff", 'i', "ffffffff", 'u', "ffffffff", 'l', "ffffffffffffffff", 'ul', "ffffffffffffffff");
%type_max2 = ('c', "80", 'uc', "80", 's', "8000", 'us', "8000", 'i', "800000000", 'u', "80000000", 'l', "8000000000000000", 'ul', "8000000000000000");

%subr_name_list = ("%", "mod", "&", "and", "|", "or", "^", "not", "<<", "lsh", ">>", "rsh", "+", "add", "-", "sub", "*", "mult", "/", "div");

&output_header;

&arith_insn("% & | ^ << >>", "c uc s us i u ul l");
&arith_insn("+ - * /",  "c uc s us i u ul l f d");
#  &arithp_insn("add sub", "+ -");

#  &arith2_insn("com not", "~ !", "i u ul l");
#  &arith2_insn("neg", "-", "i u ul l f d");

#  &branch_insn( "eq ge gt le lt ne", "== >= > <= < !=", "i u ul l p d f");
#  &branchi_insn( "eq ge gt le lt ne", "== >= > <= < !=", "i u ul l p");

#  &convert( "c d f i l s u ul us", "i u ul l");
#  &convert( "d f i l u ul", "f d");
#  &convert( "i l u ul", "c s");
#  &convert( "p ul", "ul p");

#  &load("c uc s us i u ul l f d");
#  &store("c uc s us i u ul l f d");

#  &setmovret("c uc s us i u ul l f d");

print COUT<<EOF;

int
main(int argc, char **argv)
{
    int i = 1;

    initialize_values();

    while (i < argc) {
	if (strcmp(argv[i], "-v") == 0) {
	    verbose++;
	} else {
	    printf("Only running test \\"%s\\"\\n", argv[i]);
	    only_run_test = argv[i];
	}
	i++;
    }

    $main_body
EOF
print COUT "    return failed;\n}\n";


sub arith_insn {
    local ($c_ops, $type_list) = @_;
    @c_ops = split(' ', $c_ops);
    foreach(split(' ', $c_ops)) {
      $c_op = $_;
      $print_op = $c_op;
      ($print_op eq "%") && ($print_op = "%%");
      foreach(split(' ', $type_list)) {
	  $t1 = ${_};
          $subr_name = "$subr_name_list{$c_op}_$t1";
          print COUT "void ${subr_name}_tests()\n{\n";
          $main_body .= "\n    ${subr_name}_tests();\n";
          foreach(split(' ', $type_list)) {
	    $t2 = ${_};
	    $c_type1 = $c_types{$t1};
	    $cod_type1 = $cod_types{$t1};
	    $c_type2 = $c_types{$t2};
	    $cod_type2 = $cod_types{$t2};
	    $c_pformat1 = $c_print_formats{$t1};
	    $c_pformat2 = $c_print_formats{$t2};
	    $line = __LINE__ + 3;
	    $arg_str = "$drisc_arg_formats{$t1}$drisc_arg_formats{$t2}";
	    if ($c_op eq "/") {
		$div_continue = "if (source2_$t2 == 0) {if (verbose) {printf(\"+\");fflush(stdout);} goto skip$skip_count;}";
		$div_label = "skip$skip_count: ;\n";
		if (($c_type1 eq "short") && 
		    (($t2 eq "d") || ($t2 eq "f"))) {
		    $div_continue .= "\n	if ((source1_$t1/source2_$t2) < (-32768) || ((source1_$t1/source2_$t2) > 32767)){if (verbose) {printf(\"+\");fflush(stdout);} goto skip$skip_count;}";
		} 
		$skip_count++;
	    } elsif ($c_op eq "%") {
		$div_continue = "if (source2_$t2 <= 0) {if (verbose) {printf(\"+\");fflush(stdout);} goto skip$skip_count;}";
#		$div_continue .= "if (source2_$t2  > 0) goto skip$skip_count;";
		$div_label = "skip$skip_count: ;\n";
		$skip_count++;
	    } elsif ($c_op eq ">>") {
		$div_continue = "if (source2_$_ >= sizeof(source1_$t1)) goto skip$skip_count;";
		$div_label = "skip$skip_count: ;\n";
		$skip_count++;
	    } else {
		$div_continue = "";
		$div_label = "";
	    }
	    $op1 = "${c_op}$t1";
	    $op2 = "${c_op}$t2";
	    if (defined($c_src1_values{$op1})) {
		$vals1 = $c_src1_values{$op1};
	    } else {
		$vals1 = "src1${t1}_vals";
	    }
	    if (defined($c_src2_values{$op2})) {
		$vals2 = $c_src2_values{$op2};
	    } else {
		$vals2 = "src2${_}_vals";
	    }
            $result_if = "if (expected_result != result)";
	    $expected_value = "(source1_$t1 $c_op source2_$t2)";
	    $range_decl = "";
  	    if(($t1 eq "f") || ($t2 eq "f")) {
		$result_if = "if ((result > (expected_result + range)) || (result < (expected_result - range)))";
		if (($t1 ne "f") || ($t2 ne "f")) {
		    $result_if = "range += 1.0;\n	    $result_if";
		}
		if ((${c_op} eq "+") || (${c_op} eq "-")) {
		    $range_decl ="double range = 0.000001 * (fabs((double)source1_$t1) + fabs((double)source2_$t2));";
		} elsif  ((${c_op} eq "*") || (${c_op} eq "/")) {
		    $range_decl ="double range = 0.000001 * (fabs((double)source1_$t1 $c_op (double)source2_$t2));";
		}
	    }
	    if (($t1 eq "uc") || ($t1 eq "us") || ($t1 eq "u") || ($t1 eq "ul")) {
		$result_if = "if ((expected_result == 0) || (expected_result == ($c_type1) 0x$type_max{$t1}) || (expected_result == ($c_type1) 0x$type_max2{$t1})) continue;\n	$result_if";
	    }
	    if (($t1 eq "ul") && (($t2 eq "d") || ($t2 eq "f"))) {
		$result_if = "if (source1_ul > ((ssize_t)11<<52)) continue;\n 	$result_if";
	    }
	    if (($t1 eq "uc") && ($c_op eq "/")) {
		$range_decl .= "\n	    unsigned int expect_int;";
		$div_continue .= "\n	    expect_int = (unsigned int) $expected_value;";
		$expected_value = "expect_int";
	    }
print COUT<<EOF;
  if ((only_run_test == NULL) || 
      (strcmp(only_run_test, "${t1}${c_op}${t2}") == 0)) {
     int i, j;
     static char code[] = \"{\\n\\
       $cod_type1 a = p1;\\n\\
       $cod_type2 b = p2;\\n\\
       return a $c_op b;\\n\\
     }\";

    $c_type1 (*proc)($c_type1 a, $c_type2 b);
    cod_parse_context context = new_cod_parse_context();
    int dumped = 0;

    cod_code gen_code;

    /* test for $c_type1 $c_op $c_type2 */
    if (verbose) printf(\"test for $c_type1 $print_op $c_type2\");

    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);

    cod_subroutine_declaration("$cod_type1 proc($cod_type1 p1, $cod_type2 p2)",
			       context);
    gen_code = cod_code_gen(code, context);
    proc = ($c_type1 (*)($c_type1, $c_type2))(intptr_t)gen_code->func;


    for (i=0 ; i < sizeof($vals1)/sizeof($vals1\[0\]) ; i++) {
        $c_type1 source1_$t1 = $vals1\[i\];
        for (j=0 ; j < sizeof($vals2)/sizeof($vals2\[0\]) ; j++) {
            $c_type2 source2_$t2 = $vals2\[j\];
	    $range_decl
	    $c_type1 expected_result;
	    $c_type1 result;

	    $div_continue
	    if (verbose) {printf(".");fflush(stdout);}

	    expected_result = ($c_type1) $expected_value;

	    result = proc(source1_$t1, source2_$t2);
	    $result_if {
		printf("Failed $c_type1 $print_op $c_type2 test, expected $c_pformat1, got $c_pformat1, for $c_pformat1 $print_op $c_pformat2\\n",
		       expected_result, result, source1_$t1, source2_$t2);
		if (dumped == 0) cod_dump(gen_code);
		failed++;
		dumped++;
	    }
$div_label	}
    }
    if (verbose) printf(" done\\n");
    cod_free_parse_context(context);
    cod_code_free(gen_code);
  }
EOF
	}
        print COUT "}\n";
      }
    }
}

sub arithp_insn {
    local ($dr_ops, $c_ops) = @_;
    $type_list = "p";
    @c_ops = split(' ', $c_ops);
    foreach(split(' ', $dr_ops)) {
	$dr_op = $_;
	$c_op = shift(@c_ops);
	
	foreach(split(' ', $type_list)) {
	    $dr_type = "DR_".&upperc(${_});
	    $c_type = $c_types{${_}};
	    $c_pformat = $c_print_formats{${_}};
	    $line = __LINE__ + 3;
	    $arg_str = "$drisc_arg_formats{${_}}%l";
	    $full_op = "${dr_op}${_}";
	    print COUT "\n    /* test for dr_$dr_op${_} */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$dr_op${_}\");\n";
	    if (defined($c_src1_values{$full_op})) {
		$vals = $c_src1_values{$full_op};
		print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
		print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	    } else {
		print COUT "    {\n";
		print COUT "        $c_type source1_$_ = rand1_$_;\n";
	    }
	    if (defined($c_src2_values{"addl"})) {
		$vals = $c_src2_values{'addl'};
		print COUT "        for (j=0 ; j < sizeof($vals)/sizeof($vals\[0\]) ; j++) {\n";
		print COUT "            long source2_l = $vals\[j\];\n";
	    } else {
		print COUT "        {\n";
		print COUT "            long source2_l = rand2_$_;\n";
	    }
print COUT<<EOF;

	    dr_reg_t opnd1, opnd2, dest;
	    $c_type result;
	    $c_type expected_result;
	    $c_type (*proc)($c_type a, long b);

	    if (verbose) {printf(".");fflush(stdout);}
	    $div_continue
	    dr_proc_params(c, "no name", "$arg_str");
	    opnd1 = dr_param(c, 0);
	    opnd2 = dr_param(c, 1);
	    dest = dr_getvreg(c, $dr_type);
	    dr_$dr_op${_}(c, dest, opnd1, opnd2);
	    dr_ret${_}(c, dest);
	    
	    proc = ($c_type(*)($c_type, long)) dr_end(c);
	    expected_result = ($c_type) (source1_$_ $c_op source2_l);

	    result = proc(source1_$_, source2_l);
	    if (expected_result != result) {
		printf("Failed dr_$dr_op${_} test, expected $c_pformat, got $c_pformat, for $c_pformat $c_op %l\\n",
		       expected_result, result, source1_$_, source2_l);
		dr_dump(c);
		failed++;
	    }
$div_label	}
    }
    if (verbose) printf(" done\\n");
EOF
	}
    }
}

sub arith2_insn {
    local ($dr_ops, $c_ops, $type_list) = @_;
    @c_ops = split(' ', $c_ops);
    foreach(split(' ', $dr_ops)) {
	$dr_op = $_;
	$c_op = shift(@c_ops);
	
	foreach(split(' ', $type_list)) {
	    $dr_type = "DR_".&upperc(${_});
	    $c_type = $c_types{${_}};
	    $c_pformat = $c_print_formats{${_}};
	    $arg_str = "$drisc_arg_formats{${_}}";
	    $full_op = "${dr_op}${_}";
	    print COUT "\n    /* test for dr_$dr_op${_} */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$dr_op${_}\");\n";
	    if (defined($c_src1_values{$full_op})) {
		$vals = $c_src1_values{$full_op};
		print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
		print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	    } else {
		print COUT "    {\n";
		print COUT "        $c_type source1_$_ = rand1_$_;\n";
	    }
		print COUT "        {\n";
print COUT<<EOF;

	    dr_reg_t opnd1, dest;
	    $c_type result;
	    $c_type expected_result;
	    $c_type (*proc)($c_type a);

	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "$arg_str");
	    opnd1 = dr_param(c, 0);
	    dest = dr_getvreg(c, $dr_type);
	    dr_$dr_op${_}(c, dest, opnd1);
	    dr_ret${_}(c, dest);
	    
	    proc = ($c_type(*)($c_type)) dr_end(c);
	    expected_result = $c_op source1_$_;

	    result = proc(source1_$_);
	    if (expected_result != result) {
		printf("Failed dr_$dr_op${_} test, expected $c_pformat, got $c_pformat, for $c_op $c_pformat\\n",
		       expected_result, result, source1_$_);
		dr_dump(c);
		failed++;
	    }
	}
    }
    if (verbose) printf(" done\\n");
EOF
	}
    }
}

sub arithpi_insn {
    local ($dr_ops, $c_ops, $type_list) = @_;
    @c_ops = split(' ', $c_ops);
    foreach(split(' ', $dr_ops)) {
	$dr_op = $_;
	$c_op = shift(@c_ops);
	
	foreach(split(' ', $type_list)) {
	    $dr_type = "DR_".&upperc(${_});
	    $c_type = $c_types{${_}};
	    $c_pformat = $c_print_formats{${_}};
	    $line = __LINE__ + 3;
	    $arg_str = "$drisc_arg_formats{${_}}$drisc_arg_formats{${_}}";
	    $full_op = "${dr_op}${_}i";
	    print COUT "\n    /* test for dr_$dr_op${_}i (immediate) */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$dr_op${_}i (immediate)\");\n";
	    if (defined($c_src1_values{$full_op})) {
		$vals = $c_src1_values{$full_op};
		print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
		print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	    } else {
		print COUT "    {\n";
		print COUT "        $c_type source1_$_ = rand1_$_;\n";
	    }
	    if (defined($c_src2_values{'addl'})) {
		$vals = $c_src2_values{'addl'};
		print COUT "        for (j=0 ; j < sizeof($vals)/sizeof($vals\[0\]) ; j++) {\n";
		print COUT "            long source2_l = $vals\[j\];\n";
	    } else {
		print COUT "        {\n";
		print COUT "            long source2_l = rand2_$_;\n";
	    }
print COUT<<EOF;

	    dr_reg_t opnd1, dest;
	    $c_type result;
	    $c_type expected_result;
	    $c_type (*proc)($c_type a);

	    if (verbose) {printf(".");fflush(stdout);}
	    $div_continue
	    dr_proc_params(c, "no name", "$arg_str");
	    opnd1 = dr_param(c, 0);
	    dest = dr_getvreg(c, $dr_type);
	    dr_$dr_op${_}i(c, dest, opnd1, source2_l);
	    dr_ret${_}(c, dest);
	    
	    proc = ($c_type(*)($c_type)) dr_end(c);
	    expected_result = ($c_type) (source1_$_ $c_op source2_l);

	    result = proc(source1_$_);
	    if (expected_result != result) {
		printf("Failed dr_$dr_op${_}i (immediate) test, expected $c_pformat, got $c_pformat, for $c_pformat $c_op %d\\n",
		       expected_result, result, source1_$_, source2_l);
		dr_dump(c);
		failed++;
	    }
$div_label	}
    }
    if (verbose) printf(" done\\n");
EOF
	}
    }
}

sub branch_insn {
    local ($dr_ops, $c_ops, $type_list) = @_;
    @c_ops = split(' ', $c_ops);
    foreach(split(' ', $dr_ops)) {
	$dr_op = $_;
	$c_op = shift(@c_ops);
	
	foreach(split(' ', $type_list)) {
	    $dr_type = "DR_".&upperc(${_});
	    $c_type = $c_types{${_}};
	    $c_pformat = $c_print_formats{${_}};
	    $line = __LINE__ + 3;
	    $arg_str = "$drisc_arg_formats{${_}}$drisc_arg_formats{${_}}";
	    $full_op = "b${dr_op}${_}";
	    print COUT "\n    /* test for dr_$full_op */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$full_op\");\n";
	    if (defined($c_src1_values{"br${_}"})) {
		$vals = $c_src1_values{"br${_}"};
		print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
		print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	    } else {
		print COUT "    {\n";
		print COUT "        $c_type source1_$_ = rand1_$_;\n";
	    }
	    if (defined($c_src2_values{"br${_}"})) {
		$vals = $c_src2_values{"br${_}"};
		print COUT "        for (j=0 ; j < sizeof($vals)/sizeof($vals\[0\]) ; j++) {\n";
		print COUT "            $c_type source2_$_ = $vals\[j\];\n";
	    } else {
		print COUT "        {\n";
		print COUT "            $c_type source2_$_ = rand2_$_;\n";
	    }
print COUT<<EOF;

	    dr_reg_t opnd1, opnd2;
	    int label;
	    int result;
	    int expected_result;
	    int (*proc)($c_type a, $c_type b);

	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "$arg_str");
	    opnd1 = dr_param(c, 0);
	    opnd2 = dr_param(c, 1);
	
	    label = dr_genlabel(c);
	    dr_$full_op(c, opnd1, opnd2, label);
	    dr_retii(c, 0);
	    dr_label(c, label);	    
	    dr_retii(c, 1);
	    proc = (int(*)($c_type, $c_type)) dr_end(c);
	    expected_result = source1_$_ $c_op source2_$_;

	    result = proc(source1_$_, source2_$_);
	    if (expected_result != result) {
		printf("Failed dr_$full_op test, expected %d, got %d, for $c_pformat $c_op $c_pformat\\n",
		       expected_result, result, source1_$_, source2_$_);
		dr_dump(c);
		failed++;
	    }
	}
    }
    if (verbose) printf(" done\\n");
EOF
	}
    }
}

sub branchi_insn {
    local ($dr_ops, $c_ops, $type_list) = @_;
    @c_ops = split(' ', $c_ops);
    foreach(split(' ', $dr_ops)) {
	$dr_op = $_;
	$c_op = shift(@c_ops);
	
	foreach(split(' ', $type_list)) {
	    $dr_type = "DR_".&upperc(${_});
	    $c_type = $c_types{${_}};
	    $c_pformat = $c_print_formats{${_}};
	    $line = __LINE__ + 3;
	    $arg_str = "$drisc_arg_formats{${_}}";
	    $full_op = "b${dr_op}${_}i";
	    print COUT "\n    /* test for dr_$full_op */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$full_op\");\n";
	    if (defined($c_src1_values{"br${_}"})) {
		$vals = $c_src1_values{"br${_}"};
		print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
		print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	    } else {
		print COUT "    {\n";
		print COUT "        $c_type source1_$_ = rand1_$_;\n";
	    }
	    if (defined($c_src2_values{"br${_}"})) {
		$vals = $c_src2_values{"br${_}"};
		print COUT "        for (j=0 ; j < sizeof($vals)/sizeof($vals\[0\]) ; j++) {\n";
		print COUT "            $c_type source2_$_ = $vals\[j\];\n";
	    } else {
		print COUT "        {\n";
		print COUT "            $c_type source2_$_ = rand2_$_;\n";
	    }
	    $cast = "";
	    if ($_ eq 'p') { $cast = "(long)"; }
print COUT<<EOF;

	    dr_reg_t opnd1;
	    int label;
	    int result;
	    int expected_result;
	    int (*proc)($c_type a);

	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "$arg_str");
	    opnd1 = dr_param(c, 0);
	
	    label = dr_genlabel(c);
	    dr_$full_op(c, opnd1, $cast source2_$_, label);
	    dr_retii(c, 0);
	    dr_label(c, label);	    
	    dr_retii(c, 1);
	    proc = (int(*)($c_type)) dr_end(c);
	    expected_result = source1_$_ $c_op source2_$_;

	    result = proc(source1_$_);
	    if (expected_result != result) {
		printf("Failed dr_$full_op test, expected %d, got %d, for $c_pformat $c_op $c_pformat\\n",
		       expected_result, result, source1_$_, source2_$_);
		dr_dump(c);
		failed++;
	    }
	}
    }
    if (verbose) printf(" done\\n");
EOF
	}
    }
}

sub convert {
    local($from_types, $to_types) = @_;
    foreach (split(' ', $from_types)) {
	$from = $_;
	foreach (split(' ', $to_types)) {
	    if (${from} eq ${_}) { 
	        next;
	    } 

            $to = $_;
	    $dr_to_type = "DR_".&upperc(${_});
	    $dr_from_type = "DR_".&upperc(${_});
	    $c_to_type = $c_types{${to}};
	    $c_from_type = $c_types{${from}};
	    $c_to_pformat = $c_print_formats{${to}};
	    $c_from_pformat = $c_print_formats{${from}};
	    $line = __LINE__ + 3;
	    $arg_str = "$drisc_arg_formats{${from}}";
	    $full_op = "cv${from}2${to}";
	    print COUT "\n    /* test for dr_$full_op */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$full_op\");\n";
	    $vals = "src1${from}_vals";
	    print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
	    print COUT "        $c_from_type source1_$from = $vals\[i\];\n";
print COUT<<EOF;
	{

	    dr_reg_t opnd1, dest;
	    $c_to_type result;
	    $c_to_type expected_result;
	    $c_to_type (*proc)($c_from_type a);

	    if (verbose) {printf(".");fflush(stdout);}
	    $div_continue
	    dr_proc_params(c, "no name", "$arg_str");
	    opnd1 = dr_param(c, 0);
	    dest = dr_getvreg(c, $dr_to_type);
	    dr_$full_op(c, dest, opnd1);
	    dr_ret${to}(c, dest);
	    
	    proc = ($c_to_type(*)($c_from_type)) dr_end(c);
	    expected_result = ($c_to_type) source1_${from};

	    result = proc(source1_$from);
	    if (expected_result != result) {
		printf("Failed dr_$full_op test, expected $c_to_pformat, got $c_to_pformat, for ($c_to_type) $c_from_pformat\\n",
		       expected_result, result, source1_$from);
		dr_dump(c);
		failed++;
	    }
	}
    }
    if (verbose) printf(" done\\n");
EOF
	}
    }
}

sub load {
    local ($type_list) = @_;
    foreach(("", "i")) {
	$imm = $_;
	$dr_op = ld;
	foreach(split(' ', $type_list)) {
	    $dr_type = "DR_".&upperc(${_});
	    $c_type = $c_types{${_}};
	    $c_pformat = $c_print_formats{${_}};
	    $line = __LINE__ + 3;
	    $arg_str = "%p";
	    $full_op = "${dr_op}${_}$imm";
	    print COUT "\n    /* test for dr_$full_op */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$full_op\");\n";
	    $vals = "src1${_}_vals";
	    print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
	    print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	    if ($imm eq 'i') {
		$vals = "bit_pattern_vals";
	    } else {
		$vals = "src1l_vals";
	    }
	    print COUT "        for (j=0 ; j < sizeof($vals)/sizeof($vals\[0\]) ; j++) {\n";
	    print COUT "            long offset = $vals\[j\];\n";
	    if ($imm eq 'i') {
	        print COUT<<EOF;

	    dr_reg_t opnd1, dest;
	    $c_type result;
	    $c_type expected_result;
	    $c_type (*proc)(void *p);
            void *p = ((char*)&source1_$_) - offset;
	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "%p");
	    opnd1 = dr_param(c, 0);
	    dest = dr_getvreg(c, $dr_type);
	    dr_${full_op}(c, dest, opnd1, offset);
	    dr_ret${_}(c, dest);
	    
	    proc = ($c_type(*)(void *)) dr_end(c);
	    expected_result = source1_$_;

	    result = proc(p);
	    if (expected_result != result) {
		printf("Failed dr_$dr_op${_} test, expected $c_pformat, got $c_pformat, for loading $c_pformat\\n",
		       expected_result, result, source1_$_);
		dr_dump(c);
		failed++;
	    }
$div_label	}
    }
    if (verbose) printf(" done\\n");
EOF
	    } else {
	        print COUT<<EOF;

	    dr_reg_t opnd1, opnd2, dest;
	    $c_type result;
	    $c_type expected_result;
	    $c_type (*proc)(void *p);
            void *p = ((char*)&source1_$_) - offset;
	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "%p");
	    opnd1 = dr_param(c, 0);
	    dest = dr_getvreg(c, $dr_type);
	    opnd2 = dr_getvreg(c, DR_P);
	    dr_setp(c, opnd2, offset);
	    dr_${full_op}(c, dest, opnd1, opnd2);
	    dr_ret${_}(c, dest);
	    
	    proc = ($c_type(*)(void *)) dr_end(c);
	    expected_result = source1_$_;

	    result = proc(p);
	    if (expected_result != result) {
		printf("Failed dr_$dr_op${_} test, expected $c_pformat, got $c_pformat, for loading $c_pformat\\n",
		       expected_result, result, source1_$_);
		dr_dump(c);
		failed++;
	    }
$div_label	}
    }
    if (verbose) printf(" done\\n");
EOF
	    }
	}
    }
}
sub store {
    local ($type_list) = @_;
    foreach(("", "i")) {
	$imm = $_;
	$dr_op = st;
	foreach(split(' ', $type_list)) {
	    $dr_type = "DR_".&upperc(${_});
	    $c_type = $c_types{${_}};
	    $c_pformat = $c_print_formats{${_}};
	    $line = __LINE__ + 3;
	    $arg_str = "%p";
	    $full_op = "${dr_op}${_}$imm";
	    print COUT "\n    /* test for dr_$full_op */\n";
	    print COUT "    if (verbose) printf(\"test for dr_$full_op\");\n";
	    $vals = "src1${_}_vals";
	    print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
	    print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	    if ($imm eq 'i') {
		$vals = "bit_pattern_vals";
	    } else {
		$vals = "src1l_vals";
	    }
	    print COUT "        for (j=0 ; j < sizeof($vals)/sizeof($vals\[0\]) ; j++) {\n";
	    print COUT "            long offset = $vals\[j\];\n";
	    if ($imm eq 'i') {
	        print COUT<<EOF;

	    dr_reg_t opnd1, opnd2;
	    $c_type result = 0;
	    $c_type expected_result;
	    void (*proc)(void *p, $c_type a);
            void *p = ((char*)&result) - offset;
	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "%p$drisc_arg_formats{${_}}");
	    opnd1 = dr_param(c, 0);
	    opnd2 = dr_param(c, 1);
	    dr_${full_op}(c, opnd2, opnd1, offset);
	    dr_retii(c, 0);
	    
	    proc = (void (*)(void *, $c_type)) dr_end(c);
	    expected_result = source1_$_;

	    proc(p, source1_$_);
	    if (expected_result != result) {
		printf("Failed dr_$dr_op${_}$imm test, expected $c_pformat, got $c_pformat, for storing $c_pformat\\n",
		       expected_result, result, source1_$_);
		dr_dump(c);
		failed++;
	    }
$div_label	}
    }
    if (verbose) printf(" done\\n");
EOF
	    } else {
	        print COUT<<EOF;

	    dr_reg_t opnd1, opnd2, source;
	    $c_type result;
	    $c_type expected_result;
	    void (*proc)(void *p, $c_type a);
            void *p = ((char*)&result) - offset;
	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "%p$drisc_arg_formats{${_}}");
	    opnd1 = dr_param(c, 0);
	    source = dr_param(c, 1);
	    opnd2 = dr_getvreg(c, DR_P);
	    dr_setp(c, opnd2, offset);
	    dr_${full_op}(c, source, opnd1, opnd2);
	    dr_retii(c, 0);
	    
	    proc = (void (*)(void *, $c_type)) dr_end(c);
	    expected_result = source1_$_;

	    proc(p, source1_$_);
	    if (expected_result != result) {
		printf("Failed dr_$dr_op${_} test, expected $c_pformat, got $c_pformat, for storing $c_pformat\\n",
		       expected_result, result, source1_$_);
		dr_dump(c);
		failed++;
	    }
$div_label	}
    }
    if (verbose) printf(" done\\n");
EOF
	    }
	}
    }
}

sub setmovret {
    local ($type_list) = @_;
    foreach(split(' ', $type_list)) {
	$dr_type = "DR_".&upperc(${_});
	$c_type = $c_types{${_}};
	$c_pformat = $c_print_formats{${_}};
	$arg_str = "$drisc_arg_formats{${_}}";
	print COUT "\n    /* test for set/mov/ret ${_} */\n";
	print COUT "    if (verbose) printf(\"test for set/mov/ret ${_}\");\n";
	$vals = "src1${_}_vals";
	print COUT "    for (i=0 ; i < sizeof($vals)/sizeof($vals\[0\]) ; i++) {\n";
	print COUT "        $c_type source1_$_ = $vals\[i\];\n";
	print COUT "        {\n";
print COUT<<EOF;

	    dr_reg_t opnd1, dest;
	    $c_type result;
	    $c_type expected_result;
	    $c_type (*proc)();

	    if (verbose) {printf(".");fflush(stdout);}
	    dr_proc_params(c, "no name", "");
	    opnd1 = dr_getvreg(c, $dr_type);
	    dest = dr_getvreg(c, $dr_type);
	    dr_set${_}(c, opnd1, source1_$_);
	    dr_mov${_}(c, dest, opnd1);
	    dr_ret${_}(c, dest);
	    
	    proc = ($c_type(*)()) dr_end(c);
	    expected_result = source1_$_;

	    result = proc();
	    if (expected_result != result) {
		printf("Failed set/mov/ret ${_} test, expected $c_pformat, got $c_pformat\\n",
		       expected_result, result);
		failed++;
		dumped++;
	    }
	}
    }
    if (verbose) printf(" done\\n");
EOF
    }
}

sub output_header {
    $mach = virtual;
    open(COUT, ">general.c") || die "Can't open test output";
# print COUT "# ". (__LINE__ + 1) . " \"general.ops\"\n";
print COUT<<EOF;
/* This file is generated from general.ops.  Do not edit directly. */

#include "../config.h"
#include "../cod.h"
#include "stdio.h"
#include <stdint.h>
#include "string.h"
#include "math.h"

#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define srand48(s) srand(s)
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define kill(x,y) 
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif
#else
extern double drand48();
extern long lrand48();
void srand48(long seedval);
#endif

char* only_run_test = NULL;

static char extern_string[] = "int printf(string format, ...);";

static cod_extern_entry externs[2];

EOF
    foreach $dr_type (keys %c_types) {
        $c_type = $c_types{$dr_type};
        $rand_type = $rand_types{$dr_type};
	$_ = $type_masks{$dr_type};
	if (defined $type_masks{$dr_type}) {
	    $cnt = tr/ / /;
	    $cnt = $cnt + 3;
	} else {
	    $cnt = 2;
	}
	print COUT "$c_type src1${dr_type}_vals[$cnt];\n";
	print COUT "$c_type src2${dr_type}_vals[$cnt];\n";
	print COUT "$c_type br_src${dr_type}_vals[6];\n";
    }
    print COUT "int sh_src2_vals[] = {";
    for (0..31) {
	print COUT " $_,";
    }
    print COUT " 32};\n\n";

    print COUT "int sh_src2s_vals[] = {";
    for (0..15) {
	print COUT " $_,";
    }
    print COUT " 16};\n\n";

    print COUT "int sh_src2c_vals[] = {";
    for (0..7) {
	print COUT " $_,";
    }
    print COUT " 8};\n\n";

    print COUT "unsigned long bit_pattern_vals[] = {";
    for (0..31) {
	printf COUT " 0x%x, 0x%x, ",(1<<$_), ((1<<$_) - 1);
    }
    print COUT " 0xffffffff};\n\n";

print COUT<<EOF;

int failed = 0;
int verbose = 0;

static void initialize_values()
{
EOF
    foreach $dr_type (keys %c_types) {
        $c_type = $c_types{$dr_type};
        $rand_type = $rand_types{$dr_type};
        print COUT "    $c_type rand1_$dr_type = ($c_type)${rand_type}rand48();\n";
	print COUT "    $c_type rand2_$dr_type = ($c_type)${rand_type}rand48();\n";
	$_ = $type_masks{$dr_type};
	$cnt = tr/ / /;
	if ($cnt == 0) {
	    $cnt = 2;
	} else {
	    $cnt = $cnt + 3;
	}
    }
    print COUT "    (void) rand1_p;\n    (void) rand2_p;\n";
    foreach $dr_type (split(' ', "f d c s i l uc us u ul")) {
        $rand_type = $rand_types{$dr_type};
        $c_type = $c_types{$dr_type};
	print COUT "    while (rand1_$dr_type == 0) rand1_$dr_type = ($c_type)${rand_type}rand48();\n";
	print COUT "    while (rand2_$dr_type == 0) rand2_$dr_type = ($c_type)${rand_type}rand48();\n";
    }
    foreach $dr_type (split(' ', "f d c s i l")) {
	print COUT "    src1${dr_type}_vals[0] = rand1_$dr_type;\n";
	print COUT "    src1${dr_type}_vals[1] = -rand1_$dr_type;\n";
	$index = 2;
	foreach $mask (split(' ', $type_masks{$dr_type})) {
	    print COUT "    src1${dr_type}_vals[$index] = rand1_$dr_type & 0x$mask;\n";
	    $index++;
	}
	print COUT "    src2${dr_type}_vals[0] = rand2_$dr_type;\n";
	print COUT "    src2${dr_type}_vals[1] = -rand2_$dr_type;\n";
	$index = 2;
	foreach $mask (split(' ', $type_masks{$dr_type})) {
	    print COUT "    src2${dr_type}_vals[$index] = rand2_$dr_type & 0x$mask;\n";
	    $index++;
	}
    }
    foreach $dr_type (split(' ', "uc us u ul")) {
        $c_type = $c_types{$dr_type};
	print COUT "    src1${dr_type}_vals[0] = ($c_type) rand1_$dr_type;\n";
	print COUT "    src1${dr_type}_vals[1] = ($c_type) -(long)rand1_$dr_type;\n";
	$index = 2;
	foreach $mask (split(' ', $type_masks{$dr_type})) {
	    print COUT "    src1${dr_type}_vals[$index] = ($c_type)rand1_$dr_type & 0x$mask;\n";
	    $index++;
	}
	print COUT "    src2${dr_type}_vals[0] = ($c_type) rand2_$dr_type;\n";
	print COUT "    src2${dr_type}_vals[1] = ($c_type) -(long)rand2_$dr_type;\n";
	$index = 2;
	foreach $mask (split(' ', $type_masks{$dr_type})) {
	    print COUT "    src2${dr_type}_vals[$index] = ($c_type)rand2_$dr_type & 0x$mask;\n";
	    $index++;
	}
    }
    foreach $dr_type (split(' ', "uc us u ul")) {
        $c_type = $c_types{$dr_type};
	print COUT "    br_src${dr_type}_vals[0] = ($c_type) rand1_$dr_type;\n";
	print COUT "    br_src${dr_type}_vals[1] = ($c_type) (-(long)rand1_$dr_type);\n";
	print COUT "    br_src${dr_type}_vals[2] = ($c_type) rand1_$dr_type + 1;\n";
	print COUT "    br_src${dr_type}_vals[3] = ($c_type) (-(long)rand1_$dr_type) + 1;\n";
	print COUT "    br_src${dr_type}_vals[4] = ($c_type) rand1_$dr_type - 1;\n";
	print COUT "    br_src${dr_type}_vals[5] = ($c_type) -((long)rand1_$dr_type) - 1;\n";
    }
    print COUT "    externs[0].extern_name = \"printf\";\n";
    print COUT "    externs[0].extern_value = (void*)(intptr_t)printf;\n";
    print COUT "    externs[1].extern_name = (void*)0;\n";
    print COUT "    externs[1].extern_value = (void*)0;\n";
    print COUT "}\n";
    %c_src2_values = (">>i", "sh_src2_vals", ">>l", "sh_src2_vals", 
		      ">>u", "sh_src2_vals", ">>ul", "sh_src2_vals",
		      ">>ii", "sh_src2_vals", ">>li", "sh_src2_vals", 
		      ">>ui", "sh_src2_vals", ">>uli", "sh_src2_vals",
		      ">>c", "sh_src2c_vals", ">>s", "sh_src2s_vals", 
		      ">>uc", "sh_src2c_vals", ">>us", "sh_src2s_vals",
		      ">>ci", "sh_src2c_vals", ">>si", "sh_src2s_vals", 
		      ">>uci", "sh_src2c_vals", ">>usi", "sh_src2s_vals",
		      "<<i", "sh_src2_vals", "<<l", "sh_src2_vals", 
		      "<<u", "sh_src2_vals", "<<ul", "sh_src2_vals",
		      "<<ii", "sh_src2_vals", "<<li", "sh_src2_vals", 
		      "<<ui", "sh_src2_vals", "<<uli", "sh_src2_vals",
		      "<<c", "sh_src2c_vals", "<<s", "sh_src2s_vals", 
		      "<<uc", "sh_src2c_vals", "<<us", "sh_src2s_vals",
		      "<<ci", "sh_src2c_vals", "<<si", "sh_src2s_vals", 
		      "<<uci", "sh_src2c_vals", "<<usi", "sh_src2s_vals"
);
    foreach $dr_op (split(' ', "add sub mul div mod and or xor")) {
        foreach $dr_type (split(' ', "i u l ul")) {
	    $full_op = "${dr_op}${dr_type}i";
	    $c_src2_values{"$full_op"} = "bit_pattern_vals";
	}
        foreach $dr_type (split(' ', "i u l ul")) {
	    $full_op = "${dr_op}${dr_type}";
	    $c_src2_values{"$full_op"} = "src2${dr_type}_vals";
	    $c_src1_values{"$full_op"} = "src1${dr_type}_vals";
	}
    }
    foreach $dr_type (split(' ', "i u l ul d f")) {
	$full_op = "br${dr_type}";
	$c_src2_values{"$full_op"} = "br_src${dr_type}_vals";
	$c_src1_values{"$full_op"} = "br_src${dr_type}_vals";
    }
}
