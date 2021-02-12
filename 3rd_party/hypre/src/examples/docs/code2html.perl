#!/usr/bin/perl -w
my $vernr = "0.9.1";
my $monthshort = "Jan";
my $monthlong = "Jan";
my $year = "2002";
########################################################################
#                                                                      #
# Code2HTML                                                            #
# ---------                                                            #
#                                                                      #
# Code2Html, peter@palfrader.org                                       #
#                                                                      #
# $Date$
# $Revision$
# $Id$
#                                                                      #
# AUTHOR                                                               #
#        Peter  Palfrader. Written in 1999, 2000, 2001, 2002.          #
#        A lot of other people. See CREDITS file.                      #
#                                                                      #
# DESCRIPTION                                                          #
#        code2html is a  perlscript  which  converts  a  program       #
#        source  code  to syntax highlighted HTML by applying a set    #
#        of   regular   expressions   depending   on  the  language    #
#        the source code is written.                                   #
#                                                                      #
#        see the man-page for details,                                 #
#                                                                      #
########################################################################

use strict;
use Getopt::Long;

my $FILES_DISALLOWED_IN_CGI = 1;
# you may set this to false to allow file reading from your hd in 
# cgi mode. This may be not good if your httpd runs as 'root' (yes, I've
# seen this!) and so any user could with some knowledge easily read
# your /etc/shadow for example!
my $FILES_REDIRECT_DISALLOWED = 1;
my $LANG_TEST_LENGTH = 1024;


# PP: I think Compress::Zlib could be nice for this. but it's not very widespread :(
# PP: A hash would be nicer but then it would not possible to get the keys in this very order (AFAIK)
# PP: If names contain meta characters, then those must be metaquoted (if you don't want the meta chars to be meta chars of course)
my @CGI_ENCODING = (
		    ['bzip2'     , '/usr/bin/bzip2'    , '--stdout' ],
		    ['gzip'      , '/bin/gzip'         , '--stdout' ],
		    ['compress'  , '/usr/bin/compress' , '-c' ]
		   );



# undefine the input record separator so everything gets loaded in one turn
undef $/;



my $pure_version_message = "Code2Html, version $vernr, $monthshort $year, peter\@palfrader.org";
my $version_message = "$pure_version_message\n";

my $short_short_help = "Try `code2html --help' for more information.\n";
my $short_help = 
"$pure_version_message
Usage: code2html [options] [input_file [output_file]]

Convert a program source to syntax highlighted HTML,
or any other  format for wich rules are defined.

-l, --language-mode   set language mode
    --fallback LANG   fallback language mode
-v, --verbose         prints progress information to STDER
-n, --linenumbers     print out the source code with line numbers
-P, --prefix          optional prefix to use for linenumber anchors
-N, --linknumbers     linenumbers will link to themselves
-t, --replace-tabs[=TABSTOP-WIDTH] 
                      replace <tabs> with spaces
-L, --language-file=LANGUAGE-FILE
                      specify an alternate file for definitions
-m, --modes           print all available modes
-h, --help            print this message
-V, --version         print version
-c, --content-type    prints a Content-Type header
-o, --output-format   selects the output-format
-H, --no-header       don't use the template
    --template=FILE   override template
-T, --title           set title

-w, --linewidth       max characters per line
-b, --linebreakprefix prefix of the new lines

see the man-page code2html for further help
";





my $USE_CGI_FOR_ERRORS      = 0; # is switched on in parse params if necessary
$SIG{'__DIE__'} = 
  sub {
    if ($USE_CGI_FOR_ERRORS) { print "Content-Type: text/plain\n\n", $0, ': ', $_[0], "\n"; }
    else                     { print STDERR $0, ': ', $_[0];                                };
    exit 1;
  };

$SIG{'__WARN__'} = 
  sub {
    unless ($USE_CGI_FOR_ERRORS) { print STDERR $0.': '.$_[0]; };
  };









my $DEFAULT_OUTPUTFORMAT='html';
my $DEFAULT_OUTPUTFORMAT_IN_CGI='html';
my $ENTITIES;
my %ENTITIES;


my %params = &parse_params;
if    ($params{'what_to_do'} eq 'patch_html') { &patch_html(\%params) }
elsif ($params{'what_to_do'} eq 'normal'    ) { &main(\%params)       }
else  { die("I don't know what to do :(\n") };










sub main
  {
      my %params = %{shift()};


      print STDERR "getting patterns...\n"                      if ($params{'verbose'});
      # building up the database
      # newer entries overwrite old ones
      my @CONFIG_FILES;
      push @CONFIG_FILES, "/etc/code2html.config";
      push @CONFIG_FILES, $ENV{'HOME'}."/.code2html.config"   if (defined($ENV{'HOME'}));
      push @CONFIG_FILES, split(/:/,$ENV{'CODE2HTML_CONFIG'}) if ($ENV{'CODE2HTML_CONFIG'});
      push @CONFIG_FILES, split(/:/,$params{'langfile'})      if defined($params{'langfile'});
     
      my %STYLESHEET = %{ &get_default_stylesheet } ; 
      my %LANGUAGE   = %{ &get_default_database   } ; 

      for (@CONFIG_FILES) {
	  if ( -r $_){
	      # if I use `do $_` instead of scalar eval... %LANGUAGE is not exported and imported correctly (read: at all) (PP)
	      unless (scalar eval `cat $_`) {     
		  warn "couldn't parse $_: $@" if $@;
	      };
	  };
      };




      if (defined($params{'modes'}) && $params{'modes'})
	{
	    print "Defined modes: ";
	    print join( ', ', sort keys %LANGUAGE ), ".\n" ;
	    print "Defined outputformats: ";
	    print join( ', ', sort keys %STYLESHEET ), ".\n" ;
	    exit;
	};




      
      # set outputformat
      die "Outputformat $params{'outputformat'} not defined" unless defined $STYLESHEET{$params{'outputformat'}};
      my %STYLE = % { $STYLESHEET{$params{'outputformat'}} };
      
      # load alternate template if given
      if (($params{'template'} ne "") && ( ! $params{'noheader'} )) {
	  open (FILE, $params{'template'}) || die ("Could not open template file $params{'template'}: $!");
	  $STYLE{'template'} = <FILE>;
	  close (FILE);
      };

      # set up the global ENTITIES variables ( the scalar and the hash ) from the STYLE definition
      $ENTITIES =     $ { $STYLE{'entities'} }{'listofchars'};
      %ENTITIES = % { $ { $STYLE{'entities'} }{'replace_by' } };

      # modify the header and footer so that the template variables are set correcly
      unless ($STYLE{'template'} =~ /^(.*)%%code%%(.*)$/s) {
	  die "template does not contain a %%code%% variable";
      };
      $STYLE{'header'} = $1;
      $STYLE{'footer'} = $2;
      $STYLE{'header'} =~ s/%%title%%/$params{'title'}/g;
      $STYLE{'footer'} =~ s/%%title%%/$params{'title'}/g;
      $STYLE{'header'} =~ s/%%version%%/$vernr/g;
      $STYLE{'footer'} =~ s/%%version%%/$vernr/g;



      # load the input file and set params{'langmode'} if it is not already. this is done by probing a
      # set of rules defined in %LANGUAGE
      my $code_ref;
      print STDERR "loading input file...\n"                    if ($params{'verbose'});
      $code_ref = &get_input_file(\%params, \%LANGUAGE, $params{'langmode'}, $params{'alt_langmode'});

      # select the rules for out language.
      my $language_rules_ref = $LANGUAGE{ lc($params{'langmode'}) }->{'patterns'};

      print STDERR "applying stylesheet...\n"                   if ($params{'verbose'});
      # Apply the Stylesheets
      # set 'starttag' and 'endtag' for every rule according to its 'style' value
      # the tags are defined in the stylesheet
      &apply_stylesheets_to_rules( $language_rules_ref, \%STYLE );

      print STDERR "outputting headers...\n"                    if ($params{'verbose'});
      &put_headers(\%params, \%STYLE);

      my $snippetlist_ref = [] ;
      print STDERR "creating snippet-list...\n"                 if $params{'verbose'};
      &create_snippetlist( $language_rules_ref, $$code_ref, $snippetlist_ref, \%STYLE);

      print STDERR "outputting file...\n"                       if $params{'verbose'};
      return &put_output(\%params, $snippetlist_ref, \%STYLE);
}





sub patch_html
  {
      my %params = %{shift()};
      my $code;

      open(FILEHANDLE, $params{'infile'}) || die("While opening '$params{'infile'}' for input: ".$!."\n");
      $code = <FILEHANDLE>;
      close(FILEHANDLE);
      
      $code =~ s/<!-- code2html delete start -->.*?<!-- code2html delete stop -->//gs;
      my $counter=0;
      my @chunks = split ( /(<!-- code2html add.*?-->)/s , $code);

      $code = '';
      for (@chunks)
	{
	    $code .= $_;
    	    if ($_ =~ /<!-- code2html add(.*?)(\n.*?)?-->/s)
	      {
		  my $cmdline = $1; 
		  my $input   = $2;
		  $cmdline =~ s/^[ \t]*//g;
		  $cmdline =~ s/[ \t]*$//g;
		  @ARGV = split ( / / , $cmdline);
		  my %new_params = &parse_params;


		  $new_params{'input'} = $input  if ($new_params{'infile'} eq "-");


		  undef $new_params{'outfile'};
		  ++$counter;
		  $new_params{'line_number_prefix'} = $counter unless (defined $new_params{'line_number_prefix'});

		  $new_params{'verbose'} = $params{'verbose'};

		  my $no_header = $new_params{'noheader'};
		  $new_params{'noheader'} = 1;
		  $new_params{'dont_print_output'} = 1;

		  if ($no_header)
		    {
			$code .= '<!-- code2html delete start -->'..
			  &main(\%new_params).
			    '<!-- code2html delete stop -->';
		    }
		  else
		    {
			$code .= '<!-- code2html delete start --><pre>'.
			  &main(\%new_params).
			    '</pre><!-- code2html delete stop -->';
		    };
	      };
	};


      open(FILEHANDLE, '>'.$params{'outfile'}) || die("While opening '$params{'outfile'}' for output: ".$!."\n");
      print FILEHANDLE $code;
      close(FILEHANDLE);
  };






#####################################################################
################### get_input_data ##################################
#####################################################################
# Reads the input data for the cgi script.
# in : nothing
# out: a hash with the input data
sub get_input_data
  {
    my $input_data;
    my %f;
    if($ENV{'REQUEST_METHOD'} eq 'GET')  { $input_data = $ENV{'QUERY_STRING'}; }
    else {  read(STDIN, $input_data, $ENV{'CONTENT_LENGTH'}); };
    
    
    if ($ENV{'CONTENT_TYPE'} =~ m/^multipart\/form-data; boundary=(.*)$/i)
      {
        my $boundary = quotemeta($1);
        my @blocks = split(/$boundary/, $input_data);
        
        for (@blocks)
          {
            if (my $dummy = m/name="(.*?)"/i) 
              {
                my $name = $1;
                $_ =~ s/\r\n/\n/g;
                m/\n\n(.*)\n/s;
                my $value = $1;
                $f{$name}=$value;
              };
          };
      }
    elsif ($ENV{'CONTENT_TYPE'} =~ m/^multipart\/form-data;$/i) # if the boundary is not in the enviroment variable we'll guess
      {
        my $dummy = $input_data =~ m/^(.*?)(\n|\r)/;
        my $boundary = $1;

        my @blocks = split(/$boundary/, $input_data);

        for (@blocks)
          {
            if (my $dummy = m/name="(.*?)"/i)
              {
                my $name = $1;
                $_ =~ s/\r\n/\n/g;
                m/\n\n(.*)\n/s;
                my $value = $1;
                $f{$name}=$value;
              };
          };
      }
    else
      {
        my @form_fields = split(/&/, $input_data);
        
        for (@form_fields)
          {
            my ($name, $value) = split(/=/, $_);
            $value =~ tr/+/ /;
            $value =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
            
            $f{$name} = $value;
          }
      };

    return %f;
  };

################################################################################
####################### parse_params ###########################################
################################################################################
sub parse_params
  {
      my %RESULT;

      if (defined($ENV{'GATEWAY_INTERFACE'}) && (!scalar(@ARGV))) # if there is a CGI enviroment and no parameters/options given
	{
	    $USE_CGI_FOR_ERRORS = 1;
	    $RESULT{'content-type'} = 1;
	    $RESULT{'what_to_do'} = 'normal';
	    
	    my %input = &get_input_data;

	    $input{'input-selector'} = $input{'input_selector'}        unless (defined $input{'input-selector'});
	    $input{'no-encoding'}  = $input{'no_encoding'}             unless (defined $input{'no-encoding'});
	    $input{'line-numbers'} = $input{'line_numbers'}            unless (defined $input{'line-numbers'});
	    $input{'replace-tabs'} = $input{'replace_tabs'}            unless (defined $input{'replace-tabs'});
	    $input{'language-mode'} = $input{'language_mode'}          unless (defined $input{'language-mode'});
	    $input{'cgi-input1'} = $input{'cgi_input1'}                unless (defined $input{'cgi-input1'});
	    $input{'cgi-input2'} = $input{'cgi_input2'}                unless (defined $input{'cgi-input2'});

	    if ($input{'input-selector'} =~ /^cgi[-_]input[12]$/ )
	      {
		  my $input_selector = $input{'input-selector'};
		  die("CGI parse error: $input_selector does not exist!") unless (defined $input{$input_selector});
		  $RESULT{'input'} = $input{$input_selector};
		  $RESULT{'title'} = 'code2html result of cgi input form';
	      }
	    elsif ($input{'input-selector'} eq "file")
	      {
		  die('CGI parse error: option not supported due to security reasons!') if ($FILES_DISALLOWED_IN_CGI);
		  die('CGI parse error: filename not defined!') unless (defined $input{'filename'});
		  $RESULT{'infile'} = $input{'filename'};
		  $RESULT{'title'} = $RESULT{'infile'};
	      }
	    elsif ($input{'input-selector'} eq "REDIRECT_URL")
	      {
		  die('CGI parse error: option not supported due to security reasons!') if ($FILES_REDIRECT_DISALLOWED);
		  die('CGI parse error: ENV: REDIRECT_URL not defined!') unless (defined $ENV{'REDIRECT_URL'});
		  $RESULT{'infile'} = $ENV{'DOCUMENT_ROOT'}.$ENV{'REDIRECT_URL'};
		  $RESULT{'title'} = $RESULT{'infile'};
	      }
	    else
	      {
		  die('CGI parse error: input selector not given!');
	      };
	    
	    if ((!defined ($input{'no-encoding'})) || $input{'no-encoding'})
	      {
		  for (@CGI_ENCODING)
		    {
			if ( ($ENV{'HTTP_ACCEPT_ENCODING'} =~ m/\b $_->[0] \b/x)  &&  # PP: if supported by the browser
			     (-x $_->[1]) )                                           # PP: and executable by the script
			  {
			      $RESULT{'encoding'} = $_->[0];
			      $RESULT{'encoder' } = $_->[1] .' '. $_->[2];
			      last;
			  };
		    }
	      };

	    $RESULT{'linenumbers'} = 'none';
	    if ($input{'line-numbers'} eq "yes")  { $RESULT{'linenumbers'} = 'normal'; };
	    if ($input{'line-numbers'} eq "link") { $RESULT{'linenumbers'} = 'linked'; };
	    if (defined($input{'replace_tabs'}))  { $RESULT{'replacetabs'} = $input{'replace-tabs'} };
	    if (defined($input{'fallback'}))      { $RESULT{'alt_langmode'} = $input{'fallback'} };
	    if (defined($input{'language_mode'})) { $RESULT{'langmode'} = $input{'language-mode'} };
	    if (defined($input{'title'}))         { $RESULT{'title'} = $input{'title'} };

	    $RESULT{'content_type'} = 1;
	    $RESULT{'outputformat'} = $DEFAULT_OUTPUTFORMAT_IN_CGI;
	    $RESULT{'outfile'}      = '-';
	}
      else
	{
	    my $verbose           = 0;
	    my $linenumbers       = 0;
	    my $linknumbers       = 0;
	    my $replace_tabs      = 0;
	    my $language_file     = '';
	    my $language_mode     = '';
	    my $modes             = 0;
	    my $fallback          = '';
	    my $help              = 0;
	    my $version           = 0;
	    my $content_type      = 0;
	    my $no_header         = 0;
	    my $outputformat      = $DEFAULT_OUTPUTFORMAT;
	    my $template          = '';
	    my $title             = "__NOTHING__$$"; # some magix ;(
	    my $prefix            = undef;
            my $linewidth         = undef;
            my $linebreakprefix   = undef;
            my $linebreakprefixdefault = '� ';

	    my $patch_html;


	    # Get Options does not like - as a parameters (used for STDIN and STDOUT)
	    # So we're using a stupid magix again
	    @ARGV      = map { $_ eq '-' ? "__STD__$$" : $_ } @ARGV; 

	    Getopt::Long::config('bundling');
	    unless ( GetOptions( 
				"--verbose"               , \$verbose              , 
				"-v"                      , \$verbose            ,
			      
				"--linenumbers"           , \$linenumbers          , 
				"-n"                      , \$linenumbers        , 
			      
				"--linknumbers"           , \$linknumbers          , 
				"-N"                      , \$linknumbers        , 

				"--prefix=s"              , \$prefix               ,
				"-P=s"                    , \$prefix               ,
			      
				"--replace-tabs=i"        , \$replace_tabs         , 
				"--replace_tabs=i"        , \$replace_tabs         , 
				"-t=i"                    , \$replace_tabs       , 
				
				"--language-file=s"       , \$language_file        , 
				"--language_file=s"       , \$language_file        , 
				"-L=s"                    , \$language_file        , 
				
				"--language-mode=s"       , \$language_mode        , 
				"--language_mode=s"       , \$language_mode        , 
				"-l=s"                    , \$language_mode        ,
				
				"--title=s"               , \$title                , 
				"-T=s"                    , \$title                ,
				
				"--modes"                 , \$modes                , 
				"-m"                      , \$modes                , 
				
				"--fallback=s"            , \$fallback             , 
				
				"--output=s"              , \$outputformat         , 
				"-o=s"                    , \$outputformat         , 
			      
				"--template=s"            , \$template             ,
			      
				"--help"                  , \$help                 , 
				"-h"                      , \$help                 , 

				"--version"               , \$version              , 
				"-V"                      , \$version              , 
			      
				"--content-type"          , \$content_type         , 
				"--content_type"          , \$content_type         , 
				"-c"                      , \$content_type         , 
				
				"--no-header"             , \$no_header            ,
				"--no_header"             , \$no_header            ,
				"-H"                      , \$no_header            ,

			    
				"--patch-html"            , \$patch_html           ,
				"--patch_html"            , \$patch_html           ,
				"-p"                      , \$patch_html           ,

				"--linewidth=i"           , \$linewidth            ,
				"-w=i"                    , \$linewidth            ,
				"--linebreakprefix=s"     , \$linebreakprefix      ,
				"-b=s"                    , \$linebreakprefix      ,
			       )
		   )
	      {
		  print STDERR $short_short_help;
		  exit 1;
	      }

	    #reversing magix
	    @ARGV      = map { $_ eq "__STD__$$" ? '-' : $_ } @ARGV; 

	    if ($help)              { print STDERR $short_help;    exit 0; };
	    if ($version)           { print $version_message;      exit 0; };

	    if ($patch_html)
	      {
		  $RESULT{'what_to_do'} = 'patch_html';
		  $RESULT{'verbose'}    = $verbose;

		  if (!defined ($RESULT{'infile'}        = shift(@ARGV))) { $RESULT{'infile'} = '-' };
		  if (!defined ($RESULT{'outfile'}       = shift(@ARGV))) { $RESULT{'outfile'} = $RESULT{'infile'}};
		  if (defined (shift(@ARGV))) { print STDERR  "too many parameters!\n";
						print STDERR $short_help;
						exit 1;
					    };
	      }
	    else
	      {
		  $RESULT{'what_to_do'} = 'normal';

		  $RESULT{'verbose'}        = $verbose;
		  if    ($linknumbers) { $RESULT{'linenumbers'}    = 'linked' }
		  elsif ($linenumbers) { $RESULT{'linenumbers'}    = 'normal' }
		  else                 { $RESULT{'linenumbers'}    = 'none' };
		  $RESULT{'line_number_prefix'} = $prefix;
		  $RESULT{'replacetabs'}    = $replace_tabs;
		  $RESULT{'langfile'}       = $language_file;
		  $RESULT{'modes'}          = $modes;
		  $RESULT{'alt_langmode'}   = $fallback;
		  $RESULT{'content_type'}   = $content_type;
		  $RESULT{'noheader'}       = $no_header;
		  $RESULT{'langmode'}       = $language_mode;
		  $RESULT{'template'}       = $template;
		  $RESULT{'outputformat'}   = $outputformat;
		  $RESULT{'linewidth'}      = $linewidth;
		  $RESULT{'linebreakprefix'}= $linebreakprefix;

		  if (defined ($RESULT{'linebreakprefix'}) &&
		      !defined ($RESULT{'linewidth'})) {
		      printf (STDERR "--linebreakprefix|-b does not make sense without --linewidth|-w!\n");
		      print STDERR $short_help;
		      exit 1;
		  }
		  if (defined ($RESULT{'linewidth'})) {
		      if ($RESULT{'linewidth'} <= 0) {
			  printf (STDERR "linewidth must be greater then 0!\n");
			  print STDERR $short_help;
			  exit 1;
		      }
		      if (!defined ($RESULT{'linebreakprefix'})) {
			  $RESULT{'linebreakprefix'} = $linebreakprefixdefault;
		      }
		  }

		  if (!defined ($RESULT{'infile'}          = shift(@ARGV))) { $RESULT{'infile'} = '-'};
		  if (!defined ($RESULT{'outfile'}         = shift(@ARGV))) { $RESULT{'outfile'} = '-'};
		  if (defined (shift(@ARGV))) { print STDERR  "too many parameters!\n";
						print STDERR $short_help;
						exit 1;
					    };
	      };
	      #the magix again
	    $RESULT{'title'} = $title eq  "__NOTHING__$$" ? ($RESULT{'infile'} eq '-' ? 'STDIN' : $RESULT{'infile'}) : $title;
	};


      return %RESULT;
  };


################################################################################
####################### checkTabulator #########################################
################################################################################
sub checkTabulator
{
    my ($line, $TABSTOP) = @_;
    
    while ((my $at = index($line, "\t")) != -1)
      {
	  my $cnt = ($TABSTOP - ($at % $TABSTOP));
	  my $replace_with = ' ' x $cnt if ($cnt);
	  $line =~ s/\t/$replace_with/;
      };

    return $line;
}

################################################################################
####################### splitLine ##############################################
################################################################################
sub splitLine
{
    my ($line, $linewidth, $prefix) = @_;
    
    my $length = length ($line);
    my $pos = 0;

    while ($length - $pos > $linewidth)
      {
          my $maxoff = ($pos + $linewidth > $length) ? ($length - 1)
                                                     : ($pos + $linewidth);
          my $newpos = rindex ($line, " ", $maxoff);
	  if ($newpos > $pos) {
	      $pos = $newpos;
	      $line = substr ($line, 0, $pos)."\0$prefix".substr ($line, $pos + 1, $length);
	  } else {
	      $pos = $pos + $linewidth + 1;
	      $line = substr ($line, 0, $pos)."\0$prefix".substr ($line, $pos, $length);
	  }
      };

    return $line;
}

################################################################################
####################### get_input_file #########################################
################################################################################
sub get_input_file
  {

    # in  : \%params
    # in : \%LANGUAGE;
    # in/out : $langmode;
    # in/out : $alt_langmode;
    # returns: input file
    
      my %PARAMS       = %{$_[0]};
      my %LANGUAGE     = %{$_[1]};
      my $langmode     = $_[2];
      my $alt_langmode  = $_[3];
      my $code;

      
      if (defined $PARAMS{'input'})
	{
	    $code = $PARAMS{'input'};
	    $code =~ s/\r//g;
	}
      else
	{
	    open(FILEHANDLE, $PARAMS{'infile'}) || die("While opening '$PARAMS{'infile'}' for input: ".$!."\n");
	    $code = <FILEHANDLE>;
	    close(FILEHANDLE);
	};
      
      if ($PARAMS{'replacetabs'} != 0)
	{
	    $code = join (
			  "\n",
			  map{
			      &checkTabulator($_, $PARAMS{'replacetabs'})
			  }
			  my @dummy = split(/\n/, $code)
			 );
	};



      if (defined ($PARAMS{'linewidth'}))
	{
	    $code = join (
			  "\n",
			  map{
			      &splitLine($_, $PARAMS{'linewidth'},
					 $PARAMS{'linebreakprefix'})
			  }
			  my @dummy = split(/\n/, $code)
			 );
	};
      

      
      if ((!defined($langmode)) || ($langmode eq ''))
	{
	    my $test_code = substr($code, 0, $LANG_TEST_LENGTH);
	    warn("language mode not given. guessing...\n");

	    $langmode = '';

	    for (keys %LANGUAGE)
	      {
		  if (  (($LANGUAGE{$_}->{'filename'} ne '') && ($PARAMS{'infile'} =~  m/$LANGUAGE{$_}->{filename}/))  ||
		        (($LANGUAGE{$_}->{'regex'}    ne '') && ($test_code        =~  m/$LANGUAGE{$_}->{regex}/   ))   )
		    {
			$langmode = $_;
			last;
		    };
	      };

	    if ($langmode eq '')
	      {
		  if ((defined($alt_langmode)) && ($alt_langmode ne ''))
		    {
			warn("Guessing language mode failed. Using fallback mode: '$alt_langmode'\n");
			$langmode = $alt_langmode;
			$alt_langmode = '';
		    }
		  else
		    {
			die("Guessing language mode failed.\n")
		    };
	      }
	    else
	      {
		  warn("using '$langmode'\n");
	      };
	};
      
      $_[2] = $langmode;
      $_[3] = $alt_langmode;
      return \$code;
  };


################################################################################
####################### put_headers ############################################
################################################################################
sub put_headers
{	
      my %PARAMS = %{shift()};
      my $STYLE_REF = shift();

      if (defined($PARAMS{'outfile'}))
	{
	    unless ($PARAMS{'outfile'} eq '-'){
		open(SAVEOUT, ">&STDOUT");             print SAVEOUT '';  # so perl does not typo warn
		open (STDOUT, '>'.$PARAMS{'outfile'}) || die("While redirecting STDOUT to '$PARAMS{'outfile'}' for output: ".$!."\n");
	    };

	    if (defined $PARAMS{'encoding'}) 
	      {
		  $|=1; # PP: so the header is written before the data!
		        # PP: this took me hours of debugging :(
		  print "Content-Type: $$STYLE_REF{'content-type'}\n"                                                    if ($PARAMS{'content_type'});
		  print "Content-Encoding: $PARAMS{'encoding'}\n\n";
		  open (FILEHANDLE, "|$PARAMS{'encoder'}") || die("While opening '$PARAMS{'encoder'}': ".$!."\n");
	      }
	    else
	      {
		  open( FILEHANDLE, ">&STDOUT" ) ;
		  print FILEHANDLE "Content-Type: $$STYLE_REF{'content-type'}\n\n"                                       if ($PARAMS{'content_type'});
	      };
	     
	     print FILEHANDLE $$STYLE_REF{'header'} unless $PARAMS{'noheader'};
        }
};

################################################################################
####################### apply_stylesheets_to_rules #############################
################################################################################
sub apply_stylesheets_to_rules
  {
      my ( $regexps_ref, $style_ref ) = @_;

      for ( @$regexps_ref ) {
#	  warn ("Style '".$_->{style}."' not defined in stylesheet.\n") unless defined $ { $$style_ref{'tags'} } { $_->{style} };
	  if (defined ($ { $$style_ref{'tags'} } { $_->{style} })) {
	      $_->{'starttag'} = $ { $ { $$style_ref{'tags'} } { $_->{style} } } { 'start' };
	      $_->{'endtag'}   = $ { $ { $$style_ref{'tags'} } { $_->{style} } } { 'stop' };
	  } else {
	      # no style no formating; if style == '' formating is done by childregex
	      warn ("Style '".$_->{style}."' not defined in stylesheet.\n") if ($_->{style} ne '');
	      $_->{'starttag'} = ''; #$ { $ { $$style_ref{'tags'} } { $_->{style} } } { 'start' };
	      $_->{'endtag'}   = ''; #$ { $ { $$style_ref{'tags'} } { $_->{style} } } { 'stop' };	      
	  }
	  apply_stylesheets_to_rules( $_->{childregex}, $style_ref ) if $_->{childregex};
      };
  };

################################################################################
####################### create_snippetlist #####################################
################################################################################
sub create_snippetlist
  {
    my ( $regexps_ref, $code, $snippetlist_ref, $style_ref ) = @_ ;
    my $length = length( $code );

    ## An array of regular expression sturctures, each of which is an
    ## array.  @res is kept sorted by starting position of the RExen and
    ## then by the position of the regex in the language file.  This allows
    ## us to just evaluate $res[0], and to hand write fast code that typically
    ## handles 90% of the cases without resorting to the _big_ guns.
    ##
    ## FWIW, I pronounce '@res' REEZE, as in the plural of '$re'.
    ##
    my @res ;
    
    my $pos ;
    
    for ( @$regexps_ref ) {
	pos( $code ) = 0 ;
#++$m ;
	next unless $code =~ m/($_->{regex})/gms ;

	$pos = pos( $code ) ;
#	$res[@res] = [ 
#		      $_->{regex},
#		      $ { $ { $$style_ref{'tags'} } { $_->{style} } } { 'start' },
#		      $ { $ { $$style_ref{'tags'} } { $_->{style} } } { 'stop' },
#		      $_->{childregex},
#		      $pos - length( $1 ),
#		      $pos,
#		      scalar( @res ),
#		     ] ;
	$res[@res] = [ 
		      $_->{regex},
		      $_->{starttag},
		      $_->{endtag},
		      $_->{childregex},
		      $pos - length( $1 ),
		      $pos,
		      scalar( @res ),
		     ] ;
    }
    
    ## 90% of all child regexes end up with 0 or 1 regex that needs to be
    ## worried about. Trimming out the 0's speeds things up a bit and
    ## makes the below loop simpler, since there's always at least
    ## 1 regexp.  It donsn't speed things up much by itself: the percentage 
    ## of times this fires is really small.  But it does simplify the loop
    ## below and speed it up.
    unless ( @res ) {
	$code =~ s/($ENTITIES)/$ENTITIES{$1}/ge ;
	push @$snippetlist_ref, $code ;
	return ;
    }
    
    @res = sort { $a->[4] <=> $b->[4] || $a->[6] <=> $b->[6] } @res ;
    
    ## Add a dummy at the end, which makes the logic below simpler / faster.
    $res[@res] = [
		  undef,
		  undef,
		  undef,
		  undef,
		  $length,
		  $length,
		  scalar( @res ),
		 ] ;
    
    ## These are declared here for (minor) speed improvement.
    my $re ;
    my $match_spos ;
    my $match_pos ;
    my $re_spos ;
    my $re_pos ;
    my $re_num ;
    my $prefix ;
    my $snippet ;
    my $rest ;
    my $i ;
    my $l ;
    
my @changed_res ;
my $j ;

    $pos = 0 ;
MAIN:
    while ( $pos < $length ) {
	$re = $res[0] ;
	
	$match_spos = $re->[4] ;
	$match_pos  = $re->[5] ;
	
	if ( $match_spos > $pos ) {
	    $prefix  = substr( $code, $pos, $match_spos - $pos ) ;
	    $prefix  =~ s/($ENTITIES)/$ENTITIES{$1}/ge ;
	    push @$snippetlist_ref, $prefix ;
	}
	
	if ( $match_pos > $match_spos ) {
	    $snippet = substr( $code, $match_spos, $match_pos - $match_spos ) ;
	    if ( @{$re->[3]} ) {
		push @$snippetlist_ref, $re->[1] ;
		create_snippetlist( $re->[3], $snippet, $snippetlist_ref, $style_ref ) ;
		push @$snippetlist_ref, $re->[2] ;
	    }
	    else {
		$snippet =~ s/($ENTITIES)/$ENTITIES{$1}/ge ;
		push @$snippetlist_ref, $re->[1], $snippet, $re->[2];
	    }
	}
	
	$pos = $match_pos ;
	
	##
	## Hand coded optimizations.  Luckily, the cases that arise most often
	## are the easiest to tune.
	##

# =pod

	if ( $res[1]->[4] >= $pos ) {
	    ## Only first regex needs to be moved, 2nd and later are still valid.
	    ## This is often 90% of the cases for Perl or C (others not tested,
	    ## just uncomment the $n, $o, and $p lines and try it yourself).
#++$n{1} ;
#++$m ;
	    pos( $code ) = $pos ;
	    unless ( $code =~ m/($re->[0])/gms ) {
#++$o{'0'} ;
		if ( @res == 2 ) {
		    ## If the only regexp left is the dummy, we're done.
		    $rest = substr( $code, $pos ) ;
		    $rest =~ s/($ENTITIES)/$ENTITIES{$1}/ge ;
		    push @$snippetlist_ref, $rest ;
		    last ;
		}
		shift @res ;
	    }
	    else {
		$re->[5] = $re_pos  = pos( $code ) ;
		$re->[4] = $re_spos = $re_pos - length( $1 ) ;
		
		## Walk down the array looking for $re's new home.
		## The first few loop iterations are unrolled and done manually 
		## for speed, which handles 85 to 90% of the cases where only
		## $re needs to be moved.
		##
		## Here's where that dummy regexp at the end of the array comes
		## in handy: we don't need to worry about array size here, since
		## it will always be after $re no matter what.  The unrolled
		## loop stuff is outdented to make the conditionals fit on one
		## 80 char line.
		## Element 4 in @{$res[x]} is the start position of the match.
		## Element 6 is the order in which it was declared in the lang file.
		$re_num = $re->[6] ;
		if ( ( $re_spos <=> $res[1]->[4] || $re_num <=> $res[1]->[6] ) <= 0 ) {
#++$o{'1'} ;
		    next 
		}
		$res[0] = $res[1] ;

#++$o{'2'} ;
		if ( ( $re_spos <=> $res[2]->[4] || $re_num <=> $res[2]->[6] ) <= 0 ) {
		    $res[1] = $re ;
		    next ;
		}
		$res[1] = $res[2] ;
		
		if ( ( $re_spos <=> $res[3]->[4] || $re_num <=> $res[3]->[6] ) <= 0 ) {
#++$o{'3'} ;
		    $res[2] = $re ;
		    next ;
		}
		$res[2] = $res[3] ;
		
		if ( ( $re_spos <=> $res[4]->[4] || $re_num <=> $res[4]->[6] ) <= 0 ) {
#++$o{'3'} ;
		    $res[3] = $re ;
		    next ;
		}
		$res[3] = $res[4] ;
		
		if ( ( $re_spos <=> $res[5]->[4] || $re_num <=> $res[5]->[6] ) <= 0 ) {
#++$o{'4'} ;
		    $res[4] = $re ;
		    next ;
		}
		$res[4] = $res[5] ;

#++$o{'ugh'} ;
		$i = 6 ;
		$l = $#res ;
		for ( ; $i < $l ; ++$i ) {
		    last
		      if ( 
			  ( $re_spos <=> $res[$i]->[4] || $re_num <=> $res[$i]->[6] )
			  <= 0
			 ) ;
		    $res[$i-1] = $res[$i] ;
		}
#++$p{sprintf( "%2d", $i )} ;
		$res[$i-1] = $re ;
	    }
	    
	    next ;
	}
	
# =cut

	##
	## End optimizations.  You can comment them all out and this net
	## does all the work, just more slowly.  If you do that, then
	## you also need to comment out the code below that deals with
	## the second entry in @res.
	##

#my $ni = 0 ;
	## First re always needs to be tweaked
#++$m ;
#++$ni ;
	pos( $code ) = $pos ;
	unless ( $code =~ m/($re->[0])/gms ) {
	    if ( @res == 2 ) {
		## If the only regexp left is the dummy, we're done.
		$rest = substr( $code, $pos ) ;
		$rest =~ s/($ENTITIES)/$ENTITIES{$1}/ge ;
		push @$snippetlist_ref, $rest ;
		last ;
	    }
	    shift @res ;
	    @changed_res = () ;
	    $i = 0 ;
	}
	else {
	    $re->[5] = $re_pos  = pos( $code ) ;
	    $re->[4] = $re_pos - length( $1 ) ;
	    @changed_res = ( $re ) ;
	    $i = 1 ;
	}
	
	## If the optimizations above are in, the second one always
	## needs to be tweaked, too.
	$re = $res[$i] ;
#++$m ;
#++$ni ;
	pos( $code ) = $pos ;
	unless ( $code =~ m/($re->[0])/gms ) {
	    if ( @res == 2 ) {
		## If the only regexp left is the dummy, we're done.
		$rest = substr( $code, $pos ) ;
		$rest =~ s/($ENTITIES)/$ENTITIES{$1}/ge ;
		push @$snippetlist_ref, $rest ;
		last ;
	    }
	    shift @res ;
	}
	else {
	    $re->[5] = $re_pos  = pos( $code ) ;
	    $re->[4] = $re_spos = $re_pos - length( $1 ) ;
	    if ( @changed_res &&
		 ( $changed_res[0]->[4] <=> $re_spos || 
		   $changed_res[0]->[6] <=> $re->[6]
		 ) > 0
	       ) {
		unshift @changed_res, $re ;
	    }
	    else {
		$changed_res[$i] = $re ;
	    }
	    ++$i ;
	}
	
	for ( ; ; ++$i ) {
	    local $_ = $res[$i] ;
#++$m ;
	    last if $_->[4] >= $pos ;
#++$ni ;
#++$m ;
	    pos( $code ) = $pos ;
	    unless ( $code =~ m/($_->[0])/gms ) {
		if ( @res <= 2 ) {
		    $rest = substr( $code, $pos ) ;
		    $rest =~ s/($ENTITIES)/$ENTITIES{$1}/ge ;
		    push @$snippetlist_ref, $rest ;
		    last MAIN ;
		}
		## If this regex is no longer needed, remove it by not pushing it
		## on to @changed_res.  This means we need one less slot in @res.
		shift @res ;
		redo ;
	    }

	    $_->[5] = $re_pos  = pos( $code ) ;
	    $_->[4] = $re_spos = $re_pos - length( $1 ) ;
	    
	    ## Insertion sort in to @changed_res
	    $re_num = $_->[6] ;
	    for ( $j = $#changed_res ; $j > -1 ; --$j ) {
		last
		  if ( 
		      ( $changed_res[$j]->[4] <=> $re_spos || 
			$changed_res[$j]->[6] <=> $re_num 
		      ) < 0
		     ) ;
		$changed_res[$j+1] = $changed_res[$j] ; 
	    }
	    $changed_res[$j+1] = $_ ;
	}

	## Merge sort @changed_res and @res in to @res
	$j = 0 ;
	$l = $#res ;
	for ( @changed_res ) {
	    while (
		   $i < $l &&
		   ( $_->[4] <=> $res[$i]->[4] || $_->[6] <=> $res[$i]->[6] ) > 0
		  ) {
		$res[$j++] = $res[$i++] ;
	    }
	    $res[$j++] = $_ ;
	}
# =cut
    }
};

##################################################################################
######################### create_snippetlist #####################################
##################################################################################
##sub create_snippetlist
##  {
##    my ( $regexps_ref, $code, $snippetlist_ref ) = @_ ;

##    my $length = length( $code );
##    my @regexps;
##    $regexps[scalar(@$regexps_ref)] = undef;

##    my $head_ptr    = undef;
##    my $current_ptr;
##    my $help_ptr;

##    my $index       = 0;

##    for (@$regexps_ref)
##      {
##	  $current_ptr = $regexps[$index];  #0: start_ptr 1: length 2: next_ptr, 3: regex, 4:start, 5:end, 6: child    7: index
##	  $current_ptr->[7] = $index++;
##	  $current_ptr->[6] = $$_{'childregex'}; 
##	  $current_ptr->[5] = $$_{'endtag'}; 
##	  $current_ptr->[4] = $$_{'starttag'};
##	  $current_ptr->[3] = $$_{'regex'};


##	  pos( $code ) = 0;
##	  if ( $code =~ /($current_ptr->[3])/gms ) { $current_ptr->[0] = pos ($code) - length($1); $current_ptr->[1] = length($1); } else {next};

##	  if (!defined ($head_ptr)  ||  $current_ptr->[0] < $head_ptr->[0] )
##	    {
##		$current_ptr->[2] = $head_ptr;
##		$head_ptr         = $current_ptr;
##	    }
##	  else
##	    {
##		$help_ptr = $head_ptr;
##		$help_ptr = $help_ptr->[2]
##		  while (defined ( $help_ptr->[2] ) && ($current_ptr->[0] >= $help_ptr->[2]->[0]) ); #iow: while (defined help->next && current->pos <= help->next->pos)
		
##		$current_ptr->[2] = $help_ptr->[2];
##		$help_ptr->[2]    = $current_ptr;
##	    };
##      };
    

##    my $endpos = 0;
##    my $oldhead;

##    my %entities ;
##    $entities{'&'} = '&amp;' ;
##    $entities{'<'} = '&lt;' ;
##    $entities{'>'} = '&gt;' ;
##    $entities{'"'} = '&quot;' ;
    
##    my $snippet;
##    while (defined $head_ptr)
##      {
##	  if ($head_ptr->[0] - $endpos > 0) { 
##	      $snippet = substr($code, $endpos, $head_ptr->[0] - $endpos);
##	      $snippet =~ s/($ENTITIES)/$ENTITIES{$1}/ge;          #"]);
##	      push @$snippetlist_ref, $snippet;
##	  };
##	  push  @$snippetlist_ref, $head_ptr->[4];

##	  &create_snippetlist( $head_ptr->[6], substr($code, $head_ptr->[0], $head_ptr->[1]) , $snippetlist_ref);
##	  push  @$snippetlist_ref, $head_ptr->[5];
	  
##	  $endpos = $head_ptr->[0] + $head_ptr->[1];

##	  # update & repair list :

##	  $oldhead = $head_ptr;
##	  # 1) shift now invalid matches from list

##	  $help_ptr = $head_ptr;
##	  $help_ptr = $help_ptr->[2]
##	    while (defined ( $help_ptr->[2] ) && ($endpos > $help_ptr->[2]->[0]) );
##	  $head_ptr = $help_ptr->[2];
##	  $help_ptr->[2] = undef;

##	  # 2) rematch invalid matches and insert them into the list
	  
##	  while (defined $oldhead)
##	    {
##		$current_ptr = $oldhead;
##		$oldhead = $oldhead->[2];

##		pos( $code ) = $endpos;
##		if ( $code =~ /($current_ptr->[3])/gms ) { $current_ptr->[0] = pos ($code) - length($1); $current_ptr->[1] = length($1); } else {next};
##		if (!defined ($head_ptr) ||
##		    ($current_ptr->[0] < $head_ptr->[0]) ||
##		    (
##		     ( $current_ptr->[0] == $head_ptr->[0]) &&
##		     ( $current_ptr->[7] <  $head_ptr->[7])
##		    )		     
##		   )
##		  {
##		      $current_ptr->[2] = $head_ptr;
##		      $head_ptr         = $current_ptr;
##		  }
##		else
##		  {
##		      $help_ptr = $head_ptr;
##		      $help_ptr = $help_ptr->[2]
##			while (defined ( $help_ptr->[2] ) && 
##			       (
##				($current_ptr->[0] >  $help_ptr->[2]->[0]) || 
##				(
##				 ( $current_ptr->[0] == $help_ptr->[2]->[0]) &&
##				 ( $current_ptr->[7] >  $help_ptr->[2]->[7])
##				 )
##			       )
##			      ); #iow: while (defined help->next && current->pos <= help->next->pos)   # if two patterns match at the same pos
##		                                                                                       # the one that was declared earlier is taken
		      
##		      $current_ptr->[2] = $help_ptr->[2];
##		      $help_ptr->[2]    = $current_ptr;
##		  };
##	    };
	  
##	  # 3) done
##      };
      
##    $snippet = substr($code, $endpos); $snippet =~ s/($ENTITIES)/$ENTITIES{$1}/ge;        #"   ]);
##    push  @$snippetlist_ref, $snippet;
##  };



################################################################################
####################### put_output #############################################
################################################################################
sub put_output {
    my ( $params, $snippetlist_ref, $STYLE_REF ) = @_ ;

    my $result;

    my $prefix = ''; 
    $prefix = $params->{'line_number_prefix'}.'_'  if defined $params->{'line_number_prefix'};
    $result = & { $ { $$STYLE_REF{'linenumbers'} }{$params->{'linenumbers'}} } (join ('', @$snippetlist_ref), $prefix);

    if (defined ($params{'linewidth'})) {
	$result =~ tr=\0=\n=;
    }

    print FILEHANDLE $result unless (defined $params->{'dont_print_output'} && $params->{'dont_print_output'});
    print FILEHANDLE $$STYLE_REF{'footer'}  unless $params->{'noheader'};
    
    if (defined($params->{'outfile'})) {
        unless ($params->{'outfile'} eq '-'){
	    close (FILEHANDLE);
	    close (STDOUT);
	    open (STDOUT, ">&SAVEOUT");
	};
    };
    return $result;
};




################################################################################
####################### get_default_stylesheet #################################
################################################################################
sub get_default_stylesheet
{

my %STYLESHEET;


##########
########## different color modes for html. 
# those are named html-dark, html-nobc and html-light. 
# html-light is also named html
# the only difference between html-light and html-nobc is
# that html-light defines a body background and text color.
# nobc stands for no body colors.

$STYLESHEET{'html-light'} =  { 'template'       =>
'<html>
<head>
  <title>%%title%%</title>
</head>
<body bgcolor="#ffffff" text="#000000">
<a href="%%title%%">download the original source code</a>.
<pre>
%%code%%
</pre>
<hr>
syntax highlighted by <a href="http://www.palfrader.org/code2html">Code2HTML</a>, v. %%version%%
</body>
</html>
',
			 'content-type' => 'text/html',
			 'entities'     => { 'listofchars' => '[<>&"]',   # a regex actually
					     'replace_by'  => {
							       '&' => '&amp;',
							       '<' => '&lt;',
							       '>' => '&gt;',
							       '"' => '&quot;'
							      }
					   },
			 'linenumbers'  => {
					    'none'          => sub { 
						                    return $_[0];
					                           },
					    'normal'        => sub { 
                                                                   # o as the first parameter is the joined snippetlist
                                                                   # o the second is an optional prefix, needed if more than one block
                                                                   #   in a file is highlighted. needed in patch-mode. may be empty
                                                                   # the sub should the return a scalar made up of the joined lines including linenumbers
								   my @lines = split ( /\n/, $_[0] );

                                                                   my $nr = 0;
                                                                   my $lengthofnr = length(@lines);
                                                                   my $format = qq{<a name="$_[1]line%u">%${lengthofnr}u</a> %s\n} ;
					                           join ('', map (  {$nr++; sprintf ( $format , $nr, $nr, $_ )} @lines));
					                           },
			                     'linked'       => sub { 
                                                                   # this should do the same as above only with linenumbers that link to themselves
                                                                   # If this style does not support this, use the same as above.
								   my @lines = split ( /\n/, $_[0] );

                                                                   my $nr = 0; 
                                                                   my $lengthofnr = length(@lines);
                                                                   my $format = qq{<a name="$_[1]line%u" href="#$_[1]line%u">%$ {lengthofnr}u</a> %s\n};
                                                                   join ('', map (  {$nr++; sprintf ( $format , $nr, $nr, $nr, $_ )} @lines));
					                           }
					   },
			 'tags'         => { 
					    'comment'                => { 'start' => '<font color="#444444">',
									  'stop'  => '</font>' },
					    'doc comment'            => { 'start' => '<font color="#444444"><i>',
									  'stop'  => '</i></font>' },
					    'string'                 => { 'start' => '<font color="#008000">',
									  'stop'  => '</font>' },
					    'esc string'             => { 'start' => '<font color="#77dd77">',
									  'stop'  => '</font>' },
					    'character'              => { 'start' => '<font color="#008000">',
									  'stop'  => '</font>' },
					    'esc character'          => { 'start' => '<font color="#77dd77">',
									  'stop'  => '</font>' },
					    'numeric'                => { 'start' => '<font color="#FF0000">',
									  'stop'  => '</font>' },
					    
					    'identifier'             => { 'start' => '<font color="#2040a0">',
									  'stop'  => '</font>' },
					    'predefined identifier'  => { 'start' => '<font color="#2040a0"><strong>',
									  'stop'  => '</strong></font>' },
				     
					    'type'                   => { 'start' => '<font color="#2040a0"><strong>',
									  'stop'  => '</strong></font>' },
					    'predefined type'        => { 'start' => '<font color="#2040a0"><strong>',
									  'stop'  => '</strong></font>' },
					    
					    'reserved word'          => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
					    'library function'       => { 'start' => '<font color="a52a2a"><strong>',
									  'stop'  => '</strong></font>' },
					    
					    'include'                => { 'start' => '<font color="0000ff"><strong>',
									  'stop'  => '</strong></font>' },
					    'preprocessor'           => { 'start' => '<font color="0000ff"><strong>',
									  'stop'  => '</strong></font>' },
					    
					    'braces'                 => { 'start' => '<font color="4444FF"><strong>',
									  'stop'  => '</strong></font>' },
					    'symbol'                 => { 'start' => '<font color="4444FF">',
							 	  	  'stop'  => '</font>' },

					    'function header'        => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
					    'function header name'   => { 'start' => '<font color="ff0000">',
									  'stop'  => '</font>' },
					    'function header args'   => { 'start' => '<font color="2040a0">',
									  'stop'  => '</font>' },
					    
					    'regex'                  => { 'start' => '<font color="b000d0">',
									  'stop'  => '</font>' },
					    
					    'text'                   => { 'start' => '<i>',
									  'stop'  => '</i>'},

					    # HTML
					    'entity'                 => { 'start' => '<font color="ff0000">',
									  'stop'  => '</font>' },

					    # MAKEFILE
					    'assignment'             => { 'start' => '<font color="2040a0">',
									  'stop'  => '</font>' },
					    'dependency line'        => { 'start' => '<font color="8b2252">',
									  'stop'  => '</font>' },
					    'dependency target'      => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
					    'dependency continuation'=> { 'start' => '<font color="000000"><strong>',
									  'stop'  => '</strong></font>' },
					    'continuation'           => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
					    'macro'                  => { 'start' => '<font color="2040a0">',
									  'stop'  => '</font>' },
					    'int macro'              => { 'start' => '<font color="4080ff">',
									  'stop'  => '</font>' },
					    'esc $$$'                => { 'start' => '<font color="444444">',
                                                                          'stop'  => '</font>' },

                                            # PATCH
                                            'separator'              => { 'start' => '<font color="00A040"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'line spec'              => { 'start' => '<font color="A0A000"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'deletion'               => { 'start' => '<font color="FF0000"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'insertion'              => { 'start' => '<font color="0000FF"><strong>',
                                                                          'stop'  => '</strong></font>' }

				     }
                       };
# html-light is also called html

$STYLESHEET{'html'} = $STYLESHEET{'html-light'};


# html-nobc is a modification of html-light
# in such a way, that the body tag does not define
# a background and a text color
# nobc stands for no body colors.

%{$STYLESHEET{'html-nobg'}} = %{$STYLESHEET{'html-light'}};
${ $STYLESHEET{'html-nobg'}} {'template'} = '<html>
<head>
  <title>%%title%%</title>
</head>
<body>
<pre>
%%code%%
</pre>
<hr>
syntax highlighted by <a href="http://www.palfrader.org/code2html">Code2HTML</a>, v. %%version%%
</body>
</html>
';


# html-dark is a modification of html-light
# in such a way, that the body tag does define 
# different colors and that the <font> colors are different.

%{$STYLESHEET{'html-dark'}} = %{$STYLESHEET{'html-light'}};
${ $STYLESHEET{'html-dark'}} {'template'} = '<html>
<head>
  <title>%%title%%</title>
</head>
<body bgcolor="#000000"  text="#C0C0C0" vlink="#FFFFFF" alink="#00FF00" link="#FFFFFF">
<pre>
%%code%%
</pre>
<hr>
syntax highlighted by <a href="http://www.palfrader.org/code2html">Code2HTML</a>, v. %%version%%
</body>
</html>
';
${ $STYLESHEET{'html-dark'}} {'tags'} = {
                                            'comment'                => { 'start' => '<font color="#909000">',
                                                                          'stop'  => '</font>' },
                                            'doc comment'            => { 'start' => '<font color="#909000"><i>',
                                                                          'stop'  => '</i></font>' },
                                            'string'                 => { 'start' => '<font color="yellow">',
                                                                          'stop'  => '</font>' },
                                            'esc string'             => { 'start' => '<font color="#77dd77">',
                                                                          'stop'  => '</font>' },
                                            'character'              => { 'start' => '<font color="yellow">',
                                                                          'stop'  => '</font>' },
                                            'esc character'          => { 'start' => '<font color="#77dd77">',
                                                                          'stop'  => '</font>' },
                                            'numeric'                => { 'start' => '<font color="#FF0000">',
                                                                          'stop'  => '</font>' },
                                           
                                            'identifier'             => { 'start' => '<font color="#B0B0B0">',
                                                                          'stop'  => '</font>' },
                                            'predefined identifier'  => { 'start' => '<font color="#2040a0"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                     
                                            'type'                   => { 'start' => '<font color="#2040a0"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'predefined type'        => { 'start' => '<font color="#2040a0"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            
                                            'reserved word'          => { 'start' => '<strong>',
                                                                          'stop'  => '</strong>' },
                                            'library function'       => { 'start' => '<font color="a52a2a"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            
                                            'include'                => { 'start' => '<font color="#00FF00">',
                                                                          'stop'  => '</font>' },
                                            'preprocessor'           => { 'start' => '<font color="#00FF00">',
                                                                          'stop'  => '</font>' },
                                            
                                            'braces'                 => { 'start' => '<font color="darkCyan"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'symbol'                 => { 'start' => '<font color="darkCyan">',
                                                                          'stop'  => '</font>' },

                                            'function header'        => { 'start' => '<strong>',
                                                                          'stop'  => '</strong>' },
                                            'function header name'   => { 'start' => '<font color="ff0000">',
                                                                          'stop'  => '</font>' },
                                            'function header args'   => { 'start' => '<font color="2040a0">',
                                                                          'stop'  => '</font>' },
                                            
                                            'regex'                  => { 'start' => '<font color="b000d0">',
                                                                          'stop'  => '</font>' },
                                            
                                            'text'                   => { 'start' => '<i>',
                                                                          'stop'  => '</i>'},

                                            # HTML
                                            'entity'                 => { 'start' => '<font color="ff0000">',
                                                                          'stop'  => '</font>' },

                                            # MAKEFILE
                                            'assignment'             => { 'start' => '<font color="2040a0">',
                                                                          'stop'  => '</font>' },
                                            'dependency line'        => { 'start' => '<font color="8b2252">',
                                                                          'stop'  => '</font>' },
                                            'dependency target'      => { 'start' => '<strong>',
                                                                          'stop'  => '</strong>' },
                                            'dependency continuation'=> { 'start' => '<font color="000000"><strong>',
                                                                          'stop'  => '</strong></font>' },
					    'continuation'           => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
                                            'macro'                  => { 'start' => '<font color="2040a0">',
                                                                          'stop'  => '</font>' },
                                            'int macro'              => { 'start' => '<font color="4080ff">',
                                                                          'stop'  => '</font>' },
                                            'esc $$$'                => { 'start' => '<font color="444444">',
                                                                          'stop'  => '</font>' },

                                            # PATCH
                                            'separator'              => { 'start' => '<font color="00FF00"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'line spec'              => { 'start' => '<font color="FFFF00"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'deletion'               => { 'start' => '<font color="FF0000"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'insertion'              => { 'start' => '<font color="0000FF"><strong>',
                                                                          'stop'  => '</strong></font>' }
                                     };

#####
#
# nocolor
#
%{$STYLESHEET{'html-nocolor'}} = %{$STYLESHEET{'html-nobg'}};
${ $STYLESHEET{'html-nocolor'}} {'tags'} = {
	'comment' => {
		'start' => '<I>',
		'stop'  => '</I>'
	},
	'doc comment' => {
		'start' => '',
		'stop'  => ''
	},
	'string' => {
		'start' => '<I>',
		'stop'  => '</I>'
	},
	'esc string' => {
		'start' => '',
		'stop'  => ''
	},
	'character' => {
		'start' => '',
		'stop'  => ''
	},
	'esc character' => {
		'start' => '',
		'stop'  => ''
	},
	'numeric' => {
		'start' => '',
		'stop'  => ''
	},
	'identifier' => {
		'start' => '<B>',
		'stop'  => '</B>'
	},
	'predefined identifier' => {
		'start' => '<U>',
		'stop'  => '</U>'
	},
	'type' => {
		'start' => '',
		'stop'  => ''
	},
	'predefined type' => {
		'start' => '<U>',
		'stop'  => '</U>'
	},
	'reserved word' => {
		'start' => '',
		'stop'  => ''
	},
	'library function' => {
		'start' => '',
		'stop'  => ''
	},
	'include' => {
		'start' => '',
		'stop'  => ''
	},
	'preprocessor' => {
		'start' => '<B>',
		'stop'  => '</B>'
	},
	'braces' => {
		'start' => '',
		'stop'  => ''
	},
	'symbol' => {
		'start' => '',
		'stop'  => ''
	},
	'function header' => {
		'start' => '',
		'stop'  => ''
	},
	'function header name' => {
		'start' => '',
		'stop'  => ''
	},
	'function header args' => {
		'start' => '',
		'stop'  => ''
	},
	'regex' => {
		'start' => '',
		'stop'  => ''
	},
	'text' => {
		'start' => '',
		'stop'  => ''
	},
	# HTML
	'entity' => {
		'start' => '',
		'stop'  => ''
	},
	# MAKEFILE
	'assignment' => {
		'start' => '<B>',
		'stop'  => '</B>'
	},
	'dependency line' => {
		'start' => '',
		'stop'  => ''
	},
	'dependency target' => {
		'start' => '<B>',
		'stop'  => '</B>'
	},
	'dependency continuation' => {
		'start' => '',
		'stop'  => ''
	},
	'continuation' => {
		'start' => '',
		'stop'  => ''
	},
	'macro' => {
		'start' => '<B>',
		'stop'  => '</B>'
	},
	'int macro' => {
		'start' => '<B>',
		'stop'  => '</B>'
	},
	'esc $$$' => {
		'start' => '',
		'stop'  => ''
	},
	# PATCH
	'separator'              => { 
		'start' => '<B>',
		'stop'  => '</B>' },
	'line spec'              => {
		'start' => '',
		'stop'  => ''
	},
	'deletion'               => {
		'start' => '',
		'stop'  => ''
	},
	'insertion'              => { 
		'start' => '',
		'stop'  => ''
	}
};



#####
#
# simple
#
%{$STYLESHEET{'html-simple'}} = %{$STYLESHEET{'html-nocolor'}};
${ $STYLESHEET{'html-simple'}} {'template'} = '<html>
  <head>
    <title>%%title%%</title>
  </head>

  <body>
    <h2>%%title%%</h2>
    <pre>
%%code%%
    </pre>
  </body
</html>';




# Vincent Sanders <vince@trinity.fluff.org>
# html-fntlck is a modification of html-light
# in such a way, that the body tag does define
# different colors and that the <font> colors are different.
#it is supposed to be the colours i get from emacs default font-lock mode

%{$STYLESHEET{'html-fntlck'}} = %{$STYLESHEET{'html-light'}};
${ $STYLESHEET{'html-fntlck'}} {'template'} = '<html>
<head>
  <title>%%title%%</title>
</head>
<body bgcolor="#ffffff"  text="#000000" vlink="#0000ff" alink="#00FF00" link="#000000">
<pre>
%%code%%
</pre>
<hr>
syntax highlighted by <a href="http://www.palfrader.org/code2html">Code2HTML</a>, v. %%version%%
</body>
</html>
';
${ $STYLESHEET{'html-fntlck'}} {'tags'} = {
					    'comment'                => { 'start' => '<font color="bb0000">',
									  'stop'  => '</font>' },
					    'doc comment'            => { 'start' => '<font color="bb0000"><i>',
									  'stop'  => '</i></font>' },
					    'string'                 => { 'start' => '<font color="#bb7766">',
									  'stop'  => '</font>' },
					    'esc string'             => { 'start' => '<font color="#cc8877">',
									  'stop'  => '</font>' },
					    'character'              => { 'start' => '<font color="#bb7766">',
									  'stop'  => '</font>' },
					    'esc character'          => { 'start' => '<font color="#cc8877">',
									  'stop'  => '</font>' },
					    'numeric'                => { 'start' => '<font color="#FF0000">',
									  'stop'  => '</font>' },

					    'identifier'             => { 'start' => '<font color="#0000ff">',
									  'stop'  => '</font>' },
					    'predefined identifier'  => { 'start' => '<font color="#2040a0"><strong>',
									  'stop'  => '</strong></font>' },

					    'type'                   => { 'start' => '<font color="#2040a0"><strong>',
									  'stop'  => '</strong></font>' },
					    'predefined type'        => { 'start' => '<font color="#2040a0"><strong>',
									  'stop'  => '</strong></font>' },

					    'reserved word'          => { 'start' => '<font color="b000e0">',
									  'stop'  => '</font>' },
					    'library function'       => { 'start' => '<font color="a52a2a"><strong>',
									  'stop'  => '</strong></font>' },

					    'include'                => { 'start' => '<font color="0000ff"><strong>',
									  'stop'  => '</strong></font>' },
					    'preprocessor'           => { 'start' => '<font color="0000ff"><strong>',
									  'stop'  => '</strong></font>' },

					    'braces'                 => { 'start' => '<font color="4444FF"><strong>',
									  'stop'  => '</strong></font>' },
					    'symbol'                 => { 'start' => '<font color="000000">',
							 	  	  'stop'  => '</font>' },

					    'function header'        => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
					    'function header name'   => { 'start' => '<font color="ff0000">',
									  'stop'  => '</font>' },
					    'function header args'   => { 'start' => '<font color="2040a0">',
									  'stop'  => '</font>' },

					    'regex'                  => { 'start' => '<font color="b000d0">',
									  'stop'  => '</font>' },

					    'text'                   => { 'start' => '<i>',
									  'stop'  => '</i>'},

					    # HTML
					    'entity'                 => { 'start' => '<font color="ff0000">',
									  'stop'  => '</font>' },

					    # MAKEFILE
					    'assignment'             => { 'start' => '<font color="2040a0">',
									  'stop'  => '</font>' },
					    'dependency line'        => { 'start' => '<font color="8b2252">',
									  'stop'  => '</font>' },
					    'dependency target'      => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
					    'dependency continuation'=> { 'start' => '<font color="000000"><strong>',
									  'stop'  => '</strong></font>' },
					    'continuation'           => { 'start' => '<strong>',
									  'stop'  => '</strong>' },
					    'macro'                  => { 'start' => '<font color="2040a0">',
									  'stop'  => '</font>' },
					    'int macro'              => { 'start' => '<font color="4080ff">',
									  'stop'  => '</font>' },
					    'esc $$$'                => { 'start' => '<font color="444444">',
							                  'stop'  => '</font>' },
					    
                                            # PATCH
                                            'separator'              => { 'start' => '<font color="00A040"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'line spec'              => { 'start' => '<font color="A0A000"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'deletion'               => { 'start' => '<font color="FF0000"><strong>',
                                                                          'stop'  => '</strong></font>' },
                                            'insertion'              => { 'start' => '<font color="0000FF"><strong>',
                                                                          'stop'  => '</strong></font>' }

                                     };


return \%STYLESHEET;

};



################################################################################
####################### get_default_database ###################################
################################################################################
sub get_default_database
{

my %LANGUAGE;

# written by PP
$LANGUAGE{'plain'}      = {
                            'filename'   => '',
                            'regex'      => '',
                            'patterns'   => []
                          };
 
 




# taken from nedit
# modified by PP
$LANGUAGE{'ada'}        = {
                            'filename'   => '(?i)\\.a(d[asb]?)?$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'Comments',
                                                'regex'      => '--.*?$',
                                                'style'      => 'comment',
                                                'childregex' => [],
                                              },
                                              {
                                                'name'       => 'String Literals',
                                                'regex'      => '".*?("|$)',
					        'style'      => 'string',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Character Literals',
                                                'regex'      => '\'.\'',
					        'style'      => 'character',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Ada Attributes',
                                                'regex'      => '\'[a-zA-Z][a-zA-Z_]+\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Numeric Literals',
                                                'regex'      => '(((2|8|10|16)#[_0-9a-fA-F]*#)|[0-9.]+)',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Withs Pragmas Use',
                                                'regex'      => '\\b(?i)((with|pragma|use)[ \\t\\n\\f\\r]+[a-zA-Z0-9_.]+;)+\\b',
					        'style'      => 'include',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Predefined Types',
                                                'regex'      => '\\b(?i)(boolean|character|count|duration|float|integer|long_float|long_integer|priority|short_float|short_integer|string)\\b',
					        'style'      => 'predefined type',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Predefined Subtypes',
                                                'regex'      => '\\b(?i)field|natural|number_base|positive|priority\\b',
					        'style'      => 'predefined type',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Reserved Words',
                                                'regex'      => '\\b(?i)(abort|abs|accept|access|and|array|at|begin|body|case|constant|declare|delay|delta|digits|do|else|elsif|end|entry|exception|exit|for|function|generic|goto|if|in|is|limited|loop|mod|new|not|null|of|or|others|out|package|pragma|private|procedure|raise|range|record|rem|renames|return|reverse|select|separate|subtype|task|terminate|then|type|use|when|while|with|xor)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Ada 95 Only',
                                                'regex'      => '\\b(?i)(abstract|tagged|all|protected|aliased|requeue|until)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Identifiers',
                                                'regex'      => '\\b[a-zA-Z][a-zA-Z0-9_]*\\b',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Dot All',
                                                'regex'      => '(?i)\\.all\\b',
					        'style'      => 'predefined identifier',
                                                'childregex' => []
                                              }
                                            ]
                          };
$LANGUAGE{'ada95'}      = $LANGUAGE{'ada'};















# written by JA
$LANGUAGE{'awk'}       =  {
                            'filename'   => '(?i)\\.awk$',
                            'regex'      => '^\\s*#\\s*![^\\s]*awk',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '#.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '\'\'|\'.*?([^\\\\](\\\\\\\\)*)\'|\'\\\\\\\\\'',
#                                                'regex'      => '\'\'|\'\\\\\\\\\'|\'[^\'\\\\]\'|\'[^\'].*?[^\\\\]\'',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'function header',
                                                'regex'      => 'function[\\t ]+([a-zA-Z0-9_]+)[\\t \\n]*(\\{|\\n)',
					        'style'      => 'function header',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'function coloring',
                                                                    'regex'      => '[\\t ]([a-zA-Z0-9_]+)',
					                            'style'      => 'function header name',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'regex matching I 1',
                                                'regex'      => '(\\b| )?(/)(\\\\/|[^/\\n])*(/[gimesox]*)',
					        'style'      => 'regex',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'regex matching I 2',
                                                'regex'      => '(?:\\b| )(?:(?:m|q|qq)([!"#$%&\'*+-/]))(\\\\\\2|[^\\2\\n])*(\\2[gimesox]*)',
					        'style'      => 'regex',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'regex matching II',
                                                'regex'      => '(?:\\b| )?(?:s([!"#$%&\'*+-/]))(?:\\\\\\2|[^\\2\\n])*?(\\2)[^(\\2)\\n]*?(\\2[gimesox]*)',
					        'style'      => 'regex',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'translate',
                                                'regex'      => '(?:\\b| )(?:(?:tr|y)([^\w\s]))(?:\\\\\\2|[^\\2\\n])*?(\\2)[^(\\2)\\n]*?(\\2[gimesox]*)',
					        'style'      => 'regex',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keywords',
                                                'regex'      => '\\b(BEGIN|END|ARGC|ARGIND|ARGV|CONVFMT|ENVIRON|ERRNO|FIELDWIDTHS|FILENAME|FNR|FS|IGNORECASE|NF|NR|OFMT|OFS|ORS|RS|RT|RSTART|RLENGTH|SUBSEP)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keywords 2',
                                                'regex'      => '\\b(if|while|do|for|in|break|continue|delete|exit|next|nextfile|function)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'library fns',
                                                'regex'      => '\\b(close|getline|print|printf|system|fflush|atan2|cos|exp|int|log|rand|sin|sqrt|srand|gensub|gsub|index|length|split|sprintf|sub|substr|tolower|toupper|systime|strftime)\\b',
					        'style'      => 'library function',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces and parens',
                                                'regex'      => '[\\[\\]\\{\\}\\(\\)]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => '<< stuff',
                                                'regex'      => '<<\'([^\\n]*)\';.*?^\\2$',
					        'style'      => 'text',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => '<< stuff',
                                                'regex'      => '<<([^\\n]*).*?^\\2$',
					        'style'      => 'text',
                                                'childregex' => []
                                              }
                                            ]
                           };
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# taken from nedit
# modified by PP
$LANGUAGE{'c'}          = {
                            'filename'   => '\\.[ch]$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'doc comment',
                                                'regex'      => '/\\*\\*.*?\\*/',
					        'style'      => 'doc comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '/\\*.*?\\*/',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'preprocessor line',
                                                'regex'      => '^[ \\t]*#.*?$',
					        'style'      => 'preprocessor',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'string',
                                                                    'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                                    'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					                            'style'      => 'string',
                                                                    'childregex' => [
                                                                                      {
                                                                                        'name'       => 'esc character',
                                                                                        'regex'      => '\\\\.',
					                                                'style'      => 'esc character',
                                                                                        'childregex' => []
                                                                                      }
                                                                                    ]
                                                                  },
                                                                  {
                                                                    'name'       => '<files>',
                                                                    'regex'      => '<.*?>',
					                            'style'      => 'string',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'comment',
                                                                    'regex'      => '[^/]/\\*.*?\\*/',
					                            'style'      => 'comment',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'character constant',
                                                'regex'      => '\'(\\\\)?.\'',
					        'style'      => 'character',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character', 
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b((0(x|X)[0-9a-fA-F]*)|(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))((e|E)(\\+|-)?[0-9]+)?)(L|l|UL|ul|u|U|F|f)?\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'storage keyword',
                                                'regex'      => '\\b(const|extern|auto|register|static|unsigned|signed|volatile|char|double|float|int|long|short|void|typedef|struct|union|enum)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keyword',
                                                'regex'      => '\\b(return|goto|if|else|case|default|switch|break|continue|while|do|for|sizeof)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces',
                                                'regex'      => '[\\{\\}]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'symbols',
                                                'regex'      => '([\\*\\-\\+=:;%&\\|<>\\(\\)\\[\\]!])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              },
                                              { 
                                                'name'       => 'identifiers',
                                                'regex'      => '([a-zA-Z_][a-zA-Z_0-9]*)',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              }
                                            ]
                          };
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# taken from nedit
# modified by PP
$LANGUAGE{'c++'}        = {
                            'filename'   => '\\.(c(c|pp|xx)|h(h|pp|xx)|C(C|PP|XX)?|H(H|PP|XX)?|i)$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'doc comment',
                                                'regex'      => '/\\*\\*.*?\\*/',
					        'style'      => 'doc comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '/\\*.*?\\*/',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'cplus comment',
                                                'regex'      => '//.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|"\\\\\\\\"|".*?([^\\\\](\\\\\\\\)*)"',
#                                                'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'preprocessor line',
                                                'regex'      => '^[ \\t]*#.*?$',
					        'style'      => 'preprocessor',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'string',
                                                                    'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                                    'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					                            'style'      => 'string',
                                                                    'childregex' => [
                                                                                      {
                                                                                        'name'       => 'esc character',
                                                                                        'regex'      => '\\\\.',
					                                                'style'      => 'esc character',
                                                                                        'childregex' => []
                                                                                      }
                                                                                    ]
                                                                  },
                                                                  {
                                                                    'name'       => '<files>',
                                                                    'regex'      => '<.*?>',
					                            'style'      => 'string',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'comment',
                                                                    'regex'      => '[^/]/\\*.*?\\*/',
					                            'style'      => 'comment',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'cplus comment',
                                                                    'regex'      => '//.*?$',
					                            'style'      => 'comment',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'character constant',
                                                'regex'      => '\'(\\\\)?.\'',
					        'style'      => 'character',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b((0(x|X)[0-9a-fA-F]*)|(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))((e|E)(\\+|-)?[0-9]+)?)(L|l|UL|ul|u|U|F|f)?\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'storage keyword',
                                                'regex'      => '\\b(class|typename|typeid|template|friend|virtual|inline|explicit|operator|overload|public|private|protected|const|extern|auto|register|static|mutable|unsigned|signed|volatile|char|double|float|int|long|short|bool|wchar_t|void|typedef|struct|union|enum)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => [],
                                              },
                                              {
                                                'name'       => 'keyword',
                                                'regex'      => '\\b(new|delete|this|return|goto|if|else|case|default|switch|break|continue|while|do|for|catch|throw|sizeof|true|false|namespace|using|dynamic_cast|static_cast|reinterpret_cast)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces',
                                                'regex'      => '[\\{\\}]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'symbols',
                                                'regex'      => '([\\*\\-\\+=:;%&\\|<>\\(\\)\\[\\]!])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'identifiers',
                                                'regex'      => '([a-zA-Z_][a-zA-Z_0-9]*)',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              }
                                            ]
                          };
$LANGUAGE{'cc'}         = $LANGUAGE{'c++'};
$LANGUAGE{'cpp'}        = $LANGUAGE{'c++'};
$LANGUAGE{'cxx'}        = $LANGUAGE{'c++'};










# taken from nedit
# modified by tk
$LANGUAGE{'f'}          = {
                            'filename'   => '\\.(f(f|or|77)|F(F|OR|77))$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '^[C|c].*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
					        'style'      => 'string',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b([0-9]+(\\.[0-9]*)?([DEde][-+]?[0-9]*)?|\\.[0-9]+([DEde][-+]?[0-9]*)?)\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'storage keyword',
                                                'regex'      => '\\b(BYTE|[Bb]yte|CHARACTER|[Cc]haracter|COMPLEX|[Cc]omplex|DOUBLE *COMPLEX|[Dd]ouble *[Cc]omplex|DOUBLE *PRECISION|[Dd]ouble *[Pp]recision|DOUBLE|[Dd]ouble|INTEGER|[Ii]nteger|REAL|[Rr]eal)(\\*[0-9]+)?\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keyword',
                                                'regex'      => '\\b(ACCEPT|[Aa]ccept|ASSIGN|[Aa]ssign|AUTOMATIC|[Aa]utomatic|BACKSPACE|[Bb]ackspace|BLOCK|[Bb]lock|CALL|[Cc]all|CLOSE|[Cc]lose|COMMON|[Cc]ommon|CONTINUE|[Cc]ontinue|DATA|[Dd]ata|DECODE|[Dd]ecode|DELETE|[Dd]elete|DIMENSION|[Dd]imension|DO|[Dd]o|ELSE|[Ee]lse|ELSEIF|[Ee]lseif|ENCODE|[Ee]ncode|END *FILE|[Ee]nd *[Ff]ile|ENDFILE|[Ee]ndfile|END|[Ee]nd|ENDIF|[Ee]ndif|ENTRY|[Ee]ntry|EQUIVALENCE|[Ee]quivalence|EXIT|[Ee]xit|EXTERNAL|[Ee]xternal|FORMAT|[Ff]ormat|FUNCTION|[Ff]unction|GOTO|[Gg]oto|IF|[Ii]f|IMPLICIT|[Ii]mplicit|INCLUDE|[Ii]nclude|INQUIRE|[Ii]nquire|INTRINSIC|[Ii]ntrinsic|LOGICAL|[Ll]ogical|MAP|[Mm]ap|NONE|[Nn]one|ON|[Oo]n|OPEN|[Oo]pen|PARAMETER|[Pp]arameter|PAUSE|[Pp]ause|POINTER|[Pp]ointer|PRINT|[Pp]rint|PROGRAM|[Pp]rogram|READ|[Rr]ead|RECORD|[Rr]ecord|RETURN|[Rr]eturn|REWIND|[Rr]ewind|SAVE|[Ss]ave|STATIC|[Ss]tatic|STOP|[Ss]top|STRUCTURE|[Ss]tructure|SUBROUTINE|[Ss]ubroutine|SYSTEM|[Ss]ystem|THEN|[Tt]hen|TO|[Tt]o|TYPE|[Tt]ype|UNION|[Uu]nion|UNLOCK|[Uu]nlock|VIRTUAL|[Vv]irtual|VOLATILE|[Vv]olatile|WHILE|[Ww]hile|WRITE|[Ww]rite)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'symbols',
                                                'regex'      => '([\\*\\-\\+=:;%&\\|<>\\(\\)\\[\\]!])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'identifiers',
                                                'regex'      => '([a-zA-Z_][a-zA-Z_0-9]*)',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              }
                                            ]
                          };









# written by VRS
$LANGUAGE{'gpasm'}      = {
                            'filename'   => '(?i)\\.(asm|inc)$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'args',
                                                'regex'      => '^.*$',
					        'style'      => 'symbol',
                                                'childregex' => [
						 {
						     'name'       => 'comment',
						     'regex'      => ';.*?$',
					             'style'      => 'comment',
                                                     'childregex' => []
                                                 },
                                                 {
                                                     'name'       => 'labels',
                                                     'regex'      => '^[A-Za-z_][A-Za-z_0-9]*:?',
					             'style'      => 'identifier',
                                                     'childregex' => []
                                                 },

                                                 {
                                                     'name'       => 'menonics',
                                                     'regex'      => '^[ \t]+[A-Za-z_][A-Za-z_0-9]*',
					             'style'      => 'reserved word',
                                                     'childregex' => []
                                                 },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              }


                                                                 ]
                                              }
                                            ]
                          };








# written by JA
$LANGUAGE{'groff'}      = {
                            'filename'   => '\\.groff$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '\\\\".*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              }
                                            ]
                          };
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# taken from nedit
# modified by PP
$LANGUAGE{'html'}       = {
                            'filename'   => '(?i)\\.html?$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '<!--.*?-->',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'entity',
                                                'regex'      => '\\&[-.a-zA-Z0-9#]*;?',
					        'style'      => 'entity',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'tag',
                                                'regex'      => '<(/|!)?[-.a-zA-Z0-9]*.*?>',
					        'style'      => 'predefined identifier',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'double quote string',
                                                                    'regex'      => '".*?"',
					                            'style'      => 'string',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'single quote string',
                                                                    'regex'      => '\'.*?\'',
					                            'style'      => 'string',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'brackets',
                                                                    'regex'      => '[<>]',
					                            'style'      => 'braces',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'attribute',
                                                                    'regex'      => '[^\'" ]+(?=.)',
					                            'style'      => 'identifier',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              }
                                            ]
                       };
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# taken from nedit
# modified by PP
$LANGUAGE{'java'}       = {
                            'filename'   => '\\.java$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'doc comment',
                                                'regex'      => '/\\*\\*.*?\\*/',
					        'style'      => 'doc comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '/\\*.*?\\*/',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'cplus comment',
                                                'regex'      => '//.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'single quoted',
                                                'regex'      => '\'\'|\'.*?([^\\\\](\\\\\\\\)*)\'|\'\\\\\\\\\'',
#                                                'regex'      => '\'\'|\'\\\\\\\\\'|\'[^\'\\\\]\'|\'[^\'].*?[^\\\\]\'',
					        'style'      => 'string',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b((0(x|X)[0-9a-fA-F]*)|(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))((e|E)(\\+|-)?[0-9]+)?)(L|l|UL|ul|u|U|F|f)?\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'include',
                                                'regex'      => '\\b(import|package)\\b.*?$',
					        'style'      => 'include',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\(.|\\n)',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'comment',
                                                                    'regex'      => '[^/]/\\*.*?\\*/',
					                            'style'      => 'comment',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'storage keyword',
                                                'regex'      => '\\b(abstract|boolean|byte|char|class|double|extends|final|float|int|interface|long|native|private|protected|public|short|static|transient|synchronized|void|volatile|implements)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keyword',
                                                'regex'      => '\\b(break|case|catch|continue|default|do|else|false|finally|for|if|instanceof|new|null|return|super|switch|this|throw|throws|true|try|while)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces and parens',
                                                'regex'      => '[\\{\\}\\(\\)\\[\\]]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Identifiers',
                                                'regex'      => '\\b[a-zA-Z_][a-zA-Z0-9_]*\\b',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'symbols',
                                                'regex'      => '([\\*\\-\\+=:;%&\\|<>!])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              }
                                            ]
                          };
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# taken from nedit
# modified by PP
$LANGUAGE{'javascript'} = {
                            'filename'   => '(?i)\\.js$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '/\\*.*?\\*/',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'cplus comment',
                                                'regex'      => '//.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b((0(x|X)[0-9a-fA-F]*)|(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))((e|E)(\\+|-)?[0-9]+)?)(L|l|UL|ul|u|U|F|f)?\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'events',
                                                'regex'      => '\\b(onAbort|onBlur|onClick|onChange|onDblClick|onDragDrop|onError|onFocus|onKeyDown|onKeyPress|onLoad|onMouseDown|onMouseMove|onMouseOut|onMouseOver|onMouseUp|onMove|onResize|onSelect|onSubmit|onUnload)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces',
                                                'regex'      => '[\\{\\}]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'statements',
                                                'regex'      => '\\b(break|continue|else|for|if|in|new|return|this|typeof|var|while|with)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'function',
                                                'regex'      => 'function[\\t ]+([a-zA-Z0-9_]+)[\\t \\(]+.*?[\\n{]',
					        'style'      => 'function header',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'function args',
                                                                    'regex'      => '\\(.*?\\)',
					                            'style'      => 'function header args',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'function name',
                                                                    'regex'      => '[\\t ][a-zA-Z0-9_]+',
					                            'style'      => 'function header name',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },  
                                              {
                                                'name'       => 'built in object type',
                                                'regex'      => '\\b(anchor|Applet|Area|Array|button|checkbox|Date|document|elements|FileUpload|form|frame|Function|hidden|history|Image|link|location|Math|navigator|Option|password|Plugin|radio|reset|select|string|submit|text|textarea|window)\\b',
					        'style'      => 'predefined type',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '".*?("|$)',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'colors',
                                                                    'regex'      => '(aliceblue|antiquewhite|aqua|aquamarine|azure|beige|bisque|black|blanchedalmond|blue|blueviolet|brown|burlywood|cadetblue|chartreuse|chocolate|coral|cornflowerblue|cornsilk|crimson|cyan|darkblue|darkcyan|darkgoldenrod|darkgray|darkgreen|darkkhaki|darkmagenta|darkolivegreen|darkorange|darkorchid|darkred|darksalmon|darkseagreen|darkslateblue|darkslategray|darkturquoise|darkviolet|deeppink|deepskyblue|dimgray|dodgerblue|firebrick|floralwhite|forestgreen|fuchsia|gainsboro|ghostwhite|gold|goldenrod|gray|green|greenyellow|honeydew|hotpink|indianred|indigo|ivory|khaki|lavender|lavenderblush|lawngreen|lemonchiffon|lightblue|lightcoral|lightcyan|lightgoldenrodyellow|lightgreen|lightgrey|lightpink|lightsalmon|lightseagreen|lightskyblue|lightslategray|lightsteelblue|lightyellow|lime|limegreen|linen|magenta|#008000|mediumaquamarine|mediumblue|mediumorchid|mediumpurple|mediumseagreen|mediumslateblue|mediumspringgreen|mediumturquoise|mediumvioletred|midnightblue|mintcream|mistyrose|moccasin|navajowhite|navy|oldlace|olive|olivedrab|orange|orangered|orchid|palegoldenrod|palegreen|paleturquoise|palevioletred|papayawhip|peachpuff|peru|pink|plum|powderblue|purple|red|rosybrown|royalblue|saddlebrown|salmon|sandybrown|seagreen|seashell|sienna|silver|skyblue|slateblue|slategray|snow|springgreen|steelblue|tan|teal|thistle|tomato|turquoise|violet|wheat|white|whitesmoke|yellow|yellowgreen|#[A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9])',
					                            'style'      => 'identifier',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '\'.*?(\'|$)',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'colors',
                                                                    'regex'      => '(aliceblue|antiquewhite|aqua|aquamarine|azure|beige|bisque|black|blanchedalmond|blue|blueviolet|brown|burlywood|cadetblue|chartreuse|chocolate|coral|cornflowerblue|cornsilk|crimson|cyan|darkblue|darkcyan|darkgoldenrod|darkgray|darkgreen|darkkhaki|darkmagenta|darkolivegreen|darkorange|darkorchid|darkred|darksalmon|darkseagreen|darkslateblue|darkslategray|darkturquoise|darkviolet|deeppink|deepskyblue|dimgray|dodgerblue|firebrick|floralwhite|forestgreen|fuchsia|gainsboro|ghostwhite|gold|goldenrod|gray|green|greenyellow|honeydew|hotpink|indianred|indigo|ivory|khaki|lavender|lavenderblush|lawngreen|lemonchiffon|lightblue|lightcoral|lightcyan|lightgoldenrodyellow|lightgreen|lightgrey|lightpink|lightsalmon|lightseagreen|lightskyblue|lightslategray|lightsteelblue|lightyellow|lime|limegreen|linen|magenta|#008000|mediumaquamarine|mediumblue|mediumorchid|mediumpurple|mediumseagreen|mediumslateblue|mediumspringgreen|mediumturquoise|mediumvioletred|midnightblue|mintcream|mistyrose|moccasin|navajowhite|navy|oldlace|olive|olivedrab|orange|orangered|orchid|palegoldenrod|palegreen|paleturquoise|palevioletred|papayawhip|peachpuff|peru|pink|plum|powderblue|purple|red|rosybrown|royalblue|saddlebrown|salmon|sandybrown|seagreen|seashell|sienna|silver|skyblue|slateblue|slategray|snow|springgreen|steelblue|tan|teal|thistle|tomato|turquoise|violet|wheat|white|whitesmoke|yellow|yellowgreen|#[A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9][A-Fa-f0-9])',
					                            'style'      => 'identifier',
                                                                    'childregex' => [],
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'event capturing',
                                                'regex'      => '\\b(captureEvents|releaseEvents|routeEvent|handleEvent)\\b.*?(\\)|$)',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'predefined methods',
                                                'regex'      => '\\b(abs|acos|alert|anchor|asin|atan|atan2|back|big|blink|blur|bold|ceil|charAt|clear|clearTimeout|click|close|confirm|cos|escape|eval|exp|fixed|floor|focus|fontcolor|fontsize|forward|getDate|getDay|getHours|getMinutes|getMonth|getSeconds|getTime|getTimezoneOffset|getYear|go|indexOf|isNaN|italics|javaEnabled|join|lastIndexOf|link|log|max|min|open|parse|parseFloat|parseInt|pow|prompt|random|reload|replace|reset|reverse|round|scroll|select|setDate|setHours|setMinutes|setMonth|setSeconds|setTimeout|setTime|setYear|sin|small|sort|split|sqrt|strike|sub|submit|substring|sup|taint|tan|toGMTString|toLocaleString|toLowerCase|toString|toUpperCase|unescape|untaint|UTC|write|writeln)\\b',
					        'style'      => 'library function',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'properties',
                                                'regex'      => '\\b(action|alinkColor|anchors|appCodeName|appName|appVersion|bgColor|border|checked|complete|cookie|defaultChecked|defaultSelected|defaultStatus|defaultValue|description|E|elements|enabledPlugin|encoding|fgColor|filename|forms|frames|hash|height|host|hostname|href|hspace|index|lastModified|length|linkColor|links|LN2|LN10|LOG2E|LOG10E|lowsrc|method|name|opener|options|parent|pathname|PI|port|protocol|prototype|referrer|search|selected|selectedIndex|self|SQRT1_2|SQRT2|src|status|target|text|title|top|type|URL|userAgent|value|vlinkColor|vspace|width|window)\\b',
					        'style'      => 'predefined identifier',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'operators',
                                                'regex'      => '([=;->/&|])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              }
                                            ]
                          };
$LANGUAGE{'js'}         = $LANGUAGE{'javascript'};








# written by Andreas Krennmair
# extremely incomplete

$LANGUAGE{'lisp'}       = {
                            'filename' => '\\.(lsp|l)$',
			    'regex' => '',
                            'patterns' => [
                               {
                                 'name'       => 'parens',
                                 'regex'      => '[()]',
                                 'style'      => 'braces',
                                 'childregex' => []
                               },
                               {
                                 'name'       => 'comment',
                                 'regex'      => ';.*?$',
                                 'style'      => 'comment',
                                 'childregex' => []
                               },
                               {
                                 'name'       => 'string',
                                 'regex'      => '".*?("|$)',
                                 'style'      => 'string',
                                 'childregex' => []
                               },
                               {
                                 'name'       => 'keywords',
                                 'regex'      => '\\b(defun |xyz)\\b',
                                 'style'      => 'reserved word',
                                 'childregex' => []
                               },
                               {
                                 'name'       => 'numeric constant',
                                 'regex'      => '(#\([0-9]+ [0-9]+\)|[0-9]+)',
                                 'style'      => 'numeric',
                                 'childregex' => []
                               },
                               {
                                 'name'       => 'identifiers',
                                 'regex'      => '([-a-zA-Z]+)',
                                 'style'      => 'identifier',
                                 'childregex' => []
                               }
                            ]
                          };










# written by JA
$LANGUAGE{'m4'}         = {
                            'filename'   => '\\.m4$',
                            'regex'      => '',
                            'patterns' => [
                                            {
                                              'regex'      => 'dnl.*?$',
					      'style'      => 'doc comment',
                                              'childregex' => []
                                            },
                                            {
                                              'regex'      => '#.*?$',
					      'style'      => 'comment',
                                              'childregex' => []
                                            },
                                            {
                                              'regex'      => '\\b(define|undefine|defn|pushdef|popdef|indir|builtin|changequote|changecom|changeword|m4wrap|m4exit|include|sinclude|divert|undivert|divnum|cleardiv|shift|dumpdef|traceon|traceoff|debugfile|debugmode|len|index|regexp|substr|translit|patsubst|format|incr|decr|syscmd|esyscmd|sysval|maketemp|errprint)\\b',
					      'style'      => 'reserved word',
                                              'childregex' => []
                                            },
                                            {
                                              'regex'      => '\\b(ifdef|ifelse|loops)\\b',
					      'style'      => 'reserved word',
                                              'childregex' => [
                                                                {
                                                                  'regex'      => '[$]\\$?({[^}]*}|[^a-zA-Z0-9_/\\t\\n\\.,\\\\[\\\\{\\\\(]|[0-9]+|[a-zA-Z_][a-zA-Z0-9_]*)?',
					                          'style'      => 'identifier',
                                                                  'childregex' => []
                                                                }
                                                              ]
                                            }
                                          ]
                          };















# taken from nedit
# modified by PP
$LANGUAGE{'make'}       = {
                            'filename'   => '[Mm]akefile.*',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'Comment',
                                                'regex'      => '#.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Assignment',
                                                'regex'      => '^( *| [ \\t]*)[A-Za-z0-9_+]*[ \\t]*(\\+|:)?=',
					        'style'      => 'assignment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Dependency Line',
                                                'regex'      => '^ *([A-Za-z0-9./$(){} _%+-]|\\n)*::?',
					        'style'      => 'dependency line',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'Dependency Target',
                                                                    'regex'      => '[A-Za-z0-9./$(){} _%+-]+',
								    'style'      => 'dependency target',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'Dependency Continuation',
                                                                    'regex'      => '\\\\\\n',
								    'style'      => 'dependency continuation',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'comment',
                                                                    'regex'      => '#.*?$',
								    'style'      => 'comment',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'macro',
                                                                    'regex'      => '\\$([A-Za-z0-9_]|\\([^)]*\\)|{[^}]*})',
								    'style'      => 'macro',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'int macro',
                                                                    'regex'      => '\\$([<@*?%]|\\$@)',
								    'style'      => 'int macro',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'Continuation',
                                                'regex'      => '\\\\$',
						'style'      => 'continuation',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Macro',
                                                'regex'      => '\\$([A-Za-z0-9_]|\\([^)]*\\)|{[^}]*})',
						'style'      => 'macro',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Internal Macro',
                                                'regex'      => '\\$([<@*?%]|\\$@)',
						'style'      => 'int macro',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Escaped $$$',
                                                'regex'      => '\\$\\$',
						'style'      => 'esc $$$',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'Include',
                                                'regex'      => '^include[ \\t]',
						'style'      => 'include',
                                                'childregex' => []
                                              }
                                            ]
                          };
$LANGUAGE{'makefile'} = $LANGUAGE{'make'};















# taken from nedit
# modified by PP
$LANGUAGE{'pas'}        = {
                            'filename'   => '(?i)\\.p(as)?$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment1 (*    *)',
                                                'regex'      => '\\(\\*.*?\\*\\)',
						'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'comment2 {    }',
                                                'regex'      => '\\{.*?\\}',
						'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '\'.*?(\'|$)',
						'style'      => 'string',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'preprocessor line',
                                                'regex'      => '^[ \\t]*#.*?$',
						'style'      => 'preprocessor',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'comment1 (*    *)',
                                                                    'regex'      => '\\(\\*.*?\\*\\)',
								    'style'      => 'comment',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'comment2 {    }',
                                                                    'regex'      => '\\{.*?\\}',
								    'style'      => 'comment',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'character constant',
                                                'regex'      => '\'.\'',
					        'style'      => 'character',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b((0(x|X)[0-9a-fA-F]*)|[0-9.]+((e|E)(\\+|-)?)?[0-9]*)(L|l|UL|ul|u|U|F|f)?\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'storage and ops',
                                                'regex'      => '\\b(?i)(and|array|const|div|export|file|function|import|in|label|mod|module|nil|not|only|or|packed|pow|pragma|procedure|program|protected|qualified|record|restricted|set|type|var)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keywords',
                                                'regex'      => '\\b(?i)(begin|case|do|downto|else|end|for|goto|if|of|otherwise|repeat|then|to|until|while|with)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'sumbols',
                                                'regex'      => '([\\*\\-\\+=:;<>\\(\\)\\[\\]!]|[^/]/[^/])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'identifiers',
                                                'regex'      => '([a-zA-Z_][a-zA-Z_0-9.^]*[a-zA-Z_0-9]|[a-zA-Z_][a-zA-Z_0-9]*)',
					        'style'      => 'identifier',
                                                'childregex' => [
                                                                  {
                                                                    'regex'      => '(\\.|\\^)+',
								    'style'      => 'symbol',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              }
                                            ],
                          };
$LANGUAGE{'pascal'}     = $LANGUAGE{'pas'};

 
 
 
 
 
 
 
 
 
 
 
 
 
 
# taken from nedit
# modified by PP
# modified by BS
# modified by JD
# modified by JP
$LANGUAGE{'perl'}       = {
                            'filename'   => '(?i)\\.p([lm5]|od)$',
                            'regex'      => '^\\s*#\\s*![^\\s]*perl',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '(?:#.*?(?:\r?\n\s*)+)+',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'variables',
                                                'regex'      => '[\\$@%]\\$?(?:{[^}]*}|[^a-zA-Z0-9_/\\t\\n\\.,\\\\[\\\\{\\\\(]|[0-9]+|[a-zA-Z_][a-zA-Z0-9_]*)?',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => '"" string',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
								    'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'variables',
                                                                    'regex'      => '[\\$@%]\\$?(?:{[^}]*}|[^a-zA-Z0-9_/\\t\\n\\.,\\\\[\\\\{\\\\(]|[0-9]+|[a-zA-Z_][a-zA-Z0-9_]*)?',
								    'style'      => 'identifier',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => '\'\' string',
                                                'regex'      => '\'\'|\'.*?([^\\\\](\\\\\\\\)*)\'|\'\\\\\\\\\'',
#                                                'regex'      => '\'\'|\'\\\\\\\\\'|\'[^\'\\\\]\'|\'[^\'].*?[^\\\\]\'',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
								    'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'more strings - q// qw//',
                                                'regex'      => '(?:\\b| )(?:q|qw)([^\w\s])(?:\\\\\\2|[^\\2\\n])*\\2',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
								    'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'more strings - qq// qx//',
                                                'regex'      => '(?:\\b| )(?:qq|qx)([^\w\s])(?:\\\\\\2|[^\\2\\n])*\\2',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
								    'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'variables',
                                                                    'regex'      => '[\\$@%]\\$?(?:{[^}]*}|[^a-zA-Z0-9_/\\t\\n\\.,\\\\[\\\\{\\\\(]|[0-9]+|[a-zA-Z_][a-zA-Z0-9_]*)?',
								    'style'      => 'identifier',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'subroutine header',
                                                'regex'      => 'sub[\\t ]+(?:[a-zA-Z0-9_]+)[\\t \\n]*(?:\\{|\\(|\\n)',
					        'style'      => 'function header',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'subroutine header coloring',
                                                                    'regex'      => '[\\t ][a-zA-Z0-9_]+',
								    'style'      => 'function header name',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'regex matching I',
                                                'regex'      => '(?:\\b| )?(?:/(?:\\\\/|[^/\\n])*(?:/[gimesox]*)|s([^\w\s])(?:\\\\\\2|[^\\2\\n])*?(\\2)[^(\\2)\\n]*?(\\2[gimesox]*))',
					        'style'      => 'regex',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'regex matching II',
                                                'regex'      => '(?:\\b| )(?:m|qq?|tr|y)([^\w\s])(?:\\\\\\2|[^\\2\\n])*(?:\\2[gimesox]*)',
					        'style'      => 'regex',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keywords',
                                                'regex'      => '\\b(my|local|new|if|until|while|elsif|else|eval|unless|for|foreach|continue|exit|die|last|goto|next|redo|return|local|exec|do|use|require|package|eval|BEGIN|END|eq|ne|not|\\|\\||\\&\\&|and|or)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'library functions',
                                                'regex'      => '\\b(?:a(?:bs|ccept|larm|tan2)|b(?:ind|inmode|less)|c(?:aller|hdir|hmod|homp|hop|hr|hroot|hown|losedir|lose|onnect|os|rypt)|d(?:bmclose|bmopen|efined|elete|ie|ump)|e(?:ach|nd(?:grent|hostent|netent|protoent|pwent|servent)|of|xec|xists|xp)|f(?:ctnl|ileno|lock|ork|ormat|ormline)|g(?:et(?:c|grent|grgid|grnam|hostbyaddr|hostbyname|hostent|login|netbyaddr|netbyname|netent|peername|pgrp|ppid|priority|protobyname|protobynumber|protoent|pwent|pwnam|pwuid|servbyname|servbyport|servent|sockname|sockopt)|lob|mtime|rep)|hex|i(?:mport|ndex|nt|octl)|join|keys|kill|l(?:cfirst|c|ength|ink|isten|og|ocaltime|stat)|m(?:ap|kdir|sgctl|sgget|sgrcv)|no|o(?:ct|pendir|pen|rd)|p(?:ack|ipe|op|os|rintf|rint|ush)|quotemeta|r(?:and|eaddir|ead|eadlink|ecv|ef|ename|eset|everse|ewinddir|index|mdir)|s(?:calar|eekdir|eek|elect|emctl|emget|emop|end|et(?:grent|hostent|netent|pgrp|priority|protoent|pwent|sockopt)|hift|hmctl|hmget|hmread|hmwrite|hutdown|in|leep|ocket|ocketpair|ort|plice|plit|printf|qrt|rand|tat|tudy|ubstr|ymlink|yscall|ysopen|ysread|ystem|yswrite)|t(?:elldir|ell|ie|ied|ime|imes|runcate)|u(?:c|cfirst|mask|ndef|nlink|npack|nshift|ntie|time)|values|vec|w(?:ait|aitpid|antarray|arn|rite)|qw|-[rwxoRWXOezsfdlpSbctugkTBMAC])\\b',
					        'style'      => 'library function',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces, parens and brakets',
                                                'regex'      => '[\\[\\]\\{\\}\\(\\)]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => '<< stuff',
                                                'regex'      => '<<(?:("|\')([^\\n]*)\\2|\\w*).*?^\\3$',
					        'style'      => 'text',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'POD',
                                                'regex'      => '^=.*?^(?:=cut|\\Z)',
					        'style'      => 'doc comment',
                                                'childregex' => []
                                              }
                                            ]
                          };
 














# Thanks to Matt Giwer <jull43@ij.net>
$LANGUAGE{'pov'}        = {
                            'filename'   => '(?i)\\.pov$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'doc comment',
                                                'regex'      => '/\\*\\*.*?\\*/',
					        'style'      => 'doc comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '/\\*.*?\\*/',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'cplus comment',
                                                'regex'      => '//.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'preprocessor line',
                                                'regex'      => '^[ \\t]*#.*?$',
					        'style'      => 'preprocessor',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'string',
                                                                    'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
#                                                                    'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					                            'style'      => 'string',
                                                                    'childregex' => [
                                                                                      {
                                                                                        'name'       => 'esc character',
                                                                                        'regex'      => '\\\\.',
					                                                'style'      => 'esc character',
                                                                                        'childregex' => []
                                                                                      }
                                                                                    ]
                                                                  },
                                                                  {
                                                                    'name'       => '<files>',
                                                                    'regex'      => '<.*?>',
					                            'style'      => 'string',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'comment',
                                                                    'regex'      => '[^/]/\\*.*?\\*/',
					                            'style'      => 'comment',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'cplus comment',
                                                                    'regex'      => '//.*?$',
				                    	        'style'      => 'comment',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'character constant',
                                                'regex'      => '\'(\\\\)?.\'',
					        'style'      => 'character',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character', 
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b((0(x|X)[0-9a-fA-F]*)|(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))((e|E)(\\+|-)?[0-9]+)?)(L|l|UL|ul|u|U|F|f)?\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keyword',
                                                'regex'      => '\\b(abs|absorption|acos|acosh|adaptive|adc_bailout|agate|agate_turb|all|alpha|ambient|ambient_light|angle|aperture|append|arc_angle|area_light|array|asc|asin|asinh|assumed_gamma|atan|atan2|atanh|average|background|bezier_spline|bicubic_patch|black_hole|blob|blue|blur_samples|bounded_by|box|boxed|bozo|break|brick|brick_size|brightness|brilliance|bumps|bump_map|bump_size|camera|case|caustics|ceil|checker|chr|clipped_by|clock|clock_delta|color|color_map|colour|colour_map|component|composite|concat|cone|confidence|conic_sweep|control0|control1|cos|cosh|count|crackle|crand|cube|cubic|cubic_spline|cubic_wave|cylinder|cylindrical|debug|declare|default|defined|degrees|density|density_file|density_map|dents|difference|diffuse|dimensions|dimension_size|direction|disc|distance|distance_maximum|div|eccentricity|else|emission|end|error|error_bound|exp|extinction|fade_distance|fade_power|falloff|falloff_angle|false|fclose|file_exists|filter|finish|fisheye|flatness|flip|floor|focal_point|fog|fog_alt|fog_offset|fog_type|fopen|frequency|gif|global_settings|gradient|granite|gray_threshold|green|height_field|hexagon|hf_gray_16|hierarchy|hollow|hypercomplex|if|ifdef|iff|ifndef|image_map|include|int|interior|interpolate|intersection|intervals|inverse|ior|irid|irid_wavelength|jitter|julia_fractal|lambda|lathe|leopard|light_source|linear_spline|linear_sweep|local|location|log|looks_like|look_at|low_error_factor|macro|mandel|map_type|marble|material|material_map|matrix|max|max_intersections|max_iteration|max_trace_level|media|media_attenuation|media_interaction|merge|mesh|metallic|min|minimum_reuse|mod|mortar|nearest_count|no|normal|normal_map|no_shadow|number_of_waves|object|octaves|off|offset|omega|omnimax|on|once|onion|open|orthographic|panoramic|perspective|pgm|phase|phong|phong_size|pi|pigment|pigment_map|planar|plane|png|point_at|poly|polygon|poly_wave|pot|pow|ppm|precision|prism|pwr|quadratic_spline|quadric|quartic|quaternion|quick_color|quick_colour|quilted|radial|radians|radiosity|radius|rainbow|ramp_wave|rand|range|ratio|read|reciprocal|recursion_limit|red|reflection|reflection_exponent|refraction|render|repeat|rgb|rgbf|rgbft|rgbt|right|ripples|rotate|roughness|samples|scale|scallop_wave|scattering|seed|shadowless|sin|sine_wave|sinh|sky|sky_sphere|slice|slope_map|smooth|smooth_triangle|sor|specular|sphere|spherical|spiral1|spiral2|spotlight|spotted|sqr|sqrt|statistics|str|strcmp|strength|strlen|strlwr|strupr|sturm|substr|superellipsoid|switch|sys|t|tan|tanh|text|texture|texture_map|tga|thickness|threshold|tightness|tile2|tiles|torus|track|transform|translate|transmit|triangle|triangle_wave|true|ttf|turbulence|turb_depth|type|u|ultra_wide_angle|undef|union|up|use_color|use_colour|use_index|u_steps|v|val|variance|vaxis_rotate|vcross|vdot|version|vlength|vnormalize|vrotate|v_steps|warning|warp|water_level|waves|while|width|wood|wrinkles|write|x|y|yes|z)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces',
                                                'regex'      => '[\\{\\}]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'symbols',
                                                'regex'      => '([\\*\\-\\+=:;%&\\|<>\\(\\)\\[\\]!])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              },
                                              { 
                                                'name'       => 'identifiers',
                                                'regex'      => '([a-zA-Z_][a-zA-Z_0-9]*)',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              }
                                            ]
                            };
$LANGUAGE{'povray'}     = $LANGUAGE{'pov'};
   



# by Tom Good 
$LANGUAGE{'python'}        = {
                            'filename'   => '(?i)\\.py$',
                            'regex'      => '^\\s*#\\s*![^\\s]*python',
                            'patterns'   => [
                                              {
                                                'name'       => 'python comment',
                                                'regex'      => '#.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
					      {
                                                'name'       => 'single quote string',
                                                'regex'      => '\'.*?\'',
					        'style'      => 'string',
                                                'childregex' => []
                                              },
                                                            
                                              {
                                                'name'       => 'string',
                                                'regex'      => '""|"\\\\\\\\"|".*?([^\\\\](\\\\\\\\)*)"',
                                                'regex'      => '""|".*?([^\\\\](\\\\\\\\)*)"|"\\\\\\\\"',
                                                'regex'      => '""|"\\\\\\\\"|"[^"\\\\]"|"[^"].*?[^\\\\]"',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'character constant',
                                                'regex'      => '\'(\\\\)?.\'',
					        'style'      => 'character',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '\\\\.',
					                            'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'numeric constant',
                                                'regex'      => '\\b((0(x|X)[0-9a-fA-F]*)|(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))((e|E)(\\+|-)?[0-9]+)?)(L|l|UL|ul|u|U|F|f)?\\b',
					        'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keyword',
                                                'regex'      => '\\b(and|assert|break|class|continue|del|elif|else|except|exec|finally|for|from|global|if|import|in|is|lambda|not|or|pass|print|raise|return|try|while)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces',
                                                'regex'      => '[\\{\\}]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'symbols',
                                                'regex'      => '([\\*\\-\\+=:;%&\\|<>\\(\\)\\[\\]!])',
					        'style'      => 'symbol',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'identifiers',
                                                'regex'      => '([a-zA-Z_][a-zA-Z_0-9]*)',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'function',
                                                'regex'      => '[\\t ]*def[\\t ]+([a-zA-Z0-9_]+)[\\t \\(]+.*?[\\n{]',
					        'style'      => 'function header',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'function args',
                                                                    'regex'      => '\\(.*?\\)',
					                            'style'      => 'function header args',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'function name',
                                                                    'regex'      => '[\\t ][a-zA-Z0-9_]+',
					                            'style'      => 'function header name',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },  
                                              {
                                                'name'       => 'library functions',
                                                'regex'      => '\\b(__import__|abs|apply|buffer|callable|chr|cmp|coerce|compile|complex|delatter|dir|divmod|eval|execfile|filter|float|getattr|globals|hasattr|hash|hex|id|input|int|intern|isinstance|issubclass|len|list|locals|long|map|max|min|oct|open|ord|pow|range|raw_input|reduce|reload|repr|round|setattr|slice|str|tuple|type|unichr|unicode|vars|xrange|zip)\\b',
					        'style'      => 'library function',
                                                'childregex' => []
                                              },
                                            ]
                          };

 

# by Joshua Swink <jswink@pacbell.net>
$LANGUAGE{'ruby'}       = {
                            'filename'   => '\\.rb$',
                            'regex'      => '^\\s*#\\s*![^\\s]*\\bruby\\b',
                            'patterns'   => [
                                              {
                                                'name'       => 'comment',
                                                'regex'      => '(?:#.*?(?:\r?\n\s*)+)+',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'predefined variables',
                                                'regex'      => '(?:\\$(?:[!@&`\'+\\d~=/\\\\,;.<>_*\\$?:"]|DEBUG|FILENAME|LOAD_PATH|stdin|stdout|stderr|VERBOSE|-[0adFiIlpv])|\\b(?:TRUE|FALSE|NIL|STDIN|STDOUT|STDERR|ENV|ARGF|ARGV|DATA|RUBY_VERSION|RUBY_RELEASE_DATE|RUBY_PLATFORM)\\b)',
					        'style'      => 'predefined identifier',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'variables',
                                                'regex'      => '[\\$@](?:{[^}]*}|[^\\w/\\t\\n\\.,\\\\[\\\\{\\\\(]|[0-9]+|[a-zA-Z_][\\w.]*)?',
					        'style'      => 'identifier',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => '"" string',
                                                'regex'      => '""|"(?:\\\\\\\\)+"|".*?(?:[^\\\\](?:\\\\\\\\)*)"|%[Qwx]?([^\\w\\[\\](){}<>])\\2|%[Qwx]?([^\\w\\[\\](){}<>]).*?(?:[^\\\\](?:\\\\\\\\)*)\\3|%[Qwx]?([^\\w\\[\\](){}<>])\\\\\\\\\\4|%[Qwx]?\\[\\]|%[Qwx]?\\[.*?([^\\\\](\\\\\\\\)*)\\]|%[Qwx]?\\[\\\\\\\\\\]|%[Qwx]?\\{\\}|%[Qwx]?\\{.*?([^\\\\](\\\\\\\\)*)\\}|%[Qwx]?\\{\\\\\\\\\\}|%[Qwx]?\\(\\)|%[Qwx]?\\(.*?([^\\\\](\\\\\\\\)*)\\)|%[Qwx]?\\(\\\\\\\\\\)|%[Qwx]?<>|%[Qwx]?<.*?([^\\\\](\\\\\\\\)*)>|%[Qwx]?<\\\\\\\\>',

					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex',     => '\\\\(?:x[\\da-fA-F]{2}|\d\d\d|c.|M-\\\\C-.|M-.|C-.|.)',
								    'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  },
                                                                  {
                                                                    'name'       => 'string expression',
                                                                    'regex'      => '#[\\$\\@][a-zA-Z_][\\w.]*|#\\{[\\$\\@]?[^\\}]*\\}',
								    'style'      => 'identifier',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => '\'\' string',
                                                'regex'      => '\'\'|\'(?:\\\\\\\\)+\'|\'.*?(?:[^\\\\](?:\\\\\\\\)*)\'|%q([^\\w\\[\\](){}<>])\\2|%q([^\\w\\[\\](){}<>]).*?(?:[^\\\\](?:\\\\\\\\)*)\\3|%q([^\\w\\[\\](){}<>])\\\\\\\\\\4|%q\\[\\]|%q\\[.*?([^\\\\](\\\\\\\\)*)\\]|%q\\[\\\\\\\\\\]|%q\\{\\}|%q\\{.*?([^\\\\](\\\\\\\\)*)\\}|%q\\{\\\\\\\\\\}|%q\\(\\)|%q\\(.*?([^\\\\](\\\\\\\\)*)\\)|%q\\(\\\\\\\\\\)|%q<>|%q<.*?([^\\\\](\\\\\\\\)*)>|%q<\\\\\\\\>',
					        'style'      => 'string',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'esc character',
                                                                    'regex'      => '(?:\\\\\'|\\\\\\\\)',
								    'style'      => 'esc character',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'subroutine header',
                                                'regex'      => 'def[\\t ]+\\w[\\w.]*(?:\\([^)]*\\))?',
					        'style'      => 'function header',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'arg list',
                                                                    'regex'      => '\\(.*\\)',
                                                                    'style'      => 'function header args',
                                                                    'childregex' => [
                                                                         {
                                                                         'name' => 'arg list parens',
                                                                         'regex' => '[\\(\\)]',
                                                                         'style' => 'symbol',
                                                                         'childregex' => []
                                                                         }
                                                                                    ]
                                                                  },
                                                                  {
                                                                    'name'       => 'subroutine header',
                                                                    'regex'      => '[\\t ]\w+',
								    'style'      => 'function header name',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'class header',
                                                'regex'      => 'class[\\t ]+\\w+(?:\\s*<\\s*\\w+)?',
					        'style'      => 'function header',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'class ancestor',
                                                                    'regex'      => '<\\s*\\w+',
                                                                    'style'      => 'include',
                                                                    'childregex' => [
                                                                             {
                                                                             'name' => 'inheritance doohickey',
                                                                             'regex' => '<',
                                                                             'style' => 'symbol',
                                                                             'childregex' => []
                                                                             }
                                                                                    ]
                                                                  },
                                                                  {
                                                                    'name'       => 'class main',
                                                                    'regex'      => '[\\t ]\\w+',
								    'style'      => 'type',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'regex matching 0',
                                                'regex'      => '(?:%r([^\\w\\[\\](){}<>])\\2|%r([^\\w\\[\\](){}<>]).*?(?:[^\\\\](?:\\\\\\\\)*)\\3|%r([^\\w\\[\\](){}<>])\\\\\\\\\\4|%r\\[\\]|%r\\[.*?([^\\\\](\\\\\\\\)*)\\]|%r\\[\\\\\\\\\\]|%r\\{\\}|%r\\{.*?([^\\\\](\\\\\\\\)*)\\}|%r\\{\\\\\\\\\\}|%r\\(\\)|%r\\(.*?([^\\\\](\\\\\\\\)*)\\)|%r\\(\\\\\\\\\\)|%r<>|%r<.*?([^\\\\](\\\\\\\\)*)>|%r<\\\\\\\\>)[ixpno]*',
                                                'style'      => 'regex',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'string expression',
                                                                    'regex'      => '#[\\$\\@][a-zA-Z_][\\w.]*|#\\{[\\$\\@]?[a-zA-Z_][^\\}]*\\}',
								    'style'      => 'identifier',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'regex matching I',
                                                'regex'      => '(?:\\b| )?(?:/(?:\\\\/|[^/\\n])*(?:/[ixpno]*))',
					        'style'      => 'regex',
                                                'childregex' => [
                                                                  {
                                                                    'name'       => 'string expression',
                                                                    'regex'      => '#[\\$\\@][a-zA-Z_][\\w.]*|#\\{[\\$\\@]?[a-zA-Z_][^\\}]*\\}',
								    'style'      => 'identifier',
                                                                    'childregex' => []
                                                                  }
                                                                ]
                                              },
                                              {
                                                'name'       => 'reserved words',
                                                'regex'      => '\\b(BEGIN|class|ensure|nil|self|when|END|def|false|not|super|while|alias|defined|for|or|then|yield|and|do|if|redo|true|begin|else|in|rescue|undef|break|elsif|module|retry|unless|case|end|next|return|until)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'kernel module methods',
                                                'regex',     => '\\b(Array|Float|Integer|String|at_exit|autoload|binding|caller|catch|chop|chomp|chomp!|eval|exec|exit|fail|fork|format|gets|global_variables|gsub|iterator|lambda|load|local_variables|loop|open|p|print|printf|proc|putc|puts|raise|rand|readline|readlines|require|select|sleep|split|sprintf|srand|sub|syscall|system|test|trace_var|trap|untrace_var)\\b',
					        'style'      => 'library function',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'braces, parens and brakets',
                                                'regex'      => '[\\[\\]\\{\\}\\(\\)]',
					        'style'      => 'braces',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => '<< stuff',
                                                'regex'      => '<<(?:("|\')([^\\n]*)\\2|\\w*).*?^\\3$',
					        'style'      => 'text',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'symbols',
                                                'regex'      => '(?:[:*-+<>=^!,/]+|\.\.+)',
                                                'style'      => 'symbol',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'numbers',
                                                'regex'      => '\d[\d.]*',
                                                'style'      => 'numeric',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'embedded documentation',
                                                'regex'      => '^=.*?^(?:=end|\\Z)',
					        'style'      => 'doc comment',
                                                'childregex' => []
                                              }
                                            ]
                          };
 
# taken from nedit
# modified by PP
# very inclomplete!
$LANGUAGE{'sql'}        = {
                            'filename'   => '(?i)\\.sql$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'keywords I',
                                                'regex'      => '(?i)(,|%|<|>|:=|=|\\(|\\)|\\bselect|on|from|order by|desc|where|and|or|not|null|true|false)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'comment I',
                                                'regex'      => '--.*?$',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'comment II',
                                                'regex'      => '/\\*.*?\\*/',
					        'style'      => 'comment',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'string',
                                                'regex'      => '\'\'|\'.*?([^\\\\](\\\\\\\\)*)\'|\'\\\\\\\\\'',
#                                                'regex'      => '(\'\'|\'[^\'\\\\]\'|\'[^\'].*?[^\\\\]\')',
					        'style'      => 'string',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keywords II',
                                                'regex'      => '(?i)end if;|\\b(create|replace|begin|end|function|return|fetch|open|close|into|is|in|when|others|grant|on|to|exception|show|set|out|pragma|as|package)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'keywords III',
                                                'regex'      => '(?i)\\balter\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'datatypes',
                                                'regex'      => '(?i)\\b(integer|blol|date|numeric|character|varying|varchar|char)\\b',
					        'style'      => 'predefined type',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'words',
                                                'regex'      => '(?i)\\b(constraint|key|references|primary|table|foreign|add|insert|group by)\\b',
					        'style'      => 'reserved word',
                                                'childregex' => []
                                              }
                                            ]
                            };
   
 

$LANGUAGE{'patch'}        = {
                            'filename'   => '(?i)\\.patch$|\\.diff$',
                            'regex'      => '',
                            'patterns'   => [
                                              {
                                                'name'       => 'header',
                                                'regex'      => '^Index: .*?$|^===== .*?$|^diff .*?$|^--- .*?$|^\+\+\+ .*?$',
                                                'style'      => 'separator',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'hunk',
                                                'regex'      => '^@@ .*?$',
                                                'style'      => 'line spec',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'from',
                                                'regex'      => '^-.*?$',
                                                'style'      => 'deletion',
                                                'childregex' => []
                                              },
                                              {
                                                'name'       => 'to',
                                                'regex'      => '^\+.*?$',
                                                'style'      => 'insertion',
                                                'childregex' => []
                                              }
                                            ]
                            };



#####
#
# LANGUAGE: shell script
#

$LANGUAGE{'shellscript'} = {
	'filename' => '\\.(sh|shell)$',
	'regex' => '^\\s*#\\s*![^\\s]*(sh|bash|ash|zsh|ksh)',
	'patterns' => [ {
		'name' => 'comment',
#		'regex' => '^[ \t]*[^$]?\#[^!]?.*?$',
		'regex' => '(^| )#([^\\!].)*?$',
		'style' => 'comment',
		'childregex' => []
	}, {
		'name' => 'identifier',
		'regex' => '[a-zA-Z][a-zA-Z0-9_]*=',
		'style' => '',
		'childregex' => [ {
			'name' => 'identifier',
			'regex' => '[a-zA-Z][a-zA-Z0-9_]*',
			'style' => 'identifier',
			'childregex' => []
		} ]
	}, {
		'name' => 'identifier',
		'regex' => '\\$([0-9#\\*]|[a-zA-Z][a-zA-Z0-9_]*)',
		'style' => 'identifier',
		'childregex' => []
	}, {
		'name' => 'interpreter line',
		'regex' => '^[ \t]*#!.*?$',
		'style' => 'preprocessor',
		childregex => []
	}, {
		'name' => 'string',
		'regex' => '""|"(\\\\"|[^\\"])*"',
		'style' => 'string',
		childregex => [ {
			'name' => 'identifier',
			'regex' => '\\$([0-9#\\*]|[a-zA-Z][a-zA-Z0-9_]*)',
			'style' => 'identifier',
			'childregex' => []
		} ]
	} ]
};

$LANGUAGE{'sh'} = $LANGUAGE{'shellscript'};
return \%LANGUAGE;

};
