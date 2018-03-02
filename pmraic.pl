#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  pmraic.pl
#
#        USAGE:  ./pmraic.pl --man
#                ./pmraic.pl data.dat 
#                ./pmraic.pl --bash=12 data.dat 
#                ./pmraic.pl --maui data.dat
#                ./pmraic.pl --summarize data.dat
#                ./pmraic.pl --models=JC69,JC69I,JC69G,JC69IG data.dat
#                ./pmraic.pl --modeltest data.dat
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  1) PHYML v 3. Preferrably a version downloaded from 
#                   http://code.google.com/p/phyml/. See note BUGS below.
#                
#                2) Perl modules Getopt::Long and Parallel::ForkManager. 
#                   To install:
#                      sudo perl -MCPAN -e 'install Parallel::ForkManager'
#                   or on Ubuntu
#                      sudo apt-get install libparallel-forkmanager-perl
#
#         BUGS:  A bug in phyml "v.3" for Ubuntu Linux caused alpha to be estimated
#                together with I, even if not specifically asked for. This bug
#                was fixed in later phyml versions (at least in v. 20110919)
#        NOTES:  ---
#       AUTHOR:  Johan A. A. Nylander (JN), <jnylander @ users.sourceforge.net>
#      COMPANY:  SU
#      VERSION:  1.02
#      CREATED:  11/27/2009 09:35:52 AM CET
#     REVISION:  11/16/2011 02:48:35 PM
#         TODO:  Handle --bash=6 -cpu=4 to run forked runs (6X4)
#                Check df's!
#
#===============================================================================

use strict;
use warnings;

use Getopt::Long;
use Parallel::ForkManager;
use File::Copy;
use Data::Dumper;
use Benchmark;
use Pod::Usage;


## Global variables
GLOBALS:

my $VERSION            = '1.1';
my $PROGRAM_NAME       = 'pmraic.pl';
my $CHANGES            = '11/10/2011 07:31:00 PM';
my $AUTHOR             = 'Johan Nylander';

my $infile             = q{};
my $verbose            = 1; # verbose or not
my $debug              = 0; # debug printing or not
my $summarize          = 0; # Summarize results from output unless '--nosummarize'
my $nosummarize        = 0; # Hack

my $outfile            = q{};
my $output_file_prefix = '__pMrAIC.';
my $bash_script_prefix = '__pMrAIC_script.';
my $maui_script_suffix = '.maui';
my $output_file_suffix = '.pMrAIC.txt';
my $archive_suffix     = '.pMrAIC.tgz';

my $n_models           = q{};
my @models             = ();
my @mrmodeltest_models = qw(F81 F81G F81I F81IG GTR GTRG GTRI GTRIG HKY85 HKY85G HKY85I HKY85IG JC69 JC69G JC69I JC69IG K2P K2PG K2PI K2PIG SYM SYMG SYMI SYMIG);
my $models             = q{};
my $modeltest          = 0; # Use the 24 MrBayes models by default

my $bash               = 0;
my $maui               = 0;
#my $sequential         = q{};
my %HoH                = (); 
my %run_HoH            = (); 
                            #$run_HoH{ $model_name } = {
                            #    name        => $name,
                            #    lnL         => $lnL,
                            #    df          => $df,
                            #    AIC         => $AIC,
                            #    AICc        => $AICc,
                            #    BIC         => $BIC,
                            #    mbcmd       => $mbcmd,
                            #    delta_aic   => $delta_aic,
                            #    delta_aicc  => $delta_aicc,
                            #    delta_bic   => $delta_bic,
                            #    aic_weight  => $aic_weight,
                            #    aicc_weight => $aicc_weight,
                            #    bic_weight  => $bic_weight,
                            # };
my @dimensions         = ();
my $max_num_params     = q{};
my $use_brlens         = 1;  # use branch lengths in df calculations or not
my $BAratio            = 40; # Burnham Anderson ratio nchar/nparams


## Parallel defaults
my $MAX_PROCESSES = check_processors();
my $CPU_DEFAULT   = 4;   # temporary
my $cpu           = q{}; # temporary

## Maui defaults
my $walltime       = '48:00:00'; # temporary
my $phyml_maui_bin = '/n/app/phyml'; # On Mother.bergianska.se

## Shell script defaults
my $phyml_bash_bin = q{};

## PHYML defaults
my $phyml              = q{};
my $phyml_bin          = 'phyml';
my $search             = 'NNI';
my $optimize           = 'tlr';
my $datatype           = 'nt';
my $phyml_stats_ending = '_phyml_stats.txt';  # check mac/win/linux names!
my $phyml_tree_ending  = '_phyml_tree.txt'; # check mac/win/linux names!


#-------------------------------------------------------------------------------
#  Functions
#-------------------------------------------------------------------------------
#===  FUNCTION  ================================================================
#         NAME:  calculate_aic
#      VERSION:  08/30/2011 09:34:00 PM CEST
#  DESCRIPTION:  Calculate AIC
#                AIC = -2lnL + 2K
#                L: max likelihood, K: number of free parameters
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub calculate_aic {

    my ($model) = @_;

    my $df = ($use_brlens) ? 'df' : 'dfnb'; # which df to use, with or without brlens

    my $aic = (((-2.0)*($run_HoH{$model}->{'lnl'})) + (2.0*($run_HoH{$model}->{$df})));

    return($aic);

} # end of calculate_aic


#===  FUNCTION  ================================================================
#         NAME:  calculate_aicc
#      VERSION:  08/30/2011 09:39:48 PM CEST
#  DESCRIPTION:  Calculate AICc
#                AICc = -2lnL + 2K + 2K(K+1)/n-K-1
#                L: max likelihood, K: number of free parameters, n: sample size
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub calculate_aicc {

    my ($model) = @_;

    my $nchar = $dimensions[1];
    my $df    = ($use_brlens) ? 'df' : 'dfnb'; # which df to use, with or without brlens

    my $aicc = (((-2.0)*($run_HoH{$model}->{'lnl'})) + (2.0*($run_HoH{$model}->{$df})) + ((2.0*($run_HoH{$model}->{$df})*(($run_HoH{$model}->{$df})+1.0))/($nchar-($run_HoH{$model}->{$df})-1.0)));

    return($aicc);

} # end of calculate_aicc


#===  FUNCTION  ================================================================
#         NAME:  calculate_bic
#      VERSION:  08/30/2011 09:39:41 PM CEST
#  DESCRIPTION:  Calculate BIC
#                BIC = -2lnL + Kln(n)
#                L: max likelihood, K: number of free parameters, n: sample size
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub calculate_bic {

    my ($model) = @_;

    my $nchar = $dimensions[1];
    my $df    = ($use_brlens) ? 'df' : 'dfnb'; # which df to use, with or without brlens

    my $bic = (((-2.0)*($run_HoH{$model}->{'lnl'})) + ((log($nchar)) * ($run_HoH{$model}->{$df})));

    return($bic);

} # end of calculate_bic


#===  FUNCTION  ================================================================
#         NAME:  calculate_aic_weights
#      VERSION:  08/31/2011 03:07:43 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub calculate_aic_weights {

    my $sumAICExp = 0.0;
    my $AIC_min_model = get_min_aic_model();
    my $AIC_min = $run_HoH{$AIC_min_model}->{'aic'};

    foreach my $model (keys %run_HoH) {
        $run_HoH{$model}->{'delta_aic'} = $run_HoH{$model}->{'aic'} - $AIC_min;
        $sumAICExp += exp((-0.5) * ($run_HoH{$model}->{'delta_aic'}));
    }
    foreach my $model (keys %run_HoH) {
        $run_HoH{$model}->{'aic_weight'} = exp((-0.5) * ($run_HoH{$model}->{'delta_aic'})) / $sumAICExp ;
    }

} # end of calculate_aic_weights



#===  FUNCTION  ================================================================
#         NAME:  calculate_aicc_weights
#      VERSION:  08/31/2011 03:11:58 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub calculate_aicc_weights {

    my $sumAICcExp = 0.0;
    my $AICc_min_model = get_min_aicc_model();
    my $AICc_min = $run_HoH{$AICc_min_model}->{'aicc'};

    foreach my $model (keys %run_HoH) {
        $run_HoH{$model}->{'delta_aicc'} = $run_HoH{$model}->{'aicc'} - $AICc_min;
        $sumAICcExp += exp((-0.5) * ($run_HoH{$model}->{'delta_aicc'}));
    }

    foreach my $model (keys %run_HoH) {
        $run_HoH{$model}->{'aicc_weight'} = exp((-0.5) * ($run_HoH{$model}->{'delta_aicc'})) / $sumAICcExp ;
    }

} # end of calculate_aicc_weights


#===  FUNCTION  ================================================================
#         NAME:  calculate_bic_weights
#      VERSION:  08/31/2011 03:12:06 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub calculate_bic_weights {

    my $sumBICExp = 0.0;
    my $BIC_min_model = get_min_bic_model();
    my $BIC_min = $run_HoH{$BIC_min_model}->{'bic'};

    foreach my $model (keys %run_HoH) {
        $run_HoH{$model}->{'delta_bic'} = $run_HoH{$model}->{'bic'} - $BIC_min;
        $sumBICExp += exp((-0.5) * ($run_HoH{$model}->{'delta_bic'}));
    }
    foreach my $model (keys %run_HoH) {
        $run_HoH{$model}->{'bic_weight'} = exp((-0.5) * ($run_HoH{$model}->{'delta_bic'})) / $sumBICExp ;
    }

} # end of calculate_bic_weights


#===  FUNCTION  ================================================================
#         NAME:  check_models_arg
#      VERSION:  12/08/2009 09:27:05 AM CET
#  DESCRIPTION:  Checks the models given on command line against the keys in 
#                the all_models_commands_hash.
#                Lookup code from Perl Cookbook.
#   PARAMETERS:  An array (from  --models=JC69,JC69I,GTR)
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub check_models_arg {

    my (@mods) = @_;

    @mods = split(/,/, join(',', @mods)); # From perldoc Getopt::Long
    my @keys = keys %HoH;

    my %seen;      # lookup table
    my @mods_only; # answer

    @seen{@keys} = (); # build lookup table

    foreach my $mod (@mods) {
        $mod = uc($mod);
        push(@mods_only, $mod) unless exists $seen{$mod};
    }

    if (scalar(@mods_only) > 0) {
        print STDERR  "Wrong model syntax: @mods_only\n";
        print STDERR  "Use \'$0 --help\' to see the valid choices.\n";
        exit(0);
    }
    else {
        return(@mods);
    }

} # end of check_models_arg

#
#===  FUNCTION  ================================================================
#         NAME:  check_models
#      VERSION:  12/08/2009 09:27:05 AM CET
#  DESCRIPTION:  Checks the models against the keys in 
#                the HoH hash 
#                Lookup code from Perl Cookbook.
#   PARAMETERS:  An array
#      RETURNS:  array of correct models
#         TODO:  ???
#===============================================================================
sub check_models {

    my (@mods) = @_;

    my @keys = keys %HoH;

    my %seen;      # lookup table
    my @mods_only; # answer

    @seen{@keys} = (); # build lookup table

    foreach my $mod (@mods) {
        $mod = uc($mod);
        push(@mods_only, $mod) unless exists $seen{$mod};
    }

    if (scalar(@mods_only) > 0) {
        print STDERR  "Wrong model syntax: @mods_only\n";
        print STDERR  "Use \'$0 --help\' to see the valid choices.\n";
        exit(0);
    }
    else {
        return(@mods);
    }

} # end of check_models

#===  FUNCTION  ================================================================
#         NAME:  check_optimize_arg
#      VERSION:  11/28/2009 09:22:23 PM CET
#  DESCRIPTION: Check optimize argument given to pMrAIC
#               from PHYML manual:
#               -o params
#                  This option focuses on specific parameter optimisation.
#                  params=tlr : tree topology (t), branch length (l) and rate parameters (r) are optimised.
#                  params=tl  : tree topology and branch length are optimised.
#                  params=lr  : branch length and rate parameters are optimised.
#                  params=l   : branch length are optimised.
#                  params=r   : rate parameters are optimised.
#                  params=n   : no parameter is optimised.
#
#   PARAMETERS:  $optimize_arg, string
#      RETURNS:  $optimize_arg, string
#         TODO:  ???
#===============================================================================
sub check_optimize_arg {

    my ($optimize_arg) = @_;

    my $opt_arg = lc($optimize_arg);

    if ($opt_arg eq 'tlr') {
        return($opt_arg);
    }
    elsif ($opt_arg eq 'tl') {
        return($opt_arg);
    }
    elsif ($opt_arg eq 'lr') {
        return($opt_arg);
    }
    elsif ($opt_arg eq 'l') {
        return($opt_arg);
    }
    elsif ($opt_arg eq 'r') {
        return($opt_arg);
    }
    elsif ($opt_arg eq 'n') {
        return($opt_arg);
    }
    else {
        die "Optimize argument \"$optimize_arg\" is not recognized.\n Valid options are tlr, tl, lr, l, r, or n.\n";
    }

} # end of check_optimize_arg


#===  FUNCTION  ================================================================
#         NAME:  check_processors
#      VERSION:  11/09/2011 06:18:14 PM
#  DESCRIPTION:  check number of CPUs/cores by system call
#   PARAMETERS:  Uses global $CPU_DEFAULT
#      RETURNS:  numerical (ncores or ncpus. Make sure to know which!)
#         TODO:  Make sure works on windows
#===============================================================================
sub check_processors {

    my $ncores = 0;
    my $ncpus  = 0;

    if ($^O =~ /linux/i) { # Linux: do regexp on /proc/cpuinfo
        my (@cores) = qx(grep \'^cpu cores\' /proc/cpuinfo);
        $ncpus = qx(grep -c \'^cpu cores\' /proc/cpuinfo);
        foreach my $line (@cores) {
            $line =~ /:\s*(\d+)/;
            $ncores += $1;
        }
    }
    elsif ($^O =~ /Darwin/i) {
        $ncpus = qx('/usr/bin/sysctl -n hw.ncpu'); # hw.availcpu
        $ncores = $ncpus; # temporary
    }
    else {
        $ncpus = $CPU_DEFAULT;
        $ncores = $CPU_DEFAULT;
        print STDERR "Could not search for number of CPUs/cores. Setting default cpu=$CPU_DEFAULT\n";
    }

    return $ncores; # or $ncpus

} # end of check_processors


#===  FUNCTION  ================================================================
#         NAME:  check_search_arg
#      VERSION:  11/28/2009 09:22:23 PM CET
#  DESCRIPTION:  Check search argument for phyml. From the phyml manual:
#                   -s (or --search) move
#                    Tree topology search operation option.
#                    Can be either NNI (default, fast) or SPR (a bit slower than NNI)
#                    or BEST (best of NNI and SPR search).
#   PARAMETERS:  $search_arg, string
#      RETURNS:  $search_arg
#         TODO:  ???
#===============================================================================
sub check_search_arg {

    my ($search_arg) = @_;

    my $s_arg = uc($search_arg);

    if ($s_arg eq 'NNI') {
        return($s_arg);
    }
    elsif ($s_arg eq 'SPR') {
        return($s_arg);
    }
    elsif ($s_arg eq 'BEST') {
        return($s_arg);
    }
    else {
        die "Optimize argument \"$search_arg\" is not recognized.\n Valid options are NNI, SPR, or BEST.\n";
    }

} # end of check_search_arg


#===  FUNCTION  ================================================================
#         NAME: clean_up
#      VERSION: 08/31/2011
#      PURPOSE: collect and compress output files
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#         TODO: ????
#===============================================================================
sub clean_up {

    my $glob_expression = $output_file_prefix . "*_phyml_*.txt";
    my (@model_files) = glob($glob_expression);

    ## Compress
    my $archive_name = $infile . $archive_suffix;
    my $aref = \@model_files;
    compress($aref, $archive_name);
    
    ## Remove
    $glob_expression = $output_file_prefix . "*"; # including temporary dat files
    my (@all_model_files) = glob($glob_expression);

    unlink(@all_model_files);

} # end of clean_up



#===  FUNCTION  ================================================================
#         NAME: compress
#      VERSION: 09/01/2011 07:44:11 PM
#      PURPOSE: compress files
#   PARAMETERS: array_ref
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#         TODO: THIS VERSION ONLY COMPRESSES THE FIRST FILE GIVEN A ARRAY REF! 08/31/2011 09:19:42 PM
#===============================================================================
sub compress {

    my ($array_ref, $output) = @_;

    use Archive::Tar;

    my $tar = Archive::Tar->new;
    $tar->add_files(@{$array_ref});
    $tar->write($output, COMPRESS_GZIP);

} # end of compress


#===  FUNCTION  ================================================================
#         NAME:  create_bash_scripts
#      VERSION:  08/31/2011 12:33:19 PM
#  DESCRIPTION:  Creates a number of bash scripts for running phyml
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub create_bash_scripts {

    my ($n_bash_scripts) = @_;

    ## Test if number is greater than the number of models.
    if ($n_bash_scripts > $n_models) {
        $n_bash_scripts = $n_models;
        print STDERR "Warning: number of requested bash scripts ($bash) exceeds number of models ($n_models).\nSetting bash=$n_models.\n";
    }

    ## Randomizing models to spread the work load per script
    my @keys = sort keys %run_HoH;
    my @rkeys = randomize_array(@keys);
    
    ## Set how many models to print in each script
    my $ntot = @rkeys;
    my $mods_per_bash = round_up($ntot / $n_bash_scripts); 

    ## Is file sequential?
    my $sequential = is_sequential($infile);

    ## Read the data (using the localized ARGV and file handle method, i.e., no 'open'!)
    my $infile_content = do { local( @ARGV, $/ ) = $infile ; <> } ;

    ## Print the script files
    for (my $i=1 ; $i <= $n_bash_scripts ; $i++) {
        ## Open script file
        my $bash_script_file = $bash_script_prefix . $i . ".sh";
        open my $SCRIPTFILE, ">", $bash_script_file or die "could not open script file for writing : $! \n";
        
        ## Print to script file
        print $SCRIPTFILE "#!/bin/sh\n";
        my $datestring = prettydate();
        print $SCRIPTFILE "# Shell script for running PHYML on data \'$infile\'\n";
        print $SCRIPTFILE "# Created $datestring by $PROGRAM_NAME v.$VERSION:\n\n";

        ## Print phyml binary check routine to script file
        print $SCRIPTFILE "PHYML=$phyml_bash_bin #<<<<<<< specify here the path to phyml if different from default!\n";
        print $SCRIPTFILE "if [ ! -x \"\$PHYML\" ]; then\n";
        print $SCRIPTFILE "   echo \"can not find executable phyml in the PATH. Quitting\"\n";
        print $SCRIPTFILE "   exit 1\n";
        print $SCRIPTFILE "fi\n";

        ## Print the data in the beginning of the file
        print $SCRIPTFILE "\n# Data from file \'$infile\'\n";
        print $SCRIPTFILE "DATA=\"\n";
        print $SCRIPTFILE $infile_content, "\"\n";


        ## Print models to script file
        my @mods = splice(@rkeys, 0, $mods_per_bash);
        if (0 == scalar(@mods)) {
            print STDERR "Number of requested scripts doesn't scale well with number of models. Edit output manually as desired.\n";
        }

        for my $mod (@mods) {
            my $mod_infile = $output_file_prefix . $mod . '.' . $infile;
            my $mod_treefile = $mod_infile . '_phyml_tree.txt';
            print $SCRIPTFILE "\n# Model: $mod\n";
            print $SCRIPTFILE 'echo "$DATA"', " > $mod_infile\n";
            print $SCRIPTFILE "\$PHYML --input $mod_infile --datatype $datatype $sequential --model $run_HoH{$mod}->{'phymlcmd'} --search $search -o $optimize\n";
            print $SCRIPTFILE "if [ -s \"$mod_treefile\" ] ; then\n  rm -v $mod_infile\nfi\n";
        }
        close($SCRIPTFILE) or warn "could not close script file : $! \n";
    }

} # end of create_bash_scripts


#===  FUNCTION  ================================================================
#         NAME:  create_maui_files
#      VERSION:  08/31/2011 01:08:45 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub create_maui_files {

    my $hash_ref = shift;

    my @maui_submit_files = ();

    my @maui_models = keys %{$hash_ref};

    my $datestring = prettydate();

    ## Is file sequential
    my $sequential = is_sequential($infile);

    foreach my $model (@maui_models) {
    
        ## Copy data files
        my $model_infile = $output_file_prefix . $model . '.' . $infile;
        copy($infile, $model_infile) or die "File cannot be copied.";

        ## Create maui submit files
        my $model_maui_file = $output_file_prefix . $model . "." . $infile . $maui_script_suffix;
        push @maui_submit_files, $model_maui_file;

        open my $MF , ">", $model_maui_file or die "could not open maui file : $! \n";

        print $MF "## File for submitting PHYML runs on the cluster\n";
        print $MF "## $datestring\n";
        print $MF "## Original infile: $infile\n";
        print $MF "## Edit the file where appropriate and then submit the file\n";
        print $MF "## to the queue by using the command\n";
        print $MF "## \'qsub $model_maui_file\'\n";
        print $MF "##\n";
        print $MF "#PBS -l nodes=1:ppn=1\n";
        print $MF "#PBS -l walltime=$walltime\n";
        print $MF "##\n";
        print $MF "cd \$PBS_O_WORKDIR\n";
        print $MF "echo Running on host \`hostname\`\n";
        print $MF "echo Time is \`date\`\n";
        print $MF "echo Directory is \`pwd\`\n";
        print $MF "echo This jobs runs on the following processors:\n";
        print $MF "echo \`cat \$PBS_NODEFILE\`\n";
        print $MF "echo\n";
        print $MF "## Run phyml\n";
        print $MF "$phyml_maui_bin --input $model_infile --datatype $datatype $sequential --model $hash_ref->{$model}->{'phymlcmd'} --search $search -o $optimize\n";
        print $MF "echo\n";
        print $MF "echo Stop time is \`date\`\n";
        print $MF "echo\n";
        print $MF "\n";

        close($MF) or warn "could not close maui file: $! \n";
    }


    ## Open bash script to submit the maui files
    my $maui_submit_shell_file = $output_file_prefix . 'RUN_ME_for_MAUI.' . $infile . '.sh'; 
    open my $BF, ">", $maui_submit_shell_file or die "could not open maui submit shell script : $! \n";
    print $BF "#!/bin/sh\n";
    foreach my $mf (@maui_submit_files) {
        print $BF "qsub $mf\n";
    }
    close($BF) or warn "could not close maui submit shell script : $! \n";

    ## Change mode
    my $mode = 0755;
    chmod $mode, $maui_submit_shell_file;

} # end of create_maui_files


#===  FUNCTION  ================================================================
#         NAME:  find_prog
#      VERSION:  09/01/2011 07:31:02 PM
#  DESCRIPTION:  Checks if prog (or any executable file passed as its name in a
#                string to the function) can be found.
#   PARAMETERS:  The search pattern (name of the binary)
#        USAGE:  my $prog = find_prog('prog');
#      RETURNS:  Path to prog or '0' if no success
#         TODO:  Make sure works on win
#===============================================================================
sub find_prog {

    my $PROG = shift;
    my $prog = q{};

    if ( -x $PROG ) { # if '/usr/bin/phyml'
        $prog = $PROG;
    }
    else { # if 'phyml'
        FIND_PROG:
        foreach ( split( /:/, $ENV{PATH} ) ) {
            if (-x "$_/$PROG") {
                $prog = "$_/$PROG";
                last FIND_PROG;
            }
        }
    }
    if ($prog eq '') {
        die "Can not find executable file \'$PROG\'.\nPlease specify the name of phyml v.3 in the script or using \'--phyml\'.\n";
    }
    else {
        return $prog;
    }

} # end of find_prog


#===  FUNCTION  ================================================================
#         NAME:  get_best_tree
#      VERSION:  11/28/2009 09:22:23 PM CET
#  DESCRIPTION:  get best tree from phyml tree file
#   PARAMETERS:  tree file name
#      RETURNS:  tree string
#         TODO:  ???
#===============================================================================
sub get_best_tree {

    my ($tree_file) = @_;

    my $tree = q{};

    open (my $TREEFILE, '<', $tree_file) or die "\nIn get_best_tree: Cannot open $tree_file for reading: $!\n";
    while (<$TREEFILE>) {
        chomp($tree = $_);
    }
    close($TREEFILE) or warn "could not close file handle : $!\n";
    if ($tree eq "") {
        die "\a\nERROR!\n\nCould not read the tree from PHYML.\n\n";
    }
    else {
        return($tree);
    }

} # end of get_best_tree


#===  FUNCTION  ================================================================
#         NAME:  get_dimensions
#      VERSION:  12/08/2009 01:33:55 PM CET
#  DESCRIPTION:  Get the number of taxa and characters from the infile
#   PARAMETERS:  filename
#      RETURNS:  ntax, nchar, nbranch, max_num_params
#         TODO:  ???
#===============================================================================
sub get_dimensions {

    my ($file) = shift(@_);

    my @dimensions = ();

    open my $INFILE, "<", $file or die "\a\nERROR!\n\nCannot open $file for reading: $!";
    chomp (my $firstline = <$INFILE>);
    close($INFILE);

    $_ = $firstline;

    if (/\s*(\d+)\s+(\d+)/) {
        my $ntax = $1;
        push @dimensions, $ntax;
        my $nchar = $2;
        push @dimensions, $nchar;
        my $nbranch = ((2*$ntax)-3);
        push @dimensions, $nbranch;
        my $maxNumParameters = $nbranch + 10;    # Number of branches plus parameters in GTR+I+G
        push @dimensions, $maxNumParameters;
    }
    else {
        die "\a\nERROR!\n\nCould not read number of taxa and/or characters\n\n"
    }

    return(@dimensions); 

} # end of get_dimensions


#===  FUNCTION  ================================================================
#         NAME:  get_likelihood
#      VERSION:  08/31/2011 08:17:25 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  $phyml_stat_file
#      RETURNS:  likelihood value (string)
#         TODO:  ???
#===============================================================================
sub get_likelihood {

    my ($stat_file) = @_;

    my $lnl = 999.999;

    open (my $STATFILE, '<', $stat_file) or die "\nnCannot open $stat_file for reading: $!";
    while (<$STATFILE>) {
        chomp;
        #$lnl = $1 if ($_ =~ m/likelihood\s*:[^\-]+(\-.+)/i); # A.W.
        $lnl = $1 if ($_ =~ m/Log-likelihood\s*:\s*(\-\d+\.\d+)/i); # JN
    }
    close($STATFILE) or warn "could not close file handle: $!\n";

    if ($lnl == 999.999) {
        die "\a\nERROR!\n\nCould not read the likelihood from PHYML.\n\n";
    }
    else {
        return($lnl);
    }

} # end of get_likelihood


#===  FUNCTION  ================================================================
#         NAME:  get_min_aic_model
#      VERSION:  08/31/2011 06:41:55 PM
#  DESCRIPTION:  Get the minAIC model name
#   PARAMETERS:  ???
#      RETURNS:  string, 'F81', name of AIC best model
#         TODO:  ???
#===============================================================================
sub get_min_aic_model {

    my (@array) = sort { $run_HoH{$a}->{'aic'} cmp $run_HoH{$b}->{'aic'} } keys %run_HoH;

    my $AIC_min_model = shift(@array);

    return($AIC_min_model);

} # end of get_min_aic_model


#===  FUNCTION  ================================================================
#         NAME:  get_min_aicc_model
#      VERSION:  08/31/2011 06:42:01 PM
#  DESCRIPTION:  Get the minAICc model name
#   PARAMETERS:  ???
#      RETURNS:  string, 'F81', name of AICc best model 
#         TODO:  ???
#===============================================================================
sub get_min_aicc_model {

    my (@array) = sort { $run_HoH{$a}->{'aicc'} cmp $run_HoH{$b}->{'aicc'} } keys %run_HoH;

    my $AICc_min_model = shift(@array);

    return($AICc_min_model);

} # end of get_min_aicc_model


#===  FUNCTION  ================================================================
#         NAME:  get_min_bic_model
#      VERSION:  08/31/2011 06:42:10 PM
#  DESCRIPTION:  Get the minBIC, model name, tree, and mbblock
#   PARAMETERS:  ???
#      RETURNS:  string, 'F81', name of BIC best model
#         TODO:  ???
#===============================================================================
sub get_min_bic_model {

    my (@array) = sort { $run_HoH{$a}->{'bic'} cmp $run_HoH{$b}->{'bic'} } keys %run_HoH;

    my $BIC_min_model = shift(@array);

    return($BIC_min_model);

} # end of get_min_bic_model


#===  FUNCTION  ================================================================
#         NAME:  initialize_HoH
#      VERSION:  11/10/2011 02:52:07 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  null. Populates a HoH
#         TODO:  Make sure the new phyml syntax is correct for modeltest models.
#                Especially for, e.g., using 'TN93' vs. 'custom'
#===============================================================================
sub initialize_HoH {

    ## MB commands
    my $equal      = "rates=equal";
    my $propinv    = "rates=propinv";
    my $gamma      = "rates=gamma";
    my $invgamma   = "rates=invgamma";
    my $Prset      = "Prset applyto=(1)";
    my $Lset       = "Lset applyto=(1)";
    my $nst1       = "nst=1";
    my $nst2       = "nst=2";
    my $nst6       = "nst=6";
    my $shapepr    = "shapepr=Uniform(0.1,50.0)";
    my $pinvarpr   = "pinvarpr=Uniform(0.0,1.0)";
    my $eqfreqpr   = "statefreqpr=Fixed(Equal)";
    my $uneqfreqpr = "statefreqpr=Dirichlet(1.0,1.0,1.0,1.0)";
    my $revmatpr   = "revmatpr=Dirichlet(1.0,1.0,1.0,1.0,1.0,1.0)";
    my $tratiopr   = "tratiopr=Beta(1.0,1.0)";

    ## PHYML commands
    my $pinv_e     = '--pinv e';
    my $nclasses_1 = '--nclasses 1';
    my $nclasses_4 = '--nclasses 4';
    my $alpha_e    = '--alpha e';
    my $f_m        = '-f m';
    my $f_equal    = "-f \"0.25,0.25,0.25,0.25\"";
    my $ts_tv_e    = '--ts/tv e';
    my $custom_sym = '012345';
    my $custom_k3p = '012210';
    my $custom_trn = '010020';
    my $custom_tim = '012230';
    my $custom_tvm = '012314';


    ## For df. $dimensions[2] is $nbranch
    my $nbranch = $dimensions[2];

    ## Populate the hash
    $HoH{'JC69'}   = {
        name  => 'JC69',
        mbcmd => " $Lset $nst1 $equal;\n $Prset $eqfreqpr;",
        phymlcmd => " JC69 $nclasses_1 ",
        df => 0 + $nbranch,
        dfnb => 0,
    };
    $HoH{'JC69I'}  = {
        name  => 'JC69I',
        mbcmd => " $Lset $nst1 $propinv;\n $Prset $eqfreqpr $pinvarpr",
        phymlcmd => " JC69 $nclasses_1 $pinv_e ",
        df => 1 + $nbranch,
        dfnb => 1,
    };
    $HoH{'JC69G'}  = {
        name => 'JC69G',
        mbcmd => " $Lset $nst1 $gamma;\n $Prset $eqfreqpr $shapepr;",
        phymlcmd => " JC69 $nclasses_4 $alpha_e ",
        df => 1 + $nbranch,
        dfnb => 1,
    };
    $HoH{'JC69IG'} = {
        name => 'JC69IG',
        mbcmd => " $Lset $nst1 $invgamma;\n $Prset $eqfreqpr $shapepr $pinvarpr;",
        phymlcmd => " JC69  $pinv_e $nclasses_4 $alpha_e ",
        df => 2 + $nbranch,
        dfnb => 2,
    };
    $HoH{'F81'}    = {
        name => 'F81',
        mbcmd => " $Lset $nst1 $equal;\n $Prset $uneqfreqpr;",
        phymlcmd => " F81 $f_m $nclasses_1 ",
        df => 3 + $nbranch,
        dfnb => 3,
    };
    $HoH{'F81I'}   = {
        name => 'F81I',
        mbcmd => " $Lset $nst1 $propinv;\n $Prset $uneqfreqpr $pinvarpr;",
        phymlcmd => " F81 $f_m $nclasses_1 $pinv_e ",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'F81G'}   = {
        name => 'F81G',
        mbcmd => " $Lset $nst1 $gamma;\n $Prset $uneqfreqpr $shapepr;",
        phymlcmd => " F81 $f_m $nclasses_4 $alpha_e ",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'F81IG'}  = {
        name => 'F81IG',
        mbcmd => " $Lset $nst1 $invgamma;\n $Prset $uneqfreqpr $shapepr $pinvarpr;",
        phymlcmd => " F81 $f_m $pinv_e $nclasses_4 $alpha_e ",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'K2P'}    = {
        name => 'K2P',
        mbcmd => " $Lset $nst2 $equal;\n $Prset $eqfreqpr $tratiopr;",
        phymlcmd => " K80 $ts_tv_e $nclasses_1 ",
        df => 1 + $nbranch,
        dfnb => 1,
    };
    $HoH{'K2PI'}   = {
        name => 'K2PI',
        mbcmd => " $Lset $nst2 $propinv;\n $Prset $eqfreqpr $tratiopr $pinvarpr;",
        phymlcmd => " K80 $ts_tv_e $nclasses_1 $pinv_e ",
        df => 2 + $nbranch,
        dfnb => 2,
    };
    $HoH{'K2PG'}   = {
        name => 'K2PG',
        mbcmd => " $Lset $nst2 $gamma;\n $Prset $eqfreqpr $tratiopr $shapepr;",
        phymlcmd => " K80 $ts_tv_e $nclasses_4 $alpha_e ",
        df => 2 + $nbranch,
        dfnb => 2,
    };
    $HoH{'K2PIG'}  = {
        name => 'K2PIG',
        mbcmd => " $Lset $nst2 $invgamma;\n $Prset $eqfreqpr $tratiopr $shapepr $pinvarpr;",
        phymlcmd => " K80 $ts_tv_e $pinv_e $nclasses_4 $alpha_e ",
        df => 3 + $nbranch,
        dfnb => 3,
    };
    $HoH{'HKY85'}    = {
        name => 'HKY85',
        mbcmd => " $Lset $nst2 $equal;\n $Prset $uneqfreqpr $tratiopr;",
        phymlcmd => " HKY85 $f_m $ts_tv_e $nclasses_1 ",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'HKY85I'}   = {
        name => 'HKY85I',
        mbcmd => " $Lset $nst2 $propinv;\n $Prset $uneqfreqpr $tratiopr $pinvarpr;",
        phymlcmd => " HKY85 $f_m $ts_tv_e $pinv_e $nclasses_1 ",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'HKY85G'}   = {
        name => 'HKY85G',
        mbcmd => " $Lset $nst2 $gamma;\n $Prset $uneqfreqpr $tratiopr $shapepr;",
        phymlcmd => " HKY85 $f_m $ts_tv_e $nclasses_4 $alpha_e ",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'HKY85IG'}  = {
        name => 'HKY85IG',
        mbcmd => " $Lset $nst2 $invgamma;\n $Prset $uneqfreqpr $tratiopr $shapepr $pinvarpr;",
        phymlcmd => " HKY85 $f_m $ts_tv_e $pinv_e $nclasses_4 $alpha_e ",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    # ###################################### Modeltest extra models ###############################################
    $HoH{'TRNEF'}  = {
        name => 'TRNEF',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nR\nY\n",
        phymlcmd => " TN93 $f_equal $nclasses_1 ",
        df => 2 + $nbranch,
        dfnb => 2,
    };
    $HoH{'TRNEFI'}  = {
        name => 'TRNEFI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
        phymlcmd => " TN93 $f_equal $pinv_e $nclasses_1 ",
        df => 3 + $nbranch,
        dfnb => 3,
    };
    $HoH{'TRNEFG'}  = {
        name => 'TRNEFG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nY\n",
        phymlcmd => " TN93 $f_equal $nclasses_4 $alpha_e",
        df => 3 + $nbranch,
        dfnb => 3,
    };
    $HoH{'TRNEFIG'}  = {
        name => 'TRNEFIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nV\nY\nY\n",
        phymlcmd => " TN93 $f_equal $pinv_e $nclasses_4 $alpha_e",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'TRN'}  = {
        name => 'TRN',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nF\nT\nY\nR\nY\n",
        phymlcmd => " TN93 $f_m $nclasses_1 ",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'TRNI'}  = {
        name => 'TRNI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nF\nT\nY\nR\nV\nY\nY\n",
        phymlcmd => " TN93 $f_m $pinv_e $nclasses_1 ",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'TRNG'}  = {
        name => 'TRNG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nF\nT\nY\nY\n",
        phymlcmd => " TN93 $f_m $nclasses_4 $alpha_e",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'TRNIG'}  = {
        name => 'TRNIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nF\nT\nY\nV\nY\nY\n",
        phymlcmd => " TN93 $f_m $pinv_e $nclasses_4 $alpha_e",
        df => 7 + $nbranch,
        dfnb => 7,
    };
    $HoH{'K3P'}  = {
        name => 'K3P',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nR\nY\n",
        phymlcmd => " $custom_k3p $f_equal $nclasses_1 ",
        df => 2 + $nbranch,
        dfnb => 2,
    };
    $HoH{'K3PI'}  = {
        name => 'K3PI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
        phymlcmd => " $custom_k3p $f_equal $pinv_e $nclasses_1 ",
        df => 3 + $nbranch,
        dfnb => 3,
    };
    $HoH{'K3PG'}  = {
        name => 'K3PG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nY\n",
        phymlcmd => " $custom_k3p $f_equal $nclasses_4 $alpha_e",
        df => 3 + $nbranch,
        dfnb => 3,
    };
    $HoH{'K3PIG'}  = {
        name => 'K3PIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nV\nY\nY\n",
        phymlcmd => " $custom_k3p $f_equal $pinv_e $nclasses_4 $alpha_e",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'K3PUF'}  = {
        name => 'K3PUF',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nR\nY\n",
        phymlcmd => " $custom_k3p $f_m $nclasses_1 ",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'K3PUFI'}  = {
        name => 'K3PUFI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
        phymlcmd => " $custom_k3p $f_m $pinv_e $nclasses_1 ",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'K3PUFG'}  = {
        name => 'K3PUFG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nY\n",
        phymlcmd => " $custom_k3p $f_m $nclasses_4 $alpha_e",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'K3PUFIG'}  = {
        name => 'K3PUFIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nV\nY\nY\n",
        phymlcmd => " $custom_k3p $f_m $pinv_e $nclasses_4 $alpha_e",
        df => 7 + $nbranch,
        dfnb => 7,
    };
    $HoH{'TIMEF'}  = {
        name => 'TIMEF',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nY\n",
        phymlcmd => " $custom_tim $f_equal $nclasses_1 ",
        df => 3 + $nbranch,
        dfnb => 3,
    };
    $HoH{'TIMEFI'}  = {
        name => 'TIMEFI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
        phymlcmd => " $custom_tim $f_equal $pinv_e $nclasses_1 ",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'TIMEFG'}  = {
        name => 'TIMEFG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nY\n",
        phymlcmd => " $custom_tim $f_equal $nclasses_4 $alpha_e",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'TIMEFIG'}  = {
        name => 'TIMEFIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n",
        phymlcmd => " $custom_tim $f_equal $pinv_e $nclasses_4 $alpha_e",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'TIM'}  = {
        name => 'TIM',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nY\n",
        phymlcmd => " $custom_tim $f_m $nclasses_1 ",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'TIMI'}  = {
        name => 'TIMI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
        phymlcmd => " $custom_tim $f_m $pinv_e $nclasses_1 ",
        df => 7 + $nbranch,
        dfnb => 7,
    };
    $HoH{'TIMG'}  = {
        name => 'TIMG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nY\n",
        phymlcmd => " $custom_tim $f_m $nclasses_4 $alpha_e",
        df => 7 + $nbranch,
        dfnb => 7,
    };
    $HoH{'TIMIG'}  = {
        name => 'TIMIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n",
        phymlcmd => " $custom_tim $f_m $pinv_e $nclasses_4 $alpha_e",
        df => 8 + $nbranch,
        dfnb => 8,
    };
    $HoH{'TVMEF'}  = {
        name => 'TVMEF',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n",
        phymlcmd => " $custom_tvm $f_equal $nclasses_1 ",
        df => 4 + $nbranch,
        dfnb => 4,
    };
    $HoH{'TVMEFI'}  = {
        name => 'TVMEFI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
        phymlcmd => " $custom_tvm $f_equal $pinv_e $nclasses_1 ",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'TVMEFG'}  = {
        name => 'TVMEFG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n",
        phymlcmd => " $custom_tvm $f_equal $nclasses_4 $alpha_e",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'TVMEFIG'}  = {
        name => 'TVMEFIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n",
        phymlcmd => " $custom_tvm $f_equal $pinv_e $nclasses_4 $alpha_e",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'TVM'}  = {
        name => 'TVM',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n",
        phymlcmd => " $custom_tvm $f_m $nclasses_1 ",
        df => 7 + $nbranch,
        dfnb => 7,
    };
    $HoH{'TVMI'}  = {
        name => 'TVMI',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
        phymlcmd => " $custom_tvm $f_m $pinv_e $nclasses_1 ",
        df => 8 + $nbranch,
        dfnb => 8,
    };
    $HoH{'TVMG'}  = {
        name => 'TVMG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n",
        phymlcmd => " $custom_tvm $f_m $nclasses_4 $alpha_e",
        df => 8 + $nbranch,
        dfnb => 8,
    };
    $HoH{'TVMIG'}  = {
        name => 'TVMIG',
        mbcmd => "",
        #phymlcmd => "$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n",
        phymlcmd => " $custom_tvm $f_m $pinv_e $nclasses_4 $alpha_e",
        df => 9 + $nbranch,
        dfnb => 9,
    };
    #
    ###################################### End Modeltest extra models ############################################
    $HoH{'SYM'}    = {
        name => 'SYM',
        mbcmd => " $Lset $nst6 $equal;\n $Prset $revmatpr $eqfreqpr;",
        phymlcmd => " $custom_sym $f_equal $nclasses_1 ",
        df => 5 + $nbranch,
        dfnb => 5,
    };
    $HoH{'SYMI'}   = {
        name => 'SYMI',
        mbcmd => " $Lset $nst6 $propinv;\n $Prset $revmatpr $eqfreqpr $pinvarpr;",
        phymlcmd => " $custom_sym $f_equal $pinv_e $nclasses_1 ",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'SYMG'}   = {
        name => 'SYMG',
        mbcmd => " $Lset $nst6 $gamma;\n $Prset $revmatpr $eqfreqpr;",
        phymlcmd => " $custom_sym $f_equal $nclasses_4 $alpha_e ",
        df => 6 + $nbranch,
        dfnb => 6,
    };
    $HoH{'SYMIG'}  = {
        name => 'SYMIG',
        mbcmd => " $Lset $nst6 $invgamma;\n $Prset $revmatpr $eqfreqpr $shapepr $pinvarpr;",
        phymlcmd => " $custom_sym $f_equal $pinv_e $nclasses_4 $alpha_e ",
        df => 7 + $nbranch,
        dfnb => 7,
    };
    $HoH{'GTR'}    = {
        name => 'GTR',
        mbcmd => " $Lset $nst6 $equal;\n $Prset $revmatpr $uneqfreqpr;",
        phymlcmd => " GTR $f_m $nclasses_1 ",
        df => 8 + $nbranch,
        dfnb => 8,
    };
    $HoH{'GTRI'}   = {
        name => 'GTRI',
        mbcmd => " $Lset $nst6 $propinv;\n $Prset $revmatpr $uneqfreqpr $pinvarpr;",
        phymlcmd => " GTR $f_m $pinv_e $nclasses_1 ",
        df => 9 + $nbranch,
        dfnb => 9,
    };
    $HoH{'GTRG'}     = {
        name => 'GTRG',
        mbcmd => " $Lset $nst6 $gamma;\n $Prset $revmatpr $uneqfreqpr $shapepr;",
        phymlcmd => " GTR $f_m $nclasses_4 $alpha_e ",
        df => 9 + $nbranch,
        dfnb => 9,
    };
    $HoH{'GTRIG'}  = {
        name => 'GTRIG',
        mbcmd => " $Lset $nst6 $invgamma;\n $Prset $revmatpr $uneqfreqpr $shapepr $pinvarpr;",
        phymlcmd => " GTR $f_m $pinv_e $nclasses_4 $alpha_e ",
        df => 10 + $nbranch,
        dfnb => 10,
    };

} # end of initialize_HoH





#===  FUNCTION  ================================================================
#         NAME:  is_sequential
#      VERSION:  12/08/2009 02:29:14 PM CET
#  DESCRIPTION:  Try to see if file is interleaved or not
#   PARAMETERS:  filename
#      RETURNS:  '--sequential' (=sequential), '' (=interleaved)
#         TODO:  ???
#===============================================================================
sub is_sequential {

    my ($infile_name) = @_;
    my $max_line_length = 0;
    my $nlines = 0;
    my $nchar = q{};

    ## Open file once to get the $ntax, $nchar
    open my $IN, "<", $infile_name or die "could not open infile: $! \n";
    while(<$IN>) {
        if (/\s*\d+\s+(\d+)/) {
            $nchar = $1;
            last;
        }
    }
    close($IN) or warn "could not close infile : $! \n";

    ## Open file again
    open my $FI, "<", $infile_name or die "\a\nERROR!\n\nCould not open infile for reading: $!\n";
    while (<$FI>) {
        chomp();
        my $line = $_;
        next if (/^\s+$/);
        $nlines++;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        $line =~ s/\s+/ /g;
        if(length($line) > $max_line_length) {
            $max_line_length = length($line);
        }
    }
    close($FI) or warn "Could not close file handle: $! \n";

    if ($max_line_length < $nchar) { # File is probably interleaved if
        return('');                  # no line is as long or longer than
    }                                # number of characters.
    else {                           # Else it's probably sequential
        return('--sequential');
    }

} # end of is_sequential


#===  FUNCTION  ================================================================
#         NAME:  pmraic_print
#      VERSION:  11/10/2011 06:52:00 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub pmraic_print {

    ## Set output filename
    my $pMrAIC_output = $infile . $output_file_suffix;
    open my $OUT, '>', $pMrAIC_output or die "\a\nERROR!\n\nCould not create output file $!\n";

    ## Get the date
    my $datestring = prettydate();

    ## Some dimensions
    my $ntax = $dimensions[0];
    my $nchar = $dimensions[1];
    my $max_num_params = $dimensions[3];
    my $df = 'df';

    ## Check parameters
    if ($max_num_params > $nchar) {
        $max_num_params = 10; # GTRIG
        $use_brlens = 0;
        $df = 'dfnb';

        #$estMinBASize = ($BAratio * $maxNumParametersBrLens);
        #$estMinSize = ($maxNumParametersBrLens);
        printf  $OUT "\n\nWARNING:\n";
        printf  $OUT "Number of parameters is larger than the sample size!\n";
        printf  $OUT "This causes problems in AICc calculations and you need to get\n";
        printf  $OUT "more sequence data!\n";
        printf  $OUT "\t(For %g terminal branches you probably need more\n", $ntax;
        #printf  $OUT "\tthan %g characters, and an approx. data size for\n", $estMinSize;
        #printf  $OUT "\thaving a nchar/nparam ratio over %g is about %g bp.)\n", $BAratio, $estMinBASize;
        printf  $OUT "Calculations were done anyway while not considering branch\n";
        printf  $OUT "lengths in the total number of parameters. Number of\n";
        #printf  $OUT "parameters was set to $maxNumParameters (number of free parameters\n";
        printf  $OUT "in the GTR+I+G model).\n\n"; # NOTE: max_num_params need to be set dynamically according to the current models in the run_HoH, not necessarily to 10 (GTRIG)!
    }


    ## Best models
    my $minimumAICmodel  = get_min_aic_model();
    my $minimumAICcmodel = get_min_aicc_model();
    my $minimumBICmodel  = get_min_bic_model();

    ## Print to the outputfile
    printf $OUT "\nOutput from %s version %s by %s\n",$PROGRAM_NAME, $VERSION, $AUTHOR;
    print  $OUT "-------------------------------------------------------------\n";
    printf $OUT "$datestring\n";        
    printf $OUT "\nInput data from file \"%s\" (%s taxa, %s characters)\n", $infile, $ntax, $nchar;
    printf $OUT "\nMinimum AIC  model: %s", $minimumAICmodel;
    printf $OUT "\nMinimum AICc model: %s", $minimumAICcmodel;
    printf $OUT "\nMinimum BIC  model: %s", $minimumBICmodel;
    print  $OUT "\n";

        ## See if the number of parameters are small compared to sample size
        #if($changedNparams eq TRUE) {
        #    $estMinBASize = ($BAratio * $maxNumParametersBrLens);
        #    $estMinSize = ($maxNumParametersBrLens);
        #    printf  $OUT "\n\nWARNING:\n";
        #    printf  $OUT "Number of parameters is larger than the sample size!\n";
        #    printf  $OUT "This causes problems in AICc calculations and you need to get\n";
        #    printf  $OUT "more sequence data!\n";
        #    printf  $OUT "\t(For %g terminal branches you probably need more\n", $ntax;
        #    printf  $OUT "\tthan %g characters, and an approx. data size for\n", $estMinSize;
        #    printf  $OUT "\thaving a nchar/nparam ratio over %g is about %g bp.)\n", $BAratio, $estMinBASize;
        #    printf  $OUT "Calculations were done anyway while not considering branch\n";
        #    printf  $OUT "lengths in the total number of parameters. Number of\n";
        #    printf  $OUT "parameters was set to $maxNumParameters (number of free parameters\n";
        #    printf  $OUT "in the GTR+I+G model).\n\n";
        #}

        ## Output sorted by AICc
        #elsif (($nchar/$maxNumParameters) < $BAratio) { # Esentially testing if nchar < 400
            #printf  $OUT "\n\nNote: number of parameters compared to sample size\n";
            #printf  $OUT "is low (Nchar/Nparams < %g). Consider the use of AICc.\n", $BAratio;

        ## Output sorted by AICc
        my $cumAICcWeight = 0.0;
        printf  $OUT "\n\n(Output sorted by AICc. w: Akaike weight, cw: cumulative w.)";
        printf  $OUT "\n\nModel\tdf\tlnL\tAICc\twAICc\tcwAICc";
        foreach my $model (sort { $run_HoH{$a}->{'aicc'} <=> $run_HoH{$b}->{'aicc'} } keys %run_HoH) {
            $cumAICcWeight += $run_HoH{$model}->{'aicc_weight'};
            print $model, "\n";
            print $df, "\n";
            print $run_HoH{$model}->{$df}, "\n";
            print $run_HoH{$model}->{'df'}, "\n";
            print $run_HoH{$model}->{'dfnb'}, "\n";
            warn "\n HERE (hit return to continue)\n" and getc();

            printf $OUT "\n%s\t%g\t%.4f\t%.4f\t%.4f\t%.4f", $run_HoH{$model}->{'name'}, $run_HoH{$model}->{$df}, $run_HoH{$model}->{'lnl'}, $run_HoH{$model}->{'aicc'}, $run_HoH{$model}->{'aicc_weight'}, $cumAICcWeight; # JN: Missing value here 
        }
        print  $OUT "\n\n";
        #}

        ## Output sorted by AIC
        my $cumAICWeight = 0.0;
        printf  $OUT "\n\n(Output sorted by AIC. w: Akaike weight, cw: cumulative w.)";
        printf  $OUT "\n\nModel\tdf\tlnL\tAIC\twAIC\tcwAIC";
        foreach my $model (sort { $run_HoH{$a}->{'aic'} <=> $run_HoH{$b}->{'aic'} } keys %run_HoH) {
            $cumAICWeight += $run_HoH{$model}->{'aic_weight'};
            printf $OUT "\n%s\t%g\t%.4f\t%.4f\t%.4f\t%.4f", $run_HoH{$model}->{'name'}, $run_HoH{$model}->{$df}, $run_HoH{$model}->{'lnl'}, $run_HoH{$model}->{'aic'}, $run_HoH{$model}->{'aic_weight'}, $cumAICWeight;
        }

        print  $OUT "\n\n";

        ## Output sorted by BIC
        my $cumBICWeight = 0.0;
        printf  $OUT "\n\n(Output sorted by BIC. w: Akaike weight, cw: cumulative w.)";
        printf  $OUT "\n\nModel\tdf\tlnL\tBIC\twBIC\tcwBIC";
        foreach my $model (sort { $run_HoH{$a}->{'bic'} <=> $run_HoH{$b}->{'bic'} } keys %run_HoH) {
            $cumBICWeight += $run_HoH{$model}->{'bic_weight'};
            printf $OUT "\n%s\t%g\t%.4f\t%.4f\t%.4f\t%.4f", $run_HoH{$model}->{'name'}, $run_HoH{$model}->{$df}, $run_HoH{$model}->{'lnl'}, $run_HoH{$model}->{'bic'}, $run_HoH{$model}->{'bic_weight'}, $cumBICWeight;
        }
        print  $OUT "\n\n";
        
        ## Check if the best models are among the MrBayes models
        #
        ## Print MrBayes block
        #if ($useModeltest eq FALSE) {
        my $minimumAICmbcmd = $run_HoH{$minimumAICmodel}->{'mbcmd'};
        my $minimumAICcmbcmd = $run_HoH{$minimumAICcmodel}->{'mbcmd'};
        my $minimumBICmbcmd = $run_HoH{$minimumBICmodel}->{'mbcmd'};
        my $dataSetName = "$infile";
        my $partitionName = "Dummy";

        #print  $OUT "\n\n";
        #print  $OUT "To specify default values for best AIC, AICc, or BIC models in MrBayes:\n";
        print  $OUT "\n";

        if (!$run_HoH{$minimumAICcmodel}->{'mbcmd'} eq "") {
            print  $OUT "To specify default values for best AIC model in MrBayes:\n\n";
            printf $OUT "[Mrbayes block for the best AIC model (%s)]\n", $minimumAICmodel;
            print  $OUT "BEGIN MRBAYES;\n";
            printf $OUT " Charset %s = 1 - %s;\n",$dataSetName, $nchar;
            printf $OUT " Partition %s = 1:%s;\n",$partitionName, $dataSetName;
            printf $OUT " Set partition = %s;\n",$partitionName;
            printf $OUT "%s\n", $minimumAICmbcmd;
            print  $OUT "END;\n";
            print  $OUT "\n\n";
        }
        else {
            print  $OUT "\nThe best AICc model \'$minimumAICcmodel\' is not implemented in MrBayes v.3!\n"
        }

        if (!$run_HoH{$minimumAICcmodel}->{'mbcmd'} eq "") {
            print  $OUT "To specify default values for best AICc model in MrBayes:\n\n";
            printf $OUT "[Mrbayes block for the best AICc model (%s)]\n", $minimumAICcmodel;
            print  $OUT "BEGIN MRBAYES;\n";
            printf $OUT " Charset %s = 1 - %s;\n",$dataSetName, $nchar;
            printf $OUT " Partition %s = 1:%s;\n",$partitionName, $dataSetName;
            printf $OUT " Set partition = %s;\n",$partitionName;
            printf $OUT "%s\n", $minimumAICcmbcmd;
            print  $OUT "END;\n";
            print  $OUT "\n\n";
        }
        else {
            print  $OUT "\nThe best AICc model \'$minimumAICcmodel\' is not implemented in MrBayes v.3!\n"
        }

        if (!$run_HoH{$minimumBICmodel}->{'mbcmd'} eq "") {
            print  $OUT "To specify default values for best BIC in MrBayes:\n\n";
            printf $OUT "[Mrbayes block for the best BIC model (%s)]\n", $minimumBICmodel;
            print  $OUT "BEGIN MRBAYES;\n";
            printf $OUT " Charset %s = 1 - %s;\n",$dataSetName, $nchar;
            printf $OUT " Partition %s = 1:%s;\n",$partitionName, $dataSetName;
            printf $OUT " Set partition = %s;\n",$partitionName;
            printf $OUT "%s\n", $minimumBICmbcmd;
            print  $OUT "END;\n";
        }
        else {
            print  $OUT "\nThe best BIC model \'$minimumBICmodel\' is not implemented in MrBayes v.3!\n"
        }

        #}
        print  $OUT "\n\n";
        print  $OUT "-------------------------------------------------------------\n";
        print  $OUT "End of Output\n\n";

    close($OUT);

    ## Print trees to treefiles
    open (AICTREE, ">$infile.AIC-$minimumAICmodel.tre") or die "\a\nERROR!\n\nCould not create tree file $!\n";
        print AICTREE $run_HoH{$minimumAICmodel}->{'tree'}, "\n";
    close(AICTREE);
    open (AICcTREE, ">$infile.AICc-$minimumAICcmodel.tre") or die "\a\nERROR!\n\nCould not create tree file $!\n";
        print AICcTREE $run_HoH{$minimumAICcmodel}->{'tree'}, "\n";
    close(AICcTREE);
    open (BICTREE, ">$infile.BIC-$minimumBICmodel.tre") or die "\a\nERROR!\n\nCould not create tree file $!\n";
        print BICTREE $run_HoH{$minimumBICmodel}->{'tree'}, "\n";
    close(BICTREE);


} # end of pmraic_print


#===  FUNCTION  ================================================================
#         NAME:  prettydate
#      VERSION:  12/03/2009 04:49:52 PM CET
#  DESCRIPTION:  Print time stamp
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub prettydate {

   @_ = localtime(shift || time);

   return(sprintf("%02d:%02d %02d/%02d/%04d", @_[2,1], $_[4]+1, $_[3], $_[5]+1900));

} # end of prettydate


#===  FUNCTION  ================================================================
#         NAME:  randomize_array
#      VERSION:  12/03/2009 03:01:23 PM CET
#  DESCRIPTION:  Randomises an array. Taken from http://www.perlcircus.org/arrays.shtml.
#   PARAMETERS:  array
#      RETURNS:  array
#         TODO:  ???
#===============================================================================
sub randomize_array {

    my @array = @_;

    foreach my $i (reverse 0..$#array) {
        my $r = int rand ($i+1);
        @array[$i, $r] = @array[$r, $i] unless ($i == $r);
    }

    return @array;

} # end of randomize_array


#===  FUNCTION  ================================================================
#         NAME:  round_down
#      VERSION:  12/03/2009 04:49:52 PM CET
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub round_down {

    my $n = shift;

    return(($n == int($n)) ? $n : int($n));

} # end of round_down


#===  FUNCTION  ================================================================
#         NAME:  round_up
#      VERSION:  12/03/2009 04:49:52 PM CET
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub round_up {

    my $n = shift;
    
    return(($n == int($n)) ? $n : int($n + 1));

} # end of round_up


#===  FUNCTION  ================================================================
#         NAME:  run_fork_phyml
#      VERSION:  08/30/2011 09:42:56 PM CEST
#  DESCRIPTION:  run phyml
#   PARAMETERS:  hash reference?
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
sub run_fork_phyml {

    my $hash_ref = shift;

    my $sequential = is_sequential($infile);

    ## Do the Parallel stuff
    my $pm = new Parallel::ForkManager($cpu);

    ## Calculate how many models per fork:
    #my $hash_size = (%initialize_all_model_commands_hash);
    #my $models_per_hash = $hash_size / $cpu;

    ## Divide the model hash into  $models_per_hash lists.


    #$pm->run_on_finish(
    #    sub { my ($pid, $exit_code) = @_;
    #        print "just got out of the pool ".
    #        "with PID $pid and exit code: $exit_code\n";
    #    }
    #);


    my @fork_models = keys %{$hash_ref};

    #print STDERR "Running $cpu instances of phyml.\n" if $verbose;
    #warn "Do Ctrl+c to quit now or hit return to start!\n" and getc(); # temporary


    ## Submit lists to forking:
    foreach my $model (@fork_models) {

        # Forks and returns the pid for the child:
        my $pid = $pm->start and next;

        ##... do some work with $data in the child process ...
        ## Need to copy the input file to a file with a model-specific name due to the problem of phyml writing
        ## files with the same name
        my $model_infile = $output_file_prefix . $model . '.' . $infile;
        copy($infile, $model_infile) or die "File cannot be copied.";

        ## Run PHYML
        print STDERR " model $hash_ref->{$model}->{'name'}\n" if $verbose;
        system("$phyml_bin --input $model_infile --datatype $datatype $sequential --model $hash_ref->{$model}->{'phymlcmd'} --search $search -o $optimize > /dev/null ");


        ## Parse the stat/stats file to get the likelihood
        #my $model_stat_file = $model_infile . $phyml_stats_ending;
        #my $likelihood = get_likelihood($model_stat_file);
        #print "likelihood:", $likelihood, "\n";
        #$model_likelihood_hash{$model} = $likelihood;

        #my $model_tree_file = $model_infile . $phyml_tree_ending;
        #my $best_tree = get_best_tree($model_tree_file);
        #$model_tree_hash{$model} = $best_tree;

        $pm->finish; # Terminates the child process
    }

    $pm->wait_all_children;

} # end of run_fork_phyml


#===  FUNCTION  ================================================================
#         NAME:  summarize_files 
#      VERSION:  08/21/2015 01:36:31 PM
#  DESCRIPTION:  summarize info in files
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  BUG: run_HoH is not initialized correclty if running script using --summarize
#                Need to set model name, df, and nbranches, before calling pmraic_print
#                First, make sure initialize_HoH is called, then iterate over all files to be summarized, and:
#                @models = check_models_arg(@$models);
#                @run_HoH{@models} = @HoH{@models}; # hash slice to create %run_HoH from a subset of %HoH
#===============================================================================
sub summarize_files {

    my $glob_expression = $output_file_prefix . "*";
    my (@model_files) = glob($glob_expression);
    
#    print "####### run_HoH 0: #################\n\n";
#    print Dumper(%run_HoH) and getc();

    ## Check model names in model_files, and initialze %run_HoH
    my %found_hash = (); #key: model, value:1
    foreach my $file (@model_files) {
        my (@arr) = split /\./, $file;
        my $mod = $arr[1];
        $found_hash{$mod}++;
    }
    my (@found_models) = keys (%found_hash);
    my @models = check_models_arg(@found_models);

    @run_HoH{@models} = @HoH{@models}; # hash slice to create %run_HoH from a subset of %HoH

    print "####### run_HoH 1: #################\n\n";
    print Dumper(%run_HoH) and getc();


    if ($verbose) {
        print STDERR "Summarizing output\n";
        print STDERR "Found these files:\n\n";
        foreach my $file (@model_files) {
            print STDERR $file, "\n";
        }
        warn "\n\nWill summarize the files above, hit return to continue (or Ctrl+c to quit).\nNote: to avoid interactive use, try the --noverbose flag.\n" and getc();
    }

    ## Get likelihood values and trees
    foreach my $file (@model_files) {
        if ($file =~ /$phyml_stats_ending/) {
            my (@arr) = split /\./, $file;
            my $mod = $arr[1];
            #if (exists $run_HoH{$mod}) {
                my $lnl = get_likelihood($file); 
                $run_HoH{$mod}->{'lnl'} = $lnl;
            #}
            #else {
            #    die "model name \'$mod\' from file name \'$file\' not valid.\nUse \'$PROGRAM_NAME --help\' to see the valid choices.\n"; 
            #}
        }
        elsif ($file =~ /$phyml_tree_ending/) {
            my (@arr) = split /\./, $file;
            my $mod = $arr[1];
            #if (exists $run_HoH{$mod}) {
                my $tre = get_best_tree($file); 
                $run_HoH{$mod}->{'tree'} = $tre;
            #}
            #else {
            #    die "model name \'$mod\' from file name \'$file\' not valid.\nUse \'$PROGRAM_NAME --help\' to see the valid choices.\n"; 
            #}
        }
    }

    ## Calculate AIC, AICc, BIC
    foreach my $model (keys %run_HoH) {
        $run_HoH{$model}->{'aic'}  = calculate_aic($model);
        $run_HoH{$model}->{'aicc'} = calculate_aicc($model);
        $run_HoH{$model}->{'bic'}  = calculate_bic($model);
    }

    ## Calculate AIC-, AICc-, BIC- weights
    calculate_aic_weights();
    calculate_aicc_weights();
    calculate_bic_weights();

    print Dumper(%run_HoH) if $debug;
     print Dumper(%run_HoH);warn "\n(run_HoH) (hit return to continue)\n" and getc(); 
   
    pmraic_print();

    #clean_up();

} # end of summarize_files


#===  FUNCTION  ================================================================
#         NAME:  "MAIN"
#      VERSION:  09/01/2011 11:21:00 AM
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
GETOPTIONS:

## Get options

if (@ARGV < 1) {
    print STDERR "\n Try '$0 --man' for full info\n\n";
    exit(0);
}
else {
    my $gopts = GetOptions (
        "help"         => sub { pod2usage(1); },
        "version"      => sub { print STDOUT "\n $0 version $VERSION\n  Last changes $CHANGES\n"; exit(0) },
        "man"          => sub { pod2usage(-exitstatus => 0, -verbose => 2); },
        "cpu=i"         => \$cpu,
        "infile=s"      => \$infile,
        "verbose!"      => \$verbose,
        "outfile=s"     => \$outfile,
        "phyml=s"       => \$phyml,
        "search=s"      => \$search,
        "optimize=s"    => \$optimize,
        "summarize"     => \$summarize,
        "nosummarize"   => \$nosummarize,
        "bash=i"        => \$bash,
        "maui"          => \$maui,
        "modeltest"     => \$modeltest,
        "models=s@"     => \$models,
    );
    die "wrong argument(s)?\n" unless ($gopts);
}

CHECKING:
## Find local phyml
if ($phyml) {
    $phyml_bin = find_prog($phyml);
}
else {
    $phyml_bin = find_prog($phyml_bin) unless ($maui || $bash); # No error checking for Maui/Bash users who wish to manually supply the phyml path 
}

## Get infile
if ($infile) {
    $infile = $infile;
}
elsif (@ARGV > 0) {
    $infile = shift(@ARGV);
}
## Get and check dimensions
@dimensions = get_dimensions($infile); # ntax, nchar, nbranch, max_num_params
if ($dimensions[3] > $dimensions[1]) { # if (max_num_params > nchar)
    $max_num_params = 10; # NOTE: max_num_params need to be set dynamically according to the current models in the run_HoH, not necessarily to 10 (GTRIG)!
    $use_brlens = 0;
}

#$sequential = is_sequential($infile);

## Check some args
$optimize = check_optimize_arg($optimize); # check correct syntax and lower case.
$search = check_search_arg($search); # check correct syntax and upper case.

## Initilalize some hashes
initialize_HoH(); # Initialize all 56 models

## Run a subset of the models?
if ($models) { # --model=x,x,x
    @models = check_models_arg(@$models);
    @run_HoH{@models} = @HoH{@models}; # hash slice to create %run_HoH from a subset of %HoH
    $n_models = scalar keys %run_HoH;
}
elsif ($modeltest) { # --modeltest
    %run_HoH = %HoH;
    $n_models = scalar keys %run_HoH;
}
else {
    @run_HoH{@mrmodeltest_models} = @HoH{@mrmodeltest_models}; # hash slice to create %run_HoH from a subset of %HoH
    $n_models = scalar keys %run_HoH;
}

## CPUs
if ($cpu) {
    my $requested_cpu = $cpu;
    if ($cpu > $MAX_PROCESSES) {
        $cpu = $MAX_PROCESSES;
        if ($cpu > $n_models) {
            $cpu = $n_models;
            print STDERR "Number of requested CPUs ($requested_cpu) is larger than the number of models to run ($n_models).\n";
            print STDERR "Adjusting --cpu=$cpu\n";
        }
        else {
            print STDERR "Number of requested CPUs ($requested_cpu) is larger than the available ($MAX_PROCESSES).\n";
            print STDERR "Adjusting --cpu=$cpu\n";
        }
    }
    if ($cpu > $n_models) {
        $cpu = $n_models;
        print STDERR "Number of requested CPUs ($requested_cpu) is larger than the number of models to run ($n_models).\n";
        print STDERR "Adjusting --cpu=$cpu\n";
    }
}
else {
    if ($MAX_PROCESSES > $n_models) {
        $cpu = $n_models;
    }
    else {
        $cpu = $MAX_PROCESSES;
    }
}

RUNNING:
#Get the date
my $datestring = prettydate();

## Create scripts or run directly
if ($bash > 0) {
    $summarize = 0;
    if($phyml) {
        $phyml_bash_bin = $phyml; # No error checking for Bash users who wish to manually supply the phyml path
    }
    else {
        $phyml_bash_bin = $phyml_bin;
    }
    create_bash_scripts($bash);
    print STDERR "Created bash scripts.\nCopy scripts to preferred machine(s) and submit them to the shell.\nThe output can be collected and later summarized using --summarize\n" if $verbose; 
    exit(0);
}
elsif ($maui) {
    $summarize = 0;
    if($phyml) {
        $phyml_maui_bin = $phyml; # No error checking for Maui users who wish to manually supply the phyml path
    }  
    else {
        $phyml_maui_bin = $phyml_bin;
    }
    create_maui_files(\%run_HoH);
    print STDERR "Created Maui scripts.\nIf files are OK, run the MAUI shell script.\n" if $verbose; 
    exit(0);
}
elsif ($summarize) {
    summarize_files();
    clean_up();
    print STDERR "Summarized files.\n" if $verbose; 
}
else {
    ## Run forked phyml
    print STDERR "Running $cpu instances of phyml.\n" if $verbose;
    if ($verbose) {
        warn "Hit return to start (or Ctrl+c to quit. Note: to avoid interactive use, try the --noverbose flag.)!\n" and getc();
        print STDERR "running...\n";
    }

    run_fork_phyml(\%run_HoH);

    if ( ! $nosummarize) { # hack...
        summarize_files();
        print STDERR "Summarized files.\nOutput is in $infile.*.tre, $infile.*.txt, $infile.*.tgz\n" if $verbose; 
        clean_up();
    }
}

## Print some stuff
print Dumper(%run_HoH) if $debug;
print STDERR "Done with $0\n" if $verbose;

## Exit
exit(0);


#===  POD DOCUMENTATION  =======================================================
#         NAME:  POD documentation
#      VERSION:  11/16/2011 02:48:20 PM
#  DESCRIPTION:  ???
#         TODO:  ???
#===============================================================================
=pod

=head1 NAME

pmraic.pl


=head1 VERSION

Documentation for pmraic.pl version 1.1


=head1 SYNOPSIS

pmraic.pl [options] FILE


=head1 DESCRIPTION

Script for running nucleotide substitutions using PHYML (Guindon & Gasquel, 2003) in parallel.
It will report the best AIC, AICc, BIC model(s), along with the best ML tree under that/those
models, and also the syntax for specifying the model(s) in MrBayes v.3 (Ronquist & Huelsenbeck, 2003).


=head1 OPTIONS

Mandatory arguments to long options are mandatory for short options too


=over 8

=item B<-b, --bash=>I<number>

Generate I<number> of bash scripts for running PHYML. Can be manually
distributed to several machines. Use I<--summarize> later for collecting the output.
Using this option will not start PHYML, but the analyses have to be started manually
(using the scripts).



=item B<-c, --cpu=>I<number>

Number of "CPUs" to use. Default is to run as many as the number of available cores.



=item B<-i, --infile=>I<infile>

Excplicitly specify the I<infile> to be used. The infile should have nt sequences and
be readable by phyml.



=item B<-mau, --maui>

Create submit files for the Maui cluster scheduler system. Use I<--summarize> later for collecting the output.



=item B<-models, --models=>I<arg1,arg2,...,argn>

Specify a subset of the 56 models available to be run by PHYML.
I<arg1,arg2,...,argn> should be given as a comma-separated list. No spaces.
Note that a MrBayes block will only be printed if the best model is among
the 24 nt models available in MrBayes v.3.
Valid model names are (first 24 are available in MrBayes):

    F81,F81G,F81I,F81IG,GTR,GTRG,
    GTRI,GTRIG,HKY85,HKY85G,HKY85I,HKY85IG,
    JC69,JC69G,JC69I,JC69IG,K2P,K2PG,
    K2PI,K2PIG,SYM,SYMG,SYMI,SYMIG

    K3P,K3PG,K3PI,K3PIG,K3PUF,K3PUFG,
    K3PUFI,K3PUFIG,TIM,TIMEF,TIMEFG,TIMEFI,
    TIMEFIG,TIMG,TIMI,TIMIG,TRN,TRNEF,
    TRNEFG,TRNEFI,TRNEFIG,TRNG,TRNI,TRNIG,
    TVM,TVMEF,TVMEFG,TVMEFI,TVMEFIG,TVMG,
    TVMI,TVMIG,


=item B<-modelt, --modeltest>

Run all 56 models.


=item B<-o, --optimize=>I<arg>

Specify the optimizing I<arg> for phyml. One of the following are allowed:
I<tlr>,  I<tl>,  I<lr>, I<l>,  I<r>, or  I<n>.
Default is I<tlr>, optimize all.



=item B<-p, --phyml=>I<arg>

Specify the I<PATH> to the phyml binary to be used. Good to use when
creating files for the Maui batch system or for running files on remote systems.
The default path for phyml is otherwise taken from the users current PATH.



=item B<-se, --search=>I<arg>

Specify the search I<arg> for phyml. One of the following are allowed:
I<NNI>,  I<SPR>, or  I<BEST>.
Default is I<NNI>, the default in phyml.



=item B<-su, --summarize>

Summarize the output from pmraic.pl/phyml. Can also be set to B<--nosummarize>.



=item B<-v, --verbose>

Be verbose. Can also be set to B<--noverbose>.



=item B<-h, --help>

Prints help message and exits


=item B<-v, --version>

Prints version message and exits


=item B<-man, --man>

Displays the manual page


=back


=head1 USAGE

Examples:

Run PHYML in parallel. Default is to start as many different runs as number of cores available

    pmraic.pl data.dat

Run PHYML in parallel and without asking for user input

    pmraic.pl --noverbose data.dat

Run two instances of PHYML in parallel (12 serial PHYML runs per CPU)

    pmraic.pl --cpu=2 data.dat

Generate six bash scripts for running serial (6X4) PHYML analyses.

    pmraic.pl --bash=6 data.dat 

Generate a Maui submit script for running PHYML on a cluster

    pmraic.pl --maui --phyml=/cluster/bin/phyml data.dat

Run four models in parallel

    pmraic.pl --models=JC69,JC69I,JC69G,JC69IG data.dat

Run 56 models (the "Modeltest" set) in parallel

    pmraic.pl --modeltest data.dat

Summarize output from previous run(s)

    pmraic.pl --summarize data.dat



=head1 AUTHOR

Written by Johan A. A. Nylander


=head1 REPORTING BUGS

See DEPENDENCIES below.

Please report any further bugs to I<jnylander @ users.sourceforge.net>.


=head1 DEPENDENCIES

Needs PhyML v. > 3.0 to run. Note that current precompiled versions of PhyML
("v3.0_360-500M") available from http://www.atgc-montpellier.fr/phyml/ contains a
bug preventing "I" models to be correctly analysed (alpha jointly estimated,
even when not supposed to). This might effect model evaluations, especially
if the program is used in "toggle mode". The bug seems to be fixed, however,
in latest versions (at least v."20110919") available on
http://code.google.com/p/phyml/.


Uses Perl modules Pod::Usage, Parallel::ForkManager.

Install on Ubuntu Linux:

  sudo apt-get install perl-doc libparallel-forkmanager-perl

or directly from CPAN:

  sudo perl -MCPAN -e 'install Parallel::ForkManager'




=head1 LICENSE AND COPYRIGHT

Copyright (c) 2009,2010,2011 Johan Nylander. All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. 
http://www.gnu.org/copyleft/gpl.html 


=head1 DOWNLOAD

https://github.com/nylander/pMrAIC




=cut




__END__

