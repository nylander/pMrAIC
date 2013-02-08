pMrAIC - Parallel MrAIC
======


NAME

    pmraic.pl


VERSION

    Documentation for pmraic.pl version 1.1


SYNOPSIS

    pmraic.pl [options] FILE


DESCRIPTION

    Script for running nucleotide substitutions using PHYML (Guindon &
    Gasquel, 2003) in parallel. It will report the best AIC, AICc, BIC
    model(s), along with the best ML tree under that/those models, and also
    the syntax for specifying the model(s) in MrBayes v.3 (Ronquist &
    Huelsenbeck, 2003).


OPTIONS

    Mandatory arguments to long options are mandatory for short options too

    -b, --bash=*number*
            Generate *number* of bash scripts for running PHYML. Can be
            manually distributed to several machines. Use *--summarize*
            later for collecting the output. Using this option will not
            start PHYML, but the analyses have to be started manually (using
            the scripts).

    -c, --cpu=*number*
            Number of "CPUs" to use. Default is to run as many as the number
            of available cores.

    -i, --infile=*infile*
            Excplicitly specify the *infile* to be used. The infile should
            have nt sequences and be readable by phyml.

    -mau, --maui
            Create submit files for the Maui cluster scheduler system. Use
            *--summarize* later for collecting the output.

    -models, --models=*arg1,arg2,...,argn*
            Specify a subset of the 56 models available to be run by PHYML.
            *arg1,arg2,...,argn* should be given as a comma-separated list.
            No spaces. Note that a MrBayes block will only be printed if the
            best model is among the 24 nt models available in MrBayes v.3.
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

    -modelt, --modeltest
            Run all 56 models.

    -o, --optimize=*arg*
            Specify the optimizing *arg* for phyml. One of the following are
            allowed: *tlr*, *tl*, *lr*, *l*, *r*, or *n*. Default is *tlr*,
            optimize all.

    -p, --phyml=*arg*
            Specify the *PATH* to the phyml binary to be used. Good to use
            when creating files for the Maui batch system or for running
            files on remote systems. The default path for phyml is otherwise
            taken from the users current PATH.

    -se, --search=*arg*
            Specify the search *arg* for phyml. One of the following are
            allowed: *NNI*, *SPR*, or *BEST*. Default is *NNI*, the default
            in phyml.

    -su, --summarize
            Summarize the output from pmraic.pl/phyml. Can also be set to
            --nosummarize.

    -v, --verbose
            Be verbose. Can also be set to --noverbose.

    -h, --help
            Prints help message and exits

    -v, --version
            Prints version message and exits

    -man, --man
            Displays the manual page


USAGE

    Examples:

    Run PHYML in parallel. Default is to start as many different runs as
    number of cores available

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


AUTHOR

    Written by Johan A. A. Nylander


REPORTING BUGS

    See DEPENDENCIES below.
    Please report any further bugs to *jnylander @ users.sourceforge.net*.


DEPENDENCIES

    Needs PhyML v. > 3.0 to run. Note that current precompiled versions of
    PhyML ("v3.0_360-500M") available from
    http://www.atgc-montpellier.fr/phyml/ contains a bug preventing "I"
    models to be correctly analysed (alpha jointly estimated, even when not
    supposed to). This might effect model evaluations, especially if the
    program is used in "toggle mode". The bug seems to be fixed, however, in
    latest versions (at least v."20110919") available on
    http://code.google.com/p/phyml/.

    Uses Perl modules Pod::Usage, Parallel::ForkManager.

    Install on Ubuntu Linux:

      sudo apt-get install perl-doc libparallel-forkmanager-perl

    or directly from CPAN:

      sudo perl -MCPAN -e 'install Parallel::ForkManager'


LICENSE AND COPYRIGHT

    Copyright (c) 2009,2010,2011,2012,2013 Johan Nylander. All rights reserved.

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details. http://www.gnu.org/copyleft/gpl.html


DOWNLOAD

    https://github.com/nylander/pMrAIC

