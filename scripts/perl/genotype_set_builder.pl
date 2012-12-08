#!/usr/bin/perl
# Copyright (c) 2012, Patrick Putnam
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.
#

use strict;
use Getopt::Long;
use Pod::Usage;

my $marker_count = 0;
my $individ_count = 0;
my $tplink;
my $case_control = 0;
my $help;

my $result = GetOptions( "help|?" => \$help,
                         "markers=i" => \$marker_count,
                         "individs=i" => \$individ_count,
                         "tplink=s"   => \$tplink, 
                         "case-control" => \$case_control ) or pod2usage(2);

pod2usage(2) if $help;

if( $marker_count == 0 || $individ_count == 0 ) {
    print "Please specify both --markers and --individs options for the marker and individual counts, respectively.\n";
    exit();
}

if( $tplink =~ m/^$/ ) {
    print "Please specify --tplink option.\n";
    exit();
}

if( $case_control ) {
    print "Starting to generate case/control set\n";
    make_tplink_cc ( $marker_count, $individ_count, $tplink );
} else {
    print "Starting to generate simple data set\n";
    make_tplink ( $marker_count, $individ_count, $tplink );
}

print "DONE\n";

sub make_tplink_cc {
    my ($mcount, $icount, $path) = @_;
    print $path."\n";
    open TPED, ">$path.tped" or die $!;
    open TFAM, ">$path.tfam" or die $!;
    open CASE_EXPECTED, ">$path.expected.case.dist" or die $!;
    open CTRL_EXPECTED, ">$path.expected.control.dist" or die $!;

    my @case_control = ();
    for ( my $i = 0; $i < $icount; $i++ ) {
        my $c = int( rand(2) );
        push (@case_control, $c);
        print TFAM "Fam_$i\tInd_$i\tPat_$i\tMat_$i\tx\t$c\n";
    }

    for ( my $i = 0; $i < $mcount; $i++ ) {
        my ($res, $case_dist, $ctrl_dist) = build_genotype_for_marker_cc( $icount, \@case_control );

        print TPED "0\trs$i\t0\t$i"."$res\n";
        print CASE_EXPECTED "$case_dist\n";
        print CTRL_EXPECTED "$ctrl_dist\n";
    }

    close(TPED);
    close(TFAM);
    close(CASE_EXPECTED);
    close(CTRL_EXPECTED);
}

sub make_tplink {
    my ($mcount, $icount, $path) = @_;
    open TPED, ">$path.tped" or die $!;
    open TFAM, ">$path.tfam" or die $!;
    open EXPECTED, ">$path.expected.dist" or die $!;

    for ( my $i = 0; $i < $icount; $i++ ) {
        print TFAM "Fam_$i\tInd_$i\tPat_$i\tMat_$i\tx\t1\n";
    }

    for ( my $i = 0; $i < $mcount; $i++ ) {
        my ($res, $dist) = build_genotype_for_marker( $icount );

        print TPED "0\trs$i\t0\t$i"."$res\n";
        print EXPECTED "$dist\n";
    }

    close(TPED);
    close(TFAM);
    close(EXPECTED);
}

sub build_genotype_for_marker_cc {
    my ($icount, $listptr) = @_;
    my @cc_list = @$listptr;

    my $case_xx = 0;
    my $case_aa = 0;
    my $case_ab = 0;
    my $case_bb = 0;

    my $ctrl_xx = 0;
    my $ctrl_aa = 0;
    my $ctrl_ab = 0;
    my $ctrl_bb = 0;

    my $res = "";

    for ( my $i = 0; $i < $icount; $i++ ) {
        my $v = int( rand( 100 ) );

        if( $cc_list[ $i ] ) {
            if( $v < 0 ) {
                $case_xx++;
                $res = $res."\t0\t0";
            } elsif ( $v < 34 ) {
                $case_aa++;
                $res = $res."\tA\tA";
            } elsif ( $v < 67 ) {
                $case_ab++;
                $res = $res."\tA\tC";
            } else {
                $case_bb++;
                $res = $res."\tC\tC";
            }
        } else {
            if( $v < 0) {
                $ctrl_xx++;
                $res = $res."\t0\t0";
            } elsif ( $v < 34 ) {
                $ctrl_aa++;
                $res = $res."\tA\tA";
            } elsif ( $v < 67 ) {
                $ctrl_ab++;
                $res = $res."\tA\tC";
            } else {
                $ctrl_bb++;
                $res = $res."\tC\tC";
            }
        }
    }

    my $case_dist = "$case_xx\t$case_aa\t$case_ab\t$case_bb";
    my $ctrl_dist = "$ctrl_xx\t$ctrl_aa\t$ctrl_ab\t$ctrl_bb";
    chomp($res);
    return ($res, $case_dist, $ctrl_dist);
}

sub build_genotype_for_marker {
    my $icount = $_[0];

    my $xx = 0;
    my $aa = 0;
    my $ab = 0;
    my $bb = 0;

    my $res = "";
    my $dist = "";

    for( my $i = 0; $i < $icount; $i++ ) {
        my $v = int(rand(4));

        if( $v == 0 ) { 
            $xx++;
            $res = $res."\t0\t0";
        } elsif( $v == 1 ) { 
            $aa++;
            $res = $res."\tA\tA";
        } elsif( $v == 2 ) { 
            $ab++;
            $res = $res."\tA\tC";
        } else { 
            $bb++;
            $res = $res."\tC\tC";
        }
    }

    chomp($res);
    $dist = "$xx\t$aa\t$ab\t$bb";
    return ($res, $dist);
}
