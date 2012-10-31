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

my $boost2ped = 0;
my $illumina2tped = 0;
my $ped2tped = 0;
my $illumina = "";
my $tped = "";
my $tfam = "";
my $boost = "";
my $ped = "";
my $map = "";

my $result = GetOptions( "B2P" => \$boost2ped,
                         "I2T" => \$illumina2tped,
                         "P2T" => \$ped2tped,
                         "illu|i=s" => \$illumina,
                         "tped|t=s" => \$tped,
                         "tfam|f=s" => \$tfam,
                         "boost|b=s" => \$boost,
                         "ped|p=s" => \$ped,
                         "map|m=s" => \$map);
                         

if( $boost2ped ) {
    if( $boost =~ m/^$/ ) {
        print "Expected BOOST command line option to be set\n";
        exit;
    } elsif ( $ped !~ m/\.ped$/i ) {
        print "Expected PED command line option to be set\n";
        exit;
    } elsif ( $map !~ m/\.map$/i ) {
        print "Expected MAP command line option to be set\n";
        exit;
    }

    BOOST2PED( $boost, $ped, $map );
}

if( $illumina2tped ) {
    if( $illumina =~ m/^$/ ) {
        print "Expected ILLUMINA command line option to be set\n";
        exit;
    } elsif ( $tped !~ m/\.tped$/i ) {
        print "Expected TPED command line option to be set\n";
        exit;
    } elsif ( $tfam !~ m/\.tfam$/i ) {
        print "Expected TFAM command line option to be set\n";
        exit;
    }

    Illumina2TPED( $illumina, $tped, $tfam );
}

if( $ped2tped ) {
    if( $ped !~ m/\.ped$/i ) {
        print "Expected PED command line option to be set\n";
        exit;
    } elsif( $map !~ m/\.map/i ) {
        print "Expected MAP command line option to be set\n";
        exit;
    } elsif ( $tped !~ m/\.tped/i ) {
        print "Expected TPED command line option to be set\n";
        exit;
    } elsif ( $tfam !~ m/\.tfam/i ) {
        print "Expected TFAM command line option to be set\n";
        exit;
    }

    TransposePED( @ARGV );
    print "Conversion completed\n";
    print "Verifying conversion\n";
    verifyPED2TPED( @ARGV );
}

print "Done\n";

sub BOOST2PED {
    open( BOOST, "<$_[0]") or die $!;
    open( PED, ">$_[1]") or die $!;
    open( MAP, ">$_[2]") or die $!;

    my $index = 0;
    my $i = 0;
    my $sub;
    my $snp_count = 0;
    while( my $line = <BOOST> ) {
        chomp( $line );

        print PED "Fam_$index Ind_0 Pat_$index Mat_$index x ".substr($line, 0, 1)." ";

        $sub = substr( $line, 2);
        my @j = map { ($_ == '0' ? ( "A A") : (( $_ == '1' ) ? ("A C") : ("C C"))) } split( /\s/, $sub);

        print PED join( " ", @j)."\n";

        $index++;
        #$i = length( $sub );
        if( ($#j + 1) > $snp_count) {
            $snp_count = $#j + 1;
        }
    }
    for( $i = 0; $i < $snp_count; $i++) {
        print MAP "0\trs$i\t0\t$i\n";
    }
    print "Found $index Individuals\n";
    print "Found $snp_count SNPs\n";

    close( MAP );
    close( PED );
    close( BOOST );
}

sub Illumina2TPED {
    open( ILLUMINA, "<$_[0]" ) or die $!;
    open( TPED, ">$_[1]" ) or die $!;
    open( TFAM, ">$_[2]" ) or die $!;

    my $line;

    # parse Illumina header line containing all sample IDs
    $line = <ILLUMINA>;
    chomp($line);
    my @sp = split(/\s/, $line);
    for( my $i = 0; $i < 4; $i++ ) {
        shift( @sp );
    }

    my $index = 0;
    while( $#sp >= 0 ) {
        print TFAM "$index\t$sp[0]\tFat_$index\tMot_$index\tx\t1\n";
        shift( @sp );
        $index++;
    }

    while( $line = <ILLUMINA> ) {
        chomp($line);
        @sp = split(/\s/, $line );

        print TPED "$sp[1]\t$sp[0]\t0\t0\t";

        # skip the header columns
        for( my $i = 0; $i < 4; $i++) {
            shift( @sp );
        }

        print TPED join("\t", split(//, join( '', @sp)))."\n";
    }

    close TFAM;
    close TPED;
    close ILLUMINA;
} 

sub TransposePED {
    my $tmp_file_idx = 0;

    open PED, "<$_[0]" or die $!;
    open MAP, "<$_[1]" or die $!;

    open TPED, ">$_[2].$tmp_file_idx" or die $!;
    open TFAM, ">$_[3]" or die $!;

    # Copy MAP to TPED to serve as a template and count Individuals
    my $snp_count = 0;
    my $ind_count = 0;
    my $line;
    while( $line = <MAP> ) {
        print TPED $line;
        $snp_count++;
    }

    print "Found $snp_count Markers\n";

    close MAP;
    close TPED;

    # Determine SNP count
    $line = <PED>;
    chomp($line);
    my @sp = split( /\s/, $line);

    my @buffer; # a queue buffer used for storing intermediate marker values

    my $_count = $#sp + 1;
    $_count = $_count - 6;    # first 6 columns of a PED file are Family INFO
    $_count = $_count / 2;    # b/c bi-allelic

    if( $_count != $snp_count ) {
        print "There is an error in the number of markers in the MAP and PED files\n";
        exit;
    }

    # Perform transpose on blocks of rows in PED file

    my $row_max = 0;    # Transpose 10 lines at a time 
    my $row_count = 0;
    my $i;
    my $tmp_buffer = "";

    do {
        chomp ($line);
        @sp = split( /\s/, $line );

        print TFAM "$sp[0]";
        shift( @sp );
        for( $i = 0; $i < 5; $i++ ) {
            print TFAM "\t$sp[0]";
            shift( @sp );
        }
        print TFAM "\n";

        $i = 0;
        if( $row_count == 0 ) {
            # seed an empty queue
            while( $#sp >= 0 ) {
                $tmp_buffer = " $sp[0] $sp[1]";
                shift( @sp );
                push (@buffer, $tmp_buffer );
            }
        } else {
            # if constructing multiple rows at the same time
            # pop the top element off the top of the queue
            # append the next value
            # push the value onto the end of the queue
            while( $#sp >= 0 ) {
                $tmp_buffer = $buffer[0];
                shift( @buffer );
                $tmp_buffer = $tmp_buffer." $sp[0] $sp[1]";
                shift( @sp );
                push (@buffer, $tmp_buffer );
            }
        }

        $row_count++;
        $ind_count++;
        if( $row_count >= $row_max ) {
            my $t_file = "$ARGV[2].$tmp_file_idx";
            $tmp_file_idx++;
            my $m_file = "$ARGV[2].$tmp_file_idx";
            mergeBufferAndFile( $t_file, $m_file, \@buffer );
            unlink $t_file;
            @buffer = ();
            $row_count = 0;
        }
    } while( $line = <PED> );

    my $t_file = "$ARGV[2].$tmp_file_idx";
    my $m_file = "$ARGV[2]";
    mergeBufferAndFile( $t_file, $m_file, \@buffer );
    unlink $t_file;

    print "Found $ind_count Individuals\n";

    close PED;
    close TFAM;
}

# merges a Template File with buffered data
# output is saved in a merge_file
sub mergeBufferAndFile {
    my $line;
    my ($temp_file, $merge_file, $buf) = @_;
    my @buff = @$buf;
    
    open TPED, "<$temp_file" or die $!;
    open TPED2, ">$merge_file" or die $!;
    
    while( $line = <TPED> ) {
        chomp($line);

        if( $#buff >= 0 ) {
            $line = $line.$buff[0];
            shift( @buff );
        }
        shift( @buff );
        print TPED2 "$line\n"; 
    }

    close TPED;
    close TPED2;
}

sub verifyPED2TPED {
    my ($ped_file, $map_file, $tped_file, $tfam_file) = @_;
    my $line;
    my $line2;
    my @buff = ();
    my @buff2 = ();
    my @sp;
    my @sp2;
    my @sp3;
    my $idx = 0;

    open PED, "<$ped_file" or die $!;   # row header of TFAM + column data of TPED
    open MAP, "<$map_file" or die $!;   # row header of TPED
    open TPED, "<$tped_file" or die $!;
    open TFAM, "<$tfam_file" or die $!;

    do {
        $line = <PED>;
        chomp( $line );
        #print substr($line, 0, 40)."\n";
        push @buff, [ split( /\s/, $line ) ];
        
        $line = <TPED>;
        chomp( $line );
        #print substr($line, 0, 40)."\n";
        push @buff2, [ split( /\s/, $line ) ];

        $line2 = <TFAM>;
        chomp($line2);
        #print substr($line2, 0, 20)."\n";
        @sp = split( /\s/, $line2 );

        # compare TFAM with PED
        while( $#sp >= 0 ) {
            if( $sp[0] ne $buff[$idx][0] ) {
                print "ERROR: TFAM does not match PED\n";
                print $sp[0]."\n";
                print $buff[$idx][0]."\n";
                print "$idx\n";
                exit;
            }
            shift( @sp );
            shift( @{ $buff[$idx] } );
        }

        $line2 = <MAP>;
        chomp ($line2);
        @sp = split(/\s/, $line2);

        while( $#sp >= 0 ) {
            if( $sp[0] ne $buff2[$idx][0] ) {
                print "ERROR: MAP does not match TPED\n";
                exit;
            }

            shift( @sp );
            shift( @{$buff2[$idx]} );
        }

        for( my $i = 0; $i < $idx; $i++ ) {
            for( my $j = 0; $j < 2; $j++ ) {
                #print  "[$idx][0] $buff[$idx][0] -> [$i][0] $buff2[$i][0]\n";
                if( $buff[$idx][0] ne $buff2[$i][0] ) {
                    print "ne\n";
                    exit;
                }
                shift( @{$buff[$idx]} );
                shift( @{$buff2[$i]} );
                #print  "[$i][0] $buff[$i][0] -> [$idx][0] $buff2[$idx][0]\n";
                if( $buff[$i][0] ne $buff2[$idx][0] ) {
                    print "ne\n";
                    exit;
                }

                shift( @{$buff[$i]} );
                shift( @{$buff2[$idx]} );
            }
        }

        #print  "[$idx][0] $buff[$idx][0] -> [$idx][0] $buff2[$idx][0]\n";
        if( $buff[$idx][0] ne $buff2[$idx][0] ) {
            print "ERROR: PED does not match TPED\n";
            exit;
        }
        shift( @{$buff[$idx]} );
        shift( @{$buff2[$idx]} );

        #print  "[$idx][0] $buff[$idx][0] -> [$idx][0] $buff2[$idx][0]\n";
        if( $buff[$idx][0] ne $buff2[$idx][0] ) {
            print "ERROR: PED does not match TPED\n";
            exit;
        }
        shift( @{$buff[$idx]} );
        shift( @{$buff2[$idx]} );
        $idx++; 
    } while (!eof( PED ) );

    close PED;
    close MAP;
    close TPED;
    close TFAM;
    print "Transposition verified\n";
}
