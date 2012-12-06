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

my $contin1_file = "";
my $contin2_file = "";

my $result = GetOptions( "A=s" => \$contin1_file,
                         "B=s" => \$contin2_file );

if( $contin1_file =~ m/^$/ ) {
    print "Expected a Contin1 file. Please specify C1 parameter\n";
    exit;
}

if( $contin2_file =~ m/^$/ ) {
    print "Expected a Contin2 file. Please specify C2 parameter\n";
    exit;
}

open( CON1, "<$contin1_file") or die $!;
open( CON2, "<$contin2_file") or die $!;

my $line;
my $line2;

my $table1;
my $table2;

my @case_table1;
my @control_table1;
my @case_table2;
my @control_table2;

while( !eof(CON1) and !eof(CON2) ) {
    $table1 = <CON1>;
    $table2 = <CON2>;
    chomp($table1);
    chomp($table2);

    $line = <CON1>; # Cases
    $line = <CON1>; # col-header
    $line2 = <CON2>; # Cases
    $line2 = <CON2>; # col-header

    for( my $i = 0; $i < 3; $i++ ) {
        $line = <CON1>;
        chomp($line);
        my @sp = split(/\t/, $line );
        for( my $j = 1; $j < 4; $j++ ) {
            $case_table1[$i * 3 + ($j - 1)] = $sp[$j];
        }
    }
    $line = <CON1>; # Controls
    $line = <CON1>; # col-headers
    for( my $i = 0; $i < 3; $i++ ) {
        $line = <CON1>;
        chomp($line);
        my @sp = split(/\t/, $line );
        for( my $j = 1; $j < 4; $j++ ) {
            $control_table1[$i * 3 + ($j - 1)] = $sp[$j];
        }
    }

    for( my $i = 0; $i < 3; $i++ ) {
        $line = <CON2>;
        chomp($line);
        my @sp = split(/\t/, $line );
        for( my $j = 1; $j < 4; $j++ ) {
            $case_table2[$i * 3 + ($j - 1)] = $sp[$j];
        }
    }
    $line = <CON2>; # Controls
    $line = <CON2>; # col-headers
    for( my $i = 0; $i < 3; $i++ ) {
        $line = <CON2>;
        chomp($line);
        my @sp = split(/\t/, $line );
        for( my $j = 1; $j < 4; $j++ ) {
            $control_table2[$i * 3 + ($j - 1)] = $sp[$j];
        }
    }

    if( $table1 !~ m/$table2/ ) {
        print "Unexpected table comparison: $table1 vs $table2\n";
        last;
    } elsif( $case_table1[ 4 ] !~ m/$case_table2[4]/ ) {
        print "$table1 Case Tables do not match.\n";
        last;
    } elsif( $control_table1[ 4 ] !~ m/$control_table2[4]/ ) {
        print "$table2 Control Tables do not match.\n";
        last;
    } else {
        # we know at this point that both case and control tables match in the AB_AB element
        # need to verify that the second table is a permutation of the first table
        #
        my $t1 = join(',', @case_table1 );
        my $t2 = join(',', @case_table2 );

        if( $t1 !~ m/$t2/ ) {
            $t2 = join(',', swap_columns( @case_table2 ));

            if( $t1 !~m/$t2/ ) {
                $t2 = join(',', swap_rows( @case_table2 ));
                if( $t1 !~ m/$t2/ ) {
                    $t2 = join( ',', swap_columns( swap_rows( @case_table2)));

                    if( $t1 !~ m/$t2/ ) {
                        print "$table1 Case tables do not match: $t1 v ".join(',', @case_table2)."\n";
                        last;
                    }
                }
            }
        }

        my $t1 = join(',', @control_table1 );
        my $t2 = join(',', @control_table2 );

        if( $t1 !~ m/$t2/ ) {
            $t2 = join(',', swap_columns( @control_table2 ));

            if( $t1 !~m/$t2/ ) {
                $t2 = join(',', swap_rows( @control_table2 ));
                if( $t1 !~ m/$t2/ ) {
                    $t2 = join( ',', swap_columns( swap_rows( @control_table2)));

                    if( $t1 !~ m/$t2/ ) {
                        print "$table1 Control tables do not match: $t1 v ".join(',', @control_table2)."\n";
                        last;
                    }
                }
            }
        }
    }
}

close( CON1 );
close( CON2 );

print "Done!\n";

sub rotate {
    my @t = @_;
    my @out;
    
    $out[ 4 ] = $t[4];  # center stays the same

    # row 1 becomes column 3
    $out[ 2 ] = $t[0];
    $out[ 5 ] = $t[1];
    $out[ 8 ] = $t[2];
    
    # row 3 becomes column 1
    $out[ 6 ] = $t[ 8 ];
    $out[ 3 ] = $t[ 7 ];
    $out[ 0 ] = $t[ 6 ];

    $out[ 1 ] = $t[ 3 ];
    $out[ 7 ] = $t[ 5 ];

    return @out;
}

sub swap_columns {
    my @t = @_;
    my @out;

    $out[ 1 ] = $t[1];
    $out[ 4 ] = $t[4];
    $out[ 7 ] = $t[7];

    $out[ 0 ] = $t[ 2 ];
    $out[ 3 ] = $t[ 5 ];
    $out[ 6 ] = $t[ 8 ];

    $out[ 2 ] = $t[ 0 ];
    $out[ 5 ] = $t[ 3 ];
    $out[ 8 ] = $t[ 6 ];

    return @out;
}

sub swap_rows {
    my @t = @_;
    my @out;

    $out[ 0 ] = $t[ 6 ];
    $out[ 1 ] = $t[ 7 ];
    $out[ 2 ] = $t[ 8 ];

    $out[ 3 ] = $t[ 3 ];
    $out[ 4 ] = $t[ 4 ];
    $out[ 5 ] = $t[ 5 ];

    $out[ 6 ] = $t[ 0 ];
    $out[ 7 ] = $t[ 1 ];
    $out[ 8 ] = $t[ 2 ];

    return @out;
}
