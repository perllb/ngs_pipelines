#!/usr/bin/perl -w
use strict;
use CMD::vcf2;
use Data::Dumper;

my $first = 1;

my $fn = $ARGV[0];
my $vcf = CMD::vcf2->new('file'=>$fn );

while ( my $a = $vcf->next_var() ) {

    # Print header
    if( $first ) {
	print "Chromosome\tPosition\tRef\tAlt";
	print "\t". $_->{_sample_id} foreach @{$a->{GT}};
	print "\n";
	$first = 0;
    }

    print $a->{CHROM}."\t".$a->{POS}."\t".$a->{REF}."\t".$a->{ALT};

    foreach my $gt (@{$a->{GT}}) {
	my $vaf = 0;

	# Calculate variant frequency
	if( $gt->{DP} and $gt->{AO} ne ".") {
	    $vaf = $gt->{AO} / $gt->{DP};
	}
	if( $gt->{DP} > 10 ) {
	    printf( "\t%.1f", # (%d/%d)", 
		    $vaf, 
		    #($gt->{AO} or 0), 
		    #($gt->{DP} or 0)
		);
	}
		
	else {
	    print "\t-";
	}
    }
    print "\n";
}
