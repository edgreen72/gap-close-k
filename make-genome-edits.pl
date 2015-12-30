#!/usr/bin/perl

use Getopt::Std;
use Bio::SeqIO;
use vars qw( $opt_f $opt_E );

use strict;

my ( $genome_o, ### Bio::SeqIO to fasta genome assembly to edi
     $edits_p,  ### points to ->{ ID => 
                ###               [ [start, end, before_seq, after_seq,
                ###                  K1_occur, K2_occur, spanning_reads], ...
     $edited_genome_o, ### Bio::SeqIO for writing out the edited genome
     $seq_o,
     $seq,
     $old_part,
     $edit_p,
     $total_edited,
     $total_unedited
    );
### Parse command-line options
&init();

$genome_o = Bio::SeqIO->new( -file => $opt_f, -format => 'fasta' );
$edits_p  = &parse_edits( $opt_E );
$edited_genome_o = Bio::SeqIO->new( -format => 'fasta' );

while( $seq_o = $genome_o->next_seq ) {
    if ( defined($edits_p->{ $seq_o->display_id }) ) {
	$seq = $seq_o->seq();
	# Make each edit, from the END to the beginning so coordinates
	# remain accurate
	foreach $edit_p ( sort { $b->[0] <=> $a->[0] } @{ $edits_p->{$seq_o->display_id} } ) {
	    $old_part = substr( $seq, 
				$edit_p->[0], 
				$edit_p->[1]-$edit_p->[0], 
				$edit_p->[3] );
	}
	$seq_o->seq( $seq ); # Assign the edited sequence
	$total_edited++;
    }
    else {
	$total_unedited++;
    }
    $edited_genome_o->write_seq( $seq_o );
}

print STDERR "Total sequences   edited: $total_edited\n";
print STDERR "Total sequences unedited: $total_unedited\n";

sub parse_edits {
    my $fn = shift;
    my ( $l, $id );
    my ( @a );
    my %edits;
    open( EDITS, $fn ) or die( "$!: $fn\n" );
    while( chomp( $l = <EDITS> ) ) {
	@a = split( " ", $l );
	$id = shift( @a );
	push( @{ $edits{$id} }, [@a] );
    }
    close( EDITS );
    return \%edits;
}

sub init {
    getopts( 'f:E:' );
    unless( -f $opt_f &&
	    -f $opt_E ) {
	print( "make-genome-edits.pl -f <fasta genome> -E <edits file>\n" );
	print( "Makes the edits specified in the edits file to the fasta genome\n" );
	print( "specified in the input genome file.\n" );
	exit( 0 );
    }
}
