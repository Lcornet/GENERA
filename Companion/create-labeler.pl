#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;

use Smart::Comments '###';

use Getopt::Euclid qw(:vars);

use File::Slurp;
use List::AllUtils qw(uniq);

use Bio::MUST::Core;
use aliased 'Bio::MUST::Core::Taxonomy';

# build taxonomy object
my $tax = Taxonomy->new_from_cache(tax_dir => $ARGV_taxdir);

# list of genus taxids
my @taxids = read_file($ARGV_infile, chomp => 1);

# get superkingdom rank for each taxids
my %genera_for;
$genera_for{$_} = join q{}, 
    $tax->get_taxa_from_taxid($_, 'superkingdom') for @taxids;

# filtering according to wanted superkingdom
my @wanted_taxids;
my @kingdoms = @ARGV_kingdoms ? @ARGV_kingdoms : qw(Bacteria Archaea);
for my $kingdom (@kingdoms) { 
    push @wanted_taxids, map { $genera_for{$_} eq $kingdom ? $_ : () }
        keys %genera_for;
}

# get all taxa of specified rank
my @labels = grep {$_ ne 'undef'} 
    uniq map { $tax->get_taxa_from_taxid($_, $ARGV_level) } @wanted_taxids;

say join "\n", @labels;

__END__

=pod

=head1 NAME

create-labeler.pl - Create a list of taxon at specific level

=head1 VERSION

version 1.O

=head1 USAGE

    create-labeler.pl <infile> --taxdir=<dir> [optional arguments]

=head1 REQUIRED ARGUMENTS

=over

=item <infile>

List of genera taxon ids.

=for Euclid: infile.type: readable

=item --taxdir=<dir>

Path to local mirror of the NCBI Taxonomy database.

=for Euclid: dir.type: string

=back

=head1 OPTIONS

=over

=item --level=<str>

Taxonomic level to build labeler [default: str.default].

Valid levels are: superkingdom, kingdom, subkingdom, superphylum, phylum,
subphylum, superclass, class, subclass, infraclass, superorder, order,
suborder, infraorder, parvorder, superfamily, family, subfamily, tribe,
subtribe, genus.

=for Euclid: str.type:    string
    str.default: 'phylum'

=item --kingdoms=<str>...

List of whitespace-separated superkingdom to be displayed in labeler
[default: Bacteria Archaea].
