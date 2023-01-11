#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;

use Smart::Comments '###';

use Getopt::Euclid qw(:vars);

use File::Basename;
use Path::Class 'file';
use List::AllUtils qw(sum count_by);

use Bio::MUST::Core;
use aliased 'Bio::MUST::Core::Taxonomy';

# build taxonomy object
my $tax = Taxonomy->new_from_cache(tax_dir => $ARGV_taxdir);

# TODO: fix this (make it optional?)
# build labeler
my $labeler = $tax->tax_labeler_from_list($ARGV_taxon_list);

my $method;
$method = $ARGV_auto_detect if $ARGV_auto_detect;
my @results;
for my $infile (@ARGV_infiles) {

    my ($basename) = fileparse($infile, qr{\.[^.]*}xms);
    my ($seq_n, $lca_for) = read_kraken($infile);

    # determine expected taxon from kraken infile
    my $exp_tax = auto_detect($lca_for, $method);

    format_results($lca_for, $basename, $exp_tax, $seq_n);
}

# store results
open my $out, '>', $ARGV_outfile;
say {$out} join "\n", @results;

# functions

sub read_kraken {
    my $infile = shift;
    
    open my $in, '<', $infile;
    
    my $seq_n;
    my %tree;
    my %lca_for;
    my $n_start = 0;
    my $prev_n_blank = 0;
    #~ my %test_for;
    LINE:
    while (my $line = <$in>) {
        chomp $line;
    
        # split line on tabulation
        my ($fract, $cov, $count, $rank_code, $taxon_id, $taxon) = split /\t/, $line;
    
        # compute number of total sequences
        if ($rank_code =~ m/\bU\b|\bR\b/xms) {
            $seq_n += $cov;
    
            # skip unclassified and root lines
            next LINE;
        }
    
        # catch blank in front of taxon...
        my ($blank) = ($taxon) =~ m/^(\s+)/xms;
        # and count it
        my $curr_n_blank = (split q{}, $blank) / 2;
        # remove blank from taxon
        $taxon =~ s/^\s+//;
    
    #    if ($curr_n_blank > $prev_n_blank) {
        unless ($curr_n_blank <= $prev_n_blank) {
            $tree{$curr_n_blank} = $taxon;
    
            if ($count > 0) {
                my @taxa = map { $tree{$_} } 1..$curr_n_blank;
    #~            $test_for{ join ';', @taxa } = $count;
                my $n_end = ( $n_start + $count ) - 1;
    
                $lca_for{"seq$_"} = {
                    taxon => [@taxa],
                    rel_n => 0,
                } for $n_start..$n_end;
    
                $n_start = $n_end + 1;
            }
    
            $prev_n_blank = $curr_n_blank;
            next LINE;
        }
    
        if ($curr_n_blank <= $prev_n_blank) {
    #    else {
            delete $tree{$_} 
                for $prev_n_blank..$curr_n_blank;
    
            $tree{$curr_n_blank} = $taxon;
    
            if ($count > 0) {
                my @taxa = map { $tree{$_} } 1..$curr_n_blank;
    #~            $test_for{ join ';', @taxa } = $count;
                my $n_end = ( $n_start + $count )-1;
    
                $lca_for{"seq$_"} = {
                    taxon => [@taxa],
                    rel_n => 0,
                } for $n_start..$n_end;
    
          	     $n_start = $n_end + 1;
            }
    
            $prev_n_blank = $curr_n_blank;
            next LINE;
        }
    }

    return ($seq_n, \%lca_for);
}

sub auto_detect {
    my $lca_for = shift;
    my $method  = shift
        // 'count_first';

    # TODO: implement and compare the opposite approach where we first label
    # then determine the exp_tax...

    my $exp_tax;
    if ($method eq 'count_first') {
        # count LCA lineage occurrences (i.e., sort | uniq -c)
        my %count_for = count_by { join q{;}, @{ $_->{taxon} } } values %{$lca_for};

        # compute LCA proportions
        my @taxa = sort {
            $count_for{$b} <=> $count_for{$a}
        } keys %count_for;

        $exp_tax = $labeler->classify(
            $taxa[0], { greedy => $ARGV_greedy_taxa }
        );
    }

    if ($method eq 'label_first') {
        # label LCA and count occurences
        my %count_for = count_by { 
            $labeler->classify( $_->{taxon},
            { greedy => $ARGV_greedy_taxa } ) // q{}
        } values %{$lca_for};

        my @taxa = sort {
            $count_for{$b} <=> $count_for{$a}
        } keys %count_for;

        $exp_tax = $taxa[0];
    }

    return $exp_tax;
}

sub tabulate_lcas {
    my $lca_for = shift;
    my $exp_tax = shift;
    my $seq_n   = shift;

    my $self_n    = 0;
    my $foreign_n = 0;
    my $unknown_n = 0;

    my %count_for;

    # tabulate LCAs
    while ( my ($query_id, $lca) = each %{$lca_for} ) {

        # label LCA
        my $got_tax = $labeler->classify(
            $lca->{taxon}, { greedy => $ARGV_greedy_taxa }
        );

        # update counts based on LCA label
        unless ($got_tax)                  { $unknown_n++ }
        elsif  ($got_tax eq $exp_tax)      {    $self_n++ }
        else                               { $foreign_n++;
                                             $count_for{$got_tax}++;
        }
    }

    # compute class proportions
    # TODO: decide on a data structure to include all this
    my $self_p    = sprintf "%.2f", 100 * $self_n    / $seq_n;
    my $foreign_p = sprintf "%.2f", 100 * $foreign_n / $seq_n;
    my $unknown_p = sprintf "%.2f", 100 * $unknown_n / $seq_n;
    my $unclass_p = sprintf "%.2f", 100 - ($self_p + $foreign_p + $unknown_p);

    # compute LCA proportions
    my @taxa = sort {
        $count_for{$b} <=> $count_for{$a} || $a cmp $b
    } keys %count_for;
    my @details = map {
        sprintf "%s=%.2f", $_, 100 * $count_for{$_} / $seq_n
    } @taxa;

    my $result = join "\t", $self_p, $foreign_p, $unknown_p, $unclass_p,
        join q{,}, @details;

    return $result;
}

sub format_results {
    my $lca_for  = shift;
    my $basename = shift;
    my $exp_tax  = shift;
    my $seq_n    = shift;

    # compute mean number of relatives used for LCA inference
    my @rel_counts = map { $_->{rel_n} } values %{$lca_for};
    my $rel_mean
        = @rel_counts ? sprintf "%.2f", sum(@rel_counts) / @rel_counts : 'NA';

    # format default report line
    # TODO: better specify this format
    push @results, join "\t", $basename, $exp_tax,
       tabulate_lcas($lca_for, $exp_tax, $seq_n), $rel_mean;

    return;
}


__END__

=pod

=head1 NAME

kraken-parser.pl - Parser for kraken report

=head1 VERSION

version 0.1

=head1 USAGE

    kraken.pl <infiles> --outfile=<file> --taxdir=<dir> --taxon-list=<file> \
        [optional arguments]

=head1 REQUIRED ARGUMENTS

=over

=item <infiles>

Path to input kraken report file [repeatable argument].

=for Euclid: infiles.type: readable
    repeatable

=item --outfile=<file>

Path to output file.

=for Euclid: file.type: writable

=for Euclid: file.type: readable

=item --taxdir=<dir>

Path to local mirror of the NCBI Taxonomy database.

=for Euclid: dir.type: string

=item --taxon-list=<file>

List of taxa to consider when looking for foreign sequences. This labeler file
is used throughout the program to truncate LCA lineages at specific taxonomic
levels (which can vary from one lineage to the other).

=back

=head1 OPTIONS

=over

=item --exp[ected]-tax[on]=<string>

Organism taxon [default: automatic]. Use this option when the organism does not
have an assembly accession. The specified taxon must be in the file provided
with the C<--taxon-list> option.

=for Euclid: string.type: string

=item --greedy-taxa

Enable greedy behavior when interpreting the ambiguous taxa provided in the
required argument C<--taxon-list> [default: no]

=item --auto-detect=<str>

Determine organism taxon based on kraken infile. Two mode are available:
C<count_first> and C<label_fist>. The first one count LCA lineages occurence
and then take the most abundant as main organism, while the second one
truncate LCA lineages according to the labeler file and take the most
abundant truncate LCA lineage as main organism taxon. [default: str.default].

=for Euclid: str.type: string, str eq 'count_first' || str eq 'label_first'
    str.type.error: <str> must be on count_first or label_first (not str)
    str.default:    "count_first"

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back

=head1 AUTHOR

Denis BAURAIN <denis.baurain@uliege.be>

=head1 CONTRIBUTORS

=for stopwords Valerian LUPO Mick VAN VLIERBERGHE Luc CORNET

=over 4

=item *

Valerian LUPO <valerian.lupo@doct.uliege.be>

=item *

Mick VAN VLIERBERGHE <mvanvlierberghe@doct.uliege.be>

=item *

Luc CORNET <luc.cornet@uliege.be>

=back

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2020 by University of Liege / Unit of Eukaryotic Phylogenomics / Denis BAURAIN.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
