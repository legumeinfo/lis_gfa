#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $infer_gene = 1;
my $gene_from_query_regexp;
my $evalue_threshold="1e-10";
my $transcript2genelookupfile;
my %gene2bestisoform;
my $restrict_to_best_isoform=0;
my $ignore_missing_genes=0;
GetOptions(
    "infer_gene!" => \$infer_gene,
    "gene_from_query_regexp=s", \$gene_from_query_regexp,
    "evalue_threshold=f", \$evalue_threshold,
    "transcript2genelookup=s" => \$transcript2genelookupfile,
    "restrict_to_best_isoform" => \$restrict_to_best_isoform,
    "ignore_missing_genes" => \$ignore_missing_genes,
);
#if you need a transcript2gene lookup, something like this should suffice
#zgrep '  mRNA    ' medtr.HM004.gnm1.ann1.2XTB.gene_models_main.gff3.gz | sed 's/.*ID=\([^;]*\).*Parent=\([^;]*\).*/\1\t\2/' > transcript2gene
my %transcript2genelookup;
if (defined $transcript2genelookupfile) {
	open(F,$transcript2genelookupfile) || die $!;
	%transcript2genelookup = map {chomp; my ($t, $g) = split /\t/; $t => $g;} <F>;
	close F;
	$infer_gene=0;
}
print "#gene\tfamily\tprotein\te-value\tscore\tbest_domain_score\tpct_hmm_coverage\n";
my %pep2family;
my $family;
my $model_length;
while (<>) {
    if (/^Query:\s+(\S+)\s+\[M=(\d+)\]/) {
        $family = $1;
        $model_length = $2;
    }
    elsif (/^>>\s*(\S+)(.*)/) {
        my $pep = $1;
        next unless defined $pep2family{$pep}->{family};
        #FIXME: probably no good reason to do the gene determination here, if we are no longer using desc
        my $desc = $2;
        my $gene;
        if (defined $transcript2genelookupfile) {
            $gene = $transcript2genelookup{$pep};
            if (!defined $gene) {
                my $pep2 = $pep;
                $pep2 =~ s/utrorf$//;
                $gene = $transcript2genelookup{$pep2};
                if (!defined $gene && !$ignore_missing_genes) {
                    die "$pep had no gene in $transcript2genelookupfile\n";
                }
            }
        }
        #else {
            ##FIXME: still useful?
            #the gffread phytozome ensembl conventions for fasta headers
            #(undef, $gene) = ($desc =~ /(gene|locus)[=:](\S+)/);
        #}
        if (defined $gene) {
            $pep2family{$pep}->{gene} = $gene;
        }
        #iff we have been assigned to this family, process alignment info
        if ($pep2family{$pep}->{family} eq $family) {
            <>;
            <>;
            my $hmm_coverage=0;
            while (<>) {
               last if /^$/;
               my @data = split /\s+/;
               my $score = $data[3];
               #per http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf this indicates that a hit met the inclusion thresholds
               #sometimes hits from the protein get split into separate domains, so we sum across them to get the full hmm coverage
               #(we're not checking for possible overlaps- so this could be somewhat inflated, I suppose, but probably should reflect
               #similar phenomena to cases where sequence score significantly exceeds single best domain score due to tandem gene fusions,
               #in which case we might get >100% coverage)
               if ($data[2] eq "!") {
                  $hmm_coverage += $data[8]-$data[7]+1;
                  $pep2family{$pep}->{pct_hmm_coverage} = 100*$hmm_coverage/$model_length;
               }
               #original logic; just look at coverage for the single best domain
               #if ($pep2family{$pep}->{best_domain_score} == $score) {
                  #$pep2family{$pep}->{pct_hmm_coverage} = 100*($data[8]-$data[7]+1)/$model_length;
               #}
            }
        }
        next;
    }
    elsif (/^\s+E-value/) {
        <>;
        while (<>) {
            last if /^$/;
            next if /^\s*-*\s*inclusion threshold/;
            s/^\s+//;
            my @data = split /\s+/;
            my $evalue = $data[0];
            my $score = $data[1];
            my $best_domain_score = $data[4];
            #worked fine, until we hit descriptions that had more than just the gene id
            #my $gene = $data[$#data];
            my $pep = $data[8];
            my $gene;
	    if ($infer_gene) {
		    if (defined $gene_from_query_regexp) {
			($gene) = ($pep =~ /$gene_from_query_regexp/);
		    }
		    else {
			($gene) = ($pep =~ /(.*)\.[tm]?\d+$/);
		    }
	    }
	    elsif (defined $transcript2genelookupfile) {
		    $gene = $transcript2genelookup{$pep};
	    }
	    else {
		    $gene = $pep2family{$pep}->{gene};
	    }
            if ($evalue <= $evalue_threshold) {
                if (! defined $pep2family{$pep} || ! defined $pep2family{$pep}->{score} || $score > $pep2family{$pep}->{score} || ($score == $pep2family{$pep}->{score} && $family cmp $pep2family{$pep}->{family} < 0)) {
                    $pep2family{$pep}->{evalue} = $evalue;
                    $pep2family{$pep}->{score} = $score;
                    $pep2family{$pep}->{best_domain_score} = $best_domain_score;
                    $pep2family{$pep}->{family} = $family;
                    if (defined $gene) {
                        $pep2family{$pep}->{gene} = $gene;
			if (! defined $gene2bestisoform{$gene} || $score > $gene2bestisoform{$gene}->{score}) {
				$gene2bestisoform{$gene}->{isoform} = $pep;
				$gene2bestisoform{$gene}->{score} = $score;
				$gene2bestisoform{$gene}->{best_domain_score} = $best_domain_score;
			}
                    }
                }
            }
        }
    }

}

foreach my $pep (keys %pep2family) {
    if (defined $pep2family{$pep}->{family}) {
	if (!$restrict_to_best_isoform || $gene2bestisoform{$pep2family{$pep}->{gene}}->{isoform} eq $pep) {
        if (defined $pep2family{$pep}->{gene}) {
            print $pep2family{$pep}->{gene},"\t",$pep2family{$pep}->{family},"\t", $pep,"\t", $pep2family{$pep}->{evalue}, "\t", $pep2family{$pep}->{score},"\t",$pep2family{$pep}->{best_domain_score},"\t",$pep2family{$pep}->{pct_hmm_coverage},"\n";
        }
}
    }
}
