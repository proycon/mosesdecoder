use strict;
use Getopt::Long;

my $SIZE = "2G";
GetOptions(
    "sortmem=s" => \$SIZE
);

my $fileIn = $ARGV[0];
my $fileOut = $ARGV[1];
open(SORT, "gzip -dc $fileIn |");
open(OUT, "| LC_ALL=C sort -t \$'\t' -S $SIZE | gzip -c >$fileOut");

my @buffer;
my $last = "";
while(<SORT>) {
    chomp;
    my ($s, $t, $f, @rest) = split(/\s\|\|\|\s/, $_);
    s/\s\|\|\|\s/\t/g;
    
    my @floats = split(/\s/, $f);
    
    if($s ne $last) {
        my $rank = 0;
        foreach my $pair (sort { $b->[1] <=> $a->[1] or $a->[0] cmp $b->[0] } @buffer) {
            print OUT $pair->[0], "\t", $rank, "\n";
            $rank++;
        }
        @buffer = ();
    }

    push(@buffer, [$_, $floats[2]]);
    $last = $s;
}
my $rank = 0;
foreach my $pair (sort { $b->[1] <=> $a->[1] or $a->[0] cmp $b->[0] } @buffer) {
    print OUT $pair->[0], "\t", $rank, "\n";
    $rank++;
}
