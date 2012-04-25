use strict;

my $fileIn = $ARGV[0];
my $fileOut = $ARGV[1];
open(SORT, "gzip -dc $fileIn | perl -pe 's/ \\|\\|\\| /\\0/g' | LC_ALL=C sort -t \\0 -S 8G |");
open(OUT, "| gzip -c >$fileOut");

my @buffer;
my $last = "";
while(<SORT>) {
    chomp;
    my ($s, $t, $f, @rest) = split(/\0/, $_);
    s/\0/ ||| /g;
    
    my @floats = split(/\s/, $f);
    
    if($s ne $last) {
        my $rank = 0;
        foreach my $pair (sort { $b->[1] <=> $a->[1] or $a->[0] cmp $b->[0] } @buffer) {
            print OUT $pair->[0], " ||| ", $rank, "\n";
            $rank++;
        }
        @buffer = ();
    }

    push(@buffer, [$_, $floats[2]]);
    $last = $s;
}
my $rank = 0;
foreach my $pair (sort { $b->[1] <=> $a->[1] or $a->[0] cmp $b->[0] } @buffer) {
    print OUT $pair->[0], " ||| ", $rank, "\n";
    $rank++;
}
