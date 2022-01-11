use strict; use warnings;

print "gtf_2_anno_tab.gene_level.pl <xxx.gtf>\n";
print "Out: anno.txt\n";
print "good for Equus_caballus.EquCab3.0.104.gtf\n";

die unless @ARGV == 1;

open(GTF, shift) or die;

my%id2name;
while(<GTF>){
    chomp;
    my@e = split "\t", $_;
    next if /^#/;
    next if $e[2] ne 'gene';
    next if $_ !~ m/gene_id/;

    my($id) = $_ =~ /gene_id ([^;]+);/;
    my($name);

    if ($_ !~ m/gene_name/){
        $name = $id
    }else{
        ($name) = $_ =~ /gene_name ([^;]+);/;
    }

    my($type) = $_ =~ /gene_biotype ([^;]+);/;  # change: now for mouse gencode mm10

    if (length($type) == 0){$type = "-"}

    $name =~ s/^\s+//;
    $id =~ s/^\s+//;
    $type =~ s/^\s+//;

    $id2name{$id} = "$name\t$type";
}
close GTF;

open(OUT, ">anno.txt");
print OUT "Gene\tName\tType\n";
for (sort keys %id2name){
    my$out = "$_\t$id2name{$_}\n";
    $out =~ s/\"//g;
    print OUT $out;
}
close OUT;
