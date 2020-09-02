use strict; use warnings;

while (<>){
    if (/^#/){
        $_ =~ s/##contig=<ID=/##contig=<ID=chr/;
        print
    }else{
        print("chr" . $_)
    }
}
