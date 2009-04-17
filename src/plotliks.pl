#!/usr/bin/perl -w
#

@matfn = </tmp/matlik*.txt>;
@rcfn  = </tmp/rcod*.txt>;

make_gp('/tmp/all_matlik.txt', @matfn);
make_gp('/tmp/all_rclik.txt', @rcfn);

print qq{
	set xrange [0:50]
  plot '/tmp/all_matlik.txt' w l, '/tmp/all_rclik.txt' w l
};

sub make_gp{
	my $gpfn = shift;
	my @fn   = @_;
	open FH, ">$gpfn" or die $!;
	foreach my $f (@fn){
		my $s = reformat_liks($f,300);
		print FH $s;
		print FH "\n";
	}
}

sub reformat_liks{
	my $fn     = shift;
	my $minlen = shift;
	open IN, "<$fn" or die $!;
	my $s = "";
	my $i = 0;
	my $val;
	while(<IN>){
		$val = $_;
		chomp($val);
		$val =~ s/^\s*//g;
		$s .= "$i $val\n";
		$i++;
	}
	while($i<$minlen){
		$s .= "$i $val\n";
		$i++;
	}
	$s
}
