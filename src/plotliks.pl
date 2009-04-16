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
		my $s = reformat_liks($f);
		print FH $s;
		print FH "\n";
	}
}

sub reformat_liks{
	my $fn = shift;
	open IN, "<$fn" or die $!;
	my $s = "";
	my $i = 0;
	while(<IN>){
		chomp;
		s/^\s*//g;
		#$_ = "-$_" unless(s/^-//); # invert sign
		$s .= "$i $_\n";
		$i++;
	}
	$s
}
