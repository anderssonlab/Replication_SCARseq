#!/usr/bin/perl

## usage:
## ./RFD.pl F.tab R.tab radius dradius zradius

use strict;
use List::Util qw(sum);

open(my $fin ,"$ARGV[0]") or die $!;
open(my $rin, "$ARGV[1]") or die $!;

my $radius = $ARGV[2];
my $dradius = $ARGV[3];
my $zradius = $ARGV[4];

if ($dradius > $radius)
{
	die $!;
}

if ($zradius > $radius)
{
	die $!;
}

## Keep a stack of RFD values, smoothed RFD values and window IDs for derivative calculation
my @RFD_stack = ();
my @RFD_smooth_stack = ();
my @ID_stack = ();
my $size = 0;
my $size_smooth = 0;

## initiate derivative weights as in Nonparametric Derivative Estimation (Brabanter et al.)
my @w = ();
for (my $i=0; $i<=$dradius; $i++)
{
	push(@w, $i**2);
}
my $wsum = sum(@w);
for (my $i=1; $i<=$dradius; $i++)
{
	$w[$i] /= $wsum;
}

while(my $fline = <$fin>)
{
	my $rline = <$rin>;

	my @F = split(/\t/,$fline);
	my @R = split(/\t/,$rline);

	## calculate RFD
	my $RFD = ($F[3] - $R[3])/($R[3] + $F[3]);

	## Push new item values to the stacks
	push(@RFD_stack, $RFD);
	push(@ID_stack, $F[0]);
	$size++;

	if ($size < ($radius*2+1))
	{
		if ($size < ($radius+1))
		{
			print "$ID_stack[$size-1]\t$RFD_stack[$size-1]\tNA\tNA\tNA\tNA\n";
		}
	 	next;
	}

	## calculate smoothed RFD (uniform blur)
	my $RFD_smooth = 0;
	for (my $i=0; $i<($radius*2+1); $i++)
	{
		$RFD_smooth += $RFD_stack[$i];
	}
	$RFD_smooth /= $radius*2+1;
	
	# Push new value to the stack
	push(@RFD_smooth_stack, $RFD_smooth);
	$size_smooth++;

	if ($size_smooth < ($radius*2+1))
	{
		if ($size_smooth < ($radius+1))
		{
			print "$ID_stack[$radius+$size_smooth-1]\t$RFD_stack[$radius+$size_smooth-1]\t$RFD_smooth_stack[$size_smooth-1]\tNA\tNA\tNA\n";
		}
	}
	else
	{
		## calculate RFD derivative
		my $RFD_deriv = 0;
		for (my $i=1; $i<=$dradius; $i++)
		{
			$RFD_deriv += $w[$i] * ($RFD_smooth_stack[$radius+$i]-$RFD_smooth_stack[$radius-$i]) / (2*$i+1);
		}
	
		## caclculate RFD boundary score
		my $score = $RFD_deriv*(1.0-abs($RFD_smooth_stack[$radius]));
	
		## restrict derivative reporting at zero crossings (zradius) 
		my $zero_crossing_derivative = "NA";
		my $min_RFD=1;
		my $max_RFD=-1;
		for (my $i=-1*$zradius; $i<=$zradius; $i++)
		{
			if ($RFD_smooth_stack[$radius+$i]>$max_RFD)
			{
				$max_RFD=$RFD_smooth_stack[$radius+$i];
			}
			if ($RFD_smooth_stack[$radius+$i]<$min_RFD)
			{
				$min_RFD=$RFD_smooth_stack[$radius+$i];
			}
		}
		if ($min_RFD<0 && $max_RFD>0)
		{
			$zero_crossing_derivative = $RFD_deriv;
		}
		
		## print
		print "$ID_stack[0]\t$RFD_stack[0]\t$RFD_smooth_stack[$radius]\t$RFD_deriv\t$score\t$zero_crossing_derivative\n";
		
		shift @RFD_smooth_stack;
		$size_smooth--;
	}

	## Remove first items from the stack
	shift @RFD_stack;
	shift @ID_stack;
	$size--;
}

## print remaining elements
for (my $i=$radius; $i<(2*$radius); $i++)
{
	print "$ID_stack[$i-$radius]\t$RFD_stack[$i-$radius]\t$RFD_smooth_stack[$i]\tNA\tNA\tNA\n";
}
for (my $i=$radius; $i<(2*$radius); $i++)
{
	print "$ID_stack[$i]\t$RFD_stack[$i]\tNA\tNA\tNA\tNA\n";
}

close($fin);
close($rin);
