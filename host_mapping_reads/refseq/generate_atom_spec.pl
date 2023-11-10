#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $print_usage = 0;
my $_chain = "B2";
my $_ic = "pseudo_knot";

my $usage = <<USAGE;

  This script outputs chimera atom specs corresponding to different types of interactions

  Mark Stenglein,  6/2023   

  Usage: $0 [-h] 

   [-h]                         print this message

   [-c chain]                   chain. default "$_chain"

   [-ic interaction category]   interaction category. default "$_ic"

USAGE

if ((scalar @ARGV == 0) and -t STDIN) { print $usage and exit; }

GetOptions ("h" => \$print_usage,
            "c=s" => \$_chain,
            "ic=s" => \$_ic);

my $last_ic = undef;
my $start_position = undef;
my $end_position = undef;

my %stretches = ();
my %ics = ();

my $last_chain = undef;

my $num_stretch = 0;

my $ic = undef;

while (<>)
{
   chomp;
   my @fields = split "\t";
   my $chain = $fields[0];
   my $position = $fields[1];
   $ic = $fields[2];

   if (defined $last_chain)
   {
      if ($chain ne $last_chain)
      {
         die ("ERROR: this only supports one chain\n");
      }
   }
   $last_chain = $chain;

   # keep track of observed categories
   $ics{$ic} = 1;

   # in a stretch of same category
   if ($ic eq $last_ic)
   {
      $end_position = $position;
   } 
   else 
   {
      # start of a new stretch
      warn "new stretch: $_\n";

      # record last stretch
      if (defined $start_position)
      {
         # in case of length 1 stretch, end position won't be set
         if ($end_position < $start_position) 
         { 
            $end_position = $start_position; 
         }
         push @{$stretches{$last_ic}{start}},    $start_position;
         push @{$stretches{$last_ic}{end}},      $end_position;
      }

      # keep track of this new stretch's start
      $start_position = $position;
      $num_stretch += 1;
   }

   $last_ic = $ic;
}

# last stretch
push @{$stretches{$last_ic}{start}},    $start_position;
push @{$stretches{$last_ic}{end}},      $end_position;

# create one chimera named selection for each IC

# info from chimera documentation for selecting residues
# https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/midas/frameatom_spec.html
# #3:45-83,90-98
# - residues 45-83 and 90-98 in model 3

foreach my $ic (keys %ics) 
{
   my @starts = @{$stretches{$ic}{start}};
   my @ends   = @{$stretches{$ic}{end}};

   print "~select\n";

   # construct a chimera select command to select all the stretches for this category
   # select :25-35.A5,1-10.A5
   my $cmd   = "name ". $ic . " /".$last_chain.":";

   my $num_stretches = scalar @starts;
   for (my $i = 0; $i < $num_stretches; $i++) 
   {
      my $start = $starts[$i];
      my $end   = $ends[$i];
      if ($end > $start)
      {
         $cmd .= $start."-".$end;
      }
      else
      {
         $cmd .= $start;
      }
      # remove trailing comma
      $cmd .= ",";
   }
   chop ($cmd);
   print "$cmd\n";
}


