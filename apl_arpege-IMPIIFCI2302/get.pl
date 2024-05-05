#!/usr/bin/perl -w 
#

use strict;
use lib "$ENV{HOME}/fxtran-acdc/lib";
use local::lib;
use Fxtran;
use FileHandle;
use Data::Dumper;
use File::Find;
use File::Spec;
use File::Basename;


sub getPack
{
  my $pack = '/home/gmap/mrpm/marguina/pack/49t2_main-ecrad.00.IMPIIFCI2302.x';

  my %h;

  if (-f "$pack/.find.pl")
    {
      %h = %{ do "$pack/.find.pl" };
    }
  else
    {
      my @view = do { my $fh = 'FileHandle'->new ("<$pack/.gmkview"); <$fh> };
      chomp for (@view);
     
      for my $view (reverse (@view))
        {
          my $dir = "$pack/src/$view";
          &find ({wanted => sub { my $f = $File::Find::name; $h{&basename ($f)} = $f; }, no_chdir => 1}, "$dir/.");
        }
      
      local $Data::Dumper::Terse = 1;
      'FileHandle'->new (">$pack/.find.pl")->print (&Dumper (\%h));
    }

  return \%h;
}

my $h = &getPack ();

my $f90 = shift;

system ("./strip.pl $f90");

my $d = &Fxtran::parse (location => $f90, fopts => [qw (-line-length 500)]);

my %mod = map { ($_, 1) } &F ('.//use-stmt/module-N', $d, 1);
my @mod = map { lc ($_) } sort keys (%mod);

for my $mod (@mod)
  {
    $mod =~ s/_n$/n/o;

    my $f = "$mod.F90";
    next if (-f $f);

    if (my $g = $h->{$f})
      {
        print "$g\n";
        system ("cp $g $f");
      }
    else
      {
        print "$mod not found\n";
      }
  }


