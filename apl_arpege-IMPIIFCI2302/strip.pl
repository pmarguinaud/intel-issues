#!/usr/bin/perl -w

use strict;
use lib "$ENV{HOME}/fxtran-acdc/lib";
use local::lib;
use FileHandle;
use Data::Dumper;
use File::Basename;
use FindBin qw ($Bin);

use Common;

use Fxtran;

sub stripExec
{
  my $pu = shift;

  for (&F ('.//program-unit', $pu))
    {
      $_->unbindNode ();
    }

  my $stmt = $pu->firstChild;

  (my $kind = $stmt->nodeName ()) =~ s/-stmt$//o;
  
  my ($name) = &F ('./' . $kind . '-N/N/n/text()', $stmt, 1);
  my @args = &F ('.//dummy-arg-LT//arg-N/N/n/text()', $stmt, 1);
  
  my %stmt;
  
  if (my ($result) = &F ('./result-spec/N/n', $stmt, 1))
    {
      my @en = &F ('.//EN-decl[./EN-N[translate(string(N),"abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ")="?"]]', uc ($result), $pu);
      for my $en (@en)
        {
          my $stmt = &Fxtran::stmt ($en);
          $stmt{$stmt} = $stmt;
        }
    }
  elsif (my ($name) = &F ('./function-N', $stmt, 1))
    {
      my @en = &F ('.//EN-decl[./EN-N[translate(string(N),"abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ")="?"]]', uc ($name), $pu);
      for my $en (@en)
        {
          my $stmt = &Fxtran::stmt ($en);
          $stmt{$stmt} = $stmt;
        }
    }
  
  # Keep first & last statements
  
  $stmt{$pu->firstChild} = $pu->firstChild;
  $stmt{$pu->lastChild}  = $pu->lastChild;
  
  # Keep declaration statements referencing arguments

  for my $arg (@args)
    {
# translate('some text','abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ')
      my @en = &F ('.//EN-decl[./EN-N[translate(string(N),"abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ")="?"]]', uc ($arg), $pu);
      for my $en (@en)
        {
          my $stmt = &Fxtran::stmt ($en);
          $stmt{$stmt} = $stmt;
        }
    }
  
  # Strip blocks (these may contain use statements)
  
  for (&F ('.//ANY-construct', $pu))
    {
      $_->unbindNode ();
    }
  
  # Keep use statements
  
  for (&F ('.//use-stmt', $pu))
    {
      $stmt{$_} = $_;
    }
  
  my @stmt = &F ('.//ANY-stmt', $pu);
  
  for my $stmt (@stmt)
    {
      next if ($stmt->nodeName eq 'implicit-none-stmt');
      $stmt->unbindNode () unless ($stmt{$stmt});
    }
  
  # Strip labels
  for (&F ('.//label', $pu))
    {
      $_->unbindNode ();
    }
  
  # Strip comments
  
  for (&F ('.//C', $pu))
    {
      next if ($_->textContent =~ m/^!\$acc\s+routine/o);
      $_->unbindNode ();
    }
  
  # Strip includes
  
  for (&F ('.//include', $pu))
    {
      $_->unbindNode ();
    }

  # Strip defines

  for (&F ('.//cpp[starts-with (text(),"#define ")]', $pu))
    {
      $_->unbindNode ();
    }

  my @text = &F ('.//text()[translate(.," ?","")=""]', "\n", $pu);

  for my $text (@text)
    {
      if ($text->data =~ m/\n/goms)
        {
          $text->setData ("\n");
        }
    }

}

my $f90 = shift;

my $d = &Fxtran::parse (location => $f90, fopts => [qw (-no-cpp -construct-tag -line-length 500 -no-include)]);

my @pu = &F ('./object/file/program-unit/program-unit', $d);

for my $pu (@pu)
  {
    &stripExec ($pu);
  }

#print $d->textContent;

'FileHandle'->new (">$f90")->print ($d->textContent);


