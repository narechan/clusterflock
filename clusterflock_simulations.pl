#!/usr/bin/perl -w
=head1 NAME                                                                                                 

clusterflock_simulations.pl                                                                             

=head1 SYNOPSIS                                                                                              
                                                                                                      
  clusterflock_simulations.pl --                                                                          
    Functions as a testing wrapper around clusterflock.pl. Will accept simulated 
    data and report clustering success using the Jaccard Index.
    
    Note that for the simulations script, clusterflock does not output visuals and that 
    perception (-s) is set to all, and bin lattice spatial hashing (-b) set to on.
    
Options:                                                                                         
                                                                                                          
 --help        Show brief help and exit                                                                    
 --indir       Is your directory of flocked entities (ortholog fastas)                                       
 --cfconfig    Is the configuration for your flocking simulation                                            
 --outdir      Is your output dir                                                                             
 --list        Is a list of precomputed distance stats (LD in the evolutionary case)                       
 --cfreps      Set to the number of flocking replicates you want to perform
 --jarfile     Is the complete path to the ELKI jar file installed as a dependency
 --kmeans      Is set to the number of clusters expected given the number of trees simulated
    in the data analyzed
 --finalframe  Is the number of iterations specified for each flocking replicate 

=head1 DESCRIPTION                                                                                               

Run flocking for a set of simulated data and check to see whether clustering success as 
measured using the Jaccard aligns with expectations.

=head1 AUTHOR                                                                                            
 
Apurva Narechania                                                                                                
anarechania *a|t* amnh.org

=head1 COPYRIGHT                                                                                            

This library is free software (GNU GPL);                                                                      
you can redistribute it and/or modify                                                                          
it under the same terms as Perl itself.                                                                       
                                                                                                               
=cut                                                                                               

# ----------------------------------------------------             


#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Statistics::Descriptive;

my ($help, $cfconfig, $cfreps, $procs, $outdir, $indir, $ldlist, $jarfile, $k, $finalframe);

GetOptions(
    'h|help'          => \$help,
    'c|cfconfig=s'    => \$cfconfig,
    'o|outdir=s'      => \$outdir,
    'r|cfreps=s'      => \$cfreps,
    'p|procs=s'       => \$procs,
    'i|indir=s'       => \$indir,
    'l|ldlist=s'      => \$ldlist,
    'j|jarfile=s'     => \$jarfile,
    'k|kmeans=s'      => \$k,
    'f|finalframe=s'  => \$finalframe,
	   ) or pod2usage;
pod2usage if $help;

`mkdir -p $outdir/cf`;

#####MAIN#####

# fork the clusterflock processes
my $pm = Parallel::ForkManager->new($procs);
for (my $i = 1; $i <= $cfreps; $i++){
    $pm->start and next;

    ### run clusterflock

#    print STDERR "cflocks $i\n";
    `clusterflock.pl -i $indir -c $cfconfig -l $ldlist -s all -b 1 -o $outdir/cf/$i-cf &> $outdir/cf/stderr.$i-cf`;

    ### run the elki k-means program

#    print STDERR "kmeans $i\n";
    `java -XX:+UseSerialGC -cp $jarfile de.lmu.ifi.dbs.elki.application.KDDCLIApplication -dbc.in $outdir/cf/$i-cf/logs/$finalframe.log -resulthandler ResultWriter -algorithm clustering.kmeans.KMeansLloyd -kmeans.k $k > $outdir/cf/$i-cf/flocks/$i.elki.kmeans`;

    # parse the kmeans clusters                                                        
    my $clusters = {};
    my $idseen = {};
    my $counter = -1;
    open (F, "$outdir/cf/$i-cf/flocks/$i.elki.kmeans");
    while (my $line = <F>){
	chomp $line;

	if ($line =~m/#\sCluster\:/){
	    $counter++;
	}
	elsif ($line =~m/^ID/){
	    my @line = split (/\s/, $line);
	    if (exists ($idseen->{$line[3]})){
		last;
	    }
	    else{
		$clusters->{$counter}->{$line[3]} = 1;
		$idseen->{$line[3]} = 1;
	    }
	}
	else {
	    next;
	}
    }
    close (F);
    
    # rewrite the clusters
    open (O, ">$outdir/cf/$i-cf/flocks/$i.kmeans");
    foreach my $clust (sort {$a <=> $b} keys %$clusters){
	foreach my $id (sort keys %{$clusters->{$clust}}){
	    print O "$clust\t$id\n";
	}
    }
    close (O);
    
    ### calculate external measures of clustering success
    
    my $opticsclass = {};
    open (O, "$outdir/cf/$i-cf/flocks/$i.kmeans");
    while (my $line = <O>){
	chomp $line;
	my ($cluster, $gene) = split (/\t/, $line);
	$opticsclass->{$gene} = $cluster;
    }
    close (O);

    # parse the flocking info                                                          
    opendir (I, "$indir");
    my $files = {};
    my @files = sort (readdir (I));
    shift @files;
    shift @files;
    closedir (I);
    foreach my $file (@files){
	$files->{$file} = 1;
    }

    # calculate the jaccard
    my $tp = 0;
    my $tn = 0;
    my $fp = 0;
    my $fn = 0;
    foreach my $f1 (keys %$files){
	foreach my $f2 (keys %$files){
	    next if ($f1 eq $f2);
	    
	    # get data for the true clusters                                           
	    my ($cluster1, $gene1) = split (/\_/, $f1);                                
	    my ($cluster2, $gene2) = split (/\_/, $f2);                                
#	    my ($acc1, $cluster1, $rate1) = split (/\_/, $f1);
#	    my ($acc2, $cluster2, $rate2) = split (/\_/, $f2);
	    
	    # get data for the calculated clusters                                     
	    my $bcluster1 = $opticsclass->{$f1};
	    my $bcluster2 = $opticsclass->{$f2};
	    
	    # contingencies                                                            
	    if (($cluster1 == $cluster2) and ($bcluster1 == $bcluster2)){
		$tp++;
	    }
	    elsif (($cluster1 != $cluster2) and ($bcluster1 != $bcluster2)){
		$tn++;
	    }
	    elsif (($cluster1 != $cluster2) and ($bcluster1 == $bcluster2)){
		$fp++;
	    }
	    elsif (($cluster1 == $cluster2) and ($bcluster1 != $bcluster2)){
		$fn++;
	    }
	    else{
		print STDERR "Unknown contingency\n";
		die;
	    }
#	    print STDERR "$f1\t$f2\t$tp\t$tn\t$fp\t$fn\n";
	}
	delete $files->{$f1};
    }
    
    open (T, ">$outdir/cf/$i.jaccard");
    my $jaccard   = $tp / ($tp + $fp + $fn);
    print T "$jaccard\n";
    close (T);

    $pm->finish;
}
$pm->wait_all_children;

# get averages of simulation metrics across all cf reps
my @kjaccard;
for (my $i = 1; $i <= $cfreps; $i++){
    open (K, "$outdir/cf/$i.jaccard");
    while (my $line = <K>){
	chomp $line;
	push (@kjaccard, $line);
    }
    close (K);
}

my $kjaccard = avg (\@kjaccard);
print "$kjaccard\n";

### subs ###

sub avg {
    my $array = shift;
    my $statobj = Statistics::Descriptive::Full->new();
    $statobj->add_data(@$array);
    my $mean = $statobj->mean();
    return ($mean);
}
