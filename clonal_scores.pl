use strict;
use warnings;

# clonal_scores.pl by EMB 14APR2021
# usage: perl clonal_scores.pl <sequence_data.csv> <percent AA similarity in HCDR3>
# To retrieve clonal expansion and clonal connection scores within and between B cell subsets.

open(IN, "<$ARGV[0]") || die $!;
open(OUT, ">$ARGV[0].$ARGV[1].Results.txt") || die $!;

my $similarity = $ARGV[1]/100;
my %cluster;
my %cells;
my $id = 1;
my $member = 1;
my $flag;

# Read in data and group sequences that are part of clonotypes
while(my $line = <IN>){
	if($line =~ m/Isotype/){next}
	chomp $line;
	my @data = split(',', $line);

	if(!%cluster){
		$cluster{$id}{members}++;
		$cluster{$id}{$member}{cell} = $data[0];
		$cluster{$id}{$member}{Ig} = $data[1];
		$cluster{$id}{$member}{TP} = $data[2];
		$cluster{$id}{$member}{V} = $data[3];
		$cluster{$id}{$member}{J} = $data[4];
		$cluster{$id}{$member}{HCDR3} = $data[5];
	}else{
		foreach my $seq (keys %cluster){
		$flag = 0;
			for(my $i = 1; $i < $cluster{$seq}{members} + 1; $i++){
				if($data[3] ne $cluster{$seq}{$i}{V}){next}
				if($data[4] ne $cluster{$seq}{$i}{J}){next}
				if(length($data[5]) != length($cluster{$seq}{$i}{HCDR3})){next}

				my @str1 = split(//, $data[5]);
				my @str2 = split(//, $cluster{$seq}{$i}{HCDR3});
				my $common;

				# Calculate difference in AA seq between current HCDR3 sequence and stored HCDR3
				for(my $j = 0; $j < length($data[5]); $j++){
					if($str1[$j] eq $str2[$j]){
						$common++;
					}
				}

				# Store as clonotype member if at least xx% similarity on AA level (provided by user)
				if($common && $common/length($cluster{$seq}{$i}{HCDR3}) >= $similarity){
					$cluster{$seq}{members}++;
					my $member_id = $cluster{$seq}{members};
					$cluster{$seq}{$member_id}{cell} = $data[0];
					$cluster{$seq}{$member_id}{Ig} = $data[1];
					$cluster{$seq}{$member_id}{TP} = $data[2];
					$cluster{$seq}{$member_id}{V} = $data[3];
					$cluster{$seq}{$member_id}{J} = $data[4];
					$cluster{$seq}{$member_id}{HCDR3} = $data[5];
					$flag = 1;
				}
				if($flag){last}
			}
			if($flag){last}
		}
	}

	if(!$flag){
		$id++;
		$cluster{$id}{members}++;
		$cluster{$id}{$member}{cell} = $data[0];
		$cluster{$id}{$member}{Ig} = $data[1];
		$cluster{$id}{$member}{TP} = $data[2];
		$cluster{$id}{$member}{V} = $data[3];
		$cluster{$id}{$member}{J} = $data[4];
		$cluster{$id}{$member}{HCDR3} = $data[5];
	}
}

# Tally clonal connections
foreach my $seq (sort {$a <=> $b} keys %cluster){
	if($cluster{$seq}{members} == 1){
		$cells{$cluster{$seq}{1}{TP}}{$cluster{$seq}{1}{cell}}{$cluster{$seq}{1}{Ig}}++;
	}else{
		for(my $i = 1; $i < $cluster{$seq}{members} + 1; $i++){
			$cells{$cluster{$seq}{$i}{TP}}{$cluster{$seq}{$i}{cell}}{$cluster{$seq}{$i}{Ig}}++;

			for(my $j = 1; $j < $cluster{$seq}{members} + 1; $j++){
				if($i == $j){next}
				my $cell_comparison = $cluster{$seq}{$i}{cell}."_".$cluster{$seq}{$i}{Ig}."_".$cluster{$seq}{$j}{cell}."_".$cluster{$seq}{$j}{Ig};
				my $tp1 = $cluster{$seq}{$i}{TP};
				my $tp2 = $cluster{$seq}{$j}{TP};

				# To prevent one sequence from being counted more than once in larger clonotypes,
				# each cluster member only contributes 1 count to each cell comparison with the same TP and Ig.
				# For example, if a clonotype within a B cell subset consists of 4 sequences,
				# each sequence has 3 matches, resulting in 4 x 3 = 12 counts towards the clonal score.
				# By only counting each sequence once, the total number of matches counting towards the clonal score is 4.

				if(!$cluster{$seq}{$i}{$tp1}{$tp2}{$cell_comparison}{cell1}){
					$cells{$tp1}{$tp2}{$cell_comparison}{cell1}++;
					$cluster{$seq}{$i}{$tp1}{$tp2}{$cell_comparison}{cell1} = 1;
				}
				if(!$cluster{$seq}{$j}{$tp1}{$tp2}{$cell_comparison}{cell2}){
					$cells{$tp1}{$tp2}{$cell_comparison}{cell2}++;
					$cluster{$seq}{$j}{$tp1}{$tp2}{$cell_comparison}{cell2} = 1;
				}
			}
		}
	}
}

# Print clonal expansion and clonal connection scores
print OUT "TP1", "\t", "Ig1", "\t", "Cell1", "\t", "Connections1", "\t", "Total count cells1", "\t";
print OUT "TP2", "\t", "Ig2", "\t", "Cell2", "\t", "Connections2", "\t", "Total count cells2", "\t";
print OUT "Clonal score", "\n";

my @tp = ("0", "3", "6");
my @cells = ("Naive_B_cell", "FcRL5-_classical_MBC", "FcRL5+_classical_MBC", "Atypical_MBC");
my @Ig = ("IgM", "IgG");

for(my $i1 = 0; $i1 < 3; $i1++){
	for(my $i2 = 0; $i2 < 3; $i2++){
		for(my $j1 = 0; $j1 < 4; $j1++){
			for(my $j2 = 0; $j2 < 4; $j2++){
				for(my $k1 = 0; $k1 < 2; $k1++){
					if($cells[$j1] eq "Naive_B_cell" && $Ig[$k1] eq "IgG"){next}
					for(my $k2 = 0; $k2 < 2; $k2++){
						if($cells[$j2] eq "Naive_B_cell" && $Ig[$k2] eq "IgG"){next}

						my $tp1 = $tp[$i1];
						my $tp2 = $tp[$i2];
						my $cell1 = $cells[$j1];
						my $cell2 = $cells[$j2];
						my $Ig1 = $Ig[$k1];
						my $Ig2 = $Ig[$k2];
						my $combi = $cell1."_".$Ig1."_".$cell2."_".$Ig2;

						if($cells{$tp1}{$tp2}{$combi}{cell1}){
							print OUT $tp1, "\t", $Ig1, "\t", $cell1, "\t", $cells{$tp1}{$tp2}{$combi}{cell1}, "\t", $cells{$tp1}{$cell1}{$Ig1}, "\t";
						}else{
							print OUT $tp1, "\t", $Ig1, "\t", $cell1, "\t", "0", "\t", $cells{$tp1}{$cell1}{$Ig1}, "\t";
						}
						if($cells{$tp1}{$tp2}{$combi}{cell2}){
							print OUT $tp2, "\t", $Ig2, "\t", $cell2, "\t", $cells{$tp1}{$tp2}{$combi}{cell2}, "\t", $cells{$tp2}{$cell2}{$Ig2}, "\t";
						}else{
							print OUT $tp2, "\t", $Ig2, "\t", $cell2, "\t", "0", "\t", $cells{$tp2}{$cell2}{$Ig2}, "\t";
						}

						# Print clonal connection score for the smallest of the two B cell subsets.
						if($cells{$tp1}{$tp2}{$combi}{cell1}){
							if($cells{$tp1}{$cell1}{$Ig1} < $cells{$tp2}{$cell2}{$Ig2}){
								printf OUT ("%.2f", $cells{$tp1}{$tp2}{$combi}{cell1}/$cells{$tp1}{$cell1}{$Ig1}*100);
								print OUT "\n";
							}else{
								printf OUT ("%.2f", $cells{$tp1}{$tp2}{$combi}{cell2}/$cells{$tp2}{$cell2}{$Ig2}*100);
								print OUT "\n";
							}
						}else{
							print OUT "0", "\n";
						}
					}
				}
			}
		}
	}
}

close IN;
close OUT;
