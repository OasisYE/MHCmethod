use strict;
use File::Basename;
use FindBin '$Bin';

=head1 sample

        perl MHC_sam2snp_blast.pl sam_sort_file version(hg19) > reads_snp&indel_file

=head1 done

=cut
die `pod2text $0` if(@ARGV < 1);
############################################################################################################

my$version=$ARGV[1]||"hg19";
my($sample,$gene,$exon)=(split /[\._]/,basename$ARGV[0])[0,-3,-2];
$gene=$gene."_".$exon;
my$path=dirname$ARGV[0];
my$bin=dirname$Bin;
my$mhc_data="$bin/database/MHC.$version.database";
my$mhc_type="$bin/database/database_$version/$gene.type";
my($pos_start,$pos_end,%data,%reads);

open IN,"< $mhc_data" or die "$!->$mhc_data\n";
while(<IN>){
	chomp;
	my($exon,$pos1,$pos2) = (split)[0,2,3];
	if($exon eq $gene){$pos_start=$pos1;$pos_end=$pos2;}
}

open IN,"< $mhc_type" or die "$!->$mhc_type\n";
while(<IN>){
	chomp;my@ln = split;
	my@indel;foreach my$i(3..$#ln){push @indel,$ln[$i];}
	my@snp = split /-/,$ln[2];
	@{$data{$ln[0]}{"snp"}}=@snp;
	@{$data{$ln[0]}{"indel"}}=@indel;
	$data{$ln[0]}{"base"}=$ln[1];
}
##############################################################################################################

open IN,"< $ARGV[0]";open OUT,"> $path/$sample\_$gene.fasta";
while(<IN>){
	chomp;my@ln = split;
	my($index,$len) = check($ln[10],"#");
	
	if($index == 0){$ln[9] = substr($ln[9],$len);}
	if($index > 0){$ln[9] = substr($ln[9],0,$index);}
	if(length$ln[9]<60){next;} #wj 2014-12-17

	my$flag=$ln[1] & 0x40;
	if($flag){$ln[0]=$ln[0]."/1";}
	else{$ln[0]=$ln[0]."/2";}

	print OUT ">$ln[0]\n$ln[9]\n";$reads{$ln[0]}=$ln[9];
}

if(
$gene=~/A_2/ || $gene=~/A_3/ || $gene=~/A_4/ || 
$gene=~/B_1/ || $gene=~/B_2/ || $gene=~/B_3/ || 
$gene=~/C_3/ || $gene=~/C_4/ || $gene=~/C_5/ || 
$gene=~/DRB1/ || $gene=~/DQB1/){
	open IN,"< $path/$sample\_single_unmap.sam";
	while(<IN>){
		chomp;my@ln = split;next if($ln[9]=~/N/);
		my($index,$len) = check($ln[10],"#");

		if($index == 0){$ln[9] = substr($ln[9],$len);}
		if($index > 0){$ln[9] = substr($ln[9],0,$index);}
		if(length$ln[9]<60){next;} #wj 2014-12-17

		my$flag=$ln[1] & 0x40;
		if($flag){$ln[0]=$ln[0]."/1";}
		else{$ln[0]=$ln[0]."/2";}

		print OUT ">$ln[0]\n$ln[9]\n";$reads{$ln[0]}=$ln[9];
	}
}

if($gene=~/DRB1/ || $gene=~/DQB1/){
	open IN,"< $path/$sample\_pair_unmap.sam";
	while(<IN>){
		chomp;my@ln = split;next if($ln[9]=~/N/);
		my($index,$len) = check($ln[10],"#");

		if($index == 0){$ln[9] = substr($ln[9],$len);}
		if($index > 0){$ln[9] = substr($ln[9],0,$index);}
		if(length$ln[9]<60){next;} #wj 2014-12-17

		my$flag=$ln[1] & 0x40;
		if($flag){$ln[0]=$ln[0]."/1";}
		else{$ln[0]=$ln[0]."/2";}

		print OUT ">$ln[0]\n$ln[9]\n";$reads{$ln[0]}=$ln[9];
	}
}
#########################################################################################################

my%hash;
`$bin/database/blastall -i $path/$sample\_$gene.fasta -d $bin/database/database_$version/$gene/$gene.fa -o $path/$sample\_$gene.blast -p blastn -e 10 -v 1 -b 1 -F F -m 8`;
open IN,"< $path/$sample\_$gene.blast" or die "$!->blast\n";
while(<IN>){
	chomp;my@ln=split;
	if($ln[5]!=0 || $ln[4]/$ln[3]>0.015 || $ln[3]<15){next;}
	if($ln[8]>$ln[9]){my$reverse_temp=$ln[8];$ln[8]=$ln[9];$ln[9]=$reverse_temp;}
	if($ln[8]==1 || $ln[9]==length$data{$ln[1]}{"base"} || $ln[3]==length$reads{$ln[0]}){

		my($index,$len) = check($data{$ln[1]}{"base"},"*");
		my($pos_rela_1,$pos_rela_2);
		if($index == 0){$pos_rela_1 = $pos_start+$ln[8]-1+$len;}
		else{$pos_rela_1 = $pos_start+$ln[8]-1;}
		
		my(@snp_reads,@indel_reads);
		if(@{$data{$ln[1]}{"indel"}}==0){
			$pos_rela_2 = $pos_rela_1+$ln[3]-1;
		}
		
		if(${$data{$ln[1]}{"indel"}}[0]=~/-D-/){
			my($pos_indel,$base_indel) = split /-D-/,${$data{$ln[1]}{"indel"}}[0];
			if($pos_indel<$pos_rela_1){
				$pos_rela_1=$pos_rela_1+length$base_indel;$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
			if($pos_indel>=$pos_rela_1 && $pos_indel<($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_rela_1+$ln[3]-1+length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel>=($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
		}

		if(${$data{$ln[1]}{"indel"}}[0] =~ /-I-/){
			my($pos_indel,$base_indel) = split /-I-/,${$data{$ln[1]}{"indel"}}[0];
			if(($pos_indel+length$base_indel)<$pos_rela_1){
				$pos_rela_1=$pos_rela_1-length$base_indel;$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
			if($pos_indel<=$pos_rela_1 && ($pos_indel+length$base_indel)>$pos_rela_1){
				$pos_rela_1=$pos_indel;$pos_rela_2=$pos_rela_1+$ln[3]-1-length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel>$pos_rela_1 && ($pos_indel+length$base_indel)<$pos_rela_1+$ln[3]){
				$pos_rela_2=$pos_rela_1+$ln[3]-1-length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel<$pos_rela_1 && ($pos_indel+length$base_indel)>=($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_indel+length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel>($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
			if($pos_rela_1<$pos_indel && ($pos_rela_1+$ln[3])<($pos_indel+length$base_indel)){
				$pos_rela_1=$pos_indel;$pos_rela_2=$pos_indel+1;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
		}

		foreach my$snp(@{$data{$ln[1]}{"snp"}}){
			my$pos=(split /:/,$snp)[0];
			if($pos>=$pos_rela_1 && $pos<=$pos_rela_2){push @snp_reads,$snp;}
		}

		my@new_indel_reads;
		foreach my$indel(@indel_reads){
			my$pos=(split /-/,$indel)[0];
			if($pos>=$pos_end || $pos<$pos_start){next;}
			push @new_indel_reads,$indel;
		}

		my$snp_reads=join("-",@snp_reads)||"None";
		my$indel_reads=join(":",@new_indel_reads)||"None";
		push @{$hash{$pos_rela_1."*".$pos_rela_2."*".$snp_reads."*".$indel_reads}},$ln[0];
	}
}

foreach my$kes(sort keys %hash){
	my($pos1,$pos2)=split /\*/,$kes;
	if($pos1>$pos2 || ($pos2-$pos1<15)){next;}
	my%tem;@tem{@{$hash{$kes}}}=();@{$hash{$kes}}=keys %tem;
	print "$kes\t@{$hash{$kes}}\n";
}
`rm $path/$sample\_$gene.fasta $path/$sample\_$gene.blast`;
#########################################################################################

sub check{
	my$index = index($_[0],$_[1]);
	my@base = split //,$_[0];my$len=0;
	foreach my$i($index..$#base){if($base[$i] eq $_[1]){$len++;}}
	return($index,$len);
}

