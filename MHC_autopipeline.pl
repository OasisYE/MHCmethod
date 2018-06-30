use strict;
use File::Basename;
use FindBin '$Bin';
use Getopt::Long;

=head1 Description

	MHC typing pipeline ; V1.2

=head1 Contact


=head1 Analysis steps

	HLA typing steps as follows:

		1	collect reads from each exons
		2	realignment by blast
		3	assemble reads to haps and filter haps
		5	combination all exons' haps to seclct the best two types

=head1 usage

	Options:
		-i input bam file(Alignment by bwa)
		-od output dir
		-v ref version(hg19)
		-h help

=head1 sample

	perl MHC_autopipeline_total.pl -i bam -od outdir -v version

=head1 done

=cut

my($bam,$outdir,$version);
GetOptions(
	"i:s"=>\$bam,
	"od:s"=>\$outdir,
	"v:s"=>\$version,
);
die `pod2text $0` unless($bam);
###########################################################################################

$outdir=$outdir||".";
$version=$version||"hg19";
my$bin=dirname$Bin;
my$mhc_pos="$bin/database/MHC.$version.database";
my$sample=(split /[\._]/,basename$bam)[0];
my$freq_mhc="$bin/database/MHC_freq.database";
my@gene=qw/A B C DRB1 DQB1 DPA1 DPB1 DQA1 G DMA DMB DOA DOB DRA E F H J K MICA MICB P TAP1 TAP2 V L/;
`mkdir -p -m 755 $outdir/$sample` unless(-d "$outdir/$sample");
open OUT,"> $outdir/$sample/$sample.type";

my%exon;
open IN,"< $mhc_pos" or die "$!->$mhc_pos\n";
while(<IN>){
	chomp;
	my($exon,$pos1,$pos2)=(split)[0,2,3];
	$exon{$exon}=$pos1."*".$pos2;
}
##########################################################################################################

unless(-e "$outdir/$sample/$sample\_single_unmap.sam"){
#	my$path_bam=dirname$bam;my$sam_name=(split /\./,basename$bam)[0];
#	my$unmap_bam_path="$path_bam/$sam_name.unmap.bam";
#	`$bin/database/samtools view $unmap_bam_path>$outdir/$sample/$sample\_pair_unmap.sam`;
	if ( $version eq "b38" ){ #wj change, 2014-12-16
		`$bin/database/samtools view -f 4 $bam chr6 > $outdir/$sample/$sample\_single_unmap.sam`;
	}else{
		`$bin/database/samtools view -f 4 $bam 6 > $outdir/$sample/$sample\_single_unmap.sam`;
	}
	`$bin/database/samtools view -f 4 $bam|awk '{if(\$3=="*")print}' > $outdir/$sample/$sample\_pair_unmap.sam`;
}
foreach my$kes(sort keys %exon){
	my($pos1,$pos2) = split /\*/,$exon{$kes};
	unless(-e "$outdir/$sample/$sample\_$kes.sam"){
		if ( $version eq "b38" ){ #wj change, 2014-12-16
			`$bin/database/samtools view $bam chr6:$pos1-$pos2 > $outdir/$sample/$sample\_$kes.sam`;
		}else{
			`$bin/database/samtools view $bam 6:$pos1-$pos2 > $outdir/$sample/$sample\_$kes.sam`;
		}
	}
	unless(-e "$outdir/$sample/$sample\_$kes.snp"){
		`perl $Bin/MHC_sam2snp_blast.pl $outdir/$sample/$sample\_$kes.sam $version > $outdir/$sample/$sample\_$kes.snp`;
	}
	unless(-e "$outdir/$sample/$sample\_$kes.ld"){
		`perl $Bin/MHC_snp2ld_assemble.pl $outdir/$sample/$sample\_$kes.snp $version > $outdir/$sample/$sample\_$kes.ld`;
	}
}

###########################################################################################################

foreach my$gene(@gene){
	my@type = glob "$outdir/$sample/$sample\_$gene*ld";

	my%gap;my%novel;
	my@type_database=glob "$bin/database/database_$version/*type";
	foreach my$type_database(@type_database){
		open IN,"< $type_database";
		while(<IN>){chomp;my@ln=split;$gap{$ln[0]}++;}
	}
	foreach my$kes(sort keys %gap){$gap{$kes}=@type/$gap{$kes};}

	my%type_gene;my%hash;
	foreach my$type(@type){
		my$gene_exon=(split /[\._]/,basename$type)[-2];
		open IN,"< $type" or die "$!->$type";
		while(my$ln1=<IN>){
			my$ln2=<IN>;chomp($ln1,$ln2);
			my@ln=split /\s+/,$ln1;my@reads=split /\s+/,$ln2;

			@{$hash{$ln[0]}{"reads"}}=@reads;$hash{$ln[0]}{"score"}=$ln[1];
			foreach my$i(3..$#ln){push @{$type_gene{$ln[$i]}},"$gene_exon\_$ln[0]";}
		}
	}

	############################################################################################################

	my%type_score;
	foreach my$kes(keys %type_gene){
		my$score_total=0;
		foreach my$tem_kes(@{$type_gene{$kes}}){
			my($tem_gene,$tem_hap)=split /_/,$tem_kes;
			$score_total+=$hash{$tem_hap}{"score"};
		}
		$type_score{$kes}=$score_total*$gap{$kes};
		if($score_total/@type<0.3){delete($type_gene{$kes});next;}
		$type_score{$kes}=$score_total;
	}

	my%type_comb;
	my%type_gene_sub=%type_gene;
	while(%type_gene_sub){
		my@type_gene_sub=keys %type_gene_sub;
		my$comb_type1=shift@type_gene_sub;
		my@type1_reads;
		foreach my$tem_kes(@{$type_gene_sub{$comb_type1}}){
			my($exon,$hap)=split /_/,$tem_kes;
			push @type1_reads,@{$hash{$hap}{"reads"}};
		}
		delete($type_gene_sub{$comb_type1});
		my$count_hash=0;

		foreach my$kes(keys %type_gene_sub){
			my@total_type_reads=0;

			foreach my$type_hap(@{$type_gene_sub{$kes}}){
				my($type_exon,$type_temp)=split /_/,$type_hap;
				push @total_type_reads,@{$hash{$type_temp}{"reads"}};
			}

			push @total_type_reads,@type1_reads;
			my%tem;@tem{@total_type_reads}=();@total_type_reads=keys %tem;
			$type_comb{$comb_type1."-".$kes}=@total_type_reads;
		}
		if($count_hash==0){$type_comb{$comb_type1."-".$comb_type1}++;}
	}

	my$most_reads=0;my$most_score=0;my$type_comb_total;
	foreach my$kes(keys %type_comb){if($type_comb{$kes}>$most_reads){$most_reads=$type_comb{$kes}};}
	foreach my$kes(keys %type_comb){
		if($type_comb{$kes}==$most_reads){
			my@type_all=split /-/,$kes;
			my$score_comb=$type_score{$type_all[0]}+$type_score{$type_all[1]};
			if($score_comb>$most_score){$most_score=$score_comb;$type_comb_total=$kes;}
		}
	}
	##########################################################################################################

	my@type_comb=split /-/,$type_comb_total;
	my%type_gene1;my%type_gene2;
	my%reads_stack1;my%reads_stack2;

	if($type_score{$type_comb[1]}>$type_score{$type_comb[0]}){
		my$tem_comb=$type_comb[0];$type_comb[0]=$type_comb[1];$type_comb[1]=$tem_comb;
	}

	foreach my$type_kes(@{$type_gene{$type_comb[0]}}){
		my($exon,$hap)=split /_/,$type_kes;
		foreach my$reads(@{$hash{$hap}{"reads"}}){$reads_stack1{$reads}++;}
	}
	foreach my$type_kes(@{$type_gene{$type_comb[1]}}){
		my($exon,$hap)=split /_/,$type_kes;
		foreach my$reads(@{$hash{$hap}{"reads"}}){$reads_stack2{$reads}++;}
	}

	foreach my$type_kes(keys %type_gene){
		my@uniq_stack1;my@uniq_stack2;my@own_reads;
		foreach my$type_hap(@{$type_gene{$type_kes}}){
			my($exon,$hap)=split /_/,$type_hap;
			foreach my$reads(@{$hash{$hap}{"reads"}}){
				push @own_reads,$reads;
				unless(exists $reads_stack1{$reads}){push @uniq_stack1,$reads;}
				unless(exists $reads_stack2{$reads}){push @uniq_stack2,$reads;}
			}
		}

		my$total_type_reads=@own_reads;
		if(scalar@uniq_stack1<=scalar@uniq_stack2){
			foreach my$exon_reads(@uniq_stack1){
				if(exists $reads_stack2{$exon_reads}){
					$total_type_reads=$total_type_reads-1;
				}
			}
			$type_gene1{$type_kes}=$total_type_reads*$type_score{$type_kes};
		}
		else{
			foreach my$exon_reads(@uniq_stack2){
				if(exists $reads_stack1{$exon_reads}){
					$total_type_reads=$total_type_reads-1;
				}
			}
			$type_gene2{$type_kes}=$total_type_reads*$type_score{$type_kes};
		}
	}

	my$final_type1;my$final_type2;my$final_type_score1=0;my$final_type_score2=0;
	foreach my$type_kes(keys %type_gene1){
		if($final_type_score1<$type_gene1{$type_kes}){
			$final_type1=$type_kes;$final_type_score1=$type_gene1{$type_kes};
		}
	}
	foreach my$type_kes(keys %type_gene2){
		if($final_type_score2<$type_gene2{$type_kes}){
			$final_type2=$type_kes;$final_type_score2=$type_gene2{$type_kes};
		}
	}
	###########################################################################################################

	my$type_sec_check=0;
	if($final_type2){
		my@reads_final;my$distan=0;
		foreach my$kes(@{$type_gene{$final_type1}}){
			my($exon,$hap)=split /_/,$kes;
			push @reads_final,@{$hash{$hap}{"reads"}};
		}
		my$reads_final=@reads_final;
		foreach my$kes(@{$type_gene{$final_type2}}){
			my($exon,$hap)=split /_/,$kes;
			push @reads_final,@{$hash{$hap}{"reads"}};
		}
		my%tem;@tem{@reads_final}=();@reads_final=keys %tem;
		my$reads_uniq=@reads_final-$reads_final;

		if($reads_uniq/$reads_final>0.01){$distan=1;}
		my$rate=$type_score{$final_type2}/@type;
		if($rate>=0.3 && $distan==1){$type_sec_check=1;}
	}
	#####################################################################################################

	if($type_sec_check==1){
		my@type_score_core;my@type_score_sec;

		my$second_type1;my$second_type2;my$second_type1_score=0;my$second_type2_score=0;
		foreach my$kes(sort keys %type_gene1){
			if($type_gene1{$kes}==$final_type_score1){push @type_score_core,$kes;next;}
			if($type_gene1{$kes}>$second_type1_score){$second_type1_score=$type_gene1{$kes};$second_type1=$kes;}
		}
		foreach my$kes(sort keys %type_gene2){
			if($type_gene2{$kes}==$final_type_score2){push @type_score_sec,$kes;next;}
			if($type_gene2{$kes}>$second_type2_score){$second_type2_score=$type_gene2{$kes};$second_type2=$kes;}
		}

		my$type1_score=0;my$type1_score_sub=0;my$type2_score=0;my$type2_score_sub=0;
		$type1_score=$type_score{$final_type1};
		if(exists $type_score{$second_type1}){$type1_score_sub=$type_score{$second_type1};}
		$type2_score=$type_score{$final_type2};
		if(exists $type_score{$second_type2}){$type2_score_sub=$type_score{$second_type2};}
		$second_type1||="---";$second_type2||="---";$final_type1||="---";$final_type2||="---";
		my$total_type1_score=0;my$total_type2_score=0;
		if(@type==$type1_score_sub){$total_type1_score=50;}
		else{
			$total_type1_score=sprintf("%.2f",(@type-$type1_score_sub)**2/((@type-$type1_score_sub)**2+(@type-$type1_score)**2)*100);
		}
		if(@type==$type2_score_sub){$total_type2_score=50;}
		else{
			$total_type2_score=sprintf("%.2f",(@type-$type2_score_sub)**2/((@type-$type2_score_sub)**2+(@type-$type2_score)**2)*100);
		}

		foreach my$type_check(@type_score_core){
			my@type_core=split /[\*:]/,$final_type1;my@type_check=split /[\*:]/,$type_check;
			if($type_check[1] ne $type_core[1]){$final_type1=$type_core[0]."*";last;}
			if($type_check[2] ne $type_core[2]){$final_type1=$type_core[0]."*".$type_core[1];next;}
			if($type_check[3] ne $type_core[3]){$final_type1=$type_core[0]."*".$type_core[1].":".$type_core[2];next;}
			if($type_check[4] ne $type_core[4]){$final_type1=$type_core[0]."*".$type_core[1].":".$type_core[2].":".$type_core[3];next;}
		}
		foreach my$type_check(@type_score_sec){
			my@type_sec=split /[\*:]/,$final_type2;my@type_check=split /[\*:]/,$type_check;
			if($type_check[1] ne $type_sec[1]){$final_type2=$type_sec[0]."*";last;}
			if($type_check[2] ne $type_sec[2]){$final_type2=$type_sec[0]."*".$type_sec[1];next;}
			if($type_check[3] ne $type_sec[3]){$final_type2=$type_sec[0]."*".$type_sec[1].":".$type_sec[2];next;}
			if($type_check[4] ne $type_sec[4]){$final_type2=$type_sec[0]."*".$type_sec[1].":".$type_sec[2].":".$type_sec[3];next;}
		}
		print OUT "$final_type1\t$second_type1\t$total_type1_score\t";
		print OUT "\n";

		print OUT "$final_type2\t$second_type2\t$total_type2_score\t";
		print OUT "\n";
	}
	else{
		my@type_score_all;

		my$second_type1;my$second_type1_score=0;
		foreach my$kes(sort keys %type_gene1){
			if($type_gene1{$kes}==$final_type_score1){push @type_score_all,$kes;next;}
			if($type_gene1{$kes}>$second_type1_score){$second_type1_score=$type_gene1{$kes};$second_type1=$kes;}
		}

		my$type1_score=0;my$type1_score_sub=0;
		$type1_score=$type_score{$final_type1};
		if(exists $type_score{$second_type1}){$type1_score_sub=$type_score{$second_type1};}
		$second_type1||="---";$final_type1||="---";
		my$total_type1_score=0;
		if(@type==$type1_score_sub){$total_type1_score=50;}
		else{
			$total_type1_score=sprintf("%.2f",(@type-$type1_score_sub)**2/((@type-$type1_score_sub)**2+(@type-$type1_score)**2)*100);
		}

		foreach my$type_check(@type_score_all){
			my@type_core=split /[\*:]/,$final_type1;my@type_check=split /[\*:]/,$type_check;
			if($type_check[1] ne $type_core[1]){$final_type1=$type_core[0]."*";last;}
			if($type_check[2] ne $type_core[2]){$final_type1=$type_core[0]."*".$type_core[1];next;}
			if($type_check[3] ne $type_core[3]){$final_type1=$type_core[0]."*".$type_core[1].":".$type_core[2];next;}
			if($type_check[4] ne $type_core[4]){$final_type1=$type_core[0]."*".$type_core[1].":".$type_core[2].":".$type_core[3];next;}
		}
		print OUT "$final_type1\t$second_type1\t$total_type1_score\t";
		print OUT "\n";

		print OUT "$final_type1\t$second_type1\t$total_type1_score\t";
		print OUT "\n";
	}
}
#################################################################################################################

#`rm $outdir/$sample/$sample\_*`;
`echo All done`;
