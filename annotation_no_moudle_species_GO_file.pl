#!usr/bin/perl -w

############################################################################
## 2019-7-3
## 通过GO数据库文件go.obb对非模式物种进行GO注释
############################################################################

## 首先下载物种（西瓜）的蛋白序列文件
## 利用eggnog-mapper-1.0.3软件对将蛋白序列文件与真核生物数据库进行比对

	emapper.py -m diamond -i ../../fareeha/ref/watermelon_97103_v2.protein.fa -d data/eggnog.db -o output --cpu 16	

## 根据比对后的文件提取gene名、GO号、KO号三列信息 (先进行GO注释)


#   Cla97C01G000010 GO:0016020, GO:0016021
#   Cla97C01G000030 GO:0004222, GO:0006508, GO:0008233, GO:0016020, GO:0016021, GO:0016787
#   Cla97C01G000050 GO:0005618, GO:0016787
#   Cla97C01G000070 GO:0000160, GO:0005634

## ------------------ 脚本1 ----------------------------------

	#!/usr/bin/perl -w
	use strict;

	my $file = $ARGV[0];
	open I,$file or die "$!";

	my %hash;

	while (my $line = <I>) {
		chomp $line;
		my $a = $1 if $line =~ /^(Cla\w+)/;
		my $b = $1 if $line =~ /\b(GO.*)/;
		$hash{$a} = [];
		
		my @tmp = split (/,\s/,$b);
		push @{$hash{$a}},@tmp;
	}
	close I;
	foreach my $k (keys %hash) {

		my @tmp = @{$hash{$k}};
		
		for (my $i = 0; $i < @tmp;$i++) {
			print "$tmp[$i]\t$k\n";
		}


	}


## 上述脚本处理后,得到两列,第一列为GO名,第二列为gene名,但是这不满足要求,因为富集分析需要gene名所在列不能有重复,GO号所在列重复

#   GO:0016020      Cla97C05G105170
#   GO:0016021      Cla97C05G105170
#   GO:0016301      Cla97C05G105170
#   GO:0016310      Cla97C05G105170
#   GO:0003677      Cla97C08G148310
#   GO:0003700      Cla97C08G148310
#   GO:0005634      Cla97C08G148310
#   GO:0006351      Cla97C08G148310
#   GO:0006355      Cla97C08G148310

## ------------------------- 脚本2 --------------------------


	#!/usr/bin/perl -w
	use strict;

	my ($file1,$file2) = @ARGV;

	my %hash;

	open I, $file1 or die "$!";

	while (defined(my $line = <I>)) {
		
		chomp $line;
		$hash{$line} = [];
		
	}


	open I1, $file2 or die "$!";

	while (defined(my $line = <I1>)) {
		
		chomp $line;
		my ($k,$v) = (split /\t/,$line);
		
		if (exists $hash{$k}) {
			push @{$hash{$k}},$v;
			
		}

	}


	foreach my $k (keys %hash) {
		my @tmp = @{$hash{$k}};
		for (my $i = 0; $i < @tmp; $i++) {
			print "$k\t$tmp[$i]\n";		
		}
	}

## 此脚本需要两个文件
	# 1. 非重复的GO号文件(id.txt)
		GO:0000015
		GO:0000023
		GO:0000025
		GO:0000027
		GO:0000028
		GO:0000035

	# 2. 上面的输出文件

		perl go_2_gene.pl id.txt go_gene.txt > GO2GENE.txt

	## 结果文件

		GO:0017148	Cla97C08G147460
		GO:0017148	Cla97C11G207410
		GO:0017148	Cla97C09G180320
		GO:0017148	Cla97C05G092040
		GO:0017148	Cla97C08G147410
		GO:0017148	Cla97C08G147420
		GO:0017148	Cla97C10G197160
		GO:0017148	Cla97C09G171800


## 要做富集分析还需要第二个文件--- GO号对应的name,以及namespace(BP MF CC)


## -----------------------------  脚本3 ---------------------------------

#use strict;
use Fcntl;

my ($lookup_file,$data_file) = @ARGV;

my $look_up = build_lookup($lookup_file);

#print $look_up -> {'GO:0000005'}."\n";

open I1, $data_file or die "$!";

while (defined(my $line = <I1>)) {
	chomp $line;
	my $id = $line;
	if (exists $look_up ->{$id}) {
		print $id.";".$look_up -> {$id}{'name'}.";".$look_up -> {$id}{'namespace'}."\n";
	}

}


sub build_lookup {
	my ($file) = $_[0];
	open I,$file or die "$!";
	
	require SDBM_File;
	tie (my %lookup,'SDBM_File',"lookup", O_RDWR | O_CREAT,0666) or die "$!";
	
	while (defined(my $line = <I>)) {
		chomp $line;
		if ($line =~ /^id:\s+(GO:.*)/) {
			my $k = $1;
		       #print $k;
			my $v = <I>;
			$v =~ /name:\s+(.*)/;
			$v = $1;
			my $class = <I>;
			
			$class =~ /namespace:\s?(.*)/;
			$class = $1;  	
			$lookup{$k}{'name'} = $v;
			$lookup{$k}{'namespace'} = $class;			
		}
		
	}
	close I;
	return \%lookup;
}

## 需要两个文件
	#1. GO数据库注释文件 go.obo

		http://current.geneontology.org/ontology/go.obo

	#2. GO号文件 id.txt

		perl extract_go2term.pl go.obo id.txt > term_name.txt

	
	# 结果文件(以;分隔每一列信息)


	GO:1904823;purine nucleobase transmembrane transport;biological_process
	GO:1905516;positive regulation of fertilization;biological_process
	GO:1905614;negative regulation of developmental vegetative growth;biological_process


	
######################################################################  END  ############################################################################


# 到此为止已经成功的构建了非模式物种GO注释文件,后面在 R中利用 Clusterprofile包的 enricher() 函数即可轻松的进行富集分析















