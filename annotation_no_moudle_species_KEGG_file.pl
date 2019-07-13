#!usr/bin/perl -w

####################################################################################################
## 2019-7-13
## 非模式物种KEGG pathway注释
####################################################################################################

## 先从KAAS网站进行自动注释，获取基因对应的KO号（也可以用eggnog-mapper-1.0.3软件注释的结果，两者差异很小，我自己试过）
	https://www.genome.jp/kegg/kaas/

## 对于KAAS获得的结果文件(KO_gene,txt)，可以利用GO注释的脚本转换为下面的格式
	K02469  Cla97C08G153670
	K11251  Cla97C02G041250
	K11251  Cla97C03G052950
	K11251  Cla97C03G066930
	K11251  Cla97C04G075120
	K11251  Cla97C05G096070
	K11251  Cla97C07G128860
	K11251  Cla97C09G164690
	K11251  Cla97C09G169220
	K11251  Cla97C09G169230

## 对获得KO号需要注释到对应的 Pathway路径

	# 先从KEGG网站下载对应物种的 Pathway 注释文件(若没有本物种的，可以找亲缘近的物种)
		https://www.genome.jp/kegg-bin/get_htext?cmo00001 （用甜瓜代替西瓜）

	# 写脚本获取相关信息
		
		  
		#!usr/bin/perl -w
		#use strict;
		use List::Util qw(uniq);

		my ($file) = $ARGV[0];

		open I,$file or die "$!";

		my %hash;

		#while (my $line = <I>) {
		#	chomp $line;
		#	if ($line =~ /^C\s+(\d+\s+.*)\s+\[PATH/) {
		#		my $id = $1;
		#		print $id,"\n";
		#
		#	}
		#
		#	#`last;
		#
		#
		#
		#}
		#$" = "\n";
		while (my $line = <I>) {
			next if $. == 1;
			my $id = $1 if $line =~ /^C\s+(\d+\s+.*)\s+\[PATH:cmo\d+\]/;
			$hash{$id} = [];
			#print $line;

		}
		close I;

		open I1,$file or die "$!";

		local $/ = "\nC";

		while (my $line = <I1>) {

			next if $. == 1;
			my $ko = $1 if $line =~ /\s+(\d+\s+.*)\s+\[PATH:cmo\d+\]/;
			if (exists $hash{$ko}) {
				my @tmp = $line =~ /\b(K\d+)\b/g;
				@tmp = uniq @tmp;
				#print "@tmp\n";
				push @{$hash{$ko}},@tmp;
			}

		}

		$/ = "\n";
		foreach my $k (keys %hash) {

			my @tmp = @{$hash{$k}};

			for (my $i = 0; $i < @tmp;$i++) {

				print "$k\t$tmp[$i]\n";



			}



		}

		# 此脚本输入甜瓜注释文件即可
			perl ko_id_2_pathway_id.pl cmo00001.keg > pathwaykO_2_genekn.txt
			
		# 结果文件(第一列为Pathway编号，第二列为Pathway name,第三列为KO号)
			03013   RNA transport   K13171
			03013   RNA transport   K14326
			03013   RNA transport   K14327
			03013   RNA transport   K14328

	

	# 利用KO号作为中间值，构建Pathway号对应的gene名
		
	
		# 脚本1 



		use strict;

		my ($ko,$file) = @ARGV;



		open IN, $file or die "$!";

		open IN1, $ko or die "$!";

		my %hash;
		while (<IN1>) {

			chomp;
			$hash{$_} = [];

		}



		while (my $line = <IN>) {

			chomp $line;
			my ($id, $gene) = (split /\t/,$line);
			push @{$hash{$id}},$gene;
			
		}

		$" = ",";

		foreach my $k (keys %hash) {

			print "$k\t@{$hash{$k}}\n";


		}

		# 输入文件
			# KO_gene.txt
			# KO_id 文件
				K00001
   				K00002
				K00006
				K00008

			perl pathway_2_gene.pl ko.id.txt ko_2_gene.txt > ko_gene.txt

		# 结果文件

			K12898  Cla97C01G020450,Cla97C02G049300
			K07583  Cla97C05G096360
			K14326  Cla97C10G186580
		

		# 脚本2


			#!usr/bin/perl -w
			use strict;

			my ($file1,$file2) = @ARGV;

			open I, $file1 or die "$!";

			my %hash;
			while (my $line = <I>) {

				chomp $line;
				my ($id,$gene) = (split /\t/,$line);
				$hash{$id} = $gene;
			}

			#map {print $hash{$_}."\n"} keys %hash;
			close I;


			open I2 , $file2 or die "$!";

			while (my $line = <I2>) {
				chomp $line;
				my ($id) = $line =~ /\b(K\d+)\b/;	
				#print $id."\n";
				my ($left) = $line =~ /^(.*)K\d+/;
				if (exists $hash{$id}) {
					
					print "$left$id\t$hash{$id}\n"
				}

			}



		# 输入文件
			# ko_gene.txt
			# pathwaykO_2_genekn.txt

		# 此脚本可以将 ko_gene.txt 中的gene名根据KO号加到pathwaykO_2_genekn.txt 最后一列

			perl pathway_2_gene1.pl ko_gene.txt pathwaykO_2_genekn.txt > pathwaykO_2_genekn1.txt

		# 结果文件


			03013   RNA transport   K18213  Cla97C03G061430,Cla97C06G118480,Cla97C10G189990
			03013   RNA transport   K00784  Cla97C01G015460,Cla97C09G163910,Cla97C09G182250
			03013   RNA transport   K07936  Cla97C06G120000,Cla97C06G120010,Cla97C06G120020,Cla97C08G156440

		# 脚本3
			
			
			#!/usr/bin/perl -w 
			use strict;

			my ($file) = @ARGV;

			open I, $file or die "$!";

			my %hash;
			while (<I>) {

				chomp;
				
				my ($id) = $_ =~ /^(\d+)/;
				$hash{$id} = [];
			}

			#map {print "$_\n"} keys %hash;

			close I;

			open I1, $file or die "$!";

			while (<I1>) {
				chomp;
				
				my ($id,$gene) = $_ =~ /^(\d+).*\s+(Cla.*)/;
				print $gene."\n";
				push @{$hash{$id}},$gene;
				
			}
			$" = ",";

			open O,">a.txt" or die "$!";
			foreach my $k (keys %hash) {

				my @tmp = @{$hash{$k}};
				
				for (my $i = 0; $i < @tmp; $i++) {
					print O "$k\t$tmp[$i]\n";
					
				}

			}
			
		# 输入文件
			pathwaykO_2_genekn1.txt

		# 结果文件 a.txt
			
			03015   Cla97C05G096170,Cla97C05G103240
			03015   Cla97C01G018000,Cla97C07G144230,Cla97C07G144430,Cla97C10G192490
			03015   Cla97C05G082070,Cla97C11G208000

		# 脚本4


			#!usr/bin/perl -w
			use strict;
			use List::Util qw(uniqstr);
			my ($file) = @ARGV;

			open I,$file or die "$!";

			my %hash;

			while (<I>) {
				chomp;
				my ($id,$gene) = (split /\t/);
				$hash{$id} = [];

			}

			close I;

			open I2,$file or die "$!";

			while (<I2>) {
				chomp;
				
				my ($id,$gene) = $_ =~ /^(\d+)\s+(.*)/;
				if (exists $hash{$id}) {
					push @{$hash{$id}},$gene;
				}

			}

			close I2;


			foreach my $k (keys %hash) {

				my @tmp = @{$hash{$k}};
				my $string = join (",",@tmp);
				#print $string."\n";	
				my @array = split (",",$string);
				@array = uniqstr @array;
				
				for (my $i = 0; $i < @array; $i++) {
					print "$k\t$array[$i]\n";

				}

			}


			# 输入文件

				a.tx
				
			    perl pathway_2_gene3.pl a.txt > pathway_2_gene.txt

			# 结果文件


				03040   Cla97C01G025570
				03040   Cla97C11G224500
				03040   Cla97C07G135540
				03040   Cla97C11G224800
				03040   Cla97C09G173410
				03040   Cla97C11G223280
			
			# 到此为止，我们已经获得了pathway_gene.txt文件，要做富集分析还差一个pathway_name.txt文件

			
		# 获得pathway号对应的name文件


			# 脚本5

			
				#!usr/bin/perl -w
				use strict;

				my ($file) = @ARGV;

				open IN,$file or die "$!";

				my %hash;
				while (my $line = <IN>) {

					my ($id, $name) = $line =~ /^(\d+)\s+(.*)K\w+/;
					$hash{$id} = $name;

				}

				map {print "\"$_\";$hash{$_}\n"} keys %hash;




			# 输入文件

				
				pathwaykO_2_genekn.txt


			 perl pathway_2_name.pl pathwaykO_2_genekn.txt > pathway_2_name.txt

			
			# 结果文件

				"02010";ABC transporters	
				"00430";Taurine and hypotaurine metabolism	
				"03420";Nucleotide excision repair	
				"00604";Glycosphingolipid biosynthesis - ganglio series	


		# 最终，背景文件已经准备好了，后面就可以进行KEGG富集分析	













