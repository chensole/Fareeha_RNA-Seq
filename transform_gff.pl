
#在进行西瓜转录组分析的时候，运用hisat2构建索引时，对于gff文件有一定的要求，第九列要注明gene_id transcript_id;

#然而西瓜gff文件不符合这种形式

#     gene          ID=Cla97C01G000020;Name=EVM%20prediction%20ClaScf_0011.1317
#     exon          Parent=Cla97C01G000020.1

#对于上面这种形式的gff文件，我用hisat提取剪切位点一直返回空文件

#>>> hisat2_extract_splice_sites.py watermelon_97103_v2.gene_model.gff >test

# -rw-rw-r-- 1 chenzhi chenzhi          0 7月   5 13:59 test

#于是查看了hisat2_extract_splice_sites.py的源代码发现我这个gff文件不规范，因此写了下面这个脚本来转换gff格式

# google搜索

	#1. gffread watermelon_97103_v2.gene_model.gff -o watermelon_97103_v2.gene_model.gtf


#     gene	    ID=Cla97C01G000010.1;geneID=Cla97C01G000010;gene_name=EVM%20prediction%
#     exon	    Parent=Cla97C01G000010.1



	#2. 发现有点改变，不过还是不行(需要下面这种形式才能用hisat2构建索引)

# gene   transcript_id "Cla97C01G000010.1";gene_id "Cla97C01G000010";gene_name "EVM%20prediction%20ClaScf_0011.1317";
# exon   transcript_id "Cla97C01G000010.1";gene_id "Cla97C01G000010";


# 下面这脚本基于gffread转换后的文件形式

#!/usr/bin/perl -w
use strict;

my ($file) = @ARGV;
open IN, "$file" or die "$!";
use 5.010;
while (my $line = <IN>) {
	next if ($line =~ /^#/);
	chomp $line;
	state $id;				# state声明变量很有用，下次循环的过程会记住上一次的变量值
	if ($line =~ /(geneID=[^;]+)/) {
		$id = $1;
		print $line,"\n";
		#print "$id\n";
		#last;
	}else {
		my @tmp = split(/\t/,$line);
		$tmp[8] .= ";$id";
		my $newline = join ("\t",@tmp);
		print "$newline\n";
		
	
	}



}

#上述代码运行后结果

#     gene          transcript_id=Cla97C01G000010.1;gene_id=Cla97C01G000010;gene_name=EVM%2
#     exon          transcript_id=Cla97C01G000010.1;gene_id=Cla97C01G000010



# 大致相同了，还需要将 = 转换成 "",于是又写了一个脚本


#!usr/bin/perl -w
use strict;

my ($file) = @ARGV;
open IN,"$file" or die "$!";

while (my $line = <IN>) {
	chomp $line;
	my @left = (split /\t/,$line)[ 0 .. 7];
	my $tmp = (split /\t/,$line)[8];
	my $lf = join("\t",@left);
	my @a = split (";",$tmp);
	my @b = map { { (split /=/)[0] => (split /=/)[1] } } @a;
	
	{my $new = "";
		foreach (@b) {
			my %hash = %{$_};
			my ($k) = keys %hash;
			my $v = $hash{$k};
		
			$new .= "$k \"$v\";";
		
		
			}
	print "$lf\t$new\n";
}
}

## 经过上面两步的操作终于可以了，看来写脚本还是很有用的

