package primer_blast;
use feature ':all';
use strict;
use warnings FATAL => 'all';
use autodie;
use Bio::SearchIO;
use Bio::SeqIO;

sub primer_blast{
    my ($blast_xml,$exists_probe,$target_fna) = @_;


    # 读入测序序列集，ID-序列
    my %fna_id2seq;
    # 创建一个Bio::SeqIO对象来读取FASTA文件
    my $seqio = Bio::SeqIO->new(-file =>$target_fna, -format => 'fasta');

    # 遍历文件中的每个序列
    while (my $seq = $seqio->next_seq) {
        my $now_id = $seq->id;
        my $now_seq = $seq->seq;
        $fna_id2seq{$now_id} = $now_seq;
    }

    my %statistics;
    # 引物产物判断范围
    my $product_range = 1000;
    $exists_probe ||= 1;


    # format can be 'fasta', 'blast', 'exonerate', ...
    my $searchio = Bio::SearchIO->new( -format => 'blastxml', -file   => $blast_xml );

    # hit_seq_info存储命中序列ID对应序列名字
    my (%primer_seq,%hit_seq_info,%primer_info);
    my ($next_primer_forward,$reverse_primer_forward);
    my (%all_hit_num);
    # 遍历搜索结果
    # open my $out_file,">",$out;
    while(my $result = $searchio->next_result) {
        my $algorithm_type = $result->algorithm;            # BlastOutput_program,blast软件
        my $algorithm_version = $result->algorithm_version; # BlastOutput_version,blast软件+版本

        my $query_name = $result->query_description; # Iteration_query-def,查询序列名字

        # next if ($query_name =~ /probe/);
        my $query_seq_name = $query_name =~ s/_.*$//r;

        my $primer_forward = $1 if ($query_name =~ /.*_([FR])/); # 引物的方向

        # 存储引物的本体名字
        $primer_info{$query_name} = $query_seq_name;

        # unless (defined $next_primer_forward) {
        #     # 检查第一个获得的引物方向
        #     $next_primer_forward = $primer_forward;
        #     $reverse_primer_forward = $next_primer_forward =~ tr/FR/RF/r;
        # }

        sub primer_info{
            my ($query_name,$bit_score,$evalue,$identical_length,$total_length,$identity,$gaps,$mismatch) = @_;
            my $c1 = $query_name;
            my $c2 = "Score = $bit_score bits, Expect = $evalue";
            my $c3 = "Identities = $identical_length/$total_length ($identity%), Gaps = $gaps ($mismatch mismatches)";
            return($c1,$c2,$c3);
        }

        while (my $hit = $result->next_hit) {
            # process the Bio::Search::Hit::HitI object
            while (my $hsp = $hit->next_hsp) {
                #print Dumper($hsp);
                # process the Bio::Search::HSP::HSPI object

                my $query_string = $hsp->query_string; # Hsp_qseq,查询序列
                my $hit_string = $hsp->hit_string;     # Hsp_hseq,命中序列
                my $align = $hsp->homology_string;     # Hsp_midline,比对提示

                my $hit_name = $hit->name;               # Hit_id,命中序列id
                my $hit_description = $hit->description; # Hit_def,命中序列名字


                my $identity = sprintf "%0.4f", $hsp->frac_identical; # Identities,计算相似度，保留两位小数
                $identity = $identity * 100;

                my $mismatch = () = $hsp->seq_inds('hit', 'nomatch'); # Mismatches数量
                my $gaps = $hsp->gaps;                                # Hsp_gaps，gap数量

                # my @hit_range = $hsp->range('hit');
                # say join "\t",@hit_range;
                my $query_start = $hsp->start('query'); # Hsp_query-from,查询序列起始位置
                my $query_end = $hsp->end('query');     # Hsp_query-to,查询序列结束位置
                my $hit_start = $hsp->start('hit');     # Hsp_hit-from,命中序列起始位置
                my $hit_end = $hsp->end('hit');         # Hsp_hit-to,命中序列结束位置
                # say $hit_start."\t".$hit_end;
                my $evalue = $hsp->evalue;              # Hsp_evalue,Evalue值
                my $bit_score = $hsp->bits;             # Hsp_bit-score,bit-score得分

                $hit_seq_info{$hit_name} = $hit_description;

                my $total_length = $hsp->length('total');
                my $identity_hsp = $hsp->num_identical; # 匹配碱基数

                # 猜测这个是判断序列方向的
                # Hsp_hit-frame的值表示目标序列在HSP（High-scoring Segment Pair，高分段对）中使用的阅读框。其值范围通常为-3到3，其中：
                # -3 到 -1 表示目标序列的反义链（互补链）的阅读框；
                # 0 表示目标序列与查询序列没有相位差异，通常用于非编码区域；
                # 1 到 3 表示目标序列的正义链的阅读框。
                my $query_frame = $hsp->strand('query'); # Hsp_query-frame
                my $hit_frame = $hsp->strand('hit');     # Hsp_hit-frame


                sub reverse_from_to{
                    # 反转匹配位置
                    my ($from,$to) = @_;
                    return ($to,$from);
                }

                sub reverse_seq {
                    # 获得反向互补序列
                    my $seq = shift;
                    $seq =~ tr/ACGTacgt/TGCAtgca/;
                    $seq = reverse $seq;
                    return $seq;
                }

                if ($hit_frame < 0){
                    # 由于引物的方向都是5->3，所以匹配如果是负链，则需要纠正延伸方向，原始的hit_start、hit_end有问题，总是从小到大
                    ($hit_start,$hit_end)=reverse_from_to($hit_start,$hit_end);
                    # 纠正反向引物的序列
                    $query_string = reverse_seq($query_string);
                    $hit_string = reverse_seq($hit_string);
                    $align = reverse $align;
                }
                # say $hit_name."\t".$query_name;
                $primer_seq{$hit_name}{$query_name}=
                    {
                        # query_name   => $query_name,
                        start        => $hit_start, end => $hit_end,
                        query_seq    => $query_string, hit_seq => $hit_string, align => $align,
                        total_length => $total_length, identity_length => $identity_hsp,
                        bit_score    => $bit_score, evalue => $evalue, identity => $identity, gaps => $gaps, mismatch => $mismatch
                    };
            }
        }
    }

    # 判断引物是否相向并计算扩增片段大小
    sub analyze_primers {
        my ($primer1,$primer2) = @_;

        # say $primer1->{start};
        # say $primer1->{end};
        # say $primer2->{start};
        # say $primer2->{end};

        # 判断引物是否相向
        my $are_opposite = 0;
        if (($primer1->{start} < $primer1->{end} && $primer2->{start} > $primer2->{end}) ||
            ($primer1->{start} > $primer1->{end} && $primer2->{start} < $primer2->{end})) {
            $are_opposite = 1;
        }

        # 计算扩增片段大小
        my $fragment_size;
        if ($are_opposite) {
            # 引物相向，扩增片段大小是两个引物起始位置之差的绝对值
            $fragment_size = abs($primer1->{start} - $primer2->{start});
            # say $primer1->{start}."\t".$primer2->{start};
        } else {
            # 引物相反
            $fragment_size = 0;
        }

        return ($are_opposite,$fragment_size);
    }


    my $amp_num = 0;

    my $now_hit_name = "";
    # 输出哈希测试
    while (my($hit_name, $v1) = each %primer_seq) {
        our %primers = %$v1;
        # 如果没有探针且键少于2，有探针键值少于3，直接跳过
        next if(($exists_probe <= 0 && keys %primers < 2) or ($exists_probe > 0 && keys %primers < 3));
        $now_hit_name = $hit_name;
        my $judge_primer = 0;   # 分析一个引物即可
        my ($opposite, $fragment_size);
        my ($primer_1,$primer_2,$probe);
        while (my($query_name, $value) = each %primers) {
            # say $hit_name."\t".$query_name;
            if ($query_name=~/primer/ && $judge_primer == 0){   # 只要知道一个引物知道了第二个引物
                my $primer_name = $primer_info{$query_name};
                $primer_1 = $primer_name."_F";
                $primer_2 = $primer_name."_R";
                # 分析两对引物
                ($opposite, $fragment_size) = analyze_primers($primers{$primer_1}, $primers{$primer_2});
                $judge_primer = 1;
            }
            $probe = $primer_info{$query_name} if ($query_name !~ /primer/ && $exists_probe > 0);
        }
        # if ($opposite == 1) {
        #     print "引物是相向的。";
        # } else {
        #     print "引物是相反的。";
        # }
        #
        # if ($fragment_size > 1000) {
        #     print "扩增片段大于1000bp，大小为 ${fragment_size}bp。\n";
        # } else {
        #     print "扩增片段小于或等于1000bp，大小为 ${fragment_size}bp。\n";
        # }

        # say $primers{$primer_1}->{start}."\t".$primers{$primer_2}->{start};

        # 检查探针是否在引物之间
        my %product_size;
        # say $probe."\t".$primers{$probe}->{start}."\t".$primers{$probe}->{end};
        my $primer_1_end = $primers{$primer_1}->{end};
        my $primer_2_end = $primers{$primer_2}->{end};

        my $primer_1_start = $primers{$primer_1}->{start};
        my $primer_2_start = $primers{$primer_2}->{start};

        if($primer_1_end > $primer_2_end){
            ($primer_2_end,$primer_1_end) = ($primer_1_end,$primer_2_end);
            ($primer_2_start,$primer_1_start) = ($primer_1_start,$primer_2_start);
            ($primer_2,$primer_1) = ($primer_1,$primer_2);
        }

        for ($primer_1_end..$primer_2_end){
            $product_size{$_} = 1;
        }

        # if (exists($product_size{$primers{$probe}->{start}}) && exists($product_size{$primers{$probe}->{end}})){
        #     print "探针在引物之间。";
        # } else {
        #     print "探针不在引物之间。";
        # }

        # 总结
        if ($opposite == 1 && $fragment_size < $product_range && exists($product_size{$primers{$probe}->{start}}) && exists($product_size{$primers{$probe}->{end}})){
            # say "引物是相向的，扩增片段小于${product_range}bp，探针在引物之间。";

            say ">AMP_".$amp_num." ".$now_hit_name.":".$primer_1_start."-".$primer_2_end." ".$hit_seq_info{$now_hit_name};
            $amp_num += 1;

            # 输出引物基本信息
            my ($f1,$f2,$f3)=&primer_info(
                $primer_1,$primers{$primer_1}->{bit_score},$primers{$primer_1}->{evalue},
                $primers{$primer_1}->{identity_length},$primers{$primer_1}->{total_length},
                $primers{$primer_1}->{identity},$primers{$primer_1}->{gaps},$primers{$primer_1}->{mismatch}
                            );

            my ($r1,$r2,$r3)=&primer_info(
                $primer_2,$primers{$primer_2}->{bit_score},$primers{$primer_2}->{evalue},
                $primers{$primer_2}->{identity_length},$primers{$primer_2}->{total_length},
                $primers{$primer_2}->{identity},$primers{$primer_2}->{gaps},$primers{$primer_2}->{mismatch}
                            );

            my ($p1,$p2,$p3)=&primer_info(
                $probe,$primers{$probe}->{bit_score},$primers{$probe}->{evalue},
                $primers{$probe}->{identity_length},$primers{$probe}->{total_length},
                $primers{$probe}->{identity},$primers{$probe}->{gaps},$primers{$probe}->{mismatch}
                            );

            # 打印引物信息表格
            {
                my @table = (
                    [ $f1, $p1, $r1 ],
                    [ $f2, $p2, $r2 ],
                    [ $f3, $p3, $r3 ]
                );

                # 计算每列的最大宽度
                my @col_widths;
                for my $row (@table) {
                    for my $i (0 .. $#$row) {
                        my $len = length($row->[$i]);
                        if (!defined($col_widths[$i]) || $len > $col_widths[$i]) {
                            $col_widths[$i] = $len;
                        }
                    }
                }

                # 打印表格
                print "\n";
                for my $row (@table) {
                    for my $i (0 .. $#$row) {
                        print " | " unless $i == 0;
                        printf("%-*s ", $col_widths[$i],$row->[$i]);
                    }
                    print "\n";
                }
                print "\n";
            }

            say ">> ".$primer_1;
            # 输出正向引物\正向探针（如果存在）
            my $interval_bp_1 = $primer_2_start - $primer_1_end - 1;
            my $primer_1_inter = $primers{$primer_1}->{total_length} - 2;

            # 正向引物长度信息
            say " " x 3 . "1"." " x $primer_1_inter . $primers{$primer_1}->{total_length};
            # 匹配正向引物
            my $output_primer1 = "5' ".$primers{$primer_1}->{query_seq};
            # 匹配序列匹配线及匹配序列末端位置
            my $output_primer1_align = " " x 3 . $primers{$primer_1}->{align};

            # 正向探针
            if ($primers{$probe}->{start} < $primers{$probe}->{end}){
                my $probe_primer1_inter = $primers{$probe}->{start} - $primers{$primer_1}->{end} - 1;

                $output_primer1 .= " " x $probe_primer1_inter . $primers{$probe}->{query_seq};
                my $align_inter = $probe_primer1_inter;
                $output_primer1_align .= " " x $align_inter . $primers{$probe}->{align};

                $interval_bp_1 = $primers{$primer_2}->{start} - $primers{$probe}->{end} - 1;
            }

            $output_primer1_align .= " " x $interval_bp_1 . $primer_2_start;
            say $output_primer1;
            say $output_primer1_align;


            # 输出匹配序列段
            print " " x 3;
            say substr $fna_id2seq{$now_hit_name},$primer_1_start-1,$fragment_size+1;

            # 输出反向引物\反向探针
            my ($output_primer2,$interval_bp_2_align,$interval_bp_2);
            my $output_primer2_align = " " x 3 . $primer_1_start;

            # 反向探针
            if ($primers{$probe}->{start} > $primers{$probe}->{end}){
                my $inter_probe = $primers{$probe}->{end} - $primers{$primer_1}->{start} + 3;
                my $probe_true = $primers{$probe}->{query_seq} =~ tr/ACGTacgt/TGCAtgca/r;
                $output_primer2 = " " x $inter_probe . $probe_true;
                $inter_probe -= length($output_primer2_align);
                $output_primer2_align .= " " x $inter_probe . $primers{$probe}->{align};

                $interval_bp_2_align = $primers{$primer_2}->{start} - $primers{$probe}->{end} - length($output_primer2_align) + 2;
                $interval_bp_2 = $primers{$primer_2}->{end} - $primers{$probe}->{start} - 1;
            }else{
                $interval_bp_2_align = $primer_2_end - $primer_1_start - length($primer_1_start);
                $interval_bp_2 = $primer_2_end - $primer_1_start + 3;
            }

            my $primer_2_seq = $primers{$primer_2}->{query_seq} =~ tr/ACGTacgt/TGCAtgca/r;
            $output_primer2_align .=  " " x $interval_bp_2_align . $primers{$primer_2}->{align};
            $output_primer2 .= " " x $interval_bp_2 . $primer_2_seq . " 5'";;

            say $output_primer2_align;
            say $output_primer2;
        }

        say "\n";
    }
}

1;