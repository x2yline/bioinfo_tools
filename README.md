# bioinfo_tools（python3）
一键完成某项任务的脚本(重造轮子)

### 1. kegg富集分析脚本

必备文件为geneid2symbol, hsa00001.keg, diff_gene.txt(可以为symbo或geneid)

放置同一文件夹后一键运行python kegg.py即可产生kegg_enrichment_result文件夹, 含有csv和png文件

kegg富集结果图如下:![富集结果](https://raw.githubusercontent.com/x2yline/bioinfo_tools/master/kegg/kegg_enrichment_result/enrichment.png)
参考：
http://www.biotrainee.com/thread-1347-1-1.html
http://www.cs.tau.ac.il/~rshamir/ge/09/scribe/lec14a.pdf
https://guangchuangyu.github.io/cn/2012/04/enrichment-analysis/

### 2. kegg通路图动画脚本

输入通路id与要标注的gene_id列表，对应color列表与对应的label列表后可输出通路图，注意label_list可以为None，列表长度相同

通路demo动画图如下:![通路图](https://raw.githubusercontent.com/x2yline/bioinfo_tools/master/kegg_animation/hsa03060.gif)

参考：
https://github.com/ajmazurie/kegg-animate-pathway
