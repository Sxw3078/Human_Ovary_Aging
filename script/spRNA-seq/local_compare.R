
suppressPackageStartupMessages(library(optparse))
options(bitmapType='cairo')
option_list <- list(
        make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
        make_option(c("-a", "--aggr"), help="cellranger aggr output directory, default 02_Cellranger/aggr"),
	      make_option(c("-s", "--aggr_csv"), help="csv file list sample information,also used when run cellrange aggr, default 02_Cellranger/all.csv"),
        make_option(c("-c", "--contrast"), help="compare group used when different expresion between sample/group, default pairwise "),
        make_option(c("-b", "--RMbE"), help="remove batch effect by function IntegrateData , default %default",default=TRUE),
        make_option(c("-m", "--mt"), help="mt gene pattern , default %default",default="^MT"),
        make_option(c("-r", "--mt_regress"), help="regree the mt gene in the SCTransform, default %default",default=FALSE),
        make_option(c("-f", "--mt_cutoff"),type="integer", help="mt percent filter cutfoff , default %default",default="10"),
       # make_option(c("-l", "--resolution"),type="double", help="resolution value for FindClusters, default %default",default="1.0"),
        make_option(c("-p", "--npc"),type="integer", help="number of dim used in pca  , default %default",default=30),
        make_option(c("-d", "--spatial_sep"), help="rds for spatial data of separated samples, default %default",default="suerat_Rdata/separated.sct.Rdata"),
        make_option(c("-d", "--spatial_int"), help="rds for spatial data of intergrated samples, default %default",default="suerat_Rdata/separated.sct.rds"),
        make_option(c("-z", "--singcell"), help="rds for singcell data, default %default",default="suerat_Rdata/seuratObject_celltype.rds"),      
        make_option(c("--colors"), help="color set for scatterpie, choose from NPG,NEJM,JAMA,Lancet,AAAS,report,DEF, custom, default %default",default="custom"),        
        make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default=FALSE)

)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(! is.null(opt$workdir)){
  print(paste("change R work directory to ",opt$workdir,sep=""))
  setwd(opt$workdir)
}
project_dir=getwd()

rds_dir=paste(project_dir,"suerat_Rdata",sep="/")
if(!file.exists(rds_dir)){
    dir.create(rds_dir,recursive=TRUE)
}


################################载入R包
################################

suppressPackageStartupMessages(library(SPOTlight))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressWarnings(suppressPackageStartupMessages(library(SeuratData)))
suppressMessages(suppressPackageStartupMessages(library(biomaRt)))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library("BuenColors"))
options(future.globals.maxSize = 10000 * 1024^2)

data=readRDS("../suerat_Rdata/Y18_Theca & stroma_spatial_subcell_locations.mapping.ref.rds")

region1=c('AGCCAAGCTTTGTGTC-1', 'CGAGGGTATCCAGGTG-1', 'CGTGACCAGTCCTCTG-1',
       'CGAATGACGCATAATG-1', 'CTATGCCCGAATGCAA-1', 'TCAATACGCCGTCATG-1',
       'ATACGTTATGCACGGA-1', 'CCTATATTTGTCCTGG-1', 'AAAGGCTACGGACCAT-1',
       'GATCTGCTATCTAAGG-1', 'GCTCAATGTAATACCG-1', 'CACAGGGCCGTTGTCA-1',
       'AGCGGACACTTCGTAG-1', 'GAGTAGATACTAGTTG-1', 'TGCTTCCCAAGCAGTA-1',
       'GACGGGTTGGCCCGTA-1', 'TCCACCTCTAGCCTTT-1', 'CATACACGGTTCCCAC-1',
       'CTTGGCCAAGCTGGGA-1', 'TTCGCACTGTACGACA-1', 'GAAGCCTGCACATTCC-1',
       'TCCAGGGTATATACGA-1', 'AGGAAGCTGTCCGCCG-1', 'ACCGACACATCTCCCA-1',
       'TACCTCAGTTGTCTGT-1', 'GAAACCTATACAAATG-1', 'AAGACCCAACTGAACA-1',
       'ACTACGCGTTAGAATT-1', 'TTACCCTAGGGATTGG-1', 'TAGGTGTTCCACAGAT-1',
       'GGCTTTCAATAAGGGT-1', 'CTCGCCGAATGTAGGG-1', 'GTCTCAAGGCCCGGCT-1',
       'TATACGCGTCATCACT-1', 'AGGCAGGGAGCGTACT-1', 'CGCCCAGCACGCCTAG-1',
       'TATCCGCACCGTCGGG-1', 'CAACTGCTCATCCGAT-1', 'TACCGTAGGTTAACTA-1',
       'CCAGTCTTGTCATAGA-1', 'TACCAAATAGCCCAGA-1', 'TGTGAGACTAGCCCAA-1',
       'GTCGGGAAGCAGAAAC-1', 'GAACCCTCTGTGTTCT-1', 'CCTGGTCGAATGTGGG-1',
       'CCTTTAAGGGAGCACT-1', 'AGTATTTGGCACGACC-1', 'ACCCGTGTCATCAGTA-1',
       'CCGGAGCGTACTTTCT-1', 'TGAGGAGTGCCAGCTT-1', 'GGGCAGAGCAATCGTT-1',
       'TACCTACTCCCAGTAT-1', 'ATACCTAACCAAGAAA-1', 'AATTAACGGATTTCCA-1',
       'ATATAAAGCGCTCGTG-1', 'GGGACCCGTATATCTT-1', 'AGCGCATAATGAATCG-1',
       'AATAGTCCGTCCCGAC-1', 'CCTCGCGCTGTGCGAT-1', 'TGCTCGGCGAAACCCA-1',
       'TGCCGTGGATCGTCCT-1', 'ACCCGGAAACTCCCAG-1', 'GCAAACCTTGGCCATA-1',
       'GCTGCTAAGTAGTCGA-1', 'AGCAGCCAGATGAATA-1', 'TTGTGTATGCCACCAA-1')

region2=c('GAGCACCTGTGTCCAG-1', 'GTCCCGCGACGTTATG-1', 'ATACTGCCTTACACCG-1',
       'GGCGCGGAGATCTTTC-1', 'GCTAACTGAAGTCTGA-1', 'CGTTGAATACCGCGCT-1',
       'ACGGATGGTGCGGATA-1', 'AAGCTCTTTCATGGTG-1', 'GCGCGTCATTGGTACA-1',
       'AAGTTCACTCCAAGCT-1', 'ATCAGCCTCATGCTGC-1', 'GCGCTTAAATAATTGG-1',
       'TCCGCCTGTCTACAAG-1', 'CGTCGTCCTTCGCGAA-1', 'ACCAATATGCAAGTTA-1',
       'CACCCTAACAAGATCT-1', 'CATGGGTATGCCTTAT-1', 'TCTGAAGCACGTGGTC-1',
       'TTGGACCTATAACAGT-1', 'AGTACTCTTATGCCCA-1', 'TTGTGAACCTAATCCG-1',
       'CCATTTCTACCTATTA-1', 'GTTAACATCACTTAAA-1', 'GCCTCTATACATAGCA-1',
       'ATACGTACTTAGCCAC-1', 'TATCCAATTGGTTATC-1', 'CCTATAATGAGTGCCC-1',
       'TACCAATAAAGTACCA-1', 'GGTGGACTGCTCTGGC-1', 'TTGAGAAGTTTAGCAT-1',
       'CTGTGGTCGGGAGATA-1', 'CACGTTTCGTACACAC-1', 'AACCTTTAAATACGGT-1',
       'CTGCACTCCAGTACAG-1', 'CAAGCGGCACATAATT-1', 'ATATAAATGTAGCTGC-1',
       'GTGTCGTATAGCGTTC-1', 'GCTATCATACTCATGG-1', 'TCTCCCTGGGCAGCGT-1',
       'CAATGCGAGAAGTATC-1', 'GCTAGGCACCACGGAG-1', 'TATTCCGAGCTGTTAT-1',
       'GTCCGGACCTGAAATT-1', 'CATGTAAGAGACATTT-1', 'CTAGTTACAACCCGGT-1',
       'AAGCATACTCTCCTGA-1', 'GGTAGAAGACCGCCTG-1', 'AACGATAATGCCGTAG-1',
       'CTGCTTGGCGATAGCT-1', 'TAGCGTCCGGTGTGGT-1', 'AAATAGGGTGCTATTG-1',
       'AGCACTACCTCACCAG-1', 'TCGTCTTAGGCGTTAA-1', 'CAGAACTTAGCCCTCT-1',
       'ACAGGCTTGCCCGACT-1', 'GAGCCAGCTACCTGTG-1', 'TTCGACAGAGCCCGTG-1',
       'GACGACGATCCGCGTT-1', 'GAGATGGGAGTCGACA-1', 'TGGTCCCACGCTACGG-1',
       'ATCGGAGACAGACGGC-1', 'GTTGAGTCCCGCCGGT-1', 'AAGTGTTTGGAGACGG-1',
       'TGCCTTGGCCAGGCAA-1', 'TTCTGACCGGGCTCAA-1', 'AGCTCCTTCGCACATC-1',
       'GCTATACGTCTCGGAC-1', 'TGCCAGTACGTGGAGA-1', 'CTTTAATATTGGTCGA-1',
       'GCTGCTACTGCGTAGC-1', 'ATATTCCACATAGTGA-1', 'ACCGCAATAACTGCCT-1',
       'ACTTCAGGCTGATCCC-1', 'TCTGATGTATTCTGTC-1', 'GCAGATCCATAAGACT-1',
       'TGGCAAACTAAATTAC-1', 'TCGCACCAGGAGGCAG-1')

region3=c('CGAGATTTCGCTCGGG-1', 'AAGGCAGGCTGTCTCC-1', 'ATGAGGGCAGCGGCTA-1',
       'GAGCACGGCGCCTCTT-1', 'ACCTACAGTATGTGGT-1', 'TTCCGGTATCTGTGTC-1',
       'TTCGCTAGGAAGTTGT-1', 'GTTGCGCTAACATTAC-1', 'AGCCCTGTCGCACCGT-1',
       'TGCGCCGTTAATAACG-1', 'AGTCGGTTGCGTGAGA-1', 'GAGGATAAACAGTGCT-1',
       'GGAGGCCGAAGTCGTC-1', 'TAAAGACAACCCTTTA-1', 'GCTGAACTCTCCAGGG-1',
       'AGCGCTAGAGCGATGT-1', 'TGACAACGCATGTCGC-1', 'TTCCAGACGAGATTTA-1',
       'CCTGGAAACGTTCTGC-1', 'GCCGTGGAAGAAATGT-1', 'TCCCAAAGCCCTAAAT-1',
       'TTGAGTCCCGCTGCTG-1', 'GCTAGCACCTGGGCCA-1', 'GATAGGTGTCCCGGGC-1',
       'TTGGTTCGCTCAAAGG-1', 'ATATCGGTAGGGAGAT-1', 'GACCGACGTGAAAGCA-1',
       'CCGGTAATGGCTAGTC-1', 'GTCTTGAGGAGCAGTG-1', 'TTGAGCGCCACGTGAT-1',
       'ATGGAACCTTTGCACA-1', 'CGCCGTCTACCCATCG-1', 'AGGTATGCGGACATTA-1',
       'CCTATGAAGTGGTGCC-1', 'CATACTTAGGCAATAC-1', 'CAATGTGCCAACCCTT-1',
       'CTGGGATACGCTACCC-1', 'ACTGCGGACACACCGT-1', 'GCCCAGATGCTGGAGA-1',
       'ATACGAAGGCTTTCCA-1', 'AGGGACCGGCTGCGTT-1', 'GCTTACGTAGTTAGTA-1',
       'CCTGTCCCTCACGTTA-1', 'GTTAAGTTAGAGTGGG-1', 'AACCTGTCACGGAATT-1',
       'CCGTGAGGCATTCATG-1', 'TCTGGCCGTTCAAGTT-1', 'CCGGCACGACCGTTTC-1',
       'GCAGCCTATATCACAT-1', 'AAATTCCAGGTCCAAA-1', 'TATTGACATTTCTGCC-1',
       'GGCCCGGAGCATGTCT-1', 'TTGGCGATCCGAATAT-1', 'TCTGATTGGAAATGGA-1',
       'ACCTCCGTTATTCACC-1', 'GGTATAGTGACACATA-1', 'TCTTTAGCAGGCGAAC-1',
       'TCTGATCGGGTGCTAG-1', 'GGGCGCAGCGTTACTC-1', 'CCACGTAAATTAGACT-1',
       'ACGTTCGCAATCAATT-1', 'ATTATAGCTACTTTAC-1', 'TTGGTATGGCTTGTGT-1')

region4=c('AAGCCGAAGCGGTTTA-1', 'AATATCGAATCAATGC-1', 'AGAGTAAACTTCACTA-1',
       'CCCGAGTTTCTCCGTA-1', 'GCGCCTCCCACTCCGA-1', 'CGGGAATATAGTATAC-1',
       'ATCAAACACTGTTCCA-1', 'TGAATGTCAGCCGGCC-1', 'GGTTTAGCCTTTCTTG-1',
       'TTGACTATTGTCCGGC-1', 'GCAACAGCAGTATGCG-1', 'TGGGCAATAGTTGGGT-1',
       'CGGCTCTAAAGCTGCA-1', 'AAGTCAATTGTCGTCA-1', 'CATGATGGAAGTTAGC-1',
       'ATCGCCAGTCAACATT-1', 'TCTTCTATAACCCGCC-1', 'TAGCGTCCCTCGATTG-1',
       'GACAAACATATGCAGG-1', 'AGTGAACAAACTTCTC-1', 'AAGTGCCTTGACTGTA-1',
       'ACCGCGGTGGAAGTCG-1', 'CACATTTCTTGTCAGA-1', 'CGCGGCTCAACTTGAA-1',
       'TCTTAGAGTGAACTCT-1', 'TTAAGATAGGATTGAC-1', 'CACCTTGGCGCCTTTG-1',
       'CGGCCCAGGTATATCC-1', 'GGGTTTAGGATAGGAT-1', 'AAAGTTGACTCCCGTA-1',
       'CACGCGGAACTGTTGC-1', 'AGTCCATTGGCTGATG-1', 'ACATCCTGGTAACTGT-1',
       'GCTAGACCGTCTACTG-1', 'GTCCATTACTGCTACG-1', 'GCCTCATCTGGAAATA-1',
       'ACGTATTACTCCGATC-1', 'GCTCGCTCATGTCCAA-1', 'TTGACGCTCCATGAGC-1',
       'TTACATGCCACAACTA-1', 'TATGGTTAGTGGGAGA-1', 'CAAGTGTGGTTGCAAA-1',
       'AACCCTACTGTCAATA-1', 'TCTGGGAACCTTTGAA-1', 'GCGCAAGAGCGCGCTG-1',
       'TATAGATGGTCGCAGT-1', 'ACATGGCGCCAAAGTA-1', 'GATCCGGGAATTAACA-1',
       'GTCGCGTAACCCGTTG-1', 'GTCAAGCGGACTCGGG-1', 'AATCGCCTCAGCGCCA-1',
       'GTGACCGCACACTACG-1', 'CTTGAGTTAGGGTAAT-1', 'GCTGTTGCTACCGAAC-1',
       'TTGAATTCACGTGAGG-1', 'CGCACGTGCGCTATCA-1', 'CGCTAGAGACCGCTGC-1',
       'GCAGACCCAGCACGTA-1', 'TTATATACGCTGTCAC-1', 'AGCTCTAGACGTTCCA-1',
       'CGAGGGACTGCGGTCG-1', 'CTTGTTGCTGAGTCAA-1', 'TTATGATCTTAACGAA-1',
       'CTGGGATAAATAATGG-1', 'TTCAATACTCTGAATC-1', 'AGAAGAGCGCCGTTCC-1',
       'AGTCGACGGTCTCAAG-1', 'GTATGTGGGTCTAGTT-1', 'TTAGCTGATTTGCCGT-1',
       'TATTACCATCCTGCTT-1', 'CCATCTCACCAGTGAA-1', 'ACCCGGATGACGCATC-1',
       'ATAGTTCCACCCACTC-1','TGGACTGTTCGCTCAA-1', 'CCACCCAAGGAAAGTG-1', 
       'GCGATTGTTAACGTTA-1',
       'GGTGATAAGGAGCAGT-1', 'CACGTTCGTGCTCTAG-1', 'GTACAGAGGCAAGGGT-1',
       'CCTGAACGATATATTC-1', 'TACTTGTTAGTAGTCC-1', 'CTGGCGCACAGGTCTG-1',
       'GAAGTCTCCCTAGCGA-1', 'CGAGTTTATCGGACTG-1', 'TGACGATGCACTAGAA-1',
       'AGGGTCGATGCGAACT-1', 'CATCTTACACCACCTC-1', 'CCGATCTCAACCTTAT-1',
       'GAAAGAACAGCGTTAT-1', 'TTCGGGACTAATCGCG-1', 'AACGTGCGAAAGTCTC-1',
       'CCGCACAAAGACCAAC-1', 'ACTCGTCAGTAATCCC-1', 'AAGAGGCATGGATCGC-1',
       'CTATTTGGTTACGGAT-1', 'GGGCCGGCCGAAGTAC-1', 'CCGGGCTGCTCCATAC-1',
       'CCTAGGCGTAGCGATC-1', 'ACTTATACTTACCCGG-1', 'ACCGATGGTAGCATCG-1',
       'CATAACGGACAGTCGT-1', 'TAGGGAGCTTGGGATG-1', 'TATATCCCTGGGAGGA-1',
       'GTGCGACAGGGAGTGT-1', 'ACGATCATCTTGTAAA-1', 'ATAACGCCGGAGGGTC-1',
       'TGAAAGGACCTGACTC-1', 'TGCATATGTCTGTCAC-1', 'AGTGAGACTTCCAGTA-1',
       'GCTTATGAAGCAGGAA-1', 'GGATGTCCTTACCGCA-1', 'AACTCTCAATAGAGCG-1',
       'GCGTGGTACTGGGTTA-1', 'ATATGTCTCCCTAGCC-1', 'TCATTTAAGTCTCCGA-1',
       'TGACATCGAGCGGACC-1', 'ATTGTACAACTCGGCT-1', 'TCCGCGGCAGCATCTG-1',
       'TGTAGGAGAAATTTCC-1', 'CCCAAACATGCTGCTC-1', 'TTCTAACCGAAGCTTA-1',
       'AGGGTGCTCTCGAGGG-1', 'TCTGAATTCCGTACAA-1', 'CGTCGGATAGTGTTGA-1',
       'TCTTTAAGACTATGAA-1', 'GATATTGAGATTGGCG-1', 'GCGTAAATGGCCATAA-1',
       'TACGCTATAGAAACCT-1', 'TAGCCGGCGGTCAGCG-1', 'AGTGATTCAAGCAGGA-1',
       'CACACGCGCTGTCTTA-1', 'GGTTTCAATCGGTCAG-1', 'TTACGGATGGTTCGAG-1',
       'GCTTTCAGAGGAGGTG-1', 'TACCGCGGACTTGCAG-1', 'ATTGATGAGTCCTAAC-1',
       'TAACCTACCGTCCGAG-1', 'CTAGATGTGAGTGTAA-1', 'GTTCATCGTTTGGCTG-1',
       'GTTGGATTGAGAACAC-1', 'TAGACGCCCGTACCGG-1', 'AATCTGCGTTGGGACG-1',
       'CGGCAGGGTCGGGTTG-1', 'TCTTCCCATGGGCACA-1', 'AGAATTATGGATTCGA-1',
       'TAGGTCGCCGGAACTG-1', 'CTTAGTAGGCCTACAG-1', 'ACTCCCGAATTCGTTT-1',
       'ACTTTACCCTCATGAA-1', 'GTAGCTTCCTCTTGTT-1', 'CTACCCTAAGGTCATA-1',
       'GGCTCGCGTTGAGGTA-1', 'TTAAACTCGAATTCAT-1', 'GCTTGGATCGATTAGG-1',
       'CAGACCTGTAAGTGTT-1', 'CTGACTGCGCAGCTCG-1', 'GTAATAAAGGGCTCCC-1',
       'GTGGTTACTTCTTTCG-1', 'TCCCACTCTCTTCCGG-1', 'TCGGGAACGTGCCTAG-1')

region5=c('TTGGACCATCTGGCAA-1', 'GAGCGCAAATACTCCG-1', 'CAGCAGTCCAGACTAT-1',
       'AAACAAGTATCTCCCA-1', 'ACTAGTTGCGATCGTC-1', 'CCCTCCTCGCTCGTAT-1',
       'AGACCCACCGCTGATC-1', 'CCCGCAGCGCGAACTA-1', 'TTCTGCCGCGCCTAGA-1',
       'AGACGAAGTGCCGGTC-1', 'AGGCATTGTCGTAGGG-1', 'TTCTTATCCGCTGGGT-1',
       'GCATCGGCCGTGTAGG-1', 'TCGGCGAACCCAAACC-1', 'AAGCCGAAGCGGTTTA-1',
       'CTCTATTTGGCTGCAG-1', 'ACTCCCTAGAATAGTA-1', 'TCAGACGCTATAGAAG-1',
       'CACTCGAGCTGAACAA-1', 'TAGTCGATCACGGGTT-1', 'CGCTGGTGACTACCCT-1',
       'TGGTAGAATATATGGG-1', 'AGTTTGGCACGGGTTG-1', 'ATGCGAGTCCCACCAC-1',
       'CACATATTAGCAGGAT-1', 'TCGCAAAGATGCATTT-1', 'TCGTGTACTATGGATG-1',
       'TCACAGGGAATCGCAA-1', 'TCCGCTGTCATCCCGG-1', 'TTAGTTATTCGTGGCA-1',
       'ATACGTCCACTCCTGT-1', 'TCTGAACTCGTACCCG-1', 'GTGATCACTAACGCCT-1',
       'TTGACAGGAGCTCCCG-1', 'CGTGGCCGAATATCTA-1', 'TCAAGGTTACTACACC-1',
       'ATTCATATACTGTCCA-1', 'GTCGTACCATCTCGGG-1', 'GTCGTTATTCGCTTAT-1',
       'GCGGAGAGGGAGAACG-1', 'GACTGCACCAGCCCAG-1', 'AGATACTCAAGATCGA-1',
       'TCCTACATCCACGGCC-1', 'TTAGAGTTTAGAAGGA-1', 'CTGGTTTCGAGCAAGA-1',
       'ACTTATTTATGTGCCA-1', 'AACCGCTAAGGGATGC-1', 'CGACACGCTCCGACAG-1',
       'AGGATATAGGGATTTA-1', 'ACCGGGCCTTTGTTGA-1', 'CTCTTCTATTGACTGG-1',
       'TCACCCTCTTAAGATT-1', 'CCTGACCACCGATGGT-1', 'CCGCTATCAGCACCAG-1',
       'GTGGTGGCCAAGTGAA-1', 'ACATCGTATGCAATGG-1', 'AATTGAACGCTCTGGT-1',
       'ATGGTCGCGTGGTTTC-1', 'TTCCGGTTACCCACTT-1', 'CAAGATATTATAACGT-1',
       'ATGGGCCTCGGCCTCT-1', 'TCTTACTTATGCCTCT-1', 'TGCTCCACAGTTCTTA-1',
       'TAAGGCTGAATCCCTC-1', 'TCGAAGAACCGAGCAC-1', 'AAATCGTGTACCACAA-1',
       'TGGCAGATTACGATCA-1', 'CAGGATATATCGTTGT-1', 'CTAAAGGGAAATAGGA-1',
       'CTTTAGTGCTATTATT-1', 'TACGACTGCCTCTTAG-1', 'GCGAAACTTAACTGGA-1',
       'ACAAATGGTAGTGTTT-1', 'TGTTATTGTATGTGGC-1', 'GAGTGTGCGGTACCCA-1',
       'ACACACCAGGACCAGT-1', 'AAGGTGATAAACCAGC-1', 'AAAGTGTGATTTATCT-1',
       'CTGGCTGATTCATCCT-1', 'TCTAGTGATATCGTGG-1', 'GAAATTGTCTCTATAA-1',
       'TTGCACAATTCAGAAA-1', 'ATTAAACATGCGGACC-1', 'CGAGACCCTAGAGTGT-1',
       'ATCGACCCAATACAGA-1', 'AGTTCCTATTTATGTT-1', 'GGAGACGACACCTTTG-1',
       'GCTACGACTTATTGGG-1', 'CGGTACTAGAATCAAA-1', 'AGGACAGTCGAATCCC-1',
       'GGCGAAATCTAACTTG-1', 'TAGCCATTTCAAAGTC-1', 'ATGCCATTTGCGACCA-1',
       'TGTATGGCGCAGACAG-1', 'CGACAATTTGATCTAA-1', 
       'CACAATGAGCTGCTAT-1', 'GGGCCCGTCTTAAACA-1', 'GGCGCATGAATTGATG-1',
       'TCAACGAGGAGACAAA-1', 'CATCGGACGGGTTAAT-1', 'TATCTACCACAGCGGG-1',
       'ACATCGATCGTTTACC-1', 'GAATCTGAACATTCTC-1', 'CAGTCTGTATACTGGG-1',
       'CCTAAATTAACGGTTC-1', 'CTGTGCAGGGTAGGTC-1', 'ACTTGACTCCCTCTTT-1',
       'CTAGATTTACGACGGC-1', 'TGCGTTTGTTGACACT-1', 'ACTCTCTGACTTAGGT-1',
       'GTTTGGCCCAAGTTAT-1', 'CCCAGTTAAGGCGCCG-1', 'GCTTAATGTAACTAAC-1',
       'GACAGCCAGACCTGAC-1', 'CTGGTAACACATAGAA-1', 'GGGATTTACCGCACCT-1',
       'GAATAGACGCGACCCA-1', 'GGATGAAGATCGCTGA-1', 'TACGCAGTTCTTTCCT-1',
       'GGCAAATTACTTTACT-1', 'CGGGCCTTCTTTGTAA-1', 'ATAGTGAAGCGTTCTC-1',
       'GCCGGGTTAGGGTCGC-1', 'GCCGATTGGCCAAGCT-1', 'TTATGAATGAAAGGGA-1',
       'GCGCCGTTCCACGATA-1', 'AGGTTTCACACACCTT-1', 'ACGTTAATGTCGAAGA-1',
       'GATTCGACGGTTCACG-1', 'GCTGCACGGTTTCTTA-1', 'GTGTTACTATGCGTCC-1',
       'GTGAGGACACTTAAGG-1', 'TCTGCCAGAAACTGCA-1', 'CTCGGTACCACTGCTC-1',
       'GTGCGTGTATATGAGC-1', 'CCTCGGACCGGGATAG-1', 'CTTTAGGAACACTGTT-1',
       'AGTGCTTGCACGAATA-1', 'TGAGAGATTTACCACG-1', 'AGCAACATATCTTATT-1',
       'GACCGTGCTGACGGTG-1', 'GGTACAAACATGCTAT-1', 'CGTGAAGTTAATTCAC-1',
       'TACGCCATATTCTAAT-1', 'TACATAGGCATACACC-1', 'CTGCCATGCATCACAT-1',
       'GCTGAGGCGTGAGTAT-1', 'CGCATGGTGCGATGCT-1', 'CAAGGATCGCATGTTC-1',
       'TCCAGAGCACCGGTTC-1', 'GTTTCTGCAGTCTCCC-1', 'CGTGCAGACTGGGACA-1',
       'TCCTCGGGCTGGGCTT-1', 'ATACGCCGGCGAAACC-1', 'TTCTGCGGGTTAGCGG-1',
       'GTAAGTAACAGTCTGG-1', 'ATTTGTCTTGGGAGCT-1', 'TAGGTGCTCGCCTAGC-1',
       'TCGGGCCGTCGTGGTA-1', 'TGCAGTTTCCTCCCAT-1', 'GAAACAGATGACCACC-1',
       'CTGTACTTCTTAGCAT-1', 'TAGTCCGCAGAGAATG-1', 'TACATTTCTAACGTGC-1',
       'CACTCAAGAGCTATGG-1', 'ATCATCCAATATTTGT-1', 'TGGCAGCAGTAATAGT-1',
       'CATACGGCGTCTGGGC-1', 'GCTAGTAGAGCTTGTA-1', 'CGGAGCATGGCGATCC-1',
       'GTAGCGGCTATACACT-1', 'TCTTCGAATAGACGTT-1', 'GATATGAGACACTAAC-1',
       'CGCCGCCCATGCCTGT-1', 'GTGCCCGTTCGGATTC-1', 'CGCACATGTCCACTAC-1',
       'GATAACTCGCACTGTG-1', 'ACTTATTAGGATCGGT-1', 'TTCACGAAAGGATCAC-1',
       'ACCATATCCGCAATAA-1', 'TGTACGAACAAATCCG-1', 'CGCTATTCTTAGGCTC-1',
       'TCACGTGCCCGATTCA-1', 'CACATGATTCAGCAAC-1', 'TGCTGTTGAAGAACTC-1',
       'TAGCGTTGGGTCTTAC-1', 'TAACATACAATGTGGG-1', 'CAATTTCGTATAAGGG-1',
       'CCAGCCTGGACCAATA-1')


color.ls <- list(
    'NPG' = c(
        "#E54B34", "#4CBAD4", "#009F86", "#3B5387", "#F29A7F", "#8491B3", "#91D1C1",
		"#DC0000", "#7E6047", "#CCCCCC", "#BC8B83", "#33ADAD", "#347988", "#9F7685",
		"#C1969A", "#8BB0BB", "#CE8662", "#B04929", "#A59487", "#E3907E", "#D46F5B",
		"#41B4C1", "#278C87", "#726486", "#DA988C", "#88A0B7", "#B9AC91", "#C63517",
		"#927A66", "#DBAEA4", "#97A4AB", "#21A69A", "#3A6688", "#C98882", "#A593A7",
		"#8EC0BE", "#D85935", "#985738", "#B9AFA9", "#E67059", "#E5BFB9", "#B2CED4",
		"#779F99", "#747A87", "#F2DCD5", "#A7ABB3", "#C1D1CD", "#DCA5A5", "#7E7770",
		"#CCCCCC"
    ),
    'NEJM' = c(
		"#BB3B28", "#0072B4", "#E08626", "#1F854D", "#7876B1", "#6E99AC", "#DC91FF",
		"#ED4C97", "#905B6E", "#A07C74", "#928A3C", "#587F7F", "#7487AF", "#A897D5",
		"#E970C9", "#D6435F", "#A86463", "#7C7E95", "#BA9554", "#52826E", "#858BB0",
		"#99A2C1", "#E39AE4", "#E26E95", "#747291", "#C2916E", "#6E8856", "#758298",
		"#8198AE", "#CCAAEA", "#EC82BF", "#C96265", "#BB7B72", "#5A93B4", "#E0B383",
		"#528569", "#9494B1", "#8DA3AC", "#EEC8FF", "#ED9DC2", "#90767F", "#A08E8A",
		"#928E67", "#6C7F7F", "#929BAF", "#BFB6D5", "#E9ACD9", "#D68D9B", "#A88686",
		"#898A95"
    ),
    'JAMA' = c(
		"#374D54", "#DF8E44", "#00A0D4", "#B24645", "#79AE97", "#6A6599", "#7F796B",
		"#8E6E4F", "#A7988E", "#91778A", "#A07F6C", "#748998", "#776E82", "#5C6360",
		"#655E52", "#C6946A", "#6F8CAF", "#AB6558", "#779C98", "#71698D", "#6E6D65",
		"#B67E4A", "#7B9DB1", "#A56167", "#919781", "#6F7798", "#7C7376", "#49585A",
		"#465154", "#DFB792", "#6ABAD4", "#B27C7B", "#93AEA3", "#827F99", "#7F7C75",
		"#8E7E6F", "#A7A09B", "#91848E", "#A09086", "#869198", "#7D7882", "#606362",
		"#65625C", "#C6AD98", "#8F9EAF", "#AB8882", "#8A9C9A", "#7F7B8D", "#6E6E6A",
		"#B69A80"
    ),
    'Lancet' = c(
		"#00468B", "#EC0000", "#41B43F", "#0099B3", "#925E9F", "#FDAE91", "#AC002A",
		"#ACB6B6", "#9E3B4F", "#B2811E", "#40A67D", "#6B7EA9", "#C98599", "#D7675A",
		"#B86E6B", "#6A7BA1", "#6E4D6D", "#D1793E", "#5EAD73", "#6694AE", "#AE80A1",
		"#EAA392", "#B46264", "#949CAB", "#C65356", "#8D9D4B", "#4D9F9B", "#8A7CA4",
		"#E3ABA9", "#C26161", "#B69C9A", "#596D96", "#46688B", "#EC7676", "#7BB47A",
		"#5AA6B3", "#997F9F", "#FDD6C7", "#AC566B", "#B1B6B6", "#9E6D76", "#B29A68",
		"#73A692", "#8A93A9", "#C9A7B1", "#D79F99", "#B89392", "#868EA1", "#6E5E6E",
		"#D1A588"
    )[-5],
    'AAAS' = c(
		"#3A4892", "#ED0000", "#008B45", "#631879", "#00817F", "#BA0020", "#5F549A",
		"#A10055", "#7F807F", "#A93B52", "#9F6D25", "#545A62", "#50557D", "#8A5C4E",
		"#9A3C5D", "#883A77", "#96546A", "#626489", "#7E4372", "#C85113", "#407354",
		"#5C3B7B", "#657167", "#AC293F", "#774989", "#9D3860", "#727285", "#CC2C31",
		"#6E7E35", "#5E3F6D", "#3C6C7E", "#A54137", "#83497B", "#962866", "#8D6B75",
		"#51568D", "#666D92", "#ED7777", "#468B68", "#6E4979", "#418180", "#BA5D6D",
		"#7D779A", "#A1517B", "#808080", "#A9727D", "#9F8662", "#5B5E62", "#67697D",
		"#8A736C"
    ),
    'report' = c(
		"#BF5A17", "#F0017F", "#386CB0", "#FDBF85", "#BEADD3", "#7FC97F", "#FA7F72",
		"#666666", "#B4B1B1", "#D8434F", "#A95597", "#AA949D", "#E1B6AD", "#A2BCAA",
		"#C7A978", "#B1746C", "#8C8A8A", "#C28665", "#CC5036", "#CE3F8A", "#7B7FA7",
		"#EFBB9A", "#B1B5BF", "#A7BA7B", "#D67A6F", "#787878", "#BE9B8A", "#E42F67",
		"#7D63A3", "#D4A992", "#D0B2C1", "#92C394", "#E29675", "#8C6D69", "#A09D9D",
		"#C27140", "#BF8C6B", "#F079B7", "#748EB0", "#FDDEC1", "#C9C0D3", "#A4C9A4",
		"#FABCB6", "#666666", "#B4B3B3", "#D88D93", "#A97FA0", "#AA9FA3", "#E1CCC7",
		"#AFBCB3"
    ),
    'DEF' = c(
		"#FF6600", "#FFFF66", "#009966", "#FF6666", "#666600", "#CCFFCC", "#669933",
		"#339966", "#FFB637", "#9ACB68", "#A78A65", "#B26B39", "#9BAF69", "#99CA7E",
		"#52984E", "#B0893E", "#FFAC57", "#D3E587", "#7A9371", "#D88671", "#83894D",
		"#BEE4B4", "#6C9857", "#849262", "#FFE47A", "#77B27A", "#D49281", "#8C723C",
		"#BDD6A8", "#8BB16F", "#5A986A", "#D99353", "#FFB380", "#FFFFB3", "#4D9980",
		"#FFB3B3", "#666633", "#E6FFE6", "#809966", "#669980", "#FFDB9B", "#B3CB9A",
		"#A79986", "#B28F76", "#A5AF8C", "#B2CAA4", "#759873", "#B09C77", "#FFD6AB",
		"#DCE5B6"
    ),
    'custom'=c(
    "#DC7B1E","#209EBB","#19679A","#C54733","#FDAE91","#548F40","#ACB6B6","#F0E8B3"
    )
)

##设置区域
data2=data
data@meta.data[,"Theca & stroma_1"] = data2$"Theca & stroma_0" + data2$"Theca & stroma_1" + data2$"Theca & stroma_4" + data2$"Theca & stroma_5"
data@meta.data[,"Theca & stroma_2"] = data2$"Theca & stroma_2"
data@meta.data[,"Theca & stroma_3"] = data2$"Theca & stroma_3"
data@meta.data[,"Theca & stroma_4"] = data2$"Theca & stroma_6" + data2$"Theca & stroma_8"
data@meta.data[,"Theca & stroma_5"] = data2$"Theca & stroma_7"
data@meta.data[,"Theca & stroma_6"] = NULL
data@meta.data[,"Theca & stroma_7"] = NULL
data@meta.data[,"Theca & stroma_8"] = NULL
data@meta.data[,"Theca & stroma_0"] = NULL
data@meta.data$regions=data@meta.data$orig.ident


for ( abd_thres in c(0.5,0.6,0.7,0.8,0.9,0.95)){
  data@meta.data$regions='others'
  result_dir=paste0("Theca_stroma_",abd_thres)
  if(!file.exists(result_dir)){
      dir.create(result_dir,recursive=TRUE)
  }
  TSpure=rownames(data@meta.data)[which(data@meta.data$'Theca & stroma'> abd_thres)]
  data@meta.data[intersect(TSpure,region1),"regions"]='region1'
  data@meta.data[intersect(TSpure,region2),"regions"]='region2'
  data@meta.data[intersect(TSpure,region3),"regions"]='region3'
  data@meta.data[intersect(TSpure,region4),"regions"]='region4'
  data@meta.data[intersect(TSpure,region5),"regions"]='region5'
  data@meta.data[which(data@meta.data$regions=='others'),"regions"]=NA

  cols=c(color.ls[["Lancet"]][1:5])
  names(cols)=c("region1","region2","region3","region4","region5")

  plot1=SpatialDimPlot(data,group.by = "regions", label.size = 5)+
      scale_fill_manual(values=cols)+
      labs(title="Y18")+
      theme(plot.title=element_text(color = "black", size = 20,hjust=0.5,vjust=0.5),
          legend.text=element_text(face="plain",size=12),
          legend.title = element_text(size=12,face = "bold",hjust=0.5,vjust=0.5),
          legend.key.width=unit(2,'cm'),legend.key.height=unit(2,'cm'))+
      guides(fill= guide_legend(ncol =1 , title = "Regions"))

  ggsave(plot1, width=12, height=10, filename=paste0(result_dir,"/regions.png"),limitsize = FALSE)
  ggsave(plot1, width=12, height=10, filename=paste0(result_dir,"/regions.pdf"),limitsize = FALSE)
  #提取所需要的表达量矩阵并转化
  mat <- t(data@meta.data[c(region1,region2,region3,region4,region5),c("Theca & stroma_1","Theca & stroma_2","Theca & stroma_3","Theca & stroma_4","Theca & stroma_5")])
  mat = apply(mat,2,function(x){if(sum(x)==0){x}else{scale(x)}})
  rownames(mat)=c("Theca & stroma_1","Theca & stroma_2","Theca & stroma_3","Theca & stroma_4","Theca & stroma_5")
  cluster_info=data@meta.data[c(region1,region2,region3,region4,region5),"regions"]
  names(cluster_info)=c(region1,region2,region3,region4,region5)

  cluster_info <- sort(cluster_info )
  nclusters=length(levels(as.factor(cluster_info)))
  col <- cols
  names(col) <- levels(as.factor(cluster_info))
  gene_features= rownames(mat)
  mat <- as.matrix(mat[gene_features , names(cluster_info)])
  top_anno <- HeatmapAnnotation(
      cluster = anno_block(gp = gpar(fill = col), # 设置填充色
          labels = levels(as.factor(cluster_info)), 
          labels_gp = gpar(cex = 1, col = "white"))) # 设置字体

  gene_pos <- as.factor(rownames(mat) )

  row_anno <-  rowAnnotation(mark_gene = anno_mark(at = 1:length(gene_pos), 
                                      labels = rownames(mat) ))

  library(circlize)
  col_fun = colorRamp2(c(-2,-1, 0,1, 2), c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))

  png(paste0(result_dir,"/region_heatmap.png"),
          width=9,
          height=3.5,
          units="in",
          res=600,
          bg="white")
  plot2=Heatmap(mat,name="Scaled abundance",col=col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_column_names = FALSE,
          show_row_names = FALSE,
          top_annotation = top_anno, # 在热图上边增加注释
          column_title = NULL,
          use_raster=FALSE,
          right_annotation = row_anno,
          column_split = cluster_info,
          raster_by_magick = FALSE,
          show_column_dend = FALSE
          )
  print(plot2)
  dev.off()

  pdf(paste0(result_dir,"/region_heatmap.pdf"),
          width=9,
          height=3.5,
         # units="in",
          #res=600,
          bg="white")
  plot2=Heatmap(mat,name="Scaled abundance",col=col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_column_names = FALSE,
          show_row_names = FALSE,
          top_annotation = top_anno, # 在热图上边增加注释
          column_title = NULL,
          use_raster=FALSE,
          right_annotation = row_anno,
          column_split = cluster_info,
          raster_by_magick = FALSE,
          show_column_dend = FALSE
          )
  print(plot2)
  dev.off()


  data_sub=subset(data,subset= regions!="NA")
  data_sub$regions=factor(data_sub$regions,levels=c("region1","region2","region3","region4","region5"))
  Idents(data_sub)=data_sub$regions
  cluster_markers_all <- FindAllMarkers(object = data_sub, 
                                              assay = "SCT",
                                              slot="data",
                                              verbose = opt$verbose, 
                                              only.pos = TRUE, 
                                              logfc.threshold = 0.25,
                                              min.pct = 0.1,
                                              test.use = "wilcox")   ###还有别的方法

  markers_top10=cluster_markers_all %>% arrange(cluster)  %>% filter(!grepl("^MT-",gene)) %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

  fwrite(cluster_markers_all,paste0(result_dir,"/region_markers_all.xls"),sep="\t",quote=F,row.names=T,col.names=T)

  fwrite(markers_top10,paste0(result_dir,"/region_markers_top10.xls"),sep="\t",quote=F,row.names=T,col.names=T)

  gene_features=markers_top10$gene  ## 直接用自定义的 genes 不用markers_top10

  DefaultAssay(data_sub) <- "SCT"
  data_sub <- ScaleData(data_sub, verbose = FALSE)

  #提取所需要的表达量矩阵并转化
  mat2 <- GetAssayData(data_sub ,assay="SCT", slot = "data")
  #mat2 <- log2(mat2+1)

  cluster_info <- as.factor(sort(data_sub$regions))
  nclusters=length(levels(cluster_info))
  col <-  cols[1:nclusters]
  names(col) <- levels(cluster_info)

  #genes=intersect(rownames(mat2) , gene_features)
  mat2 <- as.matrix(mat2[gene_features, names(cluster_info)])
  mat2 = t(apply(mat2,1,function(x){if(sum(x)==0){x}else{scale(x)}}))
  rownames(mat2)=gene_features
  colnames(mat2)=names(cluster_info)

  top_anno <- HeatmapAnnotation(
      cluster = anno_block(gp = gpar(fill = col), # 设置填充色
          labels = levels(cluster_info), 
          labels_gp = gpar(cex = 1, col = "white"))) # 设置字体

  gene_pos <- 1:nrow(mat2)

  row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                      labels = rownames(mat2) ))

  png(paste0(result_dir,"/markers_heatmap.png"),
          width=11,
          height=10,
          units="in",
          res=600,
          bg="white")
  plot3=Heatmap(mat2,
          cluster_rows = FALSE,name=" ",
          cluster_columns = FALSE,
          show_column_names = FALSE,
          show_row_names = FALSE,
          top_annotation = top_anno, # 在热图上边增加注释
          column_title = NULL,
          use_raster=TRUE,
          right_annotation = row_anno,
          column_split = cluster_info,
          raster_by_magick = TRUE,
          row_dend_width = unit(20, "mm"),
          )
  print(plot3)
  dev.off()

  pdf(paste0(result_dir,"/markers_heatmap.pdf"),
          width=11,
          height=10,
         # units="in",
          #res=600,
          bg="white")
  plot3=Heatmap(mat2,
          cluster_rows = FALSE,name=" ",
          cluster_columns = FALSE,
          show_column_names = FALSE,
          show_row_names = FALSE,
          top_annotation = top_anno, # 在热图上边增加注释
          column_title = NULL,
          use_raster=TRUE,
          right_annotation = row_anno,
          column_split = cluster_info,
          raster_by_magick = TRUE,
          row_dend_width = unit(20, "mm"),
          )
  print(plot3)
  dev.off()

  ############查找markers

  sample_contrast=t(combn(c("region1","region2","region3","region4","region5"),2))

  cluster_sample=as.vector(unique(Idents(data_sub)))
  #nrow=ceiling(length(sample_list[,1])/2)
  for(x in 1:nrow(sample_contrast)){
    
    diff.cluster_1=FindMarkers(data_sub, ident.1 = sample_contrast[x,1], ident.2 = sample_contrast[x,2],only.pos=TRUE, verbose = opt$verbose) %>%  filter(p_val_adj< 0.05 )
    if(nrow(diff.cluster_1)==0){
        system(paste0("echo 'no diff gene was found for this cluster' >  ",paste(result_dir,paste0(sample_contrast[x,2],"_more_then_",sample_contrast[x,1],"_diff_gene.xls"),sep="/"),sep=""))
    }else{
        fwrite(diff.cluster_1,paste(result_dir,paste0(sample_contrast[x,1],"_more_then_",sample_contrast[x,2],"_diff_gene.xls"),sep="/"),sep="\t",quote=F,row.names=T,col.names=T)	
    }

    diff.cluster_2=FindMarkers(data_sub, ident.1 = sample_contrast[x,2], ident.2 = sample_contrast[x,1],only.pos=TRUE, verbose = opt$verbose) %>%  filter(p_val_adj< 0.05 )
    if(nrow(diff.cluster_2)==0){
        system(paste0("echo 'no diff gene was found for this cluster' >  ",paste(result_dir,paste0(sample_contrast[x,2],"_more_then_",sample_contrast[x,1],"_diff_gene.xls"),sep="/"),sep=""))
    }else{
        fwrite(diff.cluster_2,paste(result_dir,paste0(sample_contrast[x,2],"_more_then_",sample_contrast[x,1],"_diff_gene.xls"),sep="/"),sep="\t",quote=F,row.names=T,col.names=T)	
    }
  }

  regions=c("region1","region2","region3","region4","region5")
  for(y in regions){

    diff.cluster = FindMarkers(data_sub, ident.1 = y, ident.2 = setdiff(regions, y), only.pos=TRUE, verbose = opt$verbose) %>%  filter(p_val_adj< 0.05 )
    
    if(nrow(diff.cluster)==0){
        system(paste0("echo 'no diff gene was found for this cluster' >  ",paste(result_dir,paste0(y,"_up_expressed_gene.xls"),sep="/"),sep=""))
    }else{
        fwrite(diff.cluster,paste(result_dir,paste0(y,"_up_expressed_gene.xls"),sep="/"),sep="\t",quote=F,row.names=T,col.names=T)	
    }
  }
}

#ls Theca_stroma_*/*_diff_gene.xls | cut -d "_" -f1-6 |xargs -P70 -i sh -c "perl ~/Git_local/Gene_enrichment/Gene_enrichment.pl -Annotation ~/Git_local/Gene_enrichment/Annotation.xls -l {}_diff_gene.xls -org hsa -o {}_Gene_enrichment > {}_Gene_enrichment.log 2>&1"

#ls Theca_stroma_*/*_diff_gene.xls | cut -d "_" -f1-6 |xargs -P70 -i sh -c "rm -r  {}_Gene_enrichment/TMP {}_Gene_enrichment/KEGG_Pathway_map"
