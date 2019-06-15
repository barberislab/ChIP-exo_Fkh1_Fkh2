library(readxl)
library(pathview)

# start in the working directory this script is located in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# move the the R_pathview figures folder
setwd('../Figures/R_pathview/')
getwd()

# List of KEGG pathways to visualize
pathways = c('sce04111','sce04113','sce04138','sce04011','sce04144','sce04139',
             'sce00240','sce00230','sce00564','sce00260',
             'sce00270','sce00480','sce00220',
             'sce00520','sce00970','sce00630',
             'sce00010','sce00020','sce00190','sce00030','sce00620',
             'sce04141','sce03013','sce03008','sce00510')
pathway_names = c('Cell cycle','Meiosis','Autophagy','MAPK signaling pathway','Endocytosis','Mitophagy',
                  'Pyrimidine metabolism', 'Purine metabolism','Glycerophospholipid metabolism','Glycine, serine and threonine metabolism',
                  'Cysteine and methionine metabolism','Glutathione metabolism','Arginine biosynthesis',
                  'Amino sugar and nucleotide sugar metabolism','Aminoacyl-tRNA biosynthesis','Glyoxylate and dicarboxylate metabolism',
                  'Glycolysis','TCA cycle','Oxidative phosphorylation','Pentose phosphate pathway','Pyruvate metabolism',
                  'Protein processing ER','RNA transport','Ribosome biogenesis','N-Glycan biosynthesis')
pathways_bottomright = c('Cell cycle','Pyruvate metabolism','Ribosome biogenesis')

# Make two sets of images: 
# 1. only our suggested targets 
# 2. Combined with literature
# (-1) Only shown with ChIP-chip
# (0) Shown by ChIP-exo alone
# (1) Shown by ChIP-exo + ChIP-chip
# The second case is done through a specific column 'KEGG_pathview' in the annotated Excel file

for (i in seq(length(pathways))) {
  pathway = pathways[i]
  pathway_name = pathway_names[i]
  
  if (is.element(pathway_name,pathways_bottomright)) {
    keypos = "bottomright"
  } else {
    keypos = "topright"
  }

  # run once for consensus view and once for our data only
  for (case in c('Mondeel_2019','consensus')) { 
    
    if (case == 'Mondeel_2019') {
      print(case)
      
      # make case dir
      dir.create('./Mondeel2019') 
      setwd("./Mondeel2019")
      print(getwd())
      
      # make folder and set as working dir
      dirname = paste(pathway_name,pathway,sep="_")
      dir.create(dirname)
      setwd(dirname)
      print(getwd())
      
      for (experiment in c('log', 'stat')) {
        print(experiment)
        
        # read from the /../pathway_name folder
        df = read_excel('../../../../Tables/pathview_targets.xlsx') # a tibble
        df = as.data.frame(df)
      
        exp_data = df[,c(paste("Target Fkh1",experiment), paste("Target Fkh2",experiment))]

        rownames(exp_data) = df[,1]

        # make image
        pv.out <- pathview(gene.data = exp_data, gene.idtype = "KEGG",
                           pathway.id = pathway, species = "sce", out.suffix = experiment, keys.align = "y", 
                           kegg.native = T, key.pos = keypos, map.symbol=T, same.layer=F,
                           discrete=list(gene=T, cpd=T), bins=list(gene=max(exp_data,na.rm=T),cpd=1), 
                           limit=list(gene=max(exp_data, na.rm=T), cpd=1), both.dirs=list(gene=F,cpd=T),
                           new.signature=F, multi.state=T, node.sum="max",
                           mid=list(gene='white',cpd='white'), high=list(gene='green',cpd='green') )
          
      } # end experiment for
      
      setwd('../../')
      print(getwd())
      
      } else { # end case if
        experiment = 'log'
        
        # make case dir
        dir.create('./Consensus')
        setwd("./Consensus")
        print(getwd())
        
        # make folder and set as working dir
        dirname = paste(pathway_name,pathway,sep="_")
        dir.create(dirname)
        setwd(dirname)
          
        # read from the /../pathway_name folder
        df = read_excel('../../../../Tables/pathview_targets.xlsx') # a tibble
        df = as.data.frame(df)
        
        exp_data = df[,c(paste("Consensus Fkh1",experiment), paste("Consensus Fkh2",experiment))]
        
        rownames(exp_data) = df[,1]

        # bidirectional
        pv.out <- pathview(gene.data = exp_data, gene.idtype = "KEGG", 
                           pathway.id = pathway, species = "sce", out.suffix = 'Fkh1,2', keys.align = "y", 
                           kegg.native = T, key.pos = keypos, map.symbol=T, same.layer=F,
                           discrete=list(gene=T, cpd=F), new.signature=F,
                           low=list(gene='yellow',cpd='yellow'), mid=list(gene='red',cpd='red'),
                           high=list(gene='green',cpd='green'), multi.state=T, both.dirs = list(gene=T,cpd=T), 
                           limit=list(gene=max(exp_data,na.rm=T),cpd=1), node.sum="max")
        
        # reset working dir
        setwd("../../")
        print(getwd())
      } # end case else 
    } # end case for
} # end pathway for

