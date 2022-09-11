# We define the following helper function in order to load or install the packages according to the condition


usePackage <- function(p) {
  if (!is.element(p, installed.packages()[, 1])) {
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

# This one below is to the the default CRAN repo

r <- getOption("repos")
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r)
rm(r)
 
 

usePackage("rotl")
usePackage("readr")
usePackage("dplyr")
usePackage("igraph")
usePackage("plyr")
usePackage("gridExtra")
usePackage("ggfortify")
usePackage("ggtree")
usePackage("plotly")
usePackage("archive")


##### Session Infos

> sessionInfo()

# R version 4.1.3 (2022-03-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur/Monterey 10.16

# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] archive_1.1.5    plotly_4.10.0    ggtree_3.2.1     ggfortify_0.4.14
#  [5] ggplot2_3.3.6    gridExtra_2.3    plyr_1.8.7       igraph_1.3.4    
#  [9] dplyr_1.0.10     readr_2.1.2      rotl_3.0.12     

# loaded via a namespace (and not attached):
#  [1] treeio_1.18.1      progress_1.2.2     tidyselect_1.1.2   purrr_0.3.4       
#  [5] lattice_0.20-45    ggfun_0.0.5        colorspace_2.0-3   vctrs_0.4.1       
#  [9] generics_0.1.3     viridisLite_0.4.1  htmltools_0.5.3    utf8_1.2.2        
# [13] XML_3.99-0.10      gridGraphics_0.5-1 rlang_1.0.5        pillar_1.8.1      
# [17] glue_1.6.2         withr_2.5.0        DBI_1.1.3          rentrez_1.2.3     
# [21] lifecycle_1.0.1    stringr_1.4.1      munsell_0.5.0      gtable_0.3.1      
# [25] htmlwidgets_1.5.4  fastmap_1.1.0      tzdb_0.2.0         parallel_4.1.3    
# [29] fansi_1.0.3        Rcpp_1.0.9         scales_1.2.1       jsonlite_1.8.0    
# [33] digest_0.6.29      hms_1.1.2          aplot_0.1.2        rncl_0.8.6        
# [37] stringi_1.7.8      grid_4.1.3         cli_3.3.0          tools_4.1.3       
# [41] yulab.utils_0.0.4  magrittr_2.0.3     lazyeval_0.2.2     patchwork_1.1.1   
# [45] tibble_3.1.8       crayon_1.5.1       ape_5.6-2          tidyr_1.2.0       
# [49] pkgconfig_2.0.3    tidytree_0.3.9     ellipsis_0.3.2     data.table_1.14.2 
# [53] ggplotify_0.1.0    prettyunits_1.1.1  assertthat_0.2.1   httr_1.4.4        
# [57] R6_2.5.1           nlme_3.1-155       compiler_4.1.3    

  ###############  

# The full taxo is downloaded from OTL 
# https://tree.opentreeoflife.org/about/taxonomy-version/ott3.3

# File can bve directly downloaded 
download.file('http://files.opentreeoflife.org/ott/ott3.3/ott3.3.tgz', destfile = "docs/data/inputs/ott3.3.tgz", method = "wget", extra = "-r -p --random-wait")

# and opened from the tar archive

taxonomy <- read_delim(archive_read("docs/data/inputs/ott3.3.tgz", file = "ott3.3/taxonomy.tsv"), col_types = cols(), delim = "|", escape_double = FALSE, trim_ws = TRUE )


PF_d <- read_delim("docs/data/inputs/pf_metadata_otoled.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)



#taxonomy <- taxonomy[2:1000, ]

taxonomy %>% 
  filter(rank =="domain")

taxonomy %>% 
  filter(name =="Asparagales")


# Here we use I graph to create a graph from the datatable of OTL tax input

sub_tax_from_to <- data.frame(from=taxonomy$parent_uid,
                                to=taxonomy$uid)


g <- graph.data.frame(sub_tax_from_to, directed=TRUE)


#plot(g)

print_all(g)

# And the subcomponent to filter for all node below Archaeplastida https://igraph.org/r/doc/subcomponent.html


sub_g <- subcomponent(g, '10210', mode = 'out') ## tracheophyta ott id

g2 <- induced_subgraph(g, sub_g)

id_sel <- as.numeric(V(g2)$name)

taxonomy_final <- taxonomy[taxonomy$uid %in% id_sel,]

taxonomy_family <- taxonomy_final %>% 
                       filter(rank == "family" & is.na(flags)) 
                       
# Ucomment the two following line if you need to launch tax_lineage() again

# taxonomy_family_lineage <- tax_lineage(taxonomy_taxon_info(taxonomy_family$uid, include_lineage = TRUE))  
# saveRDS(taxonomy_family_lineage, file="docs/data/tmp/taxonomy_family_lineage.RData")

taxonomy_family_lineage <- readRDS(file="docs/data/tmp/taxonomy_family_lineage.RData")


taxonomy_family_lineage_matt <- ldply(taxonomy_family_lineage,rbind)


taxonomy_family_lineage_wide <- reshape(taxonomy_family_lineage_matt,
                                        idvar= ".id",
                                        timevar = "rank",
                                        direction="wide"
                                        )

taxonomy_family_lineage_wide$.id <- as.integer(taxonomy_family_lineage_wide$.id)


taxonomy_family_full <-
  left_join(taxonomy_family,taxonomy_family_lineage_wide,by=c("uid" = ".id"))


length(unique(taxonomy_family_full$name.family))
#my_tree1 <- tol_induced_subtree(ott_ids = taxonomy_family_full$uid)

#sp_name <-gsub("_.*","",my_tree$tip.label)

#my_tree$tip.label <- sp_name


#species <- taxonomy_family_full$name
#g <- split(species, taxonomy_family_full$unique_name.class)
#tree_plot <- ggtree(my_tree, layout='circular') +
#  geom_tiplab(size=2, offset=5) 

#groupOTU(tree_plot, g, 'Species') + aes(color=Species) +
#  theme(legend.position="right")


###################### gtoup PF 

PF_d2 <- PF_d[!is.na(PF_d$taxon.ott_id),]

# Ucomment the two following line if you need to launch tax_lineage() again

# taxonomy_family_lineage_pf <- tax_lineage(taxonomy_taxon_info(PF_d2$taxon.ott_id, include_lineage = TRUE))  
 #saveRDS(taxonomy_family_lineage_pf, file="docs/data/tmp/taxonomy_family_lineage_pf.RData")

 taxonomy_family_lineage_pf <- readRDS(file="docs/data/tmp/taxonomy_family_lineage_pf.RData")


taxonomy_family_lineage_pf_matt <- ldply(taxonomy_family_lineage_pf,rbind)




taxonomy_family_lineage_pf_wide <- reshape(taxonomy_family_lineage_pf_matt,
                                        idvar= ".id",
                                        timevar = "rank",
                                        direction="wide"
)

family_pf <- taxonomy_family_lineage_pf_wide$unique_name.family
taxonomy_family_lineage_pf_wide$unique_name.family
pres <- rep(1,length(family_pf))

fam_pf <- data.frame(family_pf,pres)

taxonomy_pf_merge <-
  left_join(taxonomy_family_full,fam_pf,by=c("name" = "family_pf"))

length(unique(taxonomy_family_full$name))

my_tree <- tol_induced_subtree(ott_ids = taxonomy_pf_merge$uid)

sp_name <-gsub("_.*","",my_tree$tip.label)

my_tree$tip.label <- sp_name


species <- taxonomy_pf_merge$name
g <- split(species, taxonomy_pf_merge$pres)

tree_plot <- ggtree(my_tree, layout='circular') +
  geom_tiplab(size=2.5, offset=0.5) 

g1<- groupOTU(tree_plot, g, 'species') + aes(color=species) +
  theme(legend.position="right") + scale_color_manual(values = c("darkgreen","orange"))


##### add order 

species_order_tab <- data.frame(taxonomy_pf_merge$name,taxonomy_pf_merge$unique_name.order)
colnames(species_order_tab) <- c("family","order")


####detecte and exclude family error in MRCA

family <- species_order_tab$family
fam_mrca <- tapply(family,family,function(x) return(tryCatch(MRCA(tree_plot, x),error=function(e) NULL)))
fam_mrca <- sapply(fam_mrca,as.numeric)
fam_mrca_filter <-names(fam_mrca[lapply(fam_mrca, sum) > 0])
species_order_tab <- species_order_tab[species_order_tab$family %in% fam_mrca_filter,]

g2 <- split(species_order_tab$family, species_order_tab$order)
clades <- sort(sapply(g2, function(n) MRCA(tree_plot, n))) ## sort for angle position circular


##### resolver


order_fam <- taxonomy_pf_merge$unique_name.order
fami <- taxonomy_pf_merge$name
vec_select_order <- rep(0,length(order_fam))
vec_select_order[order_fam == "Icacinales"] <- 1
gx <- split(fami, vec_select_order)

tree_plot <- ggtree(my_tree, layout='circular') +
geom_tiplab(size=2.5, offset=0.5) 

#groupOTU(tree_plot, gx, 'fami') + aes(color=fami) +
#theme(legend.position="right") + scale_color_manual(values = c("gray70","black"))+
#  geom_text(aes(label=node), size=2,vjust=1,hjust=1)


################################ fix polyphyletic clades

g2 <- split(species_order_tab$family, species_order_tab$order)
clades <- sort(sapply(g2, function(n) MRCA(tree_plot, n)))

clades[names(clades) == "Asparagales"] <- 685
clades[names(clades) == "Liliales"] <- 698
clades[names(clades) == "Arecales"] <- 263
clades[names(clades) == "Dilleniales"] <- 236
clades[names(clades) == "Icacinales"] <- 150



#clades <- clades[!(clades == "Arecales")]
    
################################

clades <- sort(clades)
family_angle <- order(clades)
x1<- ggtree(my_tree,layout='circular')
x1<-x1$data
x1_fam <- x1[x1$node %in% clades,]
x1_fam$clades_node <- clades

colfunc <- colorRampPalette(c("darkseagreen4","orangered4","slategray4","palegreen3","black",
                              "navy","firebrick2","darkslategrey"))
mycolx <-colfunc(length(clades)+1) 

col_bar <- lengths(g2) 
col_bar[col_bar>0] <- "black"
col_bar[c(1:25)] <- "white" ### color in white mono clade
#col_bar[c(1,2,8,9,11)] <- "blue"

x<- ggtree(my_tree,layout='circular')+
  geom_tiplab(size=1.3, offset=1) +
  geom_cladelab(clades, names(clades),offset=14,angle=(x1_fam$angle),
                barsize=5, offset.text=4, fontsize=3,barcolour=col_bar)


tree_plot <- groupClade(x, clades, group_name='subtree') + aes(color=subtree) +
  theme(legend.position = "none") +
  scale_colour_manual(values=mycolx) +
  theme(legend.position = "None")


species <- taxonomy_pf_merge$name
g <- split(species, taxonomy_pf_merge$pres)

g1<- groupOTU(tree_plot, g, 'species') + aes(color=species) +
   scale_color_manual(values = c("slategray2","darkred"))
g1

###############  bar plot
g1$data[1]
########################


taxonomy_genus <- taxonomy_final %>% 
  filter(rank == "genus" & is.na(flags)) 

length(unique(taxonomy_genus$))

# Ucomment the two following line if you need to launch tax_lineage() again

# taxonomy_genus_lineage <- tax_lineage(taxonomy_taxon_info(taxonomy_genus$uid, include_lineage = TRUE))  
# saveRDS(taxonomy_genus_lineage, file="docs/data/tmp/taxonomy_genus_lineage.RData")

taxonomy_genus_lineage <- readRDS(file="docs/data/tmp/taxonomy_genus_lineage.RData")


taxonomy_genus_lineage_matt <- ldply(taxonomy_genus_lineage,rbind)


taxonomy_genus_lineage_wide <- reshape(taxonomy_genus_lineage_matt,
                                        idvar= ".id",
                                        timevar = "rank",
                                        direction="wide"
)

taxonomy_genus_lineage_wide$.id <- as.integer(taxonomy_genus_lineage_wide$.id)




taxonomy_genus_full <-
  left_join(taxonomy_genus,taxonomy_genus_lineage_wide,by=c("uid" = ".id"))
taxonomy_genus_full$uid <- as.character(taxonomy_genus_full$uid)

length(unique(taxonomy_genus_full$))

sp_name_matt <- subset(PF_d2,select=c(taxon.ott_id, taxon.name))
sp_name_matt$taxon.ott_id = as.character(sp_name_matt$taxon.ott_id)

taxonomy_pf_merge_sp <-
  left_join(sp_name_matt,taxonomy_family_lineage_pf_wide,by=c("taxon.ott_id" = ".id"))
pres <- rep(1,nrow(taxonomy_pf_merge_sp))

taxonomy_pf_merge_sp <- data.frame(taxonomy_pf_merge_sp,pres)



Pf_full_merge_taxo <-
  left_join(taxonomy_genus_full,taxonomy_pf_merge_sp,by=c("name" = "name.genus"))

length(unique(Pf_full_merge_taxo$name.family.x))
length(unique(Pf_full_merge_taxo$name))


Pf_full_merge_taxo <- Pf_full_merge_taxo[Pf_full_merge_taxo$name.family.x %in% taxonomy_pf_merge$name,]

ratio_phy <- tapply(Pf_full_merge_taxo$pres,Pf_full_merge_taxo$name.phylum.x,mean,na.rm=T)
ratio_phy[is.na(ratio_phy)] <- 0
ratio_phy <- mean(ratio_phy)

ratio_order <- tapply(Pf_full_merge_taxo$pres,Pf_full_merge_taxo$name.order.x,mean,na.rm=T)
ratio_order[is.na(ratio_order)] <- 0
ratio_order <- mean(ratio_order)

ratio_fam<- tapply(Pf_full_merge_taxo$pres,Pf_full_merge_taxo$name.family.x,mean,na.rm=T)
ratio_fam[is.na(ratio_fam)] <- 0
ratio_fam<- mean(ratio_fam)

ratio_genus <- tapply(Pf_full_merge_taxo$pres,Pf_full_merge_taxo$name,mean,na.rm=T)
ratio_genus[is.na(ratio_genus)] <- 0
ratio_genus <-mean(ratio_genus)

ratio_species <- length(unique(PF_d2$taxon.name))/308312 ## number of tracheophyte

Group_coverage <-c(ratio_phy,ratio_order,ratio_fam,ratio_genus,ratio_species)
group <- c("Phylum","Order","Family","Genus","Species")
matt_ratio <-  data.frame(Group_coverage,group)

matt_ratio <- within(matt_ratio, 
                     group <- factor(group, 
                                      levels=group[order(Group_coverage,decreasing=FALSE)]))

cc <- scales::seq_gradient_pal("gray50","gray20")(seq(0,1,length.out=5))

p1 <- ggplot(data=matt_ratio, aes(x=group, y=Group_coverage*100,fill=group,label = round(Group_coverage*100,1))) +
  geom_bar(stat="identity",width = 0.5)+ coord_flip() + 
  scale_fill_manual(values=cc)+ 
  geom_text(size = 3, vjust = -0.5,angle = 270)+
  theme_classic() + theme(axis.title.y=element_blank(),  
                        plot.margin = margin(1, 1, 15, 15))+
  labs(y = "Coverage in %")
p1 <-p1 +guides(fill="none")


titlex <-  ggplot()+ xlim(0,1) +ylim(0,1) + theme_void() +
          annotate(geom="text", x=0.6, y=0.5, 
                   label="Taxonom",color="black",size =7,
                   family = "serif")

 grid.arrange(p1,g1, nrow = 2,widths = c(1,2),
             layout_matrix = rbind(c(1,2), c(1,2)))

pdf(file = "docs/data/outputs/taxo_plot.pdf",   # The directory you want to save the file in
    width = 13, # The width of the plot in inches
    height = 8)

grid.arrange(p1,g1, nrow = 2,widths = c(1,2),
             layout_matrix = rbind(c(1,2), c(1,2)))

dev.off()


################### plot otl coverage 
fctx <- function (X) {length(unique(X))}

ratio_order <- tapply(Pf_full_merge_taxo$pres,Pf_full_merge_taxo$name.order.x,mean,na.rm=T)
order_size <- tapply(Pf_full_merge_taxo$name.family.x,Pf_full_merge_taxo$name.order.x,fctx)
ratio_order[is.na(ratio_order)] <- 0
matt_order <- data.frame(ratio_order,order_size) 
matt_order$name.order.x <- row.names(matt_order)

ratio_fam<- tapply(Pf_full_merge_taxo$pres,Pf_full_merge_taxo$name.family.x,mean,na.rm=T)
fam_size <- tapply(Pf_full_merge_taxo$name,Pf_full_merge_taxo$name.family.x,fctx)
ratio_fam[is.na(ratio_fam)] <- 0
matt_fam <- data.frame(ratio_fam,fam_size) 
matt_fam$name.family.x <- row.names(matt_fam)

ratio_genus <- tapply(Pf_full_merge_taxo$pres,Pf_full_merge_taxo$name,mean,na.rm=T)
ratio_genus[is.na(ratio_genus)] <- 0
matt_genus <- data.frame(ratio_genus) 
matt_genus$name <- row.names(matt_genus)

taxo_merger <- Pf_full_merge_taxo %>%
              select(name.order.x, name.family.x)
taxo_merger <- taxo_merger[!duplicated(taxo_merger),]

taxo_merger_fam <- merge(taxo_merger,matt_fam,by="name.family.x")
taxo_merger_fam <- merge(taxo_merger_fam,matt_order,by="name.order.x")
taxo_merger_fam$value <- rep(1,nrow(taxo_merger_fam))


taxo_merger <- Pf_full_merge_taxo %>%
              select(name.family.x, name)
taxo_merger <- taxo_merger[!duplicated(taxo_merger),]

taxo_merger_genus <- merge(taxo_merger,matt_fam,by="name.family.x")
taxo_merger_genus <- merge(taxo_merger_genus,matt_genus,by="name")
taxo_merger_genus$value <- rep(1,nrow(taxo_merger_genus))




pdf(file = "docs/data/outputs/family_coverage_plot.pdf",   # The directory you want to save the file in
    width = 13, # The width of the plot in inches
    height = 8)

p3 <- ggplot(taxo_merger_genus, aes(y = value, fill = as.factor(ratio_genus), x = reorder(name.family.x, -fam_size),text = name.family.x)) +
  geom_bar(stat = "identity", width = 0.5, position = "stack") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
  xlab("Family") +
  ylab("Number of genus") +
  scale_fill_manual(values = c("slategray2", "darkred"), labels = c("Absent from dataset", "Present in dataset")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(0.9, 0.9)) +
  guides(fill = guide_legend(title = ""))
mean(taxo_merger_genus$ratio_genus)
p3 

dev.off()

pdf(file = "docs/data/outputs/order_coverage_plot.pdf",   # The directory you want to save the file in
    width = 13, # The width of the plot in inches
    height = 8) 

p4 <- ggplot(taxo_merger_fam, aes(y = value, fill = as.factor(ratio_fam), x = reorder(name.order.x, -order_size),text = name.order.x)) +
  geom_bar(stat = "identity", width = 0.5, position = "stack") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  xlab("Order") +
  ylab("Number of families") +
  scale_fill_manual(values = c("slategray2", "darkred"), labels = c("Absent from dataset", "Present in dataset")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(0.9, 0.9)) +
  guides(fill = guide_legend(title = ""))
mean(taxo_merger_fam$ratio_fam)
p4

dev.off()



p3_ly <- ggplotly(p3,dynamicTicks = TRUE,tooltip = c("text"))
p3_ly %>%
        layout(legend = list(title=list(text='Presence in dataset')), hoverinfo = 'Family')

setwd("docs/data/outputs")
p3_ly %>% 
  htmlwidgets::saveWidget(file="family_coverage_plot.html", selfcontained = TRUE)
system('rm -r family_coverage_plot_files')



p4_ly <- ggplotly(p4,dynamicTicks = TRUE,tooltip = c("text"))
p4_ly %>%
        layout(legend = list(title=list(text='Presence in dataset')), hoverinfo = 'Family')

setwd("docs/data/outputs")
p4_ly %>% 
  htmlwidgets::saveWidget(file="order_coverage_plot.html", selfcontained = TRUE)
system('rm -r order_coverage_plot_files')