#http://itol.embl.de

library(rotl)
library(readr)
library(dplyr)
library(igraph)
library(plyr)
library(gridExtra)
library(ggfortify)
library(ggtree)

  ###############  

# The full taxo is downloaded from OTL 
# https://tree.opentreeoflife.org/about/taxonomy-version/ott3.3


#taxonomy <- read_delim("~/Downloads/ott3.3/taxonomy.tsv", delim = "|", escape_double = FALSE, trim_ws = TRUE)
taxonomy <- read_delim("C:/users/defossee/Downloads/ott3.3/taxonomy.tsv", delim = "|", escape_double = FALSE, trim_ws = TRUE)
PF_d <- read_delim("G:/My Drive/taf/git_repository/phylogenetic_pfdatanote_sandbox/data/inputs/pf_metadata_otoled.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)



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
                       

taxonomy_family_lineage <- tax_lineage(taxonomy_taxon_info(taxonomy_family$uid, include_lineage = TRUE))  

taxonomy_family_lineage_matt <- ldply(taxonomy_family_lineage,rbind)


taxonomy_family_lineage_wide <- reshape(taxonomy_family_lineage_matt,
                                        idvar= ".id",
                                        timevar = "rank",
                                        direction="wide"
                                        )

taxonomy_family_lineage_wide$.id <- as.integer(taxonomy_family_lineage_wide$.id)


taxonomy_family_full <-
  left_join(taxonomy_family,taxonomy_family_lineage_wide,by=c("uid" = ".id"))


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

taxonomy_family_lineage_pf <- tax_lineage(taxonomy_taxon_info(PF_d2$taxon.ott_id, include_lineage = TRUE))  

taxonomy_family_lineage_pf_matt <- ldply(taxonomy_family_lineage_pf,rbind)


taxonomy_family_lineage_pf_wide <- reshape(taxonomy_family_lineage_pf_matt,
                                        idvar= ".id",
                                        timevar = "rank",
                                        direction="wide"
)

family_pf <- taxonomy_family_lineage_pf_wide$unique_name.family
pres <- rep(1,length(family_pf))

fam_pf <- data.frame(family_pf,pres)

taxonomy_pf_merge <-
  left_join(taxonomy_family_full,fam_pf,by=c("name" = "family_pf"))



my_tree <- tol_induced_subtree(ott_ids = taxonomy_pf_merge$uid)

sp_name <-gsub("_.*","",my_tree$tip.label)

my_tree$tip.label <- sp_name


#species <- taxonomy_pf_merge$name
#g <- split(species, taxonomy_pf_merge$pres)

#tree_plot <- ggtree(my_tree, layout='circular') +
#  geom_tiplab(size=2.5, offset=0.5) 

#g1<- groupOTU(tree_plot, g, 'species') + aes(color=species) +
#  theme(legend.position="right") + scale_color_manual(values = c("darkgreen","orange"))


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

groupOTU(tree_plot, gx, 'fami') + aes(color=fami) +
theme(legend.position="right") + scale_color_manual(values = c("gray70","black"))+
  geom_text(aes(label=node), size=2,vjust=1,hjust=1)


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

########################


taxonomy_genus <- taxonomy_final %>% 
  filter(rank == "genus" & is.na(flags)) 


taxonomy_genus_lineage <- tax_lineage(taxonomy_taxon_info(taxonomy_genus$uid, include_lineage = TRUE))  

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




sp_name_matt <- subset(PF_d2,select=c(taxon.ott_id, taxon.name))
sp_name_matt$taxon.ott_id = as.character(sp_name_matt$taxon.ott_id)

taxonomy_pf_merge_sp <-
  left_join(sp_name_matt,taxonomy_family_lineage_pf_wide,by=c("taxon.ott_id" = ".id"))
pres <- rep(1,nrow(taxonomy_pf_merge_sp))

taxonomy_pf_merge_sp <- data.frame(taxonomy_pf_merge_sp,pres)



Pf_full_merge_taxo <-
  left_join(taxonomy_genus_full,taxonomy_pf_merge_sp,by=c("name" = "name.genus"))




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

pdf(file = "G:/My Drive/taf/git_repository/phylogenetic_pfdatanote_sandbox/data/outputs/taxo_plot.pdf",   # The directory you want to save the file in
    width = 13, # The width of the plot in inches
    height = 8)

grid.arrange(p1,g1, nrow = 2,widths = c(1,2),
             layout_matrix = rbind(c(1,2), c(1,2)))

dev.off()