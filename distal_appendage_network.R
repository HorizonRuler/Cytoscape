#packages
library(data.table)
library(igraph)
library(RCy3)
library(dplyr)
#protein information
keyGenes <- c('CEP164', 'CEP83', 'SCLT1', 'FBF1', 'CEP89', 'TTBK2', 'NCS1', 'C3orf14', 'KIZ', 'CCDC92', 'LRRC45', 'INPP5E', 'PDE6D', 'ANKRD26', 'MAPRE1', 'MAPRE2', 'CLASP1', 'CLASP2', 'RABL2A', 'RABL2B', 'CLUAP1', 'HSPB11', 'IFT122', 'IFT140', 'IFT172', 'IFT20', 'IFT22', 'IFT27', 'IFT43', 'IFT46', 'IFT52', 'IFT57', 'IFT74', 'IFT80', 'IFT81', 'IFT88', 'TRAF3IP1', 'TTC21B', 'TTC26', 'TTC30B', 'WDR19', 'WDR35')
humanGenes <- fread("~/Documents/2022_summer_project/Cytoscape/human_genes.txt")
#read in and reduce public data to humans and good experiments
biogridOriginal <- fread("~/Documents/2022_summer_project/Cytoscape/BIOGRID-ALL-4.4.210.tab3.txt")
biogrid <- biogridOriginal[`Organism ID Interactor A` == 9606 &
                     `Organism ID Interactor B` == 9606]
biogrid <- biogrid[`Experimental System` %in% c('Affinity Capture-MS', 
                                                'Affinity Capture-Western', 
                                                'Reconstituted Complex', 
                                                'Co-crystal Structure', 
                                                'Far Western')]
#get rid of unneeded column information
biogrid <- biogrid[, c(1, 2, 3, 8, 9)]
biogrid <- biogrid [, .(interactionID = `#BioGRID Interaction ID`, 
       geneIDA = `Entrez Gene Interactor A`, 
       geneIDB = `Entrez Gene Interactor B`, 
       geneSymbolA = `Official Symbol Interactor A`, 
       geneSymbolB = `Official Symbol Interactor B`)]
#combine information about same genes
biogrid[, geneIDA := as.integer(geneIDA)]
biogrid[, geneIDB := as.integer(geneIDB)]
biogrid[, a := ifelse(geneIDA < geneIDB, geneIDA, geneIDB)]
biogrid[, b := ifelse(geneIDA < geneIDB, geneIDB, geneIDA)]
biogrid[, x := ifelse(geneIDA < geneIDB, geneSymbolA, geneSymbolB)]
biogrid[, y := ifelse(geneIDA < geneIDB, geneSymbolB, geneSymbolA)]
bg <- biogrid[, .(pubs = .N), .(a, b, x, y)]
colnames(bg)[1:4] <- colnames(biogrid)[2:5]
#create igraph network
edges <- bg[, .(source = geneIDA, target = geneIDB, pubs)]
nodes <- unique(rbind(bg[, .(eid = geneIDA, symbol = geneSymbolA)], bg[, .(eid = geneIDB, symbol = geneSymbolB)]))
biogridGraph <- graph_from_data_frame(edges, FALSE, nodes)
V(biogridGraph)$baits <- ifelse(V(biogridGraph)$symbol %in% keyGenes, 1, 0)
E(biogridGraph)$interaction <- 'biogrid'
#make and combine close neighbor graphs based on distal appendage proteins
egoGraphs <- make_ego_graph(biogridGraph, 1, which(V(biogridGraph)$symbol == 'CEP164'))[[1]]
egoGraphs <- subgraph.edges(egoGraphs, which(E(egoGraphs)$pubs > 1), delete.vertices = FALSE)
egoGraphs <- make_ego_graph(egoGraphs, 1, which(V(egoGraphs)$symbol == 'CEP164'))[[1]]
for (i in keyGenes[2:length(keyGenes)]) {
  t <- make_ego_graph(biogridGraph, 1, which(V(biogridGraph)$symbol == i))[[1]]
  t <- subgraph.edges(t, which(E(t)$pubs > 1), delete.vertices = FALSE)
  t <- make_ego_graph(t, 1, which(V(t)$symbol == i))[[1]]
  egoGraphs <- igraph::union(egoGraphs, t)
  V(egoGraphs)$symbol <- ifelse(is.na(V(egoGraphs)$symbol_1), V(egoGraphs)$symbol_2, V(egoGraphs)$symbol_1)
  egoGraphs <- delete_vertex_attr(egoGraphs, "symbol_1")
  egoGraphs <- delete_vertex_attr(egoGraphs, "symbol_2")
}
egoGraphs <- induced_subgraph(biogridGraph, which(V(biogridGraph)$symbol %in% V(egoGraphs)$symbol))
egoGraphs <- simplify(subgraph.edges(egoGraphs, which(E(egoGraphs)$pubs > 1), delete.vertices = FALSE), edge.attr.comb = "random")
#get lab data and filter for top edges while removing duplicate experiments, adding bait ids, and ensuring each distal appendage protein has edges
labData <- readRDS("~/Documents/2022_summer_project/Cytoscape/distal_appendage_network.rds")
labData <- labData[, .SD[order(pval)][1], .(bait, eid, interaction = cline)]
labData <- labData[nsaf > 5*10^-4 & umodpept > 1 & pval < 0.01]
labData <- labData %>% left_join(labData[symbol %in% (labData[,.N, bait][, bait])][,.N, .(eid, symbol)][,.(beid = eid, symbol)], by = c('bait' = 'symbol'))
labData <- labData[!endsWith(symbol, '-bait')]
#create network of lab data edges
labGraph <- simplify(graph_from_data_frame(labData[, .(source = beid, target = eid, interaction, pval, nsaf, cov, psm, upept, umodpept, pval)], FALSE, unique(labData[, .(eid, symbol, descr)])), remove.multiple = F)
V(labGraph)$baits <- ifelse(V(labGraph)$symbol %in% labData[, bait], 1, 0)
#combine lab data with public nodes to create network of public edges
combinedNodes <- unique(c(V(egoGraphs)$name, labData[,eid]))
combinedGraph <- induced_subgraph(biogridGraph, which(V(biogridGraph)$name %in% combinedNodes))
combinedGraph <- simplify(subgraph.edges(combinedGraph, which(E(combinedGraph)$pubs > 1), delete.vertices = FALSE), edge.attr.comb = "random")
#combine public and lab data edges
finalGraph <- graph_from_data_frame(rbindlist(list(igraph::as_data_frame(combinedGraph), igraph::as_data_frame(labGraph)), fill = T), F)
V(finalGraph)$symbol <- humanGenes[match(as.integer(V(finalGraph)$name), eid), symbol]
V(finalGraph)$descr <- humanGenes[match(as.integer(V(finalGraph)$name), eid), descr]
V(finalGraph)$baits <- ifelse(V(finalGraph)$symbol %in% keyGenes, 1, 0)
igraph::as_data_frame(finalGraph, what = 'v') %>% as.data.table() %>% .[humanGenes[,.(name = as.character(eid), descr)], on = .(name), 'description' := i.descr] %>% fwrite("~/Documents/2022_summer_project/Cytoscape/finalVertices")
fwrite(as.data.table(igraph::as_data_frame(finalGraph)), "~/Documents/2022_summer_project/Cytoscape/finalEdges")
#style network
setNodeColorMapping('baits', c(0, 1), c('#CCCCFF', '#FFCCCC'))
setNodeLabelMapping('symbol')
setEdgeColorMapping('interaction', c('hek', 'rpe'), c('#00FF00', '#00FF00'), mapping.type = 'd')
setEdgeLineWidthMapping('nsaf', c(min(E(finalGraph)$nsaf, na.rm = T), mean(E(finalGraph)$nsaf, na.rm = T), max(E(finalGraph)$nsaf, na.rm = T)), c(1, 2, 30))

