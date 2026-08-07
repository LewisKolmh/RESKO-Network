# Approach A: RESKO Reproduction

## Objective

Apply the RESKO drug repositioning workflow to identify candidate
eEF1A inhibitors using side-effect similarity and biological
network enrichment.

## A1 Seed Drug Selection

Initial seed compounds:

- Plitidepsin
- Didemnin B
- Ternatin-4
- Nannocystin A
- Cytotrienin A
- Metarrestin

Status: Complete

## A2 Mapping Seed Drugs to SIDER

Status: In Progress

## A2 Results

The initial eEF1A seed compounds were mapped against the
1430 drug records contained within SIDER 4.1.

All six seed compounds were absent from the database.

Result:

0 / 6 compounds matched SIDER.

Implication:

The original RESKO side-effect similarity workflow could
not be applied directly to the initial eEF1A inhibitor set.

A seed expansion stage was therefore required.


## A3 Results

Fifteen literature-derived eEF1A-targeting compounds were
screened against SIDER 4.1.

No compounds were identified within the database.

Result:

0 / 15 compounds were represented.

This suggests that direct eEF1A inhibitors are poorly
represented within adverse-event resources used by RESKO.

# A4 Construction of an eEF1A-Centred Network

## Objective

Construct a protein interaction network centred on EEF1A1 and
EEF1A2.

## Rationale

Known eEF1A inhibitors were absent from SIDER and therefore
could not be used directly in a side-effect driven RESKO
workflow.

A network-centric approach was therefore adopted to identify
biological partners, pathways and candidate druggable proteins
associated with eEF1A.

## Data Source

STRING API

Species:
Homo sapiens (9606)

## Seed Proteins

EEF1A1
EEF1A2

## Outputs

A4_eef1a_string_interactions.csv
> summary(interactions$score)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.4020  0.9317  0.9945  0.9366  0.9990  0.9990 

 ## A5 Neo4j Network File Creation

### Objective

Convert STRING interaction data into graph nodes and edges
suitable for Neo4j import.

### Input

A4_eef1a_string_interactions.csv

### Outputs

nodes_proteins.csv

edges_interacts_with.csv

### Relationship Type

INTERACTS_WITH

### Purpose

Generate the initial eEF1A-centred protein interaction
network for downstream Neo4j Bloom visualisation.

## A6 Reactome Pathway Enrichment

### Result

46 significantly enriched Reactome pathways were identified
from the EEF1A-centred interaction network including 18 unique proteins from 84 STRING interactions.

### Interpretation

The resulting pathways provide biological context for the
EEF1A interaction neighbourhood and form the pathway layer
of the network graph.

### Output

results/A6_reactome_pathways.csv

A6 Status: Complete

STRING interaction network:
84 edges

Unique Proteins Identified:
18 Proteins

Reactome enrichment:
46 pathways

Notable pathways:

- Eukaryotic Translation Elongation
- Translation
- Peptide Chain Elongation
- Viral mRNA Translation

Significance:

Minimum adjusted p-value:
3.40e-20

 Description                                                    Count p.adjust
   <chr>                                                          <dbl>    <dbl>
 1 Eukaryotic Translation Elongation                                 12 3.40e-20
 2 Translation                                                       13 1.01e-14
 3 Selenoamino acid metabolism                                        9 9.62e-13
 4 Peptide chain elongation                                           8 5.93e-12
 5 Eukaryotic Translation Termination                                 8 7.00e-12
 6 Nonsense Mediated Decay (NMD) independent of the Exon Junctio…     8 7.04e-12
 7 Nonsense-Mediated Decay (NMD)                                      8 2.80e-11
 8 Nonsense Mediated Decay (NMD) enhanced by the Exon Junction C…     8 2.80e-11
 9 Viral mRNA Translation                                             7 2.54e-10
10 PELO:HBS1L and ABCE1 dissociate a ribosome on a non-stop mRNA      7 2.54e-10

## A6 Results

### STRING Network

84 protein interaction edges were retrieved for EEF1A1
and EEF1A2 from STRING.

### Reactome Enrichment

46 significantly enriched Reactome pathways were identified.

Top pathways:

- Eukaryotic Translation Elongation
- Translation
- Selenoamino Acid Metabolism
- Peptide Chain Elongation
- Eukaryotic Translation Termination
- Viral mRNA Translation

Minimum adjusted p-value:

3.40 × 10^-20

### Graph Construction

Protein nodes retained in pathway enrichment:
17

Protein-pathway relationships:
276

> nrow(pathway_nodes)
[1] 46
> 
> nrow(protein_pathway_edges)
[1] 276
> 
> head(pathway_nodes)
# A tibble: 6 × 3
  pathway                                                         type  p_adjust
  <chr>                                                           <chr>    <dbl>
1 Eukaryotic Translation Elongation                               Path… 3.40e-20
2 Translation                                                     Path… 1.01e-14
3 Selenoamino acid metabolism                                     Path… 9.62e-13
4 Peptide chain elongation                                        Path… 5.93e-12
5 Eukaryotic Translation Termination                              Path… 7.00e-12
6 Nonsense Mediated Decay (NMD) independent of the Exon Junction… Path… 7.04e-12
> 
> head(protein_pathway_edges)

## A7 GO Biological Process Enrichment

### Objective

Identify biological processes significantly enriched
among proteins present in the eEF1A-centred interaction network.

### Input

A6_protein_master_list.csv

### Method

Gene symbols were converted to Entrez identifiers
using org.Hs.eg.db.

GO Biological Process enrichment was performed using
clusterProfiler.

### Results

Mapped genes: 18

Enriched GO Biological Processes: 53

### Output

results/A7_GO_Biological_Process.csv

   Description                                                    Count p.adjust
   <chr>                                                          <dbl>    <dbl>
 1 cytoplasmic translation                                            9 3.37e-11
 2 translational elongation                                           5 1.12e- 6
 3 negative regulation of myoblast fusion                             3 8.96e- 4
 4 negative regulation of syncytium formation by plasma membrane…     3 8.96e- 4
 5 regulation of myoblast fusion                                      3 1.19e- 3
 6 regulation of syncytium formation by plasma membrane fusion        3 1.40e- 3
 7 myoblast fusion                                                    3 1.98e- 3
 8 syncytium formation by plasma membrane fusion                      3 1.98e- 3
 9 cell-cell fusion                                                   3 1.98e- 3
10 syncytium formation                                                3 1.98e- 3
11 chaperone-mediated autophagy                                       2 1.98e- 3
12 nucleoside bisphosphate metabolic process                          2 1.98e- 3
13 ribonucleoside bisphosphate metabolic process                      2 1.98e- 3
14 purine nucleoside bisphosphate metabolic process                   2 1.98e- 3
15 purine ribonucleoside bisphosphate metabolic process               2 1.98e- 3

## A7 Results

### GO Biological Process Enrichment

Mapped genes: 18

Enriched GO Biological Processes: 53

Unique GO Nodes: 33

Protein-GO Relationships: 93

Output Files:

- A7_GO_Biological_Process.csv
- nodes_biological_process.csv
- edges_protein_go.csv

# A8 Drug Target Integration

## Objective

Identify druggable targets within the eEF1A-centred biological network and establish relationships between network proteins and experimentally annotated ChEMBL targets.

## Rationale

Direct eEF1A inhibitors were not represented within the SIDER database and therefore could not be used directly in a side-effect driven repositioning workflow.

Following construction of the eEF1A-centred interaction network and functional enrichment analyses, proteins within the network were mapped to ChEMBL target identifiers to facilitate downstream compound discovery.

## Input

results/nodes_proteins.csv

## Protein Export

### Objective

Extract a unique list of proteins present within the eEF1A-centred interaction network.

### Output

results/A8_target_proteins.csv

### Results

Unique proteins exported: 18

Proteins identified:

- RPL18A
- RPS3
- RPL4
- EEF1G
- EEF1A1
- RPL7
- RPS2
- RPL3
- RPS3A
- EEF1B2
- ST6GALNAC1
- EEF1A2
- PAPSS1
- ETF1
- EEF1D
- HSF1
- PAPSS2
- RAN

## ChEMBL Target Mapping

### 