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