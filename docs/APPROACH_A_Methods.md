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