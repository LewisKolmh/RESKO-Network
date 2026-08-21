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

# A9 Complete ChEMBL Activity Retrieval

## Objective

Retrieve all available ChEMBL activity records associated with proteins from the eEF1A-centred interaction network that had successfully mapped to ChEMBL targets.

## Rationale

The initial ChEMBL activity-retrieval workflow was limited to the first page of results returned by the ChEMBL API. ChEMBL responses are paginated, meaning that a complete analysis required retrieval of every available results page for each target.

A paginated activity-retrieval workflow was therefore developed to ensure that all available activity records were captured.

## Input

results/A8_chembl_target_map.csv

## Data Source

ChEMBL REST API

## Method

Each mapped ChEMBL target was queried against the ChEMBL activity endpoint.

The workflow:

1. Queried each ChEMBL target identifier.
2. Retrieved results in pages of up to 1,000 activity records.
3. Continued requesting pages until the total number of records reported by the API had been reached.
4. Introduced pauses between requests to reduce API load.
5. Repeated failed requests up to three times.
6. Preserved the complete activity record returned by ChEMBL.
7. Added the original queried protein and ChEMBL target identifier to every activity record.

## Outputs

results/A9_chembl_activities_complete.csv

results/A9_chembl_retrieval_summary.csv

## Results

Mapped ChEMBL targets queried: 10

Complete activity records retrieved: 81

Unique molecules represented: 22

Protein activity-record counts:

- EEF1A1: 28
- EEF1G: 12
- RPS3: 10
- RPS3A: 8
- RPL18A: 4
- RPL3: 4
- RPL4: 4
- RPL7: 4
- RPS2: 4
- EEF1B2: 3

Activity types identified included:

- IC50
- Kd
- ED50
- Inhibition
- Ka
- Activity
- Fold change
- Ratio

The retrieved records contained both quantitative concentration measurements and contextual activity measurements. This required further evidence classification before compounds could be added to the graph.

---

# A10 ChEMBL Activity Quality Control and Evidence Classification

## Objective

Quality-control the retrieved ChEMBL records and distinguish direct inhibitory evidence from binding, functional-response, and contextual activity evidence.

## Rationale

ChEMBL activity types are biologically distinct and should not be interpreted interchangeably.

For example:

- IC50 records describe inhibitory potency.
- Kd records describe binding affinity.
- ED50 records describe functional dose-response measurements.
- Percentage inhibition records require the tested concentration for interpretation.
- Fold-change, ratio, and generic activity records do not directly describe molar potency.

A quality-control workflow was therefore required before compounds could be labelled or ranked.

## Input

results/A9_chembl_activities_complete.csv

## Method

Activity records were standardised and evaluated for:

- Presence of a ChEMBL molecule identifier.
- Presence of a numeric standard value.
- Standard activity type.
- Standard activity units.
- Activity-value relationship.
- ChEMBL data-validity warnings.
- Potential duplicate status.
- pChEMBL value.
- Assay identifier.
- Target identifier.
- Measurement type.

Quantitative concentration records were defined as IC50, Kd, or ED50 measurements reported in nanomolar units.

## Evidence Classes

Activity records were classified into the following evidence classes:

### Inhibitory Potency

IC50 records describing inhibition of the target-associated measurement.

### Binding Affinity

Kd records describing physical binding between the molecule and the target.

### Functional Dose Response

ED50 records describing a quantitative functional response.

### Percentage Inhibition

Records reporting percentage inhibition without directly comparable molar potency.

### Association Constant

Ka records describing association strength in the units reported by the original assay.

### Contextual Activity

Fold change, ratio, and generic activity measurements retained for interpretation but excluded from the initial quantitative ranking.

## Evidence Tiers

### Tier 1: Quantitative Inhibition

IC50 records in nanomolar units with an interpretable activity relation and no explicit ChEMBL validity warning.

### Tier 2: Quantitative Binding

Kd records in nanomolar units with an interpretable activity relation and no explicit ChEMBL validity warning.

### Tier 3: Quantitative Functional Response

ED50 records in nanomolar units with an interpretable activity relation and no explicit ChEMBL validity warning.

### Quantitative Record Requiring Review

Quantitative measurements with relationships or metadata requiring additional assay-level interpretation.

### Contextual Record

Measurements that could not be directly ranked using concentration-based potency.

## Potency Bands

Quantitative records were assigned broad apparent-potency bands:

- Very high apparent potency: 10 nM or lower
- High apparent potency: greater than 10 nM to 100 nM
- Moderate apparent potency: greater than 100 nM to 1,000 nM
- Lower apparent potency: greater than 1 µM to 10 µM
- Weak apparent potency: greater than 10 µM

These categories were used for visualisation and initial prioritisation only and were not treated as definitive biological classifications.

## Aggregation

Repeated measurements were aggregated separately for each:

- ChEMBL molecule
- Protein
- ChEMBL target
- Standard activity type
- Standard unit

IC50, Kd, and ED50 values were not combined.

For each compound-protein-activity combination, the following values were calculated:

- Measurement count
- Minimum activity value
- Median activity value
- Maximum activity value
- Maximum pChEMBL value
- Number of unique assays
- Number of validity warnings
- Number of potential duplicate records

## Outputs

results/A10_activity_quality_controlled.csv

results/A10_quantitative_activities.csv

results/A10_priority_inhibitory_activities.csv

results/A10_compound_protein_summary.csv

results/A10_activity_filter_summary.csv

results/A10_candidate_molecule_ids.csv

## Results

Quantitative activity records: 39

Candidate molecules with interpretable quantitative evidence: 6

The candidate molecules were:

- CHEMBL1232461
- CHEMBL1802814
- CHEMBL1802815
- CHEMBL1802973
- CHEMBL3752910
- CHEMBL5653589

The quantitative evidence included:

- IC50 inhibitory-potency records
- Kd binding-affinity records
- ED50 functional-response records

The resulting evidence classifications were retained as relationship properties rather than treating every compound as a confirmed inhibitor.

---

# A11 Candidate Molecule Annotation

## Objective

Retrieve chemical, structural, and development metadata for candidate ChEMBL molecules identified during A10.

## Input

results/A10_candidate_molecule_ids.csv

## Data Source

ChEMBL molecule API

## Method

Each candidate ChEMBL molecule identifier was queried against the ChEMBL molecule endpoint.

The following information was retrieved where available:

- Preferred compound name
- ChEMBL molecule identifier
- Molecule type
- Maximum clinical-development phase
- First approval year
- Therapeutic flag
- Natural-product flag
- Prodrug status
- Administration-route flags
- Parent ChEMBL identifier
- Canonical SMILES
- Standard InChI
- Standard InChIKey
- Molecular weight
- ALogP
- Polar surface area
- Hydrogen-bond acceptors
- Hydrogen-bond donors
- Rule-of-five violations

Where ChEMBL did not provide a preferred name, the ChEMBL molecule identifier was retained as the display name.

## Development Status Classification

Candidate molecules were classified as:

- Approved or launched
- Phase 3
- Phase 2
- Phase 1
- Preclinical or discovery
- Unknown

## Outputs

results/A11_candidate_molecule_annotations.csv

results/A11_molecule_retrieval_summary.csv

results/nodes_drugs.csv

## Results

Six candidate compound nodes were created.

One compound had a recognised preferred name:

- CHEMBL1232461: Molibresib

Other compounds were retained using their ChEMBL identifiers because no preferred compound name was available from the retrieved molecule records.

Molibresib was classified as a Phase 2 small molecule.

The remaining compounds were primarily experimental or preclinical molecules with unknown clinical-development status.

---

# A12 Compound-Protein Relationship Construction

## Objective

Convert quality-controlled ChEMBL activity evidence into graph-ready compound-protein relationships.

## Inputs

results/A10_compound_protein_summary.csv

results/A11_candidate_molecule_annotations.csv

results/nodes_proteins.csv

## Method

Activity records were converted into distinct graph relationships according to the type of experimental evidence.

## Relationship Types

### INHIBITS

Created for IC50 evidence.

This relationship indicates quantitative inhibitory-potency evidence but does not by itself establish the precise molecular mechanism of inhibition.

### BINDS_TO

Created for Kd evidence.

This relationship indicates measured physical binding affinity but does not automatically establish functional inhibition.

### HAS_FUNCTIONAL_ACTIVITY_AGAINST

Created for ED50 evidence.

This relationship indicates a measured functional dose response.

### HAS_ACTIVITY_AGAINST

Reserved for other interpretable activity records that could not be assigned to the preceding categories.

## Edge Properties

Each evidence-level compound-protein relationship retained:

- Compound ChEMBL identifier
- Compound name
- Protein symbol
- ChEMBL target identifier
- Relationship type
- Activity type
- Activity units
- Evidence class
- Evidence tier
- Evidence rank
- Potency band
- Measurement count
- Minimum activity value
- Median activity value
- Maximum activity value
- Maximum pChEMBL value
- Unique assay count
- Validity-warning count
- Potential-duplicate count
- Molecule type
- Maximum clinical phase
- Development status

## Simplified Graph Relationships

A simplified relationship table was also created containing one edge per compound-protein pair.

Where multiple activity types were associated with the same compound-protein pair, the strongest available evidence type was used as the primary display relationship.

All other activity types were retained as edge properties.

The evidence preference order was:

1. IC50 inhibitory evidence
2. Kd binding evidence
3. ED50 functional-response evidence
4. Other activity evidence

## Outputs

results/A12_edges_compound_protein_evidence.csv

results/edges_compound_protein.csv

results/A12_compound_protein_edge_summary.csv

---

# A13 Addition of the Compound Layer to the Interactive Network

## Objective

Integrate candidate compounds and quality-controlled compound-protein evidence into the interactive eEF1A-centred network.

## Inputs

### Node Files

- results/nodes_drugs.csv
- results/nodes_proteins.csv
- results/nodes_pathways.csv
- results/nodes_biological_process.csv

### Relationship Files

- results/edges_compound_protein.csv
- results/edges_interacts_with.csv
- results/edges_protein_pathway.csv
- results/edges_protein_go.csv

## Network Node Types

### Compound

Displayed as orange diamond-shaped nodes.

### Protein

Displayed as red circular nodes.

### Reactome Pathway

Displayed as green circular nodes.

### GO Biological Process

Displayed as blue circular nodes.

## Network Relationship Types

### Compound-Protein Relationships

- Red line: inhibitory-potency evidence
- Purple line: binding-affinity evidence
- Gold line: functional-response evidence
- Orange line: other compound activity

### Biological Relationships

- Grey line: protein-protein interaction
- Green line: protein-pathway association
- Blue line: protein-GO association

## Visualisation Features

The interactive graph includes:

- Connectivity-based node sizing
- STRING-score-based protein-interaction edge widths
- Search by node name
- Filtering by node type
- Hover tooltips
- Pan and zoom controls
- Reproducible layout
- Automatic layout stabilisation
- Neighbourhood highlighting
- Reduced opacity for unrelated elements following selection
- A complete visual key describing every node and relationship type

## Compound Tooltips

Compound tooltips display:

- Compound name
- ChEMBL identifier
- Molecule type
- Development status
- Associated proteins
- Activity types
- Strongest quantitative record
- Number of measurements

## Compound-Protein Edge Tooltips

Compound-protein relationship tooltips display:

- Relationship classification
- Compound name
- Protein name
- Activity types
- Minimum quantitative activity
- Number of measurements
- Development status

## Output

results/eef1a_network_A13_compounds.html

---

# A13 Candidate Compounds Identified

## Molibresib

Molibresib was the only named and clinically investigated compound identified.

The network contained IC50 evidence linking molibresib to:

- EEF1A1
- EEF1G
- RPS3
- RPS3A

The strongest apparent activity was against RPS3 at 240 nM.

The reported EEF1A1 IC50 was 10,000 nM and was classified as requiring further review.

Molibresib is therefore treated as a clinically investigated repositioning candidate with possible off-target activity in the eEF1A translation network rather than as a validated direct eEF1A inhibitor.

## CHEMBL1802814

CHEMBL1802814 had direct EEF1A1 binding evidence with a Kd of 3.74 nM.

This compound was classified as a high-priority direct eEF1A1-binding candidate.

Binding evidence does not yet establish functional inhibition.

## CHEMBL1802815

CHEMBL1802815 had direct EEF1A1 binding evidence with a Kd of 8.15 nM.

This compound was classified as a high-priority direct EEF1A1-binding candidate.

## CHEMBL1802973

CHEMBL1802973 had direct EEF1A1 binding evidence with a Kd of 57.2 nM.

This compound was classified as a moderate-to-high-affinity direct EEF1A1-binding candidate and as a useful comparator for structure-activity analysis.

## CHEMBL3752910

CHEMBL3752910 had binding and functional-response evidence across several translation-related proteins, including:

- EEF1G
- RPS3
- RPS2
- RPS3A
- RPL3
- RPL4
- RPL7
- RPL18A

The strongest records were associated with RPS3 and EEF1G.

This compound was classified as a broad translation-neighbourhood or chemogenomic probe rather than a selective eEF1A inhibitor.

## CHEMBL5653589

CHEMBL5653589 had very high-affinity binding evidence against EEF1G, with a Kd of 0.491 nM.

It also had weaker evidence against:

- RPS3
- EEF1B2
- RPS2
- RPL18A
- RPL3
- RPL4
- RPL7

This compound was classified as a high-affinity EEF1G-centred candidate with broader translation-network activity.

## Current Candidate Categories

### Direct EEF1A1-Binding Candidates

- CHEMBL1802814
- CHEMBL1802815
- CHEMBL1802973

### Clinically Investigated Repositioning Candidate

- Molibresib

### Translation-Complex or Ribosomal-Neighbourhood Candidates

- CHEMBL3752910
- CHEMBL5653589

## Interpretation

The network currently distinguishes direct target engagement from broader translation-network perturbation.

The compounds should not all be described as confirmed eEF1A inhibitors.

The current evidence supports classification as:

- Direct EEF1A1 binders
- Putative inhibitors
- Translation-complex binders
- Ribosomal-network compounds
- Clinically investigated repositioning candidates

Functional validation is required to determine whether any candidate:

- Inhibits the canonical activity of eEF1A
- Alters eEF1A conformation
- Disrupts an eEF1A protein-protein interaction
- Reduces translation elongation
- Produces indirect effects through broader translation-network perturbation

# A14 Assay and Document Provenance Curation

## Objective

Retrieve assay-level and document-level provenance for candidate compound-protein activity records and determine whether apparently duplicated Kd and ED50 measurements represented independent experimental evidence.

## Rationale

Several candidate compound-protein combinations contained numerically identical Kd and ED50 values.

Counting these values as independent evidence could artificially increase the apparent support for individual compound-protein relationships.

Assay and source-document metadata were therefore retrieved from ChEMBL before candidate prioritisation.

## Input

results/A10_activity_quality_controlled.csv

Candidate molecules were restricted using:

results/A10_candidate_molecule_ids.csv

## Data Sources

ChEMBL Assay API

ChEMBL Document API

## Assay Metadata Retrieved

For each unique ChEMBL assay identifier, the following information was retrieved:

- Assay ChEMBL identifier
- Assay type
- Assay description
- Assay organism
- Tissue
- Cell type
- Subcellular fraction
- BAO format
- BAO label
- Target-confidence score
- ChEMBL target identifier
- Target relationship type
- Source document identifier
- Source assay identifier
- Assay category
- Taxonomy identifier

## Document Metadata Retrieved

For each source document, the following information was retrieved:

- Document ChEMBL identifier
- Document type
- Title
- Abstract
- Authors
- Journal
- Publication year
- Volume
- Issue
- Page range
- DOI
- PubMed identifier
- Patent identifier

## ChEMBL Target-Confidence Interpretation

ChEMBL assay-to-target confidence scores were used to assess the reliability of each target assignment.

A confidence score of 9 represents a high-confidence assignment to a single protein target.

Lower confidence scores may represent protein families, protein complexes, homologous proteins, non-molecular targets, or incompletely curated target assignments.

## Duplicate-Evidence Assessment

Records were grouped according to:

- ChEMBL molecule identifier
- Queried protein
- Standard activity value
- Standard activity units

Within each group, the following were compared:

- Activity types
- Activity identifiers
- Assay identifiers
- Document identifiers
- Assay descriptions
- Assay types
- Target-confidence scores
- Target relationship types
- DOI
- PubMed identifier

The following flags were created:

- Same value across multiple activity types
- Likely same-assay reannotation
- Likely same-document evidence

## Outputs

results/A14_assay_annotations.csv

results/A14_document_annotations.csv

results/A14_activity_provenance_enriched.csv

results/A14_duplicate_value_review.csv

results/A14_provenance_summary.csv

## Results

Candidate activity records curated: 61

Unique assays: 35

Assays successfully retrieved: 35

Unique source documents: 3

Documents successfully retrieved: 3

Binding-assay records: 61

Functional-assay records: 0

High-confidence target assignments with confidence score 8 or higher: 61

Single-protein target assignments with confidence score 9: 61

Duplicate-value groups requiring review: 16

Identical values represented under multiple activity types: 16

Likely same-assay reannotations: 16

Likely same-document evidence groups: 16

## Duplicate Kd and ED50 Findings

Sixteen compound-protein combinations contained identical Kd and ED50 values.

Every duplicated pair shared:

- The same molecule
- The same target protein
- The same numerical activity value
- The same activity units
- The same ChEMBL assay identifier
- The same ChEMBL source document

The affected compounds were:

- CHEMBL3752910
- CHEMBL5653589

The affected proteins included:

- EEF1G
- EEF1B2
- RPS2
- RPS3
- RPS3A
- RPL3
- RPL4
- RPL7
- RPL18A

## Interpretation

The paired Kd and ED50 records were not considered independent experimental evidence.

The records were treated as alternative ChEMBL representations of the same underlying assay result.

Consequently:

- Duplicate Kd and ED50 values must contribute only once to evidence counts.
- ED50 records from these assays must not be interpreted as independent functional validation.
- Assay count and document count must be calculated after collapsing duplicate representations.
- Candidate ranking must reward independent assays and publications rather than raw ChEMBL record count.

All 61 records were associated with binding-type assays.

No candidate currently has independent functional-assay evidence within the curated ChEMBL dataset.

All 61 target assignments received a confidence score of 9, indicating high-confidence assignment to single protein targets.

## Implication for Candidate Classification

The compounds are currently supported primarily by target-binding evidence.

The following terminology should therefore be used:

- Direct binder
- Binding candidate
- Compound with inhibitory-potency annotation from a binding assay
- Translation-network binding candidate

The following terminology should not yet be used without additional functional evidence:

- Confirmed functional inhibitor
- Confirmed translation inhibitor
- Confirmed disruptor of a specific eEF1A protein-protein interaction

# A12B Provenance-Corrected Compound-Protein Relationships

## Objective

Correct the compound-protein evidence layer using the assay-level and document-level provenance identified during A14.

## Rationale

A14 demonstrated that 16 pairs of Kd and ED50 records represented alternative ChEMBL annotations of the same underlying assay results.

Each duplicate pair shared:

- The same compound
- The same protein
- The same numerical value
- The same activity units
- The same assay identifier
- The same source document

Treating these records as independent evidence would have artificially increased the apparent support for CHEMBL3752910 and CHEMBL5653589.

## Method

Quantitative activity records were grouped according to:

- ChEMBL molecule identifier
- Protein
- Corrected activity type
- Assay identifier
- Document identifier
- Activity value
- Activity units

Where Kd and ED50 records had been confirmed as duplicate representations, the records were collapsed into a single evidence record.

The Kd interpretation was retained because A14 established that all candidate activity records were derived from ChEMBL binding assays.

## Corrected Relationship Types

### BINDS_TO

Used for quantitative binding-affinity evidence.

Duplicate Kd and ED50 records from the same underlying binding assay were represented as a single BINDS_TO relationship.

### HAS_INHIBITORY_ACTIVITY_AGAINST

Used for IC50-labelled evidence obtained from binding-type assays.

This relationship replaced the previous INHIBITS relationship to avoid claiming a confirmed molecular mechanism from an activity annotation alone.

### HAS_ACTIVITY_AGAINST

Reserved for other quantitative evidence that could not be classified as direct binding or inhibitory activity.

## Relationship Properties

The corrected relationships retained:

- Original activity types
- Corrected activity type
- Corrected evidence class
- Activity identifiers
- Assay identifiers
- Document identifiers
- Minimum quantitative activity
- Median quantitative activity
- Maximum pChEMBL value
- Independent assay count
- Independent document count
- ChEMBL target-confidence score
- Assay type
- Assay description
- DOI
- PubMed identifier
- Raw record count
- Duplicate-collapse flag
- Molecule type
- Maximum clinical phase
- Development status

## Outputs

results/A12B_compound_protein_evidence_corrected.csv

results/edges_compound_protein_corrected.csv

results/A12B_correction_summary.csv

## Results

Quantitative records before provenance correction: 39

Corrected evidence records: 23

Corrected compound-protein graph relationships: 23

Duplicated Kd and ED50 groups collapsed: 16

Corrected relationship composition:

- BINDS_TO: 19
- HAS_INHIBITORY_ACTIVITY_AGAINST: 4

No independent functional-response relationships remained after provenance correction.

## Interpretation

The provenance correction reduced the apparent evidence volume without removing any unique compound-protein combinations.

The corrected network therefore represents 23 unique compound-protein evidence relationships supported by assay-level provenance.

CHEMBL3752910 and CHEMBL5653589 remain translation-network binding candidates, but the duplicated ED50 records no longer provide additional functional evidence.

Molibresib retains four inhibitory-activity relationships, although these are described cautiously because the underlying assays were classified as binding assays.

# A13B Provenance-Corrected Compound Network

## Objective

Regenerate the interactive eEF1A-centred network using the provenance-corrected compound-protein relationships.

## Input

results/edges_compound_protein_corrected.csv

## Corrections Applied

The updated graph:

- Removed duplicated Kd and ED50 representations
- Removed unsupported functional-response relationships
- Reclassified duplicate Kd and ED50 records as binding evidence
- Replaced INHIBITS with HAS_INHIBITORY_ACTIVITY_AGAINST
- Added assay and document provenance to relationship tooltips
- Added a duplicate-collapse indicator
- Updated the network key to reflect the corrected evidence classes

## Corrected Compound-Protein Relationships

- Purple line: BINDS_TO
- Red line: HAS_INHIBITORY_ACTIVITY_AGAINST
- Orange line: HAS_ACTIVITY_AGAINST

## Output

results/eef1a_network_A13B_corrected.html

# A16 Chemical-Structure Similarity Analysis 

 ## Objective 

 Compare the six network candidates with established eEF1A ligands and determine whether the candidates occupy related chemical-structure space. 

 The analysis examined: 



Whether CHEMBL1802814, CHEMBL1802815, and CHEMBL1802973 form a structural series
Which established ligand is structurally nearest to each candidate
Whether Molibresib occupies distinct chemical space
Whether CHEMBL3752910 and CHEMBL5653589 are structurally related
Whether any candidate represents a potentially distinct scaffold 

 # A16A Reference-Structure Retrieval 

 ## Input 

 results/A11_candidate_molecule_annotations.csv 

 ## Compounds 

 Network candidates: 


Molibresib, CHEMBL1232461
CHEMBL1802814
CHEMBL1802815
CHEMBL1802973
CHEMBL3752910
CHEMBL5653589 

 Reference ligands: 


Plitidepsin
Didemnin B
Ternatin-4
Nannocystin A 

 ## Method 

 Candidate structures were obtained from the A11 ChEMBL annotations. 

 Reference structures were retrieved from PubChem using PUG REST. 

 The retrieved information included: 


PubChem CID
Compound title
IUPAC name
Canonical and isomeric SMILES
InChI and InChIKey
Molecular formula
Molecular weight 

 Isomeric SMILES were used for reference compounds where available. Failed requests were retried, and the analysis was prevented from continuing if a required structure was missing. 

 ## Outputs 

 results/A16_reference_compounds.csv 

 results/A16_structures_combined.csv 

 results/A16_structure_retrieval_summary.csv 

 ## Results 

 Candidate compounds: 6 

 Candidate structures available: 6 

 Reference ligands requested: 4 

 Reference structures retrieved: 4 

 Combined structures available: 10 

 Missing structures: 0 

 # A16B Chemical-Similarity Analysis 

 ## Input 

 results/A16_structures_combined.csv 

 ## Method



 A16B was implemented in R using RDKit through the reticulate package and the resko-a16 Conda environment. 

 All structures were parsed with RDKit before analysis. 

 Morgan fingerprints were generated using: 



Radius: 2
Fingerprint size: 2,048 bits
Chirality: enabled 

 Pairwise structural similarity was calculated using the Tanimoto coefficient. 

 Every candidate was compared with all four reference ligands. References were ranked by similarity, and the nearest reference was identified for each candidate.



 Structural clustering used: 



Tanimoto distance
Average linkage
Cluster cut height: 0.60 

 A two-dimensional chemical-space projection was produced using classical multidimensional scaling. 

 RDKit was also used to calculate: 


Molecular weight
Calculated logP
Topological polar surface area 

 ## Outputs 

 results/A16_pairwise_tanimoto.csv 

 results/A16_candidate_reference_similarity.csv 

 results/A16_structural_clusters.csv 

 results/A16_similarity_heatmap.png 

 results/A16_chemical_space.png



 results/A16_structure_grid.png 

 results/A16_similarity_summary.csv 

 results/A16_chemical_similarity_report.html 

 ## Results 

 Compounds analysed: 10 

 Pairwise comparisons: 45 

 Candidate-reference comparisons: 24 

 Structural clusters identified: 6 

 Structures successfully parsed: 10 

 Required outputs created: 8 of 8 

 ## Interpretation 

 The six network candidates and four reference ligands occupied heterogeneous chemical-structure space and were separated into six clusters under the selected analysis settings. 

 The result indicates that the compounds do not represent a single uniform structural group. 

 Candidate-specific relationships, including nearest reference ligands and similarities among the direct EEF1A1 binders, are retained in the A16 output tables. 

 A candidate with low similarity to the four reference ligands may represent a distinct scaffold within this comparison set. A broader database search is required before claiming global scaffold novelty. 

 ## Limitations 

 The analysis included only ten compounds and used one fingerprint representation. 

 Similarity values depend on molecular standardisation, fingerprint settings, stereochemistry, molecular size, and the selected reference set. 

 The analysis did not compare protonation states, tautomers, three-dimensional conformations, pharmacophores, or protein-binding poses. 

 Chemical similarity does not establish: 



Shared binding to EEF1A
A shared binding site
A shared mechanism
Functional inhibition
Comparable selectivity
Comparable efficacy or toxicity 

 ## Status 

 A16A structure retrieval: Complete 

 A16B chemical-similarity analysis: Complete 

 Large-scale database expansion: Not started 

 SIDER remapping: Not started 

 HPC-scale analysis: Not started
Main improvements
Reduced length: Removed implementation details that are not essential to the main Methods record.
Tighter structure: Combined related validation and methodological points.
Scientific precision: Retained fingerprint settings, validated counts, provenance, and the principal interpretation limits.

# A17A Preparation of the SIDER Reference Library 

 ## Objective 

 Prepare a validated local copy of SIDER 4.1 for candidate mapping and later side-effect similarity analysis. 

 ## Method 

 Seven official SIDER 4.1 source files were downloaded and preserved in: 

 results/A17A_sider_4.1_raw/ 

 The files contained drug names, ATC codes, MedDRA side effects, side-effect frequencies, indications, and MedDRA mappings. 

 Side effects and indications were restricted to MedDRA preferred terms to reduce duplication from lower-level terms. 

 Drug identifiers were retained as STITCH flat and stereo identifiers. Source URLs, release information, file sizes, checksums, and retrieval status were recorded in a data manifest. 

 ## Outputs 

 results/A17A_sider_drugs.csv 

 results/A17A_sider_side_effects_pt.csv 

 results/A17A_sider_frequencies.csv 

 results/A17A_sider_indications_pt.csv 

 results/A17A_sider_data_manifest.csv 

 results/A17A_sider_summary.csv 

 ## Results 

 SIDER drugs: 1,430 

 Stereo-specific drug-side-effect relationships: 152,759 

 Unique flat-drug/preferred-term pairs: 145,321 

 Unique stereo-drug/preferred-term pairs: 152,759 

 Unique MedDRA preferred terms: 4,251 

 Preferred-term frequency records: 154,507 

 Preferred-term indication records: 15,560 

 Source files preserved: 7 

 Processed outputs validated: 6 

 ## Interpretation 

 A complete SIDER 4.1 reference library was prepared successfully. 

 The flat-drug/preferred-term representation will be used for primary side-effect similarity calculations. Stereo-specific records will be retained for identifier provenance and compound-form review. 

 Absence from SIDER will be interpreted as not represented in SIDER 4.1, not as evidence that a compound has no side effects. 

 ## Limitations 

 SIDER 4.1 is a historical resource released in 2015 and primarily represents marketed medicines. Side-effect frequency information is incomplete, and STITCH identifiers require careful reconciliation with PubChem, ChEMBL, salts, stereoisomers, and parent structures.

 # A17E Provenance Correction of Expanded SIDER Candidates 

 ## Objective 

 Validate the provisional network relationships assigned to the six SIDER chemical-neighbour candidates during A17D. 

 ## Method 

 A17D target records were read using character-only column parsing to preserve ChEMBL identifiers and missing network annotations. 

 Protein-network membership was recalculated using exact normalised protein names from the RESKO network input. 

 ChEMBL assay and document metadata were retrieved using filtered API endpoints and stable identifiers. 

 Evidence records were grouped by: 



Candidate compound
Target
Activity type
Activity value
Activity units
Assay identifier
Document identifier 

 Independent assays and documents were counted separately from raw activity records. 

 ## Outputs 

 results/A17E_relevant_activity_records.csv 

 results/A17E_assay_metadata.csv 

 results/A17E_document_metadata.csv 

 results/A17E_provenance_corrected_evidence.csv 

 results/A17E_candidate_classification_corrected.csv 

 results/A17E_correction_summary.csv 

 results/A17E_provenance_report.html 

 ## Results 

 Provisional network-related records reviewed: 3 

 Provenance-corrected evidence relationships: 3 

 Independent assays: 1 

 Independent documents: 1



 Candidate classifications corrected: 3 

 Weak 30 micromolar relationships: 3 

 Direct eEF1A candidates: 0 

 eEF1-complex candidates: 0 

 Translation-network candidates: 0 

 Chemical and side-effect comparators: 6 

 Assays with unavailable metadata: 0 

 Documents with unavailable metadata: 0 

 ## Interpretation 

 Nilotinib, imatinib, and ponatinib each contained a 30,000 nM Kd annotation for GTP-binding nuclear protein Ran. 

 The three records shared the same ChEMBL assay and source document and therefore represented one shared experimental provenance context rather than three independent pathway-level findings. 

 Exact network matching established that RAN was not a member of the defined RESKO protein network. The provisional translation-network classifications were therefore removed. 

 No direct EEF1A1, EEF1A2, eEF1-complex, or validated translation-network evidence remained for the six SIDER candidates. 

 Nilotinib, imatinib, ponatinib, alprazolam, triazolam, and temazepam were retained as chemical and side-effect comparators only. 

 ## Limitations 

 The correction evaluated database provenance and exact membership in the current RESKO protein set. 

 Absence from the defined network does not exclude broader biological relationships with translation or cellular stress pathways. 

 The 30 micromolar RAN annotations represent weak binding and do not establish cellular target engagement, functional inhibition, selectivity, clinical efficacy, or a shared eEF1A mechanism.

 # A18A Compound Detail and Supplier Manifest  

  ## Objective  

 Create a unified compound identity, scientific-evidence, clinical-annotation, and commercial-availability schema for the compounds currently represented in the RESKO workflow. 

 Prepare the compound records required for supplier searches and future display within the interactive RESKO network. 

  ## Method  

 Ten A16 compounds and six A17 SIDER-derived chemical comparators were integrated into a unified compound-detail manifest. 

 The ten A16 compounds consisted of six existing RESKO candidates and four established eEF1A reference ligands. 

 The six A17 compounds consisted of nilotinib, imatinib, ponatinib, alprazolam, triazolam, and temazepam. 

 Compound identities were reconciled using: 



 RESKO compound identifier Preferred compound name Normalised compound name ChEMBL molecule identifier PubChem compound identifier SIDER STITCH identifier Full standard InChIKey Parent connectivity InChIKey Canonical SMILES Analysis SMILES 

 Biological classifications and progression statuses were incorporated from the provenance-corrected A17E results. 

 Chemical-similarity relationships, nearest RESKO query compounds, Tanimoto similarities, SIDER representation, side-effect counts, indications, maximum clinical phases, and first-approval years were retained where available. 

 A supplier-lookup manifest was constructed using the following identifier priority: 



 Full standard InChIKey Parent connectivity InChIKey PubChem CID Canonical SMILES Compound name 

 Commercial-availability fields were created as placeholders and were kept separate from the biological and chemical evidence classifications. 

  ## Outputs  

 results/A18A_compound_detail_manifest.csv 

 results/A18A_compound_identifiers.csv 

 results/A18A_compound_evidence_summary.csv 

 results/A18A_supplier_lookup_manifest.csv 

 results/A18A_data_completeness_review.csv 

 results/A18A_summary.csv 

  ## Results  

 Compounds integrated: 16 

 A16 compounds: 10 

 A17 chemical and side-effect comparators: 6 

 Compounds with ChEMBL identifiers: 12 

 Compounds with PubChem CIDs in the input files: 10 

 Compounds with full standard InChIKeys: 16 

 Compounds with parent connectivity InChIKeys: 16 

 Compounds with analysis SMILES: 10 

 Compounds represented in SIDER 4.1: 6 

 Compounds with clinical-phase annotations: 6 

 Compounds ready for supplier lookup: 16 

 Compounds requiring critical identifier review: 0 

 Compounds ready for live-network detail integration: 16 

 Supplier checks completed during A18A: 0 

  ## Interpretation  

 A18A established a stable identity layer for the 16 compounds already present in the RESKO workflow. 

 The complete InChIKey coverage allowed all compounds to proceed to structure-based supplier lookup without relying exclusively on compound names, synonyms, or development codes. 

 The A18A compound count did not represent 16 newly discovered compounds. 

 The manifest contained ten compounds from A16 and six chemical and side-effect comparators from A17. 

 The six A17 comparators remained classified as chemical and side-effect comparators because A17E identified no validated direct eEF1A, eEF1-complex, or RESKO-network evidence. 

 Commercial availability remained a separate experimental-feasibility layer and did not alter the scientific classifications. 

  ## Limitations  

 A18A did not perform supplier searches or verify commercial products. 

 PubChem CIDs and clinical-phase annotations were incomplete in the source files for some compounds, although all compounds had sufficient structural identifiers for subsequent reconciliation. 

 Parent connectivity matching may combine different stereoisomers, salts, protonation states, or related forms and therefore cannot independently establish exact commercial-product identity. 



 # A18B PubChem Commercial-Availability Enrichment  

  ## Objective  

 Determine whether the 16 compounds in the A18A manifest had publicly accessible commercial-vendor information and create product-level records for subsequent identity reconciliation and live-network display. 

  ## Method  

 The 16 A18A compounds were mapped to PubChem using existing PubChem CIDs or full standard InChIKeys. 

 PubChem Chemical Vendors records were retrieved using the resolved PubChem compound identifiers. 

 Each compound was processed independently, and API responses were cached locally to support reproducibility and interrupted-run recovery. 

 A checkpoint was written after every compound lookup. 

 Public vendor information was converted into compound-level and product-level records containing: 



 RESKO compound identifier Compound name PubChem CID Vendor-section heading Information label Public listing text Product URL PubChem section path Date checked Commercial-data source Manual-verification status 

 Public vendor associations were classified as preliminary commercial-availability signals. 

 Exact product identity, chemical form, stereochemistry, purity, pack size, price, stock status, and lead time were not treated as verified. 

  ## Outputs  

 results/A18B_compound_commercial_summary.csv 

 results/A18B_commercial_products.csv 

 results/A18B_supplier_directory.csv 

 results/A18B_identity_review.csv 

 results/A18B_summary.csv 

 results/A18B_commercial_availability_report.html 

 results/A18B_lookup_checkpoint.csv 

 results/A18B_api_cache/ 

  ## Results  

 Compounds processed: 16 

 PubChem CIDs resolved: 16 

 Compounds with public vendor information: 11 

 Compounds without public vendor information: 5 

 Commercial product records: 11 

 Supplier-directory records: 1 

 Validated final outputs: 6 

  ## Interpretation  

 All 16 compounds were successfully reconciled to PubChem. 

 Public commercial information was identified for 11 compounds. 

 The five compounds without public vendor information remained structurally resolved and may therefore be investigated through specialist supplier databases, institutional procurement tools, or manual searches. 

 The 11 commercial records represented public PubChem vendor associations rather than independently confirmed products. 

 The supplier-directory count of one indicated that the parsed PubChem records shared a common section heading or grouping label. The count did not necessarily indicate that all products were offered by one verified supplier. 

 Product-level supplier identity and exact chemical form therefore required further reconciliation during A18C. 

  ## Limitations  

 Public vendor information may change after the recorded lookup date. 

 PubChem vendor associations do not independently confirm current stock, pricing, delivery time, purity, pack size, or suitability for biological experiments. 

 Commercial listings may correspond to salts, hydrates, solvates, formulations, stereoisomers, mixtures, isotopically labelled materials, or analytical standards rather than the exact scientific parent compound. 

 The vendor-section heading returned by PubChem may represent a database section or grouping label rather than the final supplier identity. 

 Manual verification against individual product pages remains necessary before procurement.

 # A18C Commercial Product Identity Reconciliation  

  ## Objective  

 Reconcile the public commercial listings identified during A18B with the scientific compound identities defined in A18A. 

 Determine whether the public listings represented exact scientific compounds, alternative chemical forms, analytical materials, mixtures, or records requiring manual review. 

  ## Method  

 The 11 public commercial records identified during A18B were linked to the corresponding A18A compound identifiers. 

 Product descriptions, information labels, PubChem section paths, and vendor-section headings were searched for terminology indicating different product forms. 

 Commercial records were classified using the following identity categories: 



 Public listing requiring identity review Salt form requiring review Solvate or hydrate requiring review Stereochemical identity requiring review Mixture or formulation Analytical or reference standard Isotopically labelled material 

 Full standard InChIKeys, parent connectivity InChIKeys, and canonical SMILES from A18A were retained beside the commercial records to support subsequent manual comparison. 

 No public listing was classified as an exact product match unless the exact structure and product form could be independently verified. 

 Product availability was kept separate from biological, chemical-similarity, and clinical-development evidence. 

  ## Outputs  

 results/A18C_commercial_products_reconciled.csv 

 results/A18C_commercial_identity_review.csv 

 results/A18C_compound_commercial_summary_reconciled.csv 

 results/A18C_summary.csv 

 results/A18C_commercial_identity_report.html 

  ## Results  

 Compounds reviewed: 16 

 Commercial product records classified: 11 

 Compounds with public vendor information: 11 

 Compounds requiring manual product review: 11 

 Exact product identities confirmed automatically: 0 

 Public listings requiring general identity review: 11 

 Salt-form records identified automatically: 0 

 Solvate or hydrate records identified automatically: 0 

 Stereochemical-review records identified automatically: 0 

 Mixture or formulation records identified automatically: 0 

 Analytical or reference-standard records identified automatically: 0 

 Isotopically labelled records identified automatically: 0 

  ## Interpretation  

 Public commercial information was available for 11 of the 16 compounds, but none of the associated records contained sufficient verified product-level information to establish an exact commercial identity automatically. 

 All 11 public listings were therefore retained as preliminary availability signals requiring manual product-page review. 

 The absence of automatically detected salt, solvate, stereochemical, mixture, reference-standard, or isotopically labelled records did not demonstrate that these alternative forms were absent. 

 The available PubChem records did not provide sufficient detailed product descriptions for those form classifications. 

 Commercial availability remained an experimental-feasibility annotation and did not alter the biological evidence or scientific priority of any compound. 

  ## Limitations  

 Product classifications were based on the text exposed through the public PubChem Chemical Vendors records. 

 Supplier product pages may contain additional information not present in the PubChem record. 

 Exact product identity, stereochemistry, chemical form, purity, pack size, stock, price, and lead time require manual verification before procurement. 

 A public vendor association does not guarantee that a compound is currently available or suitable for biological testing.

 # A18D Live-Network Compound Detail Data Preparation  

  ## Objective  

 Prepare unified compound, scientific-evidence, side-effect, indication, and commercial-availability records for display when a compound node is selected in the interactive RESKO network. 

 Create stable compound aliases so live-network node identifiers, compound names, ChEMBL identifiers, PubChem identifiers, and InChIKeys could be mapped to the appropriate compound-detail record. 

  ## Method  

 Compound identities and baseline annotations were retrieved from the A18A compound-detail and evidence manifests. 

 Commercial availability and product-identity classifications were retrieved from the provenance-corrected A18C outputs. 

 Optional SIDER side-effect, indication, and side-effect-similarity records from A17D were summarised at compound level. 

 The compound-detail records incorporated: 



 RESKO compound identifier Node lookup identifier Preferred compound name Compound class Compound origin ChEMBL identifier PubChem CID Full standard InChIKey Parent connectivity InChIKey Canonical SMILES Analysis SMILES Nearest RESKO query compound Tanimoto similarity Biological classification Progression status Matched network proteins Independent assay count Independent document count Maximum clinical phase First-approval year SIDER representation Side-effect record count Indication record count Closest side-effect-profile neighbour Commercial-availability status Supplier count Commercial-product record count Commercial identity classification Procurement-readiness status 

 Compound aliases were created using: 



 RESKO compound identifier Preferred compound name ChEMBL identifier PubChem CID Full standard InChIKey 

 Compound details, aliases, and commercial-product records were exported as CSV and JSON files for use by the live network. 

 All outputs were read back and validated before completion was reported. 

  ## Outputs  

 results/A18D_live_compound_details.csv 

 results/A18D_live_compound_details.json 

 results/A18D_live_commercial_products.json 

 results/A18D_live_compound_aliases.csv 

 results/A18D_live_compound_aliases.json 

 results/A18D_summary.csv 

 results/A18D_live_network_data_report.html 

  ## Results  

 Compound-detail records: 16 

 Compound-alias records: 76 

 Commercial-product records: 11 

 Compounds with public vendor information: 11 

 Compounds requiring manual commercial-product review: 11 

 Compounds ready for detail-panel display: 16 

 Validated outputs: 7 

  ## Interpretation  

 A18D created the complete data layer required to display contextual information when one of the 16 current compound nodes is selected in the live RESKO network. 

 The alias table allows network node identifiers to be reconciled with compound records even when the graph uses a ChEMBL identifier, compound name, PubChem CID, or RESKO identifier. 

 Public vendor information can be displayed for 11 compounds, but all associated commercial records remain subject to manual identity review. 

 Scientific evidence, side-effect information, clinical annotations, and commercial availability remain stored as distinct fields. 

 Commercial availability therefore does not alter a compound’s biological classification or scientific priority. 

  ## Limitations  

 A18D prepared data for the interactive interface but did not itself alter the existing visNetwork HTML widget. 

 Commercial-product records remain preliminary and do not establish exact chemical form, purity, stock, price, pack size, lead time, or suitability for biological experiments. 

 Detailed side-effect and indication information is currently available primarily for the six A17 SIDER comparators. 

 Compounds not represented in SIDER 4.1 must be displayed as not represented in that release rather than as having no side effects. 

 Live-interface integration requires the existing visNetwork source script to be extended while preserving its current network construction, styling, filtering, and ranking behaviour.

 # A19A HPC Structure and Fingerprint Pipeline Validation  

  ## Objective  

 Validate the RESKO molecular-structure and chemical-similarity workflow on the Comet high-performance computing system before screening a large external compound library. 

 Confirm that the project-local RDKit environment could parse compound structures, generate chirality-aware molecular fingerprints, calculate Tanimoto similarities, write analysis outputs, and preserve computational provenance on a Comet compute node. 

  ## Method  

 The 16 compound-detail records prepared during A18D were read using character-preserving parsing. 

 For each compound, analysis SMILES were selected where available. 

 Canonical SMILES were used when analysis SMILES were unavailable. 

 Compounds without either SMILES representation were retained in the identity-review output and were not used for fingerprint calculations. 

 Structures were parsed using RDKit 2026.03.5 in a project-local Miniforge environment running Python 3.11.15. 

 Parsed structures were converted into isomeric canonical SMILES. 

 Full InChIKeys were recalculated from the parsed structures and compared with the full InChIKeys recorded during A18A. 

 Molecular fingerprints were generated using: 



 Morgan fingerprint Radius 2 2,048 bits Chirality enabled 

 Tanimoto similarities were calculated for all unique compound pairs, including self-comparisons. 

 Results were written in CSV and Parquet formats. 

 Input and output SHA-256 checksums, software versions, compute-node identity, Slurm job identifier, CPU allocation, and fingerprint parameters were recorded in a computational-provenance file. 

 The analysis was executed through Slurm using two CPU cores and 4 GB of requested memory on the Comet short_free partition. 

  ## Outputs  

 results/A19A_structure_validation.csv 

 results/A19A_validated_structures.csv 

 results/A19A_pairwise_tanimoto.csv 

 results/A19A_pairwise_tanimoto.parquet 

 results/A19A_identifier_review.csv 

 results/A19A_summary.csv 

 results/A19A_run_provenance.json 

  ## Results  

 Input compounds: 16 

 Structures parsed successfully: 10 

 Morgan fingerprints generated: 10 

 Pairwise similarity records including self-comparisons: 55 

 Non-self pairwise comparisons: 45 

 Compounds requiring structure or identifier review: 6 

 Validated outputs: 7 

 Slurm job ID: 1946564 

 Compute node: compute016.comet.hpc.ncl.ac.uk 

 Error-log size: 0 bytes 

  ## Interpretation  

 The complete RDKit structure-processing and fingerprint workflow operated successfully on a Comet compute node. 

 All compounds with an available analysis or canonical SMILES representation were parsed and fingerprinted. 

 The 55 pairwise records represented all unique comparisons among the ten fingerprinted structures, including ten self-comparisons and 45 non-self comparisons. 

 The six compounds requiring review were not automatically considered invalid. 

 These compounds retained full InChIKeys but required enrichment of their SMILES representations or reconciliation of any detected identifier differences before inclusion in the large-scale chemical-similarity screen. 

 The successful CSV and Parquet exports confirmed that the output format required for larger HPC datasets was operational. 

  ## Limitations  

 Only ten of the 16 current compounds had a SMILES representation available to the A19A input. 

 Compound records without SMILES could not be fingerprinted even when a full InChIKey was available. 

 Morgan fingerprint similarity measures local molecular environments and does not directly establish shared target binding, mechanism, potency, selectivity, cellular activity, or clinical efficacy. 

 Large and structurally complex natural products may receive lower fingerprint similarities despite retaining biologically relevant pharmacophoric relationships.