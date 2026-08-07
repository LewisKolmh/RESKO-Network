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