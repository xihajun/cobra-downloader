# cobra-downloader
To download gene list from ensembl.org using python
## pybiomart
```
pip install pybiomart
```

### options for query
- attributes

 options | name
 :--- | --- 
 ensembl_gene_id | <biomart.Attribute name='ensembl_gene_id', display_name='Gene stable ID', description='Stable ID of the Gene'>
ensembl_gene_id_version | <biomart.Attribute name='ensembl_gene_id_version', display_name='Gene stable ID version', description='Versionned stable ID of the Gene'>
ensembl_transcript_id | <biomart.Attribute name='ensembl_transcript_id', display_name='Transcript stable ID', description='Stable ID of the Transcript'>
ensembl_transcript_id_version | <biomart.Attribute name='ensembl_transcript_id_version', display_name='Transcript stable ID version', description='Versionned stable ID of the Transcript'>
ensembl_peptide_id | <biomart.Attribute name='ensembl_peptide_id', display_name='Protein stable ID', description=''>
ensembl_peptide_id_version | <biomart.Attribute name='ensembl_peptide_id_version', display_name='Protein stable ID version', description=''>
ensembl_exon_id | <biomart.Attribute name='ensembl_exon_id', display_name='Exon stable ID', description=''>
description | <biomart.Attribute name='description', display_name='Gene description', description=''>
chromosome_name | <biomart.Attribute name='chromosome_name', display_name='Chromosome/scaffold name', description='Chromosome/scaffold name'>
start_position | <biomart.Attribute name='start_position', display_name='Gene start (bp)', description='Start Coordinate of the gene in chromosomal coordinates.'>
end_position | <biomart.Attribute name='end_position', display_name='Gene end (bp)', description='End Coordinate of the gene in chromosomal coordinates.'>
strand | <biomart.Attribute name='strand', display_name='Strand', description='Orientation of genes on chromosome/scaffold.'>
band | <biomart.Attribute name='band', display_name='Karyotype band', description='Band of chromosome/scaffold where genes occur.'>
transcript_start | <biomart.Attribute name='transcript_start', display_name='Transcript start (bp)', description='Start Coordinate of the transcript in chromosomal coordinates.'>
transcript_end | <biomart.Attribute name='transcript_end', display_name='Transcript end (bp)', description='End Coordinate of the gene in chromosomal coordinates.'>
transcription_start_site | <biomart.Attribute name='transcription_start_site', display_name='Transcription start site (TSS)', description=''>
transcript_length | <biomart.Attribute name='transcript_length', display_name='Transcript length (including UTRs and CDS)', description=''>
transcript_tsl | <biomart.Attribute name='transcript_tsl', display_name='Transcript support level (TSL)', description=''>
transcript_gencode_basic | <biomart.Attribute name='transcript_gencode_basic', display_name='GENCODE basic annotation', description=''>
transcript_appris | <biomart.Attribute name='transcript_appris', display_name='APPRIS annotation', description=''>
transcript_mane_select | <biomart.Attribute name='transcript_mane_select', display_name='RefSeq match transcript', description=''>
external_gene_name | <biomart.Attribute name='external_gene_name', display_name='Gene name', description=''>
external_gene_source | <biomart.Attribute name='external_gene_source', display_name='Source of gene name', description=''>
external_transcript_name | <biomart.Attribute name='external_transcript_name', display_name='Transcript name', description=''>
external_transcript_source_name | <biomart.Attribute name='external_transcript_source_name', display_name='Source of transcript name', description=''>
transcript_count | <biomart.Attribute name='transcript_count', display_name='Transcript count', description=''>
percentage_gene_gc_content | <biomart.Attribute name='percentage_gene_gc_content', display_name='Gene % GC content', description='Percentage gene GC content'>
gene_biotype | <biomart.Attribute name='gene_biotype', display_name='Gene type', description=''>
transcript_biotype | <biomart.Attribute name='transcript_biotype', display_name='Transcript type', description=''>
source | <biomart.Attribute name='source', display_name='Source (gene)', description=''>
transcript_source | <biomart.Attribute name='transcript_source', display_name='Source (transcript)', description=''>
version | <biomart.Attribute name='version', display_name='Version (gene)', description=''>
transcript_version | <biomart.Attribute name='transcript_version', display_name='Version (transcript)', description=''>
peptide_version | <biomart.Attribute name='peptide_version', display_name='Version (protein)', description=''>
external_synonym | <biomart.Attribute name='external_synonym', display_name='Gene Synonym', description=''>
phenotype_description | <biomart.Attribute name='phenotype_description', display_name='Phenotype description', description=''>
Source_name | <biomart.Attribute name='Source_name', display_name='Source name', description=''>
study_external_id | <biomart.Attribute name='study_external_id', display_name='Study external reference', description=''>
strain_name | <biomart.Attribute name='strain_name', display_name='Strain name', description=''>
strain_gender | <biomart.Attribute name='strain_gender', display_name='Strain gender', description=''>
p_value | <biomart.Attribute name='p_value', display_name='P value', description=''>
go_id | <biomart.Attribute name='go_id', display_name='GO term accession', description=''>
name_1006 | <biomart.Attribute name='name_1006', display_name='GO term name', description=''>
definition_1006 | <biomart.Attribute name='definition_1006', display_name='GO term definition', description=''>
go_linkage_type | <biomart.Attribute name='go_linkage_type', display_name='GO term evidence code', description=''>
namespace_1003 | <biomart.Attribute name='namespace_1003', display_name='GO domain', description=''>
goslim_goa_accession | <biomart.Attribute name='goslim_goa_accession', display_name='GOSlim GOA Accession(s)', description=''>
goslim_goa_description | <biomart.Attribute name='goslim_goa_description', display_name='GOSlim GOA Description', description=''>
ccds | <biomart.Attribute name='ccds', display_name='CCDS ID', description=''>
chembl | <biomart.Attribute name='chembl', display_name='ChEMBL ID', description=''>
clone_based_ensembl_gene | <biomart.Attribute name='clone_based_ensembl_gene', display_name='Clone-based (Ensembl) gene ID', description=''>
clone_based_ensembl_transcript | <biomart.Attribute name='clone_based_ensembl_transcript', display_name='Clone-based (Ensembl) transcript ID', description=''>
dbass3_name | <biomart.Attribute name='dbass3_name', display_name="DataBase of Aberrant 3' Splice Sites name", description=''>
dbass3_id | <biomart.Attribute name='dbass3_id', display_name="DataBase of Aberrant 3' Splice Sites ID", description=''>
dbass5_name | <biomart.Attribute name='dbass5_name', display_name="DataBase of Aberrant 5' Splice Sites name", description=''>
dbass5_id | <biomart.Attribute name='dbass5_id', display_name="DataBase of Aberrant 5' Splice Sites ID", description=''>
ens_hs_transcript | <biomart.Attribute name='ens_hs_transcript', display_name='Ensembl Human Transcript ID', description=''>
ens_hs_translation | <biomart.Attribute name='ens_hs_translation', display_name='Ensembl Human Translation ID', description=''>
entrezgene_trans_name | <biomart.Attribute name='entrezgene_trans_name', display_name='EntrezGene transcript name ID', description=''>
embl | <biomart.Attribute name='embl', display_name='European Nucleotide Archive ID', description=''>
arrayexpress | <biomart.Attribute name='arrayexpress', display_name='Expression Atlas ID', description=''>
genedb | <biomart.Attribute name='genedb', display_name='GeneDB ID', description=''>
hgnc_id | <biomart.Attribute name='hgnc_id', display_name='HGNC ID', description=''>
hgnc_symbol | <biomart.Attribute name='hgnc_symbol', display_name='HGNC symbol', description=''>
hpa_accession | <biomart.Attribute name='hpa_accession', display_name='Human Protein Atlas accession', description=''>
hpa_id | <biomart.Attribute name='hpa_id', display_name='Human Protein Atlas ID', description=''>
protein_id | <biomart.Attribute name='protein_id', display_name='INSDC protein ID', description=''>
kegg_enzyme | <biomart.Attribute name='kegg_enzyme', display_name='KEGG Pathway and Enzyme ID', description=''>
ens_lrg_gene | <biomart.Attribute name='ens_lrg_gene', display_name='LRG display in Ensembl gene ID', description=''>
ens_lrg_transcript | <biomart.Attribute name='ens_lrg_transcript', display_name='LRG display in Ensembl transcript ID', description=''>
merops | <biomart.Attribute name='merops', display_name='MEROPS - the Peptidase Database ID', description=''>
mim_gene_description | <biomart.Attribute name='mim_gene_description', display_name='MIM gene description', description=''>
mim_gene_accession | <biomart.Attribute name='mim_gene_accession', display_name='MIM gene accession', description=''>
mim_morbid_description | <biomart.Attribute name='mim_morbid_description', display_name='MIM morbid description', description=''>
mim_morbid_accession | <biomart.Attribute name='mim_morbid_accession', display_name='MIM morbid accession', description=''>
mirbase_accession | <biomart.Attribute name='mirbase_accession', display_name='miRBase accession', description=''>
mirbase_id | <biomart.Attribute name='mirbase_id', display_name='miRBase ID', description=''>
mirbase_trans_name | <biomart.Attribute name='mirbase_trans_name', display_name='miRBase transcript name ID', description=''>
entrezgene_description | <biomart.Attribute name='entrezgene_description', display_name='NCBI gene (formerly Entrezgene) description', description=''>
entrezgene_accession | <biomart.Attribute name='entrezgene_accession', display_name='NCBI gene (formerly Entrezgene) accession', description=''>
entrezgene_id | <biomart.Attribute name='entrezgene_id', display_name='NCBI gene (formerly Entrezgene) ID', description=''>
pdb | <biomart.Attribute name='pdb', display_name='PDB ID', description=''>
reactome | <biomart.Attribute name='reactome', display_name='Reactome ID', description=''>
reactome_gene | <biomart.Attribute name='reactome_gene', display_name='Reactome gene ID', description=''>
reactome_transcript | <biomart.Attribute name='reactome_transcript', display_name='Reactome transcript ID', description=''>
refseq_mrna | <biomart.Attribute name='refseq_mrna', display_name='RefSeq mRNA ID', description=''>
refseq_mrna_predicted | <biomart.Attribute name='refseq_mrna_predicted', display_name='RefSeq mRNA predicted ID', description=''>
refseq_ncrna | <biomart.Attribute name='refseq_ncrna', display_name='RefSeq ncRNA ID', description=''>
refseq_ncrna_predicted | <biomart.Attribute name='refseq_ncrna_predicted', display_name='RefSeq ncRNA predicted ID', description=''>
refseq_peptide | <biomart.Attribute name='refseq_peptide', display_name='RefSeq peptide ID', description=''>
refseq_peptide_predicted | <biomart.Attribute name='refseq_peptide_predicted', display_name='RefSeq peptide predicted ID', description=''>
rnacentral | <biomart.Attribute name='rnacentral', display_name='RNAcentral ID', description=''>
hgnc_trans_name | <biomart.Attribute name='hgnc_trans_name', display_name='Transcript name ID', description=''>
ucsc | <biomart.Attribute name='ucsc', display_name='UCSC Stable ID', description=''>
uniparc | <biomart.Attribute name='uniparc', display_name='UniParc ID', description=''>
uniprot_gn_symbol | <biomart.Attribute name='uniprot_gn_symbol', display_name='UniProtKB Gene Name symbol', description=''>
uniprot_gn_id | <biomart.Attribute name='uniprot_gn_id', display_name='UniProtKB Gene Name ID', description=''>
uniprotswissprot | <biomart.Attribute name='uniprotswissprot', display_name='UniProtKB/Swiss-Prot ID', description=''>
uniprotsptrembl | <biomart.Attribute name='uniprotsptrembl', display_name='UniProtKB/TrEMBL ID', description=''>
wikigene_description | <biomart.Attribute name='wikigene_description', display_name='WikiGene description', description=''>
wikigene_name | <biomart.Attribute name='wikigene_name', display_name='WikiGene name', description=''>
wikigene_id | <biomart.Attribute name='wikigene_id', display_name='WikiGene ID', description=''>
affy_hc_g110 | <biomart.Attribute name='affy_hc_g110', display_name='AFFY HC G110 probe', description=''>
affy_hg_focus | <biomart.Attribute name='affy_hg_focus', display_name='AFFY HG Focus probe', description=''>
affy_hg_u133a | <biomart.Attribute name='affy_hg_u133a', display_name='AFFY HG U133A probe', description=''>
affy_hg_u133a_2 | <biomart.Attribute name='affy_hg_u133a_2', display_name='AFFY HG U133A 2 probe', description=''>
affy_hg_u133b | <biomart.Attribute name='affy_hg_u133b', display_name='AFFY HG U133B probe', description=''>
affy_hg_u133_plus_2 | <biomart.Attribute name='affy_hg_u133_plus_2', display_name='AFFY HG U133 Plus 2 probe', description=''>
affy_hg_u95a | <biomart.Attribute name='affy_hg_u95a', display_name='AFFY HG U95A probe', description=''>
affy_hg_u95av2 | <biomart.Attribute name='affy_hg_u95av2', display_name='AFFY HG U95Av2 probe', description=''>
affy_hg_u95b | <biomart.Attribute name='affy_hg_u95b', display_name='AFFY HG U95B probe', description=''>
affy_hg_u95c | <biomart.Attribute name='affy_hg_u95c', display_name='AFFY HG U95C probe', description=''>
affy_hg_u95d | <biomart.Attribute name='affy_hg_u95d', display_name='AFFY HG U95D probe', description=''>
affy_hg_u95e | <biomart.Attribute name='affy_hg_u95e', display_name='AFFY HG U95E probe', description=''>
affy_hta_2_0 | <biomart.Attribute name='affy_hta_2_0', display_name='AFFY HTA 2 0 probe', description=''>
affy_huex_1_0_st_v2 | <biomart.Attribute name='affy_huex_1_0_st_v2', display_name='AFFY HuEx 1 0 st v2 probe', description=''>
affy_hugenefl | <biomart.Attribute name='affy_hugenefl', display_name='AFFY HuGeneFL probe', description=''>
affy_hugene_1_0_st_v1 | <biomart.Attribute name='affy_hugene_1_0_st_v1', display_name='AFFY HuGene 1 0 st v1 probe', description=''>
affy_hugene_2_0_st_v1 | <biomart.Attribute name='affy_hugene_2_0_st_v1', display_name='AFFY HuGene 2 0 st v1 probe', description=''>
affy_primeview | <biomart.Attribute name='affy_primeview', display_name='AFFY PrimeView probe', description=''>
affy_u133_x3p | <biomart.Attribute name='affy_u133_x3p', display_name='AFFY U133 X3P probe', description=''>
agilent_cgh_44b | <biomart.Attribute name='agilent_cgh_44b', display_name='AGILENT CGH 44b probe', description=''>
agilent_gpl6848 | <biomart.Attribute name='agilent_gpl6848', display_name='AGILENT GPL6848 probe', description=''>
agilent_sureprint_g3_ge_8x60k | <biomart.Attribute name='agilent_sureprint_g3_ge_8x60k', display_name='AGILENT SurePrint G3 GE 8x60k probe', description=''>
agilent_sureprint_g3_ge_8x60k_v2 | <biomart.Attribute name='agilent_sureprint_g3_ge_8x60k_v2', display_name='AGILENT SurePrint G3 GE 8x60k v2 probe', description=''>
agilent_wholegenome | <biomart.Attribute name='agilent_wholegenome', display_name='AGILENT WholeGenome probe', description=''>
agilent_wholegenome_4x44k_v1 | <biomart.Attribute name='agilent_wholegenome_4x44k_v1', display_name='AGILENT WholeGenome 4x44k v1 probe', description=''>
agilent_wholegenome_4x44k_v2 | <biomart.Attribute name='agilent_wholegenome_4x44k_v2', display_name='AGILENT WholeGenome 4x44k v2 probe', description=''>
codelink_codelink | <biomart.Attribute name='codelink_codelink', display_name='CODELINK CODELINK probe', description=''>
illumina_humanht_12_v3 | <biomart.Attribute name='illumina_humanht_12_v3', display_name='ILLUMINA HumanHT 12 V3 probe', description=''>
illumina_humanht_12_v4 | <biomart.Attribute name='illumina_humanht_12_v4', display_name='ILLUMINA HumanHT 12 V4 probe', description=''>
illumina_humanref_8_v3 | <biomart.Attribute name='illumina_humanref_8_v3', display_name='ILLUMINA HumanRef 8 V3 probe', description=''>
illumina_humanwg_6_v1 | <biomart.Attribute name='illumina_humanwg_6_v1', display_name='ILLUMINA HumanWG 6 V1 probe', description=''>
illumina_humanwg_6_v2 | <biomart.Attribute name='illumina_humanwg_6_v2', display_name='ILLUMINA HumanWG 6 V2 probe', description=''>
illumina_humanwg_6_v3 | <biomart.Attribute name='illumina_humanwg_6_v3', display_name='ILLUMINA HumanWG 6 V3 probe', description=''>
phalanx_onearray | <biomart.Attribute name='phalanx_onearray', display_name='PHALANX OneArray probe', description=''>
family | <biomart.Attribute name='family', display_name='Ensembl Protein Family ID(s)', description=''>
family_description | <biomart.Attribute name='family_description', display_name='Ensembl Family Description', description=''>
cdd | <biomart.Attribute name='cdd', display_name='CDD ID', description=''>
cdd_start | <biomart.Attribute name='cdd_start', display_name='CDD start', description=''>
cdd_end | <biomart.Attribute name='cdd_end', display_name='CDD end', description=''>
gene3d | <biomart.Attribute name='gene3d', display_name='Gene3D ID', description=''>
gene3d_start | <biomart.Attribute name='gene3d_start', display_name='Gene3D start', description=''>
gene3d_end | <biomart.Attribute name='gene3d_end', display_name='Gene3D end', description=''>
hamap | <biomart.Attribute name='hamap', display_name='HAMAP ID', description=''>
hamap_start | <biomart.Attribute name='hamap_start', display_name='HAMAP start', description=''>
hamap_end | <biomart.Attribute name='hamap_end', display_name='HAMAP end', description=''>
hmmpanther | <biomart.Attribute name='hmmpanther', display_name='PANTHER ID', description=''>
hmmpanther_start | <biomart.Attribute name='hmmpanther_start', display_name='PANTHER start', description=''>
hmmpanther_end | <biomart.Attribute name='hmmpanther_end', display_name='PANTHER end', description=''>
pfam | <biomart.Attribute name='pfam', display_name='Pfam ID', description=''>
pfam_start | <biomart.Attribute name='pfam_start', display_name='Pfam start', description=''>
pfam_end | <biomart.Attribute name='pfam_end', display_name='Pfam end', description=''>
pirsf | <biomart.Attribute name='pirsf', display_name='PIRSF ID', description=''>
pirsf_start | <biomart.Attribute name='pirsf_start', display_name='PIRSF start', description=''>
pirsf_end | <biomart.Attribute name='pirsf_end', display_name='PIRSF end', description=''>
prints | <biomart.Attribute name='prints', display_name='Prints ID', description=''>
prints_start | <biomart.Attribute name='prints_start', display_name='Prints start', description=''>
prints_end | <biomart.Attribute name='prints_end', display_name='Prints end', description=''>
scanprosite | <biomart.Attribute name='scanprosite', display_name='PROSITE patterns ID', description=''>
scanprosite_start | <biomart.Attribute name='scanprosite_start', display_name='PROSITE patterns start', description=''>
scanprosite_end | <biomart.Attribute name='scanprosite_end', display_name='PROSITE patterns end', description=''>
pfscan | <biomart.Attribute name='pfscan', display_name='PROSITE profiles ID', description=''>
pfscan_start | <biomart.Attribute name='pfscan_start', display_name='PROSITE profiles start', description=''>
pfscan_end | <biomart.Attribute name='pfscan_end', display_name='PROSITE profiles end', description=''>
sfld | <biomart.Attribute name='sfld', display_name='SFLD ID', description=''>
sfld_start | <biomart.Attribute name='sfld_start', display_name='SFLD start', description=''>
sfld_end | <biomart.Attribute name='sfld_end', display_name='SFLD end', description=''>
smart | <biomart.Attribute name='smart', display_name='SMART ID', description=''>
smart_start | <biomart.Attribute name='smart_start', display_name='SMART start', description=''>
smart_end | <biomart.Attribute name='smart_end', display_name='SMART end', description=''>
superfamily | <biomart.Attribute name='superfamily', display_name='Superfamily ID', description=''>
superfamily_start | <biomart.Attribute name='superfamily_start', display_name='Superfamily start', description=''>
superfamily_end | <biomart.Attribute name='superfamily_end', display_name='Superfamily end', description=''>
tigrfam | <biomart.Attribute name='tigrfam', display_name='TIGRFAM ID', description=''>
tigrfam_start | <biomart.Attribute name='tigrfam_start', display_name='TIGRFAM start', description=''>
tigrfam_end | <biomart.Attribute name='tigrfam_end', display_name='TIGRFAM end', description=''>
interpro | <biomart.Attribute name='interpro', display_name='Interpro ID', description='InterPro Accession ID'>
interpro_short_description | <biomart.Attribute name='interpro_short_description', display_name='Interpro Short Description', description='InterPro Short Description'>
interpro_description | <biomart.Attribute name='interpro_description', display_name='Interpro Description', description='InterPro Description'>
interpro_start | <biomart.Attribute name='interpro_start', display_name='Interpro start', description=''>
interpro_end | <biomart.Attribute name='interpro_end', display_name='Interpro end', description=''>
mobidblite | <biomart.Attribute name='mobidblite', display_name='MobiDB lite', description=''>
mobidblite_start | <biomart.Attribute name='mobidblite_start', display_name='MobiDB lite start', description=''>
mobidblite_end | <biomart.Attribute name='mobidblite_end', display_name='MobiDB lite end', description=''>
ncoils | <biomart.Attribute name='ncoils', display_name='Coiled-coils (Ncoils)', description=''>
ncoils_start | <biomart.Attribute name='ncoils_start', display_name='Coiled-coils (Ncoils) start', description=''>
ncoils_end | <biomart.Attribute name='ncoils_end', display_name='Coiled-coils (Ncoils) end', description=''>
seg | <biomart.Attribute name='seg', display_name='Low complexity (Seg)', description=''>
seg_start | <biomart.Attribute name='seg_start', display_name='Low complexity (Seg) start', description=''>
seg_end | <biomart.Attribute name='seg_end', display_name='Low complexity (Seg) end', description=''>
sifts_import | <biomart.Attribute name='sifts_import', display_name='PDB-ENSP mappings', description=''>
sifts_import_start | <biomart.Attribute name='sifts_import_start', display_name='PDB-ENSP mappings start', description=''>
sifts_import_end | <biomart.Attribute name='sifts_import_end', display_name='PDB-ENSP mappings end', description=''>
signalp | <biomart.Attribute name='signalp', display_name='Cleavage site (Signalp)', description=''>
signalp_start | <biomart.Attribute name='signalp_start', display_name='Cleavage site (Signalp) start', description=''>
signalp_end | <biomart.Attribute name='signalp_end', display_name='Cleavage site (Signalp) end', description=''>
tmhmm | <biomart.Attribute name='tmhmm', display_name='Transmembrane helices', description=''>
tmhmm_start | <biomart.Attribute name='tmhmm_start', display_name='Transmembrane helices start', description=''>
tmhmm_end | <biomart.Attribute name='tmhmm_end', display_name='Transmembrane helices end', description=''>
structure_gene_stable_id | <biomart.Attribute name='structure_gene_stable_id', display_name='', description=''>
structure_gene_stable_id_version | <biomart.Attribute name='structure_gene_stable_id_version', display_name='', description=''>
structure_gene_version | <biomart.Attribute name='structure_gene_version', display_name='', description=''>
structure_transcript_stable_id | <biomart.Attribute name='structure_transcript_stable_id', display_name='', description=''>
structure_transcript_stable_id_version | <biomart.Attribute name='structure_transcript_stable_id_version', display_name='', description=''>
structure_transcript_version | <biomart.Attribute name='structure_transcript_version', display_name='', description=''>
structure_translation_stable_id | <biomart.Attribute name='structure_translation_stable_id', display_name='', description=''>
structure_translation_stable_id_version | <biomart.Attribute name='structure_translation_stable_id_version', display_name='', description=''>
structure_peptide_version | <biomart.Attribute name='structure_peptide_version', display_name='', description=''>
structure_canonical_transcript_id | <biomart.Attribute name='structure_canonical_transcript_id', display_name='', description=''>
structure_chrom_name | <biomart.Attribute name='structure_chrom_name', display_name='', description=''>
structure_gene_chrom_start | <biomart.Attribute name='structure_gene_chrom_start', display_name='', description=''>
structure_gene_chrom_end | <biomart.Attribute name='structure_gene_chrom_end', display_name='', description=''>
structure_transcript_chrom_start | <biomart.Attribute name='structure_transcript_chrom_start', display_name='', description=''>
structure_transcript_chrom_end | <biomart.Attribute name='structure_transcript_chrom_end', display_name='', description=''>
structure_transcription_start_site | <biomart.Attribute name='structure_transcription_start_site', display_name='', description=''>
structure_transcript_length | <biomart.Attribute name='structure_transcript_length', display_name='', description=''>
structure_transcript_chrom_strand | <biomart.Attribute name='structure_transcript_chrom_strand', display_name='', description=''>
structure_external_gene_name | <biomart.Attribute name='structure_external_gene_name', display_name='', description=''>
structure_external_source_name | <biomart.Attribute name='structure_external_source_name', display_name='', description=''>
structure_5_utr_start | <biomart.Attribute name='structure_5_utr_start', display_name='', description=''>
structure_5_utr_end | <biomart.Attribute name='structure_5_utr_end', display_name='', description=''>
structure_3_utr_start | <biomart.Attribute name='structure_3_utr_start', display_name='', description=''>
structure_3_utr_end | <biomart.Attribute name='structure_3_utr_end', display_name='', description=''>
structure_cds_length | <biomart.Attribute name='structure_cds_length', display_name='', description=''>
structure_cdna_length | <biomart.Attribute name='structure_cdna_length', display_name='', description=''>
structure_peptide_length | <biomart.Attribute name='structure_peptide_length', display_name='', description=''>
struct_transcript_count | <biomart.Attribute name='struct_transcript_count', display_name='', description=''>
structure_translation_count | <biomart.Attribute name='structure_translation_count', display_name='', description=''>
structure_description | <biomart.Attribute name='structure_description', display_name='', description=''>
structure_biotype | <biomart.Attribute name='structure_biotype', display_name='', description=''>
structure_ensembl_exon_id | <biomart.Attribute name='structure_ensembl_exon_id', display_name='', description=''>
exon_chrom_start | <biomart.Attribute name='exon_chrom_start', display_name='Exon region start (bp)', description=''>
exon_chrom_end | <biomart.Attribute name='exon_chrom_end', display_name='Exon region end (bp)', description=''>
is_constitutive | <biomart.Attribute name='is_constitutive', display_name='Constitutive exon', description=''>
rank | <biomart.Attribute name='rank', display_name='Exon rank in transcript', description=''>
phase | <biomart.Attribute name='phase', display_name='Start phase', description=''>
end_phase | <biomart.Attribute name='end_phase', display_name='End phase', description=''>
cdna_coding_start | <biomart.Attribute name='cdna_coding_start', display_name='CDS start (within cDNA)', description=''>
cdna_coding_end | <biomart.Attribute name='cdna_coding_end', display_name='CDS end (within cDNA)', description=''>
genomic_coding_start | <biomart.Attribute name='genomic_coding_start', display_name='Genomic coding start', description=''>
genomic_coding_end | <biomart.Attribute name='genomic_coding_end', display_name='Genomic coding end', description=''>
structure_cds_start | <biomart.Attribute name='structure_cds_start', display_name='', description=''>
structure_cds_end | <biomart.Attribute name='structure_cds_end', display_name='', description=''>
homologs_ensembl_gene_id | <biomart.Attribute name='homologs_ensembl_gene_id', display_name='', description=''>
homologs_gene_stable_id_version | <biomart.Attribute name='homologs_gene_stable_id_version', display_name='', description=''>
homologs_gene_version | <biomart.Attribute name='homologs_gene_version', display_name='', description=''>
homologs_ensembl_transcript_id | <biomart.Attribute name='homologs_ensembl_transcript_id', display_name='', description=''>
homologs_transcript_stable_id_version | <biomart.Attribute name='homologs_transcript_stable_id_version', display_name='', description=''>
homologs_transcript_version | <biomart.Attribute name='homologs_transcript_version', display_name='', description=''>
homologs_ensembl_peptide_id | <biomart.Attribute name='homologs_ensembl_peptide_id', display_name='', description=''>
homologs_translation_stable_id_version | <biomart.Attribute name='homologs_translation_stable_id_version', display_name='', description=''>
homologs_peptide_version | <biomart.Attribute name='homologs_peptide_version', display_name='', description=''>
homologs_canonical_transcript_id | <biomart.Attribute name='homologs_canonical_transcript_id', display_name='', description=''>
homologs_chromosome_name | <biomart.Attribute name='homologs_chromosome_name', display_name='', description=''>
homologs_start_position | <biomart.Attribute name='homologs_start_position', display_name='', description=''>
homologs_end_position | <biomart.Attribute name='homologs_end_position', display_name='', description=''>
homologs_strand | <biomart.Attribute name='homologs_strand', display_name='', description=''>
homologs_band | <biomart.Attribute name='homologs_band', display_name='', description=''>
homologs_external_gene_name | <biomart.Attribute name='homologs_external_gene_name', display_name='', description=''>
homologs_external_gene_source | <biomart.Attribute name='homologs_external_gene_source', display_name='', description=''>
homologs_ensembl_CDS_length | <biomart.Attribute name='homologs_ensembl_CDS_length', display_name='', description=''>
homologs_ensembl_cDNA_length | <biomart.Attribute name='homologs_ensembl_cDNA_length', display_name='', description=''>
homologs_ensembl_peptide_length | <biomart.Attribute name='homologs_ensembl_peptide_length', display_name='', description=''>
homologs_transcript_count | <biomart.Attribute name='homologs_transcript_count', display_name='', description=''>
homologs_percentage_gc_content | <biomart.Attribute name='homologs_percentage_gc_content', display_name='', description=''>
homologs_description | <biomart.Attribute name='homologs_description', display_name='', description=''>
cabingdonii_homolog_ensembl_gene | <biomart.Attribute name='cabingdonii_homolog_ensembl_gene', display_name='Abingdon island giant tortoise gene stable ID', description=''>
cabingdonii_homolog_associated_gene_name | <biomart.Attribute name='cabingdonii_homolog_associated_gene_name', display_name='Abingdon island giant tortoise gene name', description=''>
cabingdonii_homolog_ensembl_peptide | <biomart.Attribute name='cabingdonii_homolog_ensembl_peptide', display_name='Abingdon island giant tortoise protein or transcript stable ID', description=''>
cabingdonii_homolog_chromosome | <biomart.Attribute name='cabingdonii_homolog_chromosome', display_name='Abingdon island giant tortoise chromosome/scaffold name', description=''>
cabingdonii_homolog_chrom_start | <biomart.Attribute name='cabingdonii_homolog_chrom_start', display_name='Abingdon island giant tortoise chromosome/scaffold start (bp)', description=''>
cabingdonii_homolog_chrom_end | <biomart.Attribute name='cabingdonii_homolog_chrom_end', display_name='Abingdon island giant tortoise chromosome/scaffold end (bp)', description=''>
cabingdonii_homolog_canonical_transcript_protein | <biomart.Attribute name='cabingdonii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cabingdonii_homolog_subtype | <biomart.Attribute name='cabingdonii_homolog_subtype', display_name='Last common ancestor with Abingdon island giant tortoise', description=''>
cabingdonii_homolog_orthology_type | <biomart.Attribute name='cabingdonii_homolog_orthology_type', display_name='Abingdon island giant tortoise homology type', description=''>
cabingdonii_homolog_perc_id | <biomart.Attribute name='cabingdonii_homolog_perc_id', display_name='%id. target Abingdon island giant tortoise gene identical to query gene', description=''>
cabingdonii_homolog_perc_id_r1 | <biomart.Attribute name='cabingdonii_homolog_perc_id_r1', display_name='%id. query gene identical to target Abingdon island giant tortoise gene', description=''>
cabingdonii_homolog_goc_score | <biomart.Attribute name='cabingdonii_homolog_goc_score', display_name='Abingdon island giant tortoise Gene-order conservation score', description=''>
cabingdonii_homolog_wga_coverage | <biomart.Attribute name='cabingdonii_homolog_wga_coverage', display_name='Abingdon island giant tortoise Whole-genome alignment coverage', description=''>
cabingdonii_homolog_orthology_confidence | <biomart.Attribute name='cabingdonii_homolog_orthology_confidence', display_name='Abingdon island giant tortoise orthology confidence [0 low, 1 high]', description=''>
gagassizii_homolog_ensembl_gene | <biomart.Attribute name='gagassizii_homolog_ensembl_gene', display_name="Agassiz's desert tortoise gene stable ID", description=''>
gagassizii_homolog_associated_gene_name | <biomart.Attribute name='gagassizii_homolog_associated_gene_name', display_name="Agassiz's desert tortoise gene name", description=''>
gagassizii_homolog_ensembl_peptide | <biomart.Attribute name='gagassizii_homolog_ensembl_peptide', display_name="Agassiz's desert tortoise protein or transcript stable ID", description=''>
gagassizii_homolog_chromosome | <biomart.Attribute name='gagassizii_homolog_chromosome', display_name="Agassiz's desert tortoise chromosome/scaffold name", description=''>
gagassizii_homolog_chrom_start | <biomart.Attribute name='gagassizii_homolog_chrom_start', display_name="Agassiz's desert tortoise chromosome/scaffold start (bp)", description=''>
gagassizii_homolog_chrom_end | <biomart.Attribute name='gagassizii_homolog_chrom_end', display_name="Agassiz's desert tortoise chromosome/scaffold end (bp)", description=''>
gagassizii_homolog_canonical_transcript_protein | <biomart.Attribute name='gagassizii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
gagassizii_homolog_subtype | <biomart.Attribute name='gagassizii_homolog_subtype', display_name="Last common ancestor with Agassiz's desert tortoise", description=''>
gagassizii_homolog_orthology_type | <biomart.Attribute name='gagassizii_homolog_orthology_type', display_name="Agassiz's desert tortoise homology type", description=''>
gagassizii_homolog_perc_id | <biomart.Attribute name='gagassizii_homolog_perc_id', display_name="%id. target Agassiz's desert tortoise gene identical to query gene", description=''>
gagassizii_homolog_perc_id_r1 | <biomart.Attribute name='gagassizii_homolog_perc_id_r1', display_name="%id. query gene identical to target Agassiz's desert tortoise gene", description=''>
gagassizii_homolog_goc_score | <biomart.Attribute name='gagassizii_homolog_goc_score', display_name="Agassiz's desert tortoise Gene-order conservation score", description=''>
gagassizii_homolog_wga_coverage | <biomart.Attribute name='gagassizii_homolog_wga_coverage', display_name="Agassiz's desert tortoise Whole-genome alignment coverage", description=''>
gagassizii_homolog_orthology_confidence | <biomart.Attribute name='gagassizii_homolog_orthology_confidence', display_name="Agassiz's desert tortoise orthology confidence [0 low, 1 high]", description=''>
mspretus_homolog_ensembl_gene | <biomart.Attribute name='mspretus_homolog_ensembl_gene', display_name='Algerian mouse gene stable ID', description=''>
mspretus_homolog_associated_gene_name | <biomart.Attribute name='mspretus_homolog_associated_gene_name', display_name='Algerian mouse gene name', description=''>
mspretus_homolog_ensembl_peptide | <biomart.Attribute name='mspretus_homolog_ensembl_peptide', display_name='Algerian mouse protein or transcript stable ID', description=''>
mspretus_homolog_chromosome | <biomart.Attribute name='mspretus_homolog_chromosome', display_name='Algerian mouse chromosome/scaffold name', description=''>
mspretus_homolog_chrom_start | <biomart.Attribute name='mspretus_homolog_chrom_start', display_name='Algerian mouse chromosome/scaffold start (bp)', description=''>
mspretus_homolog_chrom_end | <biomart.Attribute name='mspretus_homolog_chrom_end', display_name='Algerian mouse chromosome/scaffold end (bp)', description=''>
mspretus_homolog_canonical_transcript_protein | <biomart.Attribute name='mspretus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mspretus_homolog_subtype | <biomart.Attribute name='mspretus_homolog_subtype', display_name='Last common ancestor with Algerian mouse', description=''>
mspretus_homolog_orthology_type | <biomart.Attribute name='mspretus_homolog_orthology_type', display_name='Algerian mouse homology type', description=''>
mspretus_homolog_perc_id | <biomart.Attribute name='mspretus_homolog_perc_id', display_name='%id. target Algerian mouse gene identical to query gene', description=''>
mspretus_homolog_perc_id_r1 | <biomart.Attribute name='mspretus_homolog_perc_id_r1', display_name='%id. query gene identical to target Algerian mouse gene', description=''>
mspretus_homolog_goc_score | <biomart.Attribute name='mspretus_homolog_goc_score', display_name='Algerian mouse Gene-order conservation score', description=''>
mspretus_homolog_wga_coverage | <biomart.Attribute name='mspretus_homolog_wga_coverage', display_name='Algerian mouse Whole-genome alignment coverage', description=''>
mspretus_homolog_orthology_confidence | <biomart.Attribute name='mspretus_homolog_orthology_confidence', display_name='Algerian mouse orthology confidence [0 low, 1 high]', description=''>
vpacos_homolog_ensembl_gene | <biomart.Attribute name='vpacos_homolog_ensembl_gene', display_name='Alpaca gene stable ID', description=''>
vpacos_homolog_associated_gene_name | <biomart.Attribute name='vpacos_homolog_associated_gene_name', display_name='Alpaca gene name', description=''>
vpacos_homolog_ensembl_peptide | <biomart.Attribute name='vpacos_homolog_ensembl_peptide', display_name='Alpaca protein or transcript stable ID', description=''>
vpacos_homolog_chromosome | <biomart.Attribute name='vpacos_homolog_chromosome', display_name='Alpaca chromosome/scaffold name', description=''>
vpacos_homolog_chrom_start | <biomart.Attribute name='vpacos_homolog_chrom_start', display_name='Alpaca chromosome/scaffold start (bp)', description=''>
vpacos_homolog_chrom_end | <biomart.Attribute name='vpacos_homolog_chrom_end', display_name='Alpaca chromosome/scaffold end (bp)', description=''>
vpacos_homolog_canonical_transcript_protein | <biomart.Attribute name='vpacos_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
vpacos_homolog_subtype | <biomart.Attribute name='vpacos_homolog_subtype', display_name='Last common ancestor with Alpaca', description=''>
vpacos_homolog_orthology_type | <biomart.Attribute name='vpacos_homolog_orthology_type', display_name='Alpaca homology type', description=''>
vpacos_homolog_perc_id | <biomart.Attribute name='vpacos_homolog_perc_id', display_name='%id. target Alpaca gene identical to query gene', description=''>
vpacos_homolog_perc_id_r1 | <biomart.Attribute name='vpacos_homolog_perc_id_r1', display_name='%id. query gene identical to target Alpaca gene', description=''>
vpacos_homolog_goc_score | <biomart.Attribute name='vpacos_homolog_goc_score', display_name='Alpaca Gene-order conservation score', description=''>
vpacos_homolog_wga_coverage | <biomart.Attribute name='vpacos_homolog_wga_coverage', display_name='Alpaca Whole-genome alignment coverage', description=''>
vpacos_homolog_orthology_confidence | <biomart.Attribute name='vpacos_homolog_orthology_confidence', display_name='Alpaca orthology confidence [0 low, 1 high]', description=''>
mmmarmota_homolog_ensembl_gene | <biomart.Attribute name='mmmarmota_homolog_ensembl_gene', display_name='Alpine marmot gene stable ID', description=''>
mmmarmota_homolog_associated_gene_name | <biomart.Attribute name='mmmarmota_homolog_associated_gene_name', display_name='Alpine marmot gene name', description=''>
mmmarmota_homolog_ensembl_peptide | <biomart.Attribute name='mmmarmota_homolog_ensembl_peptide', display_name='Alpine marmot protein or transcript stable ID', description=''>
mmmarmota_homolog_chromosome | <biomart.Attribute name='mmmarmota_homolog_chromosome', display_name='Alpine marmot chromosome/scaffold name', description=''>
mmmarmota_homolog_chrom_start | <biomart.Attribute name='mmmarmota_homolog_chrom_start', display_name='Alpine marmot chromosome/scaffold start (bp)', description=''>
mmmarmota_homolog_chrom_end | <biomart.Attribute name='mmmarmota_homolog_chrom_end', display_name='Alpine marmot chromosome/scaffold end (bp)', description=''>
mmmarmota_homolog_canonical_transcript_protein | <biomart.Attribute name='mmmarmota_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mmmarmota_homolog_subtype | <biomart.Attribute name='mmmarmota_homolog_subtype', display_name='Last common ancestor with Alpine marmot', description=''>
mmmarmota_homolog_orthology_type | <biomart.Attribute name='mmmarmota_homolog_orthology_type', display_name='Alpine marmot homology type', description=''>
mmmarmota_homolog_perc_id | <biomart.Attribute name='mmmarmota_homolog_perc_id', display_name='%id. target Alpine marmot gene identical to query gene', description=''>
mmmarmota_homolog_perc_id_r1 | <biomart.Attribute name='mmmarmota_homolog_perc_id_r1', display_name='%id. query gene identical to target Alpine marmot gene', description=''>
mmmarmota_homolog_goc_score | <biomart.Attribute name='mmmarmota_homolog_goc_score', display_name='Alpine marmot Gene-order conservation score', description=''>
mmmarmota_homolog_wga_coverage | <biomart.Attribute name='mmmarmota_homolog_wga_coverage', display_name='Alpine marmot Whole-genome alignment coverage', description=''>
mmmarmota_homolog_orthology_confidence | <biomart.Attribute name='mmmarmota_homolog_orthology_confidence', display_name='Alpine marmot orthology confidence [0 low, 1 high]', description=''>
pformosa_homolog_ensembl_gene | <biomart.Attribute name='pformosa_homolog_ensembl_gene', display_name='Amazon molly gene stable ID', description=''>
pformosa_homolog_associated_gene_name | <biomart.Attribute name='pformosa_homolog_associated_gene_name', display_name='Amazon molly gene name', description=''>
pformosa_homolog_ensembl_peptide | <biomart.Attribute name='pformosa_homolog_ensembl_peptide', display_name='Amazon molly protein or transcript stable ID', description=''>
pformosa_homolog_chromosome | <biomart.Attribute name='pformosa_homolog_chromosome', display_name='Amazon molly chromosome/scaffold name', description=''>
pformosa_homolog_chrom_start | <biomart.Attribute name='pformosa_homolog_chrom_start', display_name='Amazon molly chromosome/scaffold start (bp)', description=''>
pformosa_homolog_chrom_end | <biomart.Attribute name='pformosa_homolog_chrom_end', display_name='Amazon molly chromosome/scaffold end (bp)', description=''>
pformosa_homolog_canonical_transcript_protein | <biomart.Attribute name='pformosa_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pformosa_homolog_subtype | <biomart.Attribute name='pformosa_homolog_subtype', display_name='Last common ancestor with Amazon molly', description=''>
pformosa_homolog_orthology_type | <biomart.Attribute name='pformosa_homolog_orthology_type', display_name='Amazon molly homology type', description=''>
pformosa_homolog_perc_id | <biomart.Attribute name='pformosa_homolog_perc_id', display_name='%id. target Amazon molly gene identical to query gene', description=''>
pformosa_homolog_perc_id_r1 | <biomart.Attribute name='pformosa_homolog_perc_id_r1', display_name='%id. query gene identical to target Amazon molly gene', description=''>
pformosa_homolog_goc_score | <biomart.Attribute name='pformosa_homolog_goc_score', display_name='Amazon molly Gene-order conservation score', description=''>
pformosa_homolog_wga_coverage | <biomart.Attribute name='pformosa_homolog_wga_coverage', display_name='Amazon molly Whole-genome alignment coverage', description=''>
pformosa_homolog_orthology_confidence | <biomart.Attribute name='pformosa_homolog_orthology_confidence', display_name='Amazon molly orthology confidence [0 low, 1 high]', description=''>
ccanadensis_homolog_ensembl_gene | <biomart.Attribute name='ccanadensis_homolog_ensembl_gene', display_name='American beaver gene stable ID', description=''>
ccanadensis_homolog_associated_gene_name | <biomart.Attribute name='ccanadensis_homolog_associated_gene_name', display_name='American beaver gene name', description=''>
ccanadensis_homolog_ensembl_peptide | <biomart.Attribute name='ccanadensis_homolog_ensembl_peptide', display_name='American beaver protein or transcript stable ID', description=''>
ccanadensis_homolog_chromosome | <biomart.Attribute name='ccanadensis_homolog_chromosome', display_name='American beaver chromosome/scaffold name', description=''>
ccanadensis_homolog_chrom_start | <biomart.Attribute name='ccanadensis_homolog_chrom_start', display_name='American beaver chromosome/scaffold start (bp)', description=''>
ccanadensis_homolog_chrom_end | <biomart.Attribute name='ccanadensis_homolog_chrom_end', display_name='American beaver chromosome/scaffold end (bp)', description=''>
ccanadensis_homolog_canonical_transcript_protein | <biomart.Attribute name='ccanadensis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ccanadensis_homolog_subtype | <biomart.Attribute name='ccanadensis_homolog_subtype', display_name='Last common ancestor with American beaver', description=''>
ccanadensis_homolog_orthology_type | <biomart.Attribute name='ccanadensis_homolog_orthology_type', display_name='American beaver homology type', description=''>
ccanadensis_homolog_perc_id | <biomart.Attribute name='ccanadensis_homolog_perc_id', display_name='%id. target American beaver gene identical to query gene', description=''>
ccanadensis_homolog_perc_id_r1 | <biomart.Attribute name='ccanadensis_homolog_perc_id_r1', display_name='%id. query gene identical to target American beaver gene', description=''>
ccanadensis_homolog_goc_score | <biomart.Attribute name='ccanadensis_homolog_goc_score', display_name='American beaver Gene-order conservation score', description=''>
ccanadensis_homolog_wga_coverage | <biomart.Attribute name='ccanadensis_homolog_wga_coverage', display_name='American beaver Whole-genome alignment coverage', description=''>
ccanadensis_homolog_orthology_confidence | <biomart.Attribute name='ccanadensis_homolog_orthology_confidence', display_name='American beaver orthology confidence [0 low, 1 high]', description=''>
bbbison_homolog_ensembl_gene | <biomart.Attribute name='bbbison_homolog_ensembl_gene', display_name='American bison gene stable ID', description=''>
bbbison_homolog_associated_gene_name | <biomart.Attribute name='bbbison_homolog_associated_gene_name', display_name='American bison gene name', description=''>
bbbison_homolog_ensembl_peptide | <biomart.Attribute name='bbbison_homolog_ensembl_peptide', display_name='American bison protein or transcript stable ID', description=''>
bbbison_homolog_chromosome | <biomart.Attribute name='bbbison_homolog_chromosome', display_name='American bison chromosome/scaffold name', description=''>
bbbison_homolog_chrom_start | <biomart.Attribute name='bbbison_homolog_chrom_start', display_name='American bison chromosome/scaffold start (bp)', description=''>
bbbison_homolog_chrom_end | <biomart.Attribute name='bbbison_homolog_chrom_end', display_name='American bison chromosome/scaffold end (bp)', description=''>
bbbison_homolog_canonical_transcript_protein | <biomart.Attribute name='bbbison_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
bbbison_homolog_subtype | <biomart.Attribute name='bbbison_homolog_subtype', display_name='Last common ancestor with American bison', description=''>
bbbison_homolog_orthology_type | <biomart.Attribute name='bbbison_homolog_orthology_type', display_name='American bison homology type', description=''>
bbbison_homolog_perc_id | <biomart.Attribute name='bbbison_homolog_perc_id', display_name='%id. target American bison gene identical to query gene', description=''>
bbbison_homolog_perc_id_r1 | <biomart.Attribute name='bbbison_homolog_perc_id_r1', display_name='%id. query gene identical to target American bison gene', description=''>
bbbison_homolog_goc_score | <biomart.Attribute name='bbbison_homolog_goc_score', display_name='American bison Gene-order conservation score', description=''>
bbbison_homolog_wga_coverage | <biomart.Attribute name='bbbison_homolog_wga_coverage', display_name='American bison Whole-genome alignment coverage', description=''>
bbbison_homolog_orthology_confidence | <biomart.Attribute name='bbbison_homolog_orthology_confidence', display_name='American bison orthology confidence [0 low, 1 high]', description=''>
uamericanus_homolog_ensembl_gene | <biomart.Attribute name='uamericanus_homolog_ensembl_gene', display_name='American black bear gene stable ID', description=''>
uamericanus_homolog_associated_gene_name | <biomart.Attribute name='uamericanus_homolog_associated_gene_name', display_name='American black bear gene name', description=''>
uamericanus_homolog_ensembl_peptide | <biomart.Attribute name='uamericanus_homolog_ensembl_peptide', display_name='American black bear protein or transcript stable ID', description=''>
uamericanus_homolog_chromosome | <biomart.Attribute name='uamericanus_homolog_chromosome', display_name='American black bear chromosome/scaffold name', description=''>
uamericanus_homolog_chrom_start | <biomart.Attribute name='uamericanus_homolog_chrom_start', display_name='American black bear chromosome/scaffold start (bp)', description=''>
uamericanus_homolog_chrom_end | <biomart.Attribute name='uamericanus_homolog_chrom_end', display_name='American black bear chromosome/scaffold end (bp)', description=''>
uamericanus_homolog_canonical_transcript_protein | <biomart.Attribute name='uamericanus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
uamericanus_homolog_subtype | <biomart.Attribute name='uamericanus_homolog_subtype', display_name='Last common ancestor with American black bear', description=''>
uamericanus_homolog_orthology_type | <biomart.Attribute name='uamericanus_homolog_orthology_type', display_name='American black bear homology type', description=''>
uamericanus_homolog_perc_id | <biomart.Attribute name='uamericanus_homolog_perc_id', display_name='%id. target American black bear gene identical to query gene', description=''>
uamericanus_homolog_perc_id_r1 | <biomart.Attribute name='uamericanus_homolog_perc_id_r1', display_name='%id. query gene identical to target American black bear gene', description=''>
uamericanus_homolog_goc_score | <biomart.Attribute name='uamericanus_homolog_goc_score', display_name='American black bear Gene-order conservation score', description=''>
uamericanus_homolog_wga_coverage | <biomart.Attribute name='uamericanus_homolog_wga_coverage', display_name='American black bear Whole-genome alignment coverage', description=''>
uamericanus_homolog_orthology_confidence | <biomart.Attribute name='uamericanus_homolog_orthology_confidence', display_name='American black bear orthology confidence [0 low, 1 high]', description=''>
nvison_homolog_ensembl_gene | <biomart.Attribute name='nvison_homolog_ensembl_gene', display_name='American mink gene stable ID', description=''>
nvison_homolog_associated_gene_name | <biomart.Attribute name='nvison_homolog_associated_gene_name', display_name='American mink gene name', description=''>
nvison_homolog_ensembl_peptide | <biomart.Attribute name='nvison_homolog_ensembl_peptide', display_name='American mink protein or transcript stable ID', description=''>
nvison_homolog_chromosome | <biomart.Attribute name='nvison_homolog_chromosome', display_name='American mink chromosome/scaffold name', description=''>
nvison_homolog_chrom_start | <biomart.Attribute name='nvison_homolog_chrom_start', display_name='American mink chromosome/scaffold start (bp)', description=''>
nvison_homolog_chrom_end | <biomart.Attribute name='nvison_homolog_chrom_end', display_name='American mink chromosome/scaffold end (bp)', description=''>
nvison_homolog_canonical_transcript_protein | <biomart.Attribute name='nvison_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
nvison_homolog_subtype | <biomart.Attribute name='nvison_homolog_subtype', display_name='Last common ancestor with American mink', description=''>
nvison_homolog_orthology_type | <biomart.Attribute name='nvison_homolog_orthology_type', display_name='American mink homology type', description=''>
nvison_homolog_perc_id | <biomart.Attribute name='nvison_homolog_perc_id', display_name='%id. target American mink gene identical to query gene', description=''>
nvison_homolog_perc_id_r1 | <biomart.Attribute name='nvison_homolog_perc_id_r1', display_name='%id. query gene identical to target American mink gene', description=''>
nvison_homolog_goc_score | <biomart.Attribute name='nvison_homolog_goc_score', display_name='American mink Gene-order conservation score', description=''>
nvison_homolog_wga_coverage | <biomart.Attribute name='nvison_homolog_wga_coverage', display_name='American mink Whole-genome alignment coverage', description=''>
nvison_homolog_orthology_confidence | <biomart.Attribute name='nvison_homolog_orthology_confidence', display_name='American mink orthology confidence [0 low, 1 high]', description=''>
capalliatus_homolog_ensembl_gene | <biomart.Attribute name='capalliatus_homolog_ensembl_gene', display_name='Angola colobus gene stable ID', description=''>
capalliatus_homolog_associated_gene_name | <biomart.Attribute name='capalliatus_homolog_associated_gene_name', display_name='Angola colobus gene name', description=''>
capalliatus_homolog_ensembl_peptide | <biomart.Attribute name='capalliatus_homolog_ensembl_peptide', display_name='Angola colobus protein or transcript stable ID', description=''>
capalliatus_homolog_chromosome | <biomart.Attribute name='capalliatus_homolog_chromosome', display_name='Angola colobus chromosome/scaffold name', description=''>
capalliatus_homolog_chrom_start | <biomart.Attribute name='capalliatus_homolog_chrom_start', display_name='Angola colobus chromosome/scaffold start (bp)', description=''>
capalliatus_homolog_chrom_end | <biomart.Attribute name='capalliatus_homolog_chrom_end', display_name='Angola colobus chromosome/scaffold end (bp)', description=''>
capalliatus_homolog_canonical_transcript_protein | <biomart.Attribute name='capalliatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
capalliatus_homolog_subtype | <biomart.Attribute name='capalliatus_homolog_subtype', display_name='Last common ancestor with Angola colobus', description=''>
capalliatus_homolog_orthology_type | <biomart.Attribute name='capalliatus_homolog_orthology_type', display_name='Angola colobus homology type', description=''>
capalliatus_homolog_perc_id | <biomart.Attribute name='capalliatus_homolog_perc_id', display_name='%id. target Angola colobus gene identical to query gene', description=''>
capalliatus_homolog_perc_id_r1 | <biomart.Attribute name='capalliatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Angola colobus gene', description=''>
capalliatus_homolog_goc_score | <biomart.Attribute name='capalliatus_homolog_goc_score', display_name='Angola colobus Gene-order conservation score', description=''>
capalliatus_homolog_wga_coverage | <biomart.Attribute name='capalliatus_homolog_wga_coverage', display_name='Angola colobus Whole-genome alignment coverage', description=''>
capalliatus_homolog_orthology_confidence | <biomart.Attribute name='capalliatus_homolog_orthology_confidence', display_name='Angola colobus orthology confidence [0 low, 1 high]', description=''>
acarolinensis_homolog_ensembl_gene | <biomart.Attribute name='acarolinensis_homolog_ensembl_gene', display_name='Anole lizard gene stable ID', description=''>
acarolinensis_homolog_associated_gene_name | <biomart.Attribute name='acarolinensis_homolog_associated_gene_name', display_name='Anole lizard gene name', description=''>
acarolinensis_homolog_ensembl_peptide | <biomart.Attribute name='acarolinensis_homolog_ensembl_peptide', display_name='Anole lizard protein or transcript stable ID', description=''>
acarolinensis_homolog_chromosome | <biomart.Attribute name='acarolinensis_homolog_chromosome', display_name='Anole lizard chromosome/scaffold name', description=''>
acarolinensis_homolog_chrom_start | <biomart.Attribute name='acarolinensis_homolog_chrom_start', display_name='Anole lizard chromosome/scaffold start (bp)', description=''>
acarolinensis_homolog_chrom_end | <biomart.Attribute name='acarolinensis_homolog_chrom_end', display_name='Anole lizard chromosome/scaffold end (bp)', description=''>
acarolinensis_homolog_canonical_transcript_protein | <biomart.Attribute name='acarolinensis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
acarolinensis_homolog_subtype | <biomart.Attribute name='acarolinensis_homolog_subtype', display_name='Last common ancestor with Anole lizard', description=''>
acarolinensis_homolog_orthology_type | <biomart.Attribute name='acarolinensis_homolog_orthology_type', display_name='Anole lizard homology type', description=''>
acarolinensis_homolog_perc_id | <biomart.Attribute name='acarolinensis_homolog_perc_id', display_name='%id. target Anole lizard gene identical to query gene', description=''>
acarolinensis_homolog_perc_id_r1 | <biomart.Attribute name='acarolinensis_homolog_perc_id_r1', display_name='%id. query gene identical to target Anole lizard gene', description=''>
acarolinensis_homolog_goc_score | <biomart.Attribute name='acarolinensis_homolog_goc_score', display_name='Anole lizard Gene-order conservation score', description=''>
acarolinensis_homolog_wga_coverage | <biomart.Attribute name='acarolinensis_homolog_wga_coverage', display_name='Anole lizard Whole-genome alignment coverage', description=''>
acarolinensis_homolog_orthology_confidence | <biomart.Attribute name='acarolinensis_homolog_orthology_confidence', display_name='Anole lizard orthology confidence [0 low, 1 high]', description=''>
cdromedarius_homolog_ensembl_gene | <biomart.Attribute name='cdromedarius_homolog_ensembl_gene', display_name='Arabian camel gene stable ID', description=''>
cdromedarius_homolog_associated_gene_name | <biomart.Attribute name='cdromedarius_homolog_associated_gene_name', display_name='Arabian camel gene name', description=''>
cdromedarius_homolog_ensembl_peptide | <biomart.Attribute name='cdromedarius_homolog_ensembl_peptide', display_name='Arabian camel protein or transcript stable ID', description=''>
cdromedarius_homolog_chromosome | <biomart.Attribute name='cdromedarius_homolog_chromosome', display_name='Arabian camel chromosome/scaffold name', description=''>
cdromedarius_homolog_chrom_start | <biomart.Attribute name='cdromedarius_homolog_chrom_start', display_name='Arabian camel chromosome/scaffold start (bp)', description=''>
cdromedarius_homolog_chrom_end | <biomart.Attribute name='cdromedarius_homolog_chrom_end', display_name='Arabian camel chromosome/scaffold end (bp)', description=''>
cdromedarius_homolog_canonical_transcript_protein | <biomart.Attribute name='cdromedarius_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cdromedarius_homolog_subtype | <biomart.Attribute name='cdromedarius_homolog_subtype', display_name='Last common ancestor with Arabian camel', description=''>
cdromedarius_homolog_orthology_type | <biomart.Attribute name='cdromedarius_homolog_orthology_type', display_name='Arabian camel homology type', description=''>
cdromedarius_homolog_perc_id | <biomart.Attribute name='cdromedarius_homolog_perc_id', display_name='%id. target Arabian camel gene identical to query gene', description=''>
cdromedarius_homolog_perc_id_r1 | <biomart.Attribute name='cdromedarius_homolog_perc_id_r1', display_name='%id. query gene identical to target Arabian camel gene', description=''>
cdromedarius_homolog_goc_score | <biomart.Attribute name='cdromedarius_homolog_goc_score', display_name='Arabian camel Gene-order conservation score', description=''>
cdromedarius_homolog_wga_coverage | <biomart.Attribute name='cdromedarius_homolog_wga_coverage', display_name='Arabian camel Whole-genome alignment coverage', description=''>
cdromedarius_homolog_orthology_confidence | <biomart.Attribute name='cdromedarius_homolog_orthology_confidence', display_name='Arabian camel orthology confidence [0 low, 1 high]', description=''>
smerianae_homolog_ensembl_gene | <biomart.Attribute name='smerianae_homolog_ensembl_gene', display_name='Argentine black and white tegu gene stable ID', description=''>
smerianae_homolog_associated_gene_name | <biomart.Attribute name='smerianae_homolog_associated_gene_name', display_name='Argentine black and white tegu gene name', description=''>
smerianae_homolog_ensembl_peptide | <biomart.Attribute name='smerianae_homolog_ensembl_peptide', display_name='Argentine black and white tegu protein or transcript stable ID', description=''>
smerianae_homolog_chromosome | <biomart.Attribute name='smerianae_homolog_chromosome', display_name='Argentine black and white tegu chromosome/scaffold name', description=''>
smerianae_homolog_chrom_start | <biomart.Attribute name='smerianae_homolog_chrom_start', display_name='Argentine black and white tegu chromosome/scaffold start (bp)', description=''>
smerianae_homolog_chrom_end | <biomart.Attribute name='smerianae_homolog_chrom_end', display_name='Argentine black and white tegu chromosome/scaffold end (bp)', description=''>
smerianae_homolog_canonical_transcript_protein | <biomart.Attribute name='smerianae_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
smerianae_homolog_subtype | <biomart.Attribute name='smerianae_homolog_subtype', display_name='Last common ancestor with Argentine black and white tegu', description=''>
smerianae_homolog_orthology_type | <biomart.Attribute name='smerianae_homolog_orthology_type', display_name='Argentine black and white tegu homology type', description=''>
smerianae_homolog_perc_id | <biomart.Attribute name='smerianae_homolog_perc_id', display_name='%id. target Argentine black and white tegu gene identical to query gene', description=''>
smerianae_homolog_perc_id_r1 | <biomart.Attribute name='smerianae_homolog_perc_id_r1', display_name='%id. query gene identical to target Argentine black and white tegu gene', description=''>
smerianae_homolog_goc_score | <biomart.Attribute name='smerianae_homolog_goc_score', display_name='Argentine black and white tegu Gene-order conservation score', description=''>
smerianae_homolog_wga_coverage | <biomart.Attribute name='smerianae_homolog_wga_coverage', display_name='Argentine black and white tegu Whole-genome alignment coverage', description=''>
smerianae_homolog_orthology_confidence | <biomart.Attribute name='smerianae_homolog_orthology_confidence', display_name='Argentine black and white tegu orthology confidence [0 low, 1 high]', description=''>
dnovemcinctus_homolog_ensembl_gene | <biomart.Attribute name='dnovemcinctus_homolog_ensembl_gene', display_name='Armadillo gene stable ID', description=''>
dnovemcinctus_homolog_associated_gene_name | <biomart.Attribute name='dnovemcinctus_homolog_associated_gene_name', display_name='Armadillo gene name', description=''>
dnovemcinctus_homolog_ensembl_peptide | <biomart.Attribute name='dnovemcinctus_homolog_ensembl_peptide', display_name='Armadillo protein or transcript stable ID', description=''>
dnovemcinctus_homolog_chromosome | <biomart.Attribute name='dnovemcinctus_homolog_chromosome', display_name='Armadillo chromosome/scaffold name', description=''>
dnovemcinctus_homolog_chrom_start | <biomart.Attribute name='dnovemcinctus_homolog_chrom_start', display_name='Armadillo chromosome/scaffold start (bp)', description=''>
dnovemcinctus_homolog_chrom_end | <biomart.Attribute name='dnovemcinctus_homolog_chrom_end', display_name='Armadillo chromosome/scaffold end (bp)', description=''>
dnovemcinctus_homolog_canonical_transcript_protein | <biomart.Attribute name='dnovemcinctus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
dnovemcinctus_homolog_subtype | <biomart.Attribute name='dnovemcinctus_homolog_subtype', display_name='Last common ancestor with Armadillo', description=''>
dnovemcinctus_homolog_orthology_type | <biomart.Attribute name='dnovemcinctus_homolog_orthology_type', display_name='Armadillo homology type', description=''>
dnovemcinctus_homolog_perc_id | <biomart.Attribute name='dnovemcinctus_homolog_perc_id', display_name='%id. target Armadillo gene identical to query gene', description=''>
dnovemcinctus_homolog_perc_id_r1 | <biomart.Attribute name='dnovemcinctus_homolog_perc_id_r1', display_name='%id. query gene identical to target Armadillo gene', description=''>
dnovemcinctus_homolog_goc_score | <biomart.Attribute name='dnovemcinctus_homolog_goc_score', display_name='Armadillo Gene-order conservation score', description=''>
dnovemcinctus_homolog_wga_coverage | <biomart.Attribute name='dnovemcinctus_homolog_wga_coverage', display_name='Armadillo Whole-genome alignment coverage', description=''>
dnovemcinctus_homolog_orthology_confidence | <biomart.Attribute name='dnovemcinctus_homolog_orthology_confidence', display_name='Armadillo orthology confidence [0 low, 1 high]', description=''>
sformosus_homolog_ensembl_gene | <biomart.Attribute name='sformosus_homolog_ensembl_gene', display_name='Asian bonytongue gene stable ID', description=''>
sformosus_homolog_associated_gene_name | <biomart.Attribute name='sformosus_homolog_associated_gene_name', display_name='Asian bonytongue gene name', description=''>
sformosus_homolog_ensembl_peptide | <biomart.Attribute name='sformosus_homolog_ensembl_peptide', display_name='Asian bonytongue protein or transcript stable ID', description=''>
sformosus_homolog_chromosome | <biomart.Attribute name='sformosus_homolog_chromosome', display_name='Asian bonytongue chromosome/scaffold name', description=''>
sformosus_homolog_chrom_start | <biomart.Attribute name='sformosus_homolog_chrom_start', display_name='Asian bonytongue chromosome/scaffold start (bp)', description=''>
sformosus_homolog_chrom_end | <biomart.Attribute name='sformosus_homolog_chrom_end', display_name='Asian bonytongue chromosome/scaffold end (bp)', description=''>
sformosus_homolog_canonical_transcript_protein | <biomart.Attribute name='sformosus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sformosus_homolog_subtype | <biomart.Attribute name='sformosus_homolog_subtype', display_name='Last common ancestor with Asian bonytongue', description=''>
sformosus_homolog_orthology_type | <biomart.Attribute name='sformosus_homolog_orthology_type', display_name='Asian bonytongue homology type', description=''>
sformosus_homolog_perc_id | <biomart.Attribute name='sformosus_homolog_perc_id', display_name='%id. target Asian bonytongue gene identical to query gene', description=''>
sformosus_homolog_perc_id_r1 | <biomart.Attribute name='sformosus_homolog_perc_id_r1', display_name='%id. query gene identical to target Asian bonytongue gene', description=''>
sformosus_homolog_goc_score | <biomart.Attribute name='sformosus_homolog_goc_score', display_name='Asian bonytongue Gene-order conservation score', description=''>
sformosus_homolog_wga_coverage | <biomart.Attribute name='sformosus_homolog_wga_coverage', display_name='Asian bonytongue Whole-genome alignment coverage', description=''>
sformosus_homolog_orthology_confidence | <biomart.Attribute name='sformosus_homolog_orthology_confidence', display_name='Asian bonytongue orthology confidence [0 low, 1 high]', description=''>
charengus_homolog_ensembl_gene | <biomart.Attribute name='charengus_homolog_ensembl_gene', display_name='Atlantic herring gene stable ID', description=''>
charengus_homolog_associated_gene_name | <biomart.Attribute name='charengus_homolog_associated_gene_name', display_name='Atlantic herring gene name', description=''>
charengus_homolog_ensembl_peptide | <biomart.Attribute name='charengus_homolog_ensembl_peptide', display_name='Atlantic herring protein or transcript stable ID', description=''>
charengus_homolog_chromosome | <biomart.Attribute name='charengus_homolog_chromosome', display_name='Atlantic herring chromosome/scaffold name', description=''>
charengus_homolog_chrom_start | <biomart.Attribute name='charengus_homolog_chrom_start', display_name='Atlantic herring chromosome/scaffold start (bp)', description=''>
charengus_homolog_chrom_end | <biomart.Attribute name='charengus_homolog_chrom_end', display_name='Atlantic herring chromosome/scaffold end (bp)', description=''>
charengus_homolog_canonical_transcript_protein | <biomart.Attribute name='charengus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
charengus_homolog_subtype | <biomart.Attribute name='charengus_homolog_subtype', display_name='Last common ancestor with Atlantic herring', description=''>
charengus_homolog_orthology_type | <biomart.Attribute name='charengus_homolog_orthology_type', display_name='Atlantic herring homology type', description=''>
charengus_homolog_perc_id | <biomart.Attribute name='charengus_homolog_perc_id', display_name='%id. target Atlantic herring gene identical to query gene', description=''>
charengus_homolog_perc_id_r1 | <biomart.Attribute name='charengus_homolog_perc_id_r1', display_name='%id. query gene identical to target Atlantic herring gene', description=''>
charengus_homolog_goc_score | <biomart.Attribute name='charengus_homolog_goc_score', display_name='Atlantic herring Gene-order conservation score', description=''>
charengus_homolog_wga_coverage | <biomart.Attribute name='charengus_homolog_wga_coverage', display_name='Atlantic herring Whole-genome alignment coverage', description=''>
charengus_homolog_orthology_confidence | <biomart.Attribute name='charengus_homolog_orthology_confidence', display_name='Atlantic herring orthology confidence [0 low, 1 high]', description=''>
ssalar_homolog_ensembl_gene | <biomart.Attribute name='ssalar_homolog_ensembl_gene', display_name='Atlantic salmon gene stable ID', description=''>
ssalar_homolog_associated_gene_name | <biomart.Attribute name='ssalar_homolog_associated_gene_name', display_name='Atlantic salmon gene name', description=''>
ssalar_homolog_ensembl_peptide | <biomart.Attribute name='ssalar_homolog_ensembl_peptide', display_name='Atlantic salmon protein or transcript stable ID', description=''>
ssalar_homolog_chromosome | <biomart.Attribute name='ssalar_homolog_chromosome', display_name='Atlantic salmon chromosome/scaffold name', description=''>
ssalar_homolog_chrom_start | <biomart.Attribute name='ssalar_homolog_chrom_start', display_name='Atlantic salmon chromosome/scaffold start (bp)', description=''>
ssalar_homolog_chrom_end | <biomart.Attribute name='ssalar_homolog_chrom_end', display_name='Atlantic salmon chromosome/scaffold end (bp)', description=''>
ssalar_homolog_canonical_transcript_protein | <biomart.Attribute name='ssalar_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ssalar_homolog_subtype | <biomart.Attribute name='ssalar_homolog_subtype', display_name='Last common ancestor with Atlantic salmon', description=''>
ssalar_homolog_orthology_type | <biomart.Attribute name='ssalar_homolog_orthology_type', display_name='Atlantic salmon homology type', description=''>
ssalar_homolog_perc_id | <biomart.Attribute name='ssalar_homolog_perc_id', display_name='%id. target Atlantic salmon gene identical to query gene', description=''>
ssalar_homolog_perc_id_r1 | <biomart.Attribute name='ssalar_homolog_perc_id_r1', display_name='%id. query gene identical to target Atlantic salmon gene', description=''>
ssalar_homolog_goc_score | <biomart.Attribute name='ssalar_homolog_goc_score', display_name='Atlantic salmon Gene-order conservation score', description=''>
ssalar_homolog_wga_coverage | <biomart.Attribute name='ssalar_homolog_wga_coverage', display_name='Atlantic salmon Whole-genome alignment coverage', description=''>
ssalar_homolog_orthology_confidence | <biomart.Attribute name='ssalar_homolog_orthology_confidence', display_name='Atlantic salmon orthology confidence [0 low, 1 high]', description=''>
cporosus_homolog_ensembl_gene | <biomart.Attribute name='cporosus_homolog_ensembl_gene', display_name='Australian saltwater crocodile gene stable ID', description=''>
cporosus_homolog_associated_gene_name | <biomart.Attribute name='cporosus_homolog_associated_gene_name', display_name='Australian saltwater crocodile gene name', description=''>
cporosus_homolog_ensembl_peptide | <biomart.Attribute name='cporosus_homolog_ensembl_peptide', display_name='Australian saltwater crocodile protein or transcript stable ID', description=''>
cporosus_homolog_chromosome | <biomart.Attribute name='cporosus_homolog_chromosome', display_name='Australian saltwater crocodile chromosome/scaffold name', description=''>
cporosus_homolog_chrom_start | <biomart.Attribute name='cporosus_homolog_chrom_start', display_name='Australian saltwater crocodile chromosome/scaffold start (bp)', description=''>
cporosus_homolog_chrom_end | <biomart.Attribute name='cporosus_homolog_chrom_end', display_name='Australian saltwater crocodile chromosome/scaffold end (bp)', description=''>
cporosus_homolog_canonical_transcript_protein | <biomart.Attribute name='cporosus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cporosus_homolog_subtype | <biomart.Attribute name='cporosus_homolog_subtype', display_name='Last common ancestor with Australian saltwater crocodile', description=''>
cporosus_homolog_orthology_type | <biomart.Attribute name='cporosus_homolog_orthology_type', display_name='Australian saltwater crocodile homology type', description=''>
cporosus_homolog_perc_id | <biomart.Attribute name='cporosus_homolog_perc_id', display_name='%id. target Australian saltwater crocodile gene identical to query gene', description=''>
cporosus_homolog_perc_id_r1 | <biomart.Attribute name='cporosus_homolog_perc_id_r1', display_name='%id. query gene identical to target Australian saltwater crocodile gene', description=''>
cporosus_homolog_goc_score | <biomart.Attribute name='cporosus_homolog_goc_score', display_name='Australian saltwater crocodile Gene-order conservation score', description=''>
cporosus_homolog_wga_coverage | <biomart.Attribute name='cporosus_homolog_wga_coverage', display_name='Australian saltwater crocodile Whole-genome alignment coverage', description=''>
cporosus_homolog_orthology_confidence | <biomart.Attribute name='cporosus_homolog_orthology_confidence', display_name='Australian saltwater crocodile orthology confidence [0 low, 1 high]', description=''>
lbergylta_homolog_ensembl_gene | <biomart.Attribute name='lbergylta_homolog_ensembl_gene', display_name='Ballan wrasse gene stable ID', description=''>
lbergylta_homolog_associated_gene_name | <biomart.Attribute name='lbergylta_homolog_associated_gene_name', display_name='Ballan wrasse gene name', description=''>
lbergylta_homolog_ensembl_peptide | <biomart.Attribute name='lbergylta_homolog_ensembl_peptide', display_name='Ballan wrasse protein or transcript stable ID', description=''>
lbergylta_homolog_chromosome | <biomart.Attribute name='lbergylta_homolog_chromosome', display_name='Ballan wrasse chromosome/scaffold name', description=''>
lbergylta_homolog_chrom_start | <biomart.Attribute name='lbergylta_homolog_chrom_start', display_name='Ballan wrasse chromosome/scaffold start (bp)', description=''>
lbergylta_homolog_chrom_end | <biomart.Attribute name='lbergylta_homolog_chrom_end', display_name='Ballan wrasse chromosome/scaffold end (bp)', description=''>
lbergylta_homolog_canonical_transcript_protein | <biomart.Attribute name='lbergylta_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lbergylta_homolog_subtype | <biomart.Attribute name='lbergylta_homolog_subtype', display_name='Last common ancestor with Ballan wrasse', description=''>
lbergylta_homolog_orthology_type | <biomart.Attribute name='lbergylta_homolog_orthology_type', display_name='Ballan wrasse homology type', description=''>
lbergylta_homolog_perc_id | <biomart.Attribute name='lbergylta_homolog_perc_id', display_name='%id. target Ballan wrasse gene identical to query gene', description=''>
lbergylta_homolog_perc_id_r1 | <biomart.Attribute name='lbergylta_homolog_perc_id_r1', display_name='%id. query gene identical to target Ballan wrasse gene', description=''>
lbergylta_homolog_goc_score | <biomart.Attribute name='lbergylta_homolog_goc_score', display_name='Ballan wrasse Gene-order conservation score', description=''>
lbergylta_homolog_wga_coverage | <biomart.Attribute name='lbergylta_homolog_wga_coverage', display_name='Ballan wrasse Whole-genome alignment coverage', description=''>
lbergylta_homolog_orthology_confidence | <biomart.Attribute name='lbergylta_homolog_orthology_confidence', display_name='Ballan wrasse orthology confidence [0 low, 1 high]', description=''>
lcalcarifer_homolog_ensembl_gene | <biomart.Attribute name='lcalcarifer_homolog_ensembl_gene', display_name='Barramundi perch gene stable ID', description=''>
lcalcarifer_homolog_associated_gene_name | <biomart.Attribute name='lcalcarifer_homolog_associated_gene_name', display_name='Barramundi perch gene name', description=''>
lcalcarifer_homolog_ensembl_peptide | <biomart.Attribute name='lcalcarifer_homolog_ensembl_peptide', display_name='Barramundi perch protein or transcript stable ID', description=''>
lcalcarifer_homolog_chromosome | <biomart.Attribute name='lcalcarifer_homolog_chromosome', display_name='Barramundi perch chromosome/scaffold name', description=''>
lcalcarifer_homolog_chrom_start | <biomart.Attribute name='lcalcarifer_homolog_chrom_start', display_name='Barramundi perch chromosome/scaffold start (bp)', description=''>
lcalcarifer_homolog_chrom_end | <biomart.Attribute name='lcalcarifer_homolog_chrom_end', display_name='Barramundi perch chromosome/scaffold end (bp)', description=''>
lcalcarifer_homolog_canonical_transcript_protein | <biomart.Attribute name='lcalcarifer_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lcalcarifer_homolog_subtype | <biomart.Attribute name='lcalcarifer_homolog_subtype', display_name='Last common ancestor with Barramundi perch', description=''>
lcalcarifer_homolog_orthology_type | <biomart.Attribute name='lcalcarifer_homolog_orthology_type', display_name='Barramundi perch homology type', description=''>
lcalcarifer_homolog_perc_id | <biomart.Attribute name='lcalcarifer_homolog_perc_id', display_name='%id. target Barramundi perch gene identical to query gene', description=''>
lcalcarifer_homolog_perc_id_r1 | <biomart.Attribute name='lcalcarifer_homolog_perc_id_r1', display_name='%id. query gene identical to target Barramundi perch gene', description=''>
lcalcarifer_homolog_goc_score | <biomart.Attribute name='lcalcarifer_homolog_goc_score', display_name='Barramundi perch Gene-order conservation score', description=''>
lcalcarifer_homolog_wga_coverage | <biomart.Attribute name='lcalcarifer_homolog_wga_coverage', display_name='Barramundi perch Whole-genome alignment coverage', description=''>
lcalcarifer_homolog_orthology_confidence | <biomart.Attribute name='lcalcarifer_homolog_orthology_confidence', display_name='Barramundi perch orthology confidence [0 low, 1 high]', description=''>
lsdomestica_homolog_ensembl_gene | <biomart.Attribute name='lsdomestica_homolog_ensembl_gene', display_name='Bengalese finch gene stable ID', description=''>
lsdomestica_homolog_associated_gene_name | <biomart.Attribute name='lsdomestica_homolog_associated_gene_name', display_name='Bengalese finch gene name', description=''>
lsdomestica_homolog_ensembl_peptide | <biomart.Attribute name='lsdomestica_homolog_ensembl_peptide', display_name='Bengalese finch protein or transcript stable ID', description=''>
lsdomestica_homolog_chromosome | <biomart.Attribute name='lsdomestica_homolog_chromosome', display_name='Bengalese finch chromosome/scaffold name', description=''>
lsdomestica_homolog_chrom_start | <biomart.Attribute name='lsdomestica_homolog_chrom_start', display_name='Bengalese finch chromosome/scaffold start (bp)', description=''>
lsdomestica_homolog_chrom_end | <biomart.Attribute name='lsdomestica_homolog_chrom_end', display_name='Bengalese finch chromosome/scaffold end (bp)', description=''>
lsdomestica_homolog_canonical_transcript_protein | <biomart.Attribute name='lsdomestica_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lsdomestica_homolog_subtype | <biomart.Attribute name='lsdomestica_homolog_subtype', display_name='Last common ancestor with Bengalese finch', description=''>
lsdomestica_homolog_orthology_type | <biomart.Attribute name='lsdomestica_homolog_orthology_type', display_name='Bengalese finch homology type', description=''>
lsdomestica_homolog_perc_id | <biomart.Attribute name='lsdomestica_homolog_perc_id', display_name='%id. target Bengalese finch gene identical to query gene', description=''>
lsdomestica_homolog_perc_id_r1 | <biomart.Attribute name='lsdomestica_homolog_perc_id_r1', display_name='%id. query gene identical to target Bengalese finch gene', description=''>
lsdomestica_homolog_goc_score | <biomart.Attribute name='lsdomestica_homolog_goc_score', display_name='Bengalese finch Gene-order conservation score', description=''>
lsdomestica_homolog_wga_coverage | <biomart.Attribute name='lsdomestica_homolog_wga_coverage', display_name='Bengalese finch Whole-genome alignment coverage', description=''>
lsdomestica_homolog_orthology_confidence | <biomart.Attribute name='lsdomestica_homolog_orthology_confidence', display_name='Bengalese finch orthology confidence [0 low, 1 high]', description=''>
rbieti_homolog_ensembl_gene | <biomart.Attribute name='rbieti_homolog_ensembl_gene', display_name='Black snub-nosed monkey gene stable ID', description=''>
rbieti_homolog_associated_gene_name | <biomart.Attribute name='rbieti_homolog_associated_gene_name', display_name='Black snub-nosed monkey gene name', description=''>
rbieti_homolog_ensembl_peptide | <biomart.Attribute name='rbieti_homolog_ensembl_peptide', display_name='Black snub-nosed monkey protein or transcript stable ID', description=''>
rbieti_homolog_chromosome | <biomart.Attribute name='rbieti_homolog_chromosome', display_name='Black snub-nosed monkey chromosome/scaffold name', description=''>
rbieti_homolog_chrom_start | <biomart.Attribute name='rbieti_homolog_chrom_start', display_name='Black snub-nosed monkey chromosome/scaffold start (bp)', description=''>
rbieti_homolog_chrom_end | <biomart.Attribute name='rbieti_homolog_chrom_end', display_name='Black snub-nosed monkey chromosome/scaffold end (bp)', description=''>
rbieti_homolog_canonical_transcript_protein | <biomart.Attribute name='rbieti_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
rbieti_homolog_subtype | <biomart.Attribute name='rbieti_homolog_subtype', display_name='Last common ancestor with Black snub-nosed monkey', description=''>
rbieti_homolog_orthology_type | <biomart.Attribute name='rbieti_homolog_orthology_type', display_name='Black snub-nosed monkey homology type', description=''>
rbieti_homolog_perc_id | <biomart.Attribute name='rbieti_homolog_perc_id', display_name='%id. target Black snub-nosed monkey gene identical to query gene', description=''>
rbieti_homolog_perc_id_r1 | <biomart.Attribute name='rbieti_homolog_perc_id_r1', display_name='%id. query gene identical to target Black snub-nosed monkey gene', description=''>
rbieti_homolog_goc_score | <biomart.Attribute name='rbieti_homolog_goc_score', display_name='Black snub-nosed monkey Gene-order conservation score', description=''>
rbieti_homolog_wga_coverage | <biomart.Attribute name='rbieti_homolog_wga_coverage', display_name='Black snub-nosed monkey Whole-genome alignment coverage', description=''>
rbieti_homolog_orthology_confidence | <biomart.Attribute name='rbieti_homolog_orthology_confidence', display_name='Black snub-nosed monkey orthology confidence [0 low, 1 high]', description=''>
oaureus_homolog_ensembl_gene | <biomart.Attribute name='oaureus_homolog_ensembl_gene', display_name='Blue tilapia gene stable ID', description=''>
oaureus_homolog_associated_gene_name | <biomart.Attribute name='oaureus_homolog_associated_gene_name', display_name='Blue tilapia gene name', description=''>
oaureus_homolog_ensembl_peptide | <biomart.Attribute name='oaureus_homolog_ensembl_peptide', display_name='Blue tilapia protein or transcript stable ID', description=''>
oaureus_homolog_chromosome | <biomart.Attribute name='oaureus_homolog_chromosome', display_name='Blue tilapia chromosome/scaffold name', description=''>
oaureus_homolog_chrom_start | <biomart.Attribute name='oaureus_homolog_chrom_start', display_name='Blue tilapia chromosome/scaffold start (bp)', description=''>
oaureus_homolog_chrom_end | <biomart.Attribute name='oaureus_homolog_chrom_end', display_name='Blue tilapia chromosome/scaffold end (bp)', description=''>
oaureus_homolog_canonical_transcript_protein | <biomart.Attribute name='oaureus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
oaureus_homolog_subtype | <biomart.Attribute name='oaureus_homolog_subtype', display_name='Last common ancestor with Blue tilapia', description=''>
oaureus_homolog_orthology_type | <biomart.Attribute name='oaureus_homolog_orthology_type', display_name='Blue tilapia homology type', description=''>
oaureus_homolog_perc_id | <biomart.Attribute name='oaureus_homolog_perc_id', display_name='%id. target Blue tilapia gene identical to query gene', description=''>
oaureus_homolog_perc_id_r1 | <biomart.Attribute name='oaureus_homolog_perc_id_r1', display_name='%id. query gene identical to target Blue tilapia gene', description=''>
oaureus_homolog_goc_score | <biomart.Attribute name='oaureus_homolog_goc_score', display_name='Blue tilapia Gene-order conservation score', description=''>
oaureus_homolog_wga_coverage | <biomart.Attribute name='oaureus_homolog_wga_coverage', display_name='Blue tilapia Whole-genome alignment coverage', description=''>
oaureus_homolog_orthology_confidence | <biomart.Attribute name='oaureus_homolog_orthology_confidence', display_name='Blue tilapia orthology confidence [0 low, 1 high]', description=''>
lcoronata_homolog_ensembl_gene | <biomart.Attribute name='lcoronata_homolog_ensembl_gene', display_name='Blue-crowned manakin gene stable ID', description=''>
lcoronata_homolog_associated_gene_name | <biomart.Attribute name='lcoronata_homolog_associated_gene_name', display_name='Blue-crowned manakin gene name', description=''>
lcoronata_homolog_ensembl_peptide | <biomart.Attribute name='lcoronata_homolog_ensembl_peptide', display_name='Blue-crowned manakin protein or transcript stable ID', description=''>
lcoronata_homolog_chromosome | <biomart.Attribute name='lcoronata_homolog_chromosome', display_name='Blue-crowned manakin chromosome/scaffold name', description=''>
lcoronata_homolog_chrom_start | <biomart.Attribute name='lcoronata_homolog_chrom_start', display_name='Blue-crowned manakin chromosome/scaffold start (bp)', description=''>
lcoronata_homolog_chrom_end | <biomart.Attribute name='lcoronata_homolog_chrom_end', display_name='Blue-crowned manakin chromosome/scaffold end (bp)', description=''>
lcoronata_homolog_canonical_transcript_protein | <biomart.Attribute name='lcoronata_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lcoronata_homolog_subtype | <biomart.Attribute name='lcoronata_homolog_subtype', display_name='Last common ancestor with Blue-crowned manakin', description=''>
lcoronata_homolog_orthology_type | <biomart.Attribute name='lcoronata_homolog_orthology_type', display_name='Blue-crowned manakin homology type', description=''>
lcoronata_homolog_perc_id | <biomart.Attribute name='lcoronata_homolog_perc_id', display_name='%id. target Blue-crowned manakin gene identical to query gene', description=''>
lcoronata_homolog_perc_id_r1 | <biomart.Attribute name='lcoronata_homolog_perc_id_r1', display_name='%id. query gene identical to target Blue-crowned manakin gene', description=''>
lcoronata_homolog_goc_score | <biomart.Attribute name='lcoronata_homolog_goc_score', display_name='Blue-crowned manakin Gene-order conservation score', description=''>
lcoronata_homolog_wga_coverage | <biomart.Attribute name='lcoronata_homolog_wga_coverage', display_name='Blue-crowned manakin Whole-genome alignment coverage', description=''>
lcoronata_homolog_orthology_confidence | <biomart.Attribute name='lcoronata_homolog_orthology_confidence', display_name='Blue-crowned manakin orthology confidence [0 low, 1 high]', description=''>
gwilldenowi_homolog_ensembl_gene | <biomart.Attribute name='gwilldenowi_homolog_ensembl_gene', display_name='Blunt-snouted clingfish gene stable ID', description=''>
gwilldenowi_homolog_associated_gene_name | <biomart.Attribute name='gwilldenowi_homolog_associated_gene_name', display_name='Blunt-snouted clingfish gene name', description=''>
gwilldenowi_homolog_ensembl_peptide | <biomart.Attribute name='gwilldenowi_homolog_ensembl_peptide', display_name='Blunt-snouted clingfish protein or transcript stable ID', description=''>
gwilldenowi_homolog_chromosome | <biomart.Attribute name='gwilldenowi_homolog_chromosome', display_name='Blunt-snouted clingfish chromosome/scaffold name', description=''>
gwilldenowi_homolog_chrom_start | <biomart.Attribute name='gwilldenowi_homolog_chrom_start', display_name='Blunt-snouted clingfish chromosome/scaffold start (bp)', description=''>
gwilldenowi_homolog_chrom_end | <biomart.Attribute name='gwilldenowi_homolog_chrom_end', display_name='Blunt-snouted clingfish chromosome/scaffold end (bp)', description=''>
gwilldenowi_homolog_canonical_transcript_protein | <biomart.Attribute name='gwilldenowi_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
gwilldenowi_homolog_subtype | <biomart.Attribute name='gwilldenowi_homolog_subtype', display_name='Last common ancestor with Blunt-snouted clingfish', description=''>
gwilldenowi_homolog_orthology_type | <biomart.Attribute name='gwilldenowi_homolog_orthology_type', display_name='Blunt-snouted clingfish homology type', description=''>
gwilldenowi_homolog_perc_id | <biomart.Attribute name='gwilldenowi_homolog_perc_id', display_name='%id. target Blunt-snouted clingfish gene identical to query gene', description=''>
gwilldenowi_homolog_perc_id_r1 | <biomart.Attribute name='gwilldenowi_homolog_perc_id_r1', display_name='%id. query gene identical to target Blunt-snouted clingfish gene', description=''>
gwilldenowi_homolog_goc_score | <biomart.Attribute name='gwilldenowi_homolog_goc_score', display_name='Blunt-snouted clingfish Gene-order conservation score', description=''>
gwilldenowi_homolog_wga_coverage | <biomart.Attribute name='gwilldenowi_homolog_wga_coverage', display_name='Blunt-snouted clingfish Whole-genome alignment coverage', description=''>
gwilldenowi_homolog_orthology_confidence | <biomart.Attribute name='gwilldenowi_homolog_orthology_confidence', display_name='Blunt-snouted clingfish orthology confidence [0 low, 1 high]', description=''>
sbboliviensis_homolog_ensembl_gene | <biomart.Attribute name='sbboliviensis_homolog_ensembl_gene', display_name='Bolivian squirrel monkey gene stable ID', description=''>
sbboliviensis_homolog_associated_gene_name | <biomart.Attribute name='sbboliviensis_homolog_associated_gene_name', display_name='Bolivian squirrel monkey gene name', description=''>
sbboliviensis_homolog_ensembl_peptide | <biomart.Attribute name='sbboliviensis_homolog_ensembl_peptide', display_name='Bolivian squirrel monkey protein or transcript stable ID', description=''>
sbboliviensis_homolog_chromosome | <biomart.Attribute name='sbboliviensis_homolog_chromosome', display_name='Bolivian squirrel monkey chromosome/scaffold name', description=''>
sbboliviensis_homolog_chrom_start | <biomart.Attribute name='sbboliviensis_homolog_chrom_start', display_name='Bolivian squirrel monkey chromosome/scaffold start (bp)', description=''>
sbboliviensis_homolog_chrom_end | <biomart.Attribute name='sbboliviensis_homolog_chrom_end', display_name='Bolivian squirrel monkey chromosome/scaffold end (bp)', description=''>
sbboliviensis_homolog_canonical_transcript_protein | <biomart.Attribute name='sbboliviensis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sbboliviensis_homolog_subtype | <biomart.Attribute name='sbboliviensis_homolog_subtype', display_name='Last common ancestor with Bolivian squirrel monkey', description=''>
sbboliviensis_homolog_orthology_type | <biomart.Attribute name='sbboliviensis_homolog_orthology_type', display_name='Bolivian squirrel monkey homology type', description=''>
sbboliviensis_homolog_perc_id | <biomart.Attribute name='sbboliviensis_homolog_perc_id', display_name='%id. target Bolivian squirrel monkey gene identical to query gene', description=''>
sbboliviensis_homolog_perc_id_r1 | <biomart.Attribute name='sbboliviensis_homolog_perc_id_r1', display_name='%id. query gene identical to target Bolivian squirrel monkey gene', description=''>
sbboliviensis_homolog_goc_score | <biomart.Attribute name='sbboliviensis_homolog_goc_score', display_name='Bolivian squirrel monkey Gene-order conservation score', description=''>
sbboliviensis_homolog_wga_coverage | <biomart.Attribute name='sbboliviensis_homolog_wga_coverage', display_name='Bolivian squirrel monkey Whole-genome alignment coverage', description=''>
sbboliviensis_homolog_orthology_confidence | <biomart.Attribute name='sbboliviensis_homolog_orthology_confidence', display_name='Bolivian squirrel monkey orthology confidence [0 low, 1 high]', description=''>
ppaniscus_homolog_ensembl_gene | <biomart.Attribute name='ppaniscus_homolog_ensembl_gene', display_name='Bonobo gene stable ID', description=''>
ppaniscus_homolog_associated_gene_name | <biomart.Attribute name='ppaniscus_homolog_associated_gene_name', display_name='Bonobo gene name', description=''>
ppaniscus_homolog_ensembl_peptide | <biomart.Attribute name='ppaniscus_homolog_ensembl_peptide', display_name='Bonobo protein or transcript stable ID', description=''>
ppaniscus_homolog_chromosome | <biomart.Attribute name='ppaniscus_homolog_chromosome', display_name='Bonobo chromosome/scaffold name', description=''>
ppaniscus_homolog_chrom_start | <biomart.Attribute name='ppaniscus_homolog_chrom_start', display_name='Bonobo chromosome/scaffold start (bp)', description=''>
ppaniscus_homolog_chrom_end | <biomart.Attribute name='ppaniscus_homolog_chrom_end', display_name='Bonobo chromosome/scaffold end (bp)', description=''>
ppaniscus_homolog_canonical_transcript_protein | <biomart.Attribute name='ppaniscus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ppaniscus_homolog_subtype | <biomart.Attribute name='ppaniscus_homolog_subtype', display_name='Last common ancestor with Bonobo', description=''>
ppaniscus_homolog_orthology_type | <biomart.Attribute name='ppaniscus_homolog_orthology_type', display_name='Bonobo homology type', description=''>
ppaniscus_homolog_perc_id | <biomart.Attribute name='ppaniscus_homolog_perc_id', display_name='%id. target Bonobo gene identical to query gene', description=''>
ppaniscus_homolog_perc_id_r1 | <biomart.Attribute name='ppaniscus_homolog_perc_id_r1', display_name='%id. query gene identical to target Bonobo gene', description=''>
ppaniscus_homolog_goc_score | <biomart.Attribute name='ppaniscus_homolog_goc_score', display_name='Bonobo Gene-order conservation score', description=''>
ppaniscus_homolog_wga_coverage | <biomart.Attribute name='ppaniscus_homolog_wga_coverage', display_name='Bonobo Whole-genome alignment coverage', description=''>
ppaniscus_homolog_orthology_confidence | <biomart.Attribute name='ppaniscus_homolog_orthology_confidence', display_name='Bonobo orthology confidence [0 low, 1 high]', description=''>
caperea_homolog_ensembl_gene | <biomart.Attribute name='caperea_homolog_ensembl_gene', display_name='Brazilian guinea pig gene stable ID', description=''>
caperea_homolog_associated_gene_name | <biomart.Attribute name='caperea_homolog_associated_gene_name', display_name='Brazilian guinea pig gene name', description=''>
caperea_homolog_ensembl_peptide | <biomart.Attribute name='caperea_homolog_ensembl_peptide', display_name='Brazilian guinea pig protein or transcript stable ID', description=''>
caperea_homolog_chromosome | <biomart.Attribute name='caperea_homolog_chromosome', display_name='Brazilian guinea pig chromosome/scaffold name', description=''>
caperea_homolog_chrom_start | <biomart.Attribute name='caperea_homolog_chrom_start', display_name='Brazilian guinea pig chromosome/scaffold start (bp)', description=''>
caperea_homolog_chrom_end | <biomart.Attribute name='caperea_homolog_chrom_end', display_name='Brazilian guinea pig chromosome/scaffold end (bp)', description=''>
caperea_homolog_canonical_transcript_protein | <biomart.Attribute name='caperea_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
caperea_homolog_subtype | <biomart.Attribute name='caperea_homolog_subtype', display_name='Last common ancestor with Brazilian guinea pig', description=''>
caperea_homolog_orthology_type | <biomart.Attribute name='caperea_homolog_orthology_type', display_name='Brazilian guinea pig homology type', description=''>
caperea_homolog_perc_id | <biomart.Attribute name='caperea_homolog_perc_id', display_name='%id. target Brazilian guinea pig gene identical to query gene', description=''>
caperea_homolog_perc_id_r1 | <biomart.Attribute name='caperea_homolog_perc_id_r1', display_name='%id. query gene identical to target Brazilian guinea pig gene', description=''>
caperea_homolog_goc_score | <biomart.Attribute name='caperea_homolog_goc_score', display_name='Brazilian guinea pig Gene-order conservation score', description=''>
caperea_homolog_wga_coverage | <biomart.Attribute name='caperea_homolog_wga_coverage', display_name='Brazilian guinea pig Whole-genome alignment coverage', description=''>
caperea_homolog_orthology_confidence | <biomart.Attribute name='caperea_homolog_orthology_confidence', display_name='Brazilian guinea pig orthology confidence [0 low, 1 high]', description=''>
strutta_homolog_ensembl_gene | <biomart.Attribute name='strutta_homolog_ensembl_gene', display_name='Brown trout gene stable ID', description=''>
strutta_homolog_associated_gene_name | <biomart.Attribute name='strutta_homolog_associated_gene_name', display_name='Brown trout gene name', description=''>
strutta_homolog_ensembl_peptide | <biomart.Attribute name='strutta_homolog_ensembl_peptide', display_name='Brown trout protein or transcript stable ID', description=''>
strutta_homolog_chromosome | <biomart.Attribute name='strutta_homolog_chromosome', display_name='Brown trout chromosome/scaffold name', description=''>
strutta_homolog_chrom_start | <biomart.Attribute name='strutta_homolog_chrom_start', display_name='Brown trout chromosome/scaffold start (bp)', description=''>
strutta_homolog_chrom_end | <biomart.Attribute name='strutta_homolog_chrom_end', display_name='Brown trout chromosome/scaffold end (bp)', description=''>
strutta_homolog_canonical_transcript_protein | <biomart.Attribute name='strutta_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
strutta_homolog_subtype | <biomart.Attribute name='strutta_homolog_subtype', display_name='Last common ancestor with Brown trout', description=''>
strutta_homolog_orthology_type | <biomart.Attribute name='strutta_homolog_orthology_type', display_name='Brown trout homology type', description=''>
strutta_homolog_perc_id | <biomart.Attribute name='strutta_homolog_perc_id', display_name='%id. target Brown trout gene identical to query gene', description=''>
strutta_homolog_perc_id_r1 | <biomart.Attribute name='strutta_homolog_perc_id_r1', display_name='%id. query gene identical to target Brown trout gene', description=''>
strutta_homolog_goc_score | <biomart.Attribute name='strutta_homolog_goc_score', display_name='Brown trout Gene-order conservation score', description=''>
strutta_homolog_wga_coverage | <biomart.Attribute name='strutta_homolog_wga_coverage', display_name='Brown trout Whole-genome alignment coverage', description=''>
strutta_homolog_orthology_confidence | <biomart.Attribute name='strutta_homolog_orthology_confidence', display_name='Brown trout orthology confidence [0 low, 1 high]', description=''>
mundulatus_homolog_ensembl_gene | <biomart.Attribute name='mundulatus_homolog_ensembl_gene', display_name='Budgerigar gene stable ID', description=''>
mundulatus_homolog_associated_gene_name | <biomart.Attribute name='mundulatus_homolog_associated_gene_name', display_name='Budgerigar gene name', description=''>
mundulatus_homolog_ensembl_peptide | <biomart.Attribute name='mundulatus_homolog_ensembl_peptide', display_name='Budgerigar protein or transcript stable ID', description=''>
mundulatus_homolog_chromosome | <biomart.Attribute name='mundulatus_homolog_chromosome', display_name='Budgerigar chromosome/scaffold name', description=''>
mundulatus_homolog_chrom_start | <biomart.Attribute name='mundulatus_homolog_chrom_start', display_name='Budgerigar chromosome/scaffold start (bp)', description=''>
mundulatus_homolog_chrom_end | <biomart.Attribute name='mundulatus_homolog_chrom_end', display_name='Budgerigar chromosome/scaffold end (bp)', description=''>
mundulatus_homolog_canonical_transcript_protein | <biomart.Attribute name='mundulatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mundulatus_homolog_subtype | <biomart.Attribute name='mundulatus_homolog_subtype', display_name='Last common ancestor with Budgerigar', description=''>
mundulatus_homolog_orthology_type | <biomart.Attribute name='mundulatus_homolog_orthology_type', display_name='Budgerigar homology type', description=''>
mundulatus_homolog_perc_id | <biomart.Attribute name='mundulatus_homolog_perc_id', display_name='%id. target Budgerigar gene identical to query gene', description=''>
mundulatus_homolog_perc_id_r1 | <biomart.Attribute name='mundulatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Budgerigar gene', description=''>
mundulatus_homolog_goc_score | <biomart.Attribute name='mundulatus_homolog_goc_score', display_name='Budgerigar Gene-order conservation score', description=''>
mundulatus_homolog_wga_coverage | <biomart.Attribute name='mundulatus_homolog_wga_coverage', display_name='Budgerigar Whole-genome alignment coverage', description=''>
mundulatus_homolog_orthology_confidence | <biomart.Attribute name='mundulatus_homolog_orthology_confidence', display_name='Budgerigar orthology confidence [0 low, 1 high]', description=''>
hburtoni_homolog_ensembl_gene | <biomart.Attribute name='hburtoni_homolog_ensembl_gene', display_name="Burton's mouthbrooder gene stable ID", description=''>
hburtoni_homolog_associated_gene_name | <biomart.Attribute name='hburtoni_homolog_associated_gene_name', display_name="Burton's mouthbrooder gene name", description=''>
hburtoni_homolog_ensembl_peptide | <biomart.Attribute name='hburtoni_homolog_ensembl_peptide', display_name="Burton's mouthbrooder protein or transcript stable ID", description=''>
hburtoni_homolog_chromosome | <biomart.Attribute name='hburtoni_homolog_chromosome', display_name="Burton's mouthbrooder chromosome/scaffold name", description=''>
hburtoni_homolog_chrom_start | <biomart.Attribute name='hburtoni_homolog_chrom_start', display_name="Burton's mouthbrooder chromosome/scaffold start (bp)", description=''>
hburtoni_homolog_chrom_end | <biomart.Attribute name='hburtoni_homolog_chrom_end', display_name="Burton's mouthbrooder chromosome/scaffold end (bp)", description=''>
hburtoni_homolog_canonical_transcript_protein | <biomart.Attribute name='hburtoni_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
hburtoni_homolog_subtype | <biomart.Attribute name='hburtoni_homolog_subtype', display_name="Last common ancestor with Burton's mouthbrooder", description=''>
hburtoni_homolog_orthology_type | <biomart.Attribute name='hburtoni_homolog_orthology_type', display_name="Burton's mouthbrooder homology type", description=''>
hburtoni_homolog_perc_id | <biomart.Attribute name='hburtoni_homolog_perc_id', display_name="%id. target Burton's mouthbrooder gene identical to query gene", description=''>
hburtoni_homolog_perc_id_r1 | <biomart.Attribute name='hburtoni_homolog_perc_id_r1', display_name="%id. query gene identical to target Burton's mouthbrooder gene", description=''>
hburtoni_homolog_goc_score | <biomart.Attribute name='hburtoni_homolog_goc_score', display_name="Burton's mouthbrooder Gene-order conservation score", description=''>
hburtoni_homolog_wga_coverage | <biomart.Attribute name='hburtoni_homolog_wga_coverage', display_name="Burton's mouthbrooder Whole-genome alignment coverage", description=''>
hburtoni_homolog_orthology_confidence | <biomart.Attribute name='hburtoni_homolog_orthology_confidence', display_name="Burton's mouthbrooder orthology confidence [0 low, 1 high]", description=''>
ogarnettii_homolog_ensembl_gene | <biomart.Attribute name='ogarnettii_homolog_ensembl_gene', display_name='Bushbaby gene stable ID', description=''>
ogarnettii_homolog_associated_gene_name | <biomart.Attribute name='ogarnettii_homolog_associated_gene_name', display_name='Bushbaby gene name', description=''>
ogarnettii_homolog_ensembl_peptide | <biomart.Attribute name='ogarnettii_homolog_ensembl_peptide', display_name='Bushbaby protein or transcript stable ID', description=''>
ogarnettii_homolog_chromosome | <biomart.Attribute name='ogarnettii_homolog_chromosome', display_name='Bushbaby chromosome/scaffold name', description=''>
ogarnettii_homolog_chrom_start | <biomart.Attribute name='ogarnettii_homolog_chrom_start', display_name='Bushbaby chromosome/scaffold start (bp)', description=''>
ogarnettii_homolog_chrom_end | <biomart.Attribute name='ogarnettii_homolog_chrom_end', display_name='Bushbaby chromosome/scaffold end (bp)', description=''>
ogarnettii_homolog_canonical_transcript_protein | <biomart.Attribute name='ogarnettii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ogarnettii_homolog_subtype | <biomart.Attribute name='ogarnettii_homolog_subtype', display_name='Last common ancestor with Bushbaby', description=''>
ogarnettii_homolog_orthology_type | <biomart.Attribute name='ogarnettii_homolog_orthology_type', display_name='Bushbaby homology type', description=''>
ogarnettii_homolog_perc_id | <biomart.Attribute name='ogarnettii_homolog_perc_id', display_name='%id. target Bushbaby gene identical to query gene', description=''>
ogarnettii_homolog_perc_id_r1 | <biomart.Attribute name='ogarnettii_homolog_perc_id_r1', display_name='%id. query gene identical to target Bushbaby gene', description=''>
ogarnettii_homolog_goc_score | <biomart.Attribute name='ogarnettii_homolog_goc_score', display_name='Bushbaby Gene-order conservation score', description=''>
ogarnettii_homolog_wga_coverage | <biomart.Attribute name='ogarnettii_homolog_wga_coverage', display_name='Bushbaby Whole-genome alignment coverage', description=''>
ogarnettii_homolog_orthology_confidence | <biomart.Attribute name='ogarnettii_homolog_orthology_confidence', display_name='Bushbaby orthology confidence [0 low, 1 high]', description=''>
cintestinalis_homolog_ensembl_gene | <biomart.Attribute name='cintestinalis_homolog_ensembl_gene', display_name='C.intestinalis gene stable ID', description=''>
cintestinalis_homolog_associated_gene_name | <biomart.Attribute name='cintestinalis_homolog_associated_gene_name', display_name='C.intestinalis gene name', description=''>
cintestinalis_homolog_ensembl_peptide | <biomart.Attribute name='cintestinalis_homolog_ensembl_peptide', display_name='C.intestinalis protein or transcript stable ID', description=''>
cintestinalis_homolog_chromosome | <biomart.Attribute name='cintestinalis_homolog_chromosome', display_name='C.intestinalis chromosome/scaffold name', description=''>
cintestinalis_homolog_chrom_start | <biomart.Attribute name='cintestinalis_homolog_chrom_start', display_name='C.intestinalis chromosome/scaffold start (bp)', description=''>
cintestinalis_homolog_chrom_end | <biomart.Attribute name='cintestinalis_homolog_chrom_end', display_name='C.intestinalis chromosome/scaffold end (bp)', description=''>
cintestinalis_homolog_canonical_transcript_protein | <biomart.Attribute name='cintestinalis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cintestinalis_homolog_subtype | <biomart.Attribute name='cintestinalis_homolog_subtype', display_name='Last common ancestor with C.intestinalis', description=''>
cintestinalis_homolog_orthology_type | <biomart.Attribute name='cintestinalis_homolog_orthology_type', display_name='C.intestinalis homology type', description=''>
cintestinalis_homolog_perc_id | <biomart.Attribute name='cintestinalis_homolog_perc_id', display_name='%id. target C.intestinalis gene identical to query gene', description=''>
cintestinalis_homolog_perc_id_r1 | <biomart.Attribute name='cintestinalis_homolog_perc_id_r1', display_name='%id. query gene identical to target C.intestinalis gene', description=''>
cintestinalis_homolog_wga_coverage | <biomart.Attribute name='cintestinalis_homolog_wga_coverage', display_name='C.intestinalis Whole-genome alignment coverage', description=''>
cintestinalis_homolog_orthology_confidence | <biomart.Attribute name='cintestinalis_homolog_orthology_confidence', display_name='C.intestinalis orthology confidence [0 low, 1 high]', description=''>
csavignyi_homolog_ensembl_gene | <biomart.Attribute name='csavignyi_homolog_ensembl_gene', display_name='C.savignyi gene stable ID', description=''>
csavignyi_homolog_associated_gene_name | <biomart.Attribute name='csavignyi_homolog_associated_gene_name', display_name='C.savignyi gene name', description=''>
csavignyi_homolog_ensembl_peptide | <biomart.Attribute name='csavignyi_homolog_ensembl_peptide', display_name='C.savignyi protein or transcript stable ID', description=''>
csavignyi_homolog_chromosome | <biomart.Attribute name='csavignyi_homolog_chromosome', display_name='C.savignyi chromosome/scaffold name', description=''>
csavignyi_homolog_chrom_start | <biomart.Attribute name='csavignyi_homolog_chrom_start', display_name='C.savignyi chromosome/scaffold start (bp)', description=''>
csavignyi_homolog_chrom_end | <biomart.Attribute name='csavignyi_homolog_chrom_end', display_name='C.savignyi chromosome/scaffold end (bp)', description=''>
csavignyi_homolog_canonical_transcript_protein | <biomart.Attribute name='csavignyi_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
csavignyi_homolog_subtype | <biomart.Attribute name='csavignyi_homolog_subtype', display_name='Last common ancestor with C.savignyi', description=''>
csavignyi_homolog_orthology_type | <biomart.Attribute name='csavignyi_homolog_orthology_type', display_name='C.savignyi homology type', description=''>
csavignyi_homolog_perc_id | <biomart.Attribute name='csavignyi_homolog_perc_id', display_name='%id. target C.savignyi gene identical to query gene', description=''>
csavignyi_homolog_perc_id_r1 | <biomart.Attribute name='csavignyi_homolog_perc_id_r1', display_name='%id. query gene identical to target C.savignyi gene', description=''>
csavignyi_homolog_wga_coverage | <biomart.Attribute name='csavignyi_homolog_wga_coverage', display_name='C.savignyi Whole-genome alignment coverage', description=''>
csavignyi_homolog_orthology_confidence | <biomart.Attribute name='csavignyi_homolog_orthology_confidence', display_name='C.savignyi orthology confidence [0 low, 1 high]', description=''>
celegans_homolog_ensembl_gene | <biomart.Attribute name='celegans_homolog_ensembl_gene', display_name='Caenorhabditis elegans gene stable ID', description=''>
celegans_homolog_associated_gene_name | <biomart.Attribute name='celegans_homolog_associated_gene_name', display_name='Caenorhabditis elegans gene name', description=''>
celegans_homolog_ensembl_peptide | <biomart.Attribute name='celegans_homolog_ensembl_peptide', display_name='Caenorhabditis elegans protein or transcript stable ID', description=''>
celegans_homolog_chromosome | <biomart.Attribute name='celegans_homolog_chromosome', display_name='Caenorhabditis elegans chromosome/scaffold name', description=''>
celegans_homolog_chrom_start | <biomart.Attribute name='celegans_homolog_chrom_start', display_name='Caenorhabditis elegans chromosome/scaffold start (bp)', description=''>
celegans_homolog_chrom_end | <biomart.Attribute name='celegans_homolog_chrom_end', display_name='Caenorhabditis elegans chromosome/scaffold end (bp)', description=''>
celegans_homolog_canonical_transcript_protein | <biomart.Attribute name='celegans_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
celegans_homolog_subtype | <biomart.Attribute name='celegans_homolog_subtype', display_name='Last common ancestor with Caenorhabditis elegans', description=''>
celegans_homolog_orthology_type | <biomart.Attribute name='celegans_homolog_orthology_type', display_name='Caenorhabditis elegans homology type', description=''>
celegans_homolog_perc_id | <biomart.Attribute name='celegans_homolog_perc_id', display_name='%id. target Caenorhabditis elegans gene identical to query gene', description=''>
celegans_homolog_perc_id_r1 | <biomart.Attribute name='celegans_homolog_perc_id_r1', display_name='%id. query gene identical to target Caenorhabditis elegans gene', description=''>
celegans_homolog_orthology_confidence | <biomart.Attribute name='celegans_homolog_orthology_confidence', display_name='Caenorhabditis elegans orthology confidence [0 low, 1 high]', description=''>
lcanadensis_homolog_ensembl_gene | <biomart.Attribute name='lcanadensis_homolog_ensembl_gene', display_name='Canada lynx gene stable ID', description=''>
lcanadensis_homolog_associated_gene_name | <biomart.Attribute name='lcanadensis_homolog_associated_gene_name', display_name='Canada lynx gene name', description=''>
lcanadensis_homolog_ensembl_peptide | <biomart.Attribute name='lcanadensis_homolog_ensembl_peptide', display_name='Canada lynx protein or transcript stable ID', description=''>
lcanadensis_homolog_chromosome | <biomart.Attribute name='lcanadensis_homolog_chromosome', display_name='Canada lynx chromosome/scaffold name', description=''>
lcanadensis_homolog_chrom_start | <biomart.Attribute name='lcanadensis_homolog_chrom_start', display_name='Canada lynx chromosome/scaffold start (bp)', description=''>
lcanadensis_homolog_chrom_end | <biomart.Attribute name='lcanadensis_homolog_chrom_end', display_name='Canada lynx chromosome/scaffold end (bp)', description=''>
lcanadensis_homolog_canonical_transcript_protein | <biomart.Attribute name='lcanadensis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lcanadensis_homolog_subtype | <biomart.Attribute name='lcanadensis_homolog_subtype', display_name='Last common ancestor with Canada lynx', description=''>
lcanadensis_homolog_orthology_type | <biomart.Attribute name='lcanadensis_homolog_orthology_type', display_name='Canada lynx homology type', description=''>
lcanadensis_homolog_perc_id | <biomart.Attribute name='lcanadensis_homolog_perc_id', display_name='%id. target Canada lynx gene identical to query gene', description=''>
lcanadensis_homolog_perc_id_r1 | <biomart.Attribute name='lcanadensis_homolog_perc_id_r1', display_name='%id. query gene identical to target Canada lynx gene', description=''>
lcanadensis_homolog_goc_score | <biomart.Attribute name='lcanadensis_homolog_goc_score', display_name='Canada lynx Gene-order conservation score', description=''>
lcanadensis_homolog_wga_coverage | <biomart.Attribute name='lcanadensis_homolog_wga_coverage', display_name='Canada lynx Whole-genome alignment coverage', description=''>
lcanadensis_homolog_orthology_confidence | <biomart.Attribute name='lcanadensis_homolog_orthology_confidence', display_name='Canada lynx orthology confidence [0 low, 1 high]', description=''>
ccapucinus_homolog_ensembl_gene | <biomart.Attribute name='ccapucinus_homolog_ensembl_gene', display_name='Capuchin gene stable ID', description=''>
ccapucinus_homolog_associated_gene_name | <biomart.Attribute name='ccapucinus_homolog_associated_gene_name', display_name='Capuchin gene name', description=''>
ccapucinus_homolog_ensembl_peptide | <biomart.Attribute name='ccapucinus_homolog_ensembl_peptide', display_name='Capuchin protein or transcript stable ID', description=''>
ccapucinus_homolog_chromosome | <biomart.Attribute name='ccapucinus_homolog_chromosome', display_name='Capuchin chromosome/scaffold name', description=''>
ccapucinus_homolog_chrom_start | <biomart.Attribute name='ccapucinus_homolog_chrom_start', display_name='Capuchin chromosome/scaffold start (bp)', description=''>
ccapucinus_homolog_chrom_end | <biomart.Attribute name='ccapucinus_homolog_chrom_end', display_name='Capuchin chromosome/scaffold end (bp)', description=''>
ccapucinus_homolog_canonical_transcript_protein | <biomart.Attribute name='ccapucinus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ccapucinus_homolog_subtype | <biomart.Attribute name='ccapucinus_homolog_subtype', display_name='Last common ancestor with Capuchin', description=''>
ccapucinus_homolog_orthology_type | <biomart.Attribute name='ccapucinus_homolog_orthology_type', display_name='Capuchin homology type', description=''>
ccapucinus_homolog_perc_id | <biomart.Attribute name='ccapucinus_homolog_perc_id', display_name='%id. target Capuchin gene identical to query gene', description=''>
ccapucinus_homolog_perc_id_r1 | <biomart.Attribute name='ccapucinus_homolog_perc_id_r1', display_name='%id. query gene identical to target Capuchin gene', description=''>
ccapucinus_homolog_goc_score | <biomart.Attribute name='ccapucinus_homolog_goc_score', display_name='Capuchin Gene-order conservation score', description=''>
ccapucinus_homolog_wga_coverage | <biomart.Attribute name='ccapucinus_homolog_wga_coverage', display_name='Capuchin Whole-genome alignment coverage', description=''>
ccapucinus_homolog_orthology_confidence | <biomart.Attribute name='ccapucinus_homolog_orthology_confidence', display_name='Capuchin orthology confidence [0 low, 1 high]', description=''>
fcatus_homolog_ensembl_gene | <biomart.Attribute name='fcatus_homolog_ensembl_gene', display_name='Cat gene stable ID', description=''>
fcatus_homolog_associated_gene_name | <biomart.Attribute name='fcatus_homolog_associated_gene_name', display_name='Cat gene name', description=''>
fcatus_homolog_ensembl_peptide | <biomart.Attribute name='fcatus_homolog_ensembl_peptide', display_name='Cat protein or transcript stable ID', description=''>
fcatus_homolog_chromosome | <biomart.Attribute name='fcatus_homolog_chromosome', display_name='Cat chromosome/scaffold name', description=''>
fcatus_homolog_chrom_start | <biomart.Attribute name='fcatus_homolog_chrom_start', display_name='Cat chromosome/scaffold start (bp)', description=''>
fcatus_homolog_chrom_end | <biomart.Attribute name='fcatus_homolog_chrom_end', display_name='Cat chromosome/scaffold end (bp)', description=''>
fcatus_homolog_canonical_transcript_protein | <biomart.Attribute name='fcatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
fcatus_homolog_subtype | <biomart.Attribute name='fcatus_homolog_subtype', display_name='Last common ancestor with Cat', description=''>
fcatus_homolog_orthology_type | <biomart.Attribute name='fcatus_homolog_orthology_type', display_name='Cat homology type', description=''>
fcatus_homolog_perc_id | <biomart.Attribute name='fcatus_homolog_perc_id', display_name='%id. target Cat gene identical to query gene', description=''>
fcatus_homolog_perc_id_r1 | <biomart.Attribute name='fcatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Cat gene', description=''>
fcatus_homolog_goc_score | <biomart.Attribute name='fcatus_homolog_goc_score', display_name='Cat Gene-order conservation score', description=''>
fcatus_homolog_wga_coverage | <biomart.Attribute name='fcatus_homolog_wga_coverage', display_name='Cat Whole-genome alignment coverage', description=''>
fcatus_homolog_orthology_confidence | <biomart.Attribute name='fcatus_homolog_orthology_confidence', display_name='Cat orthology confidence [0 low, 1 high]', description=''>
pvitticeps_homolog_ensembl_gene | <biomart.Attribute name='pvitticeps_homolog_ensembl_gene', display_name='Central bearded dragon gene stable ID', description=''>
pvitticeps_homolog_associated_gene_name | <biomart.Attribute name='pvitticeps_homolog_associated_gene_name', display_name='Central bearded dragon gene name', description=''>
pvitticeps_homolog_ensembl_peptide | <biomart.Attribute name='pvitticeps_homolog_ensembl_peptide', display_name='Central bearded dragon protein or transcript stable ID', description=''>
pvitticeps_homolog_chromosome | <biomart.Attribute name='pvitticeps_homolog_chromosome', display_name='Central bearded dragon chromosome/scaffold name', description=''>
pvitticeps_homolog_chrom_start | <biomart.Attribute name='pvitticeps_homolog_chrom_start', display_name='Central bearded dragon chromosome/scaffold start (bp)', description=''>
pvitticeps_homolog_chrom_end | <biomart.Attribute name='pvitticeps_homolog_chrom_end', display_name='Central bearded dragon chromosome/scaffold end (bp)', description=''>
pvitticeps_homolog_canonical_transcript_protein | <biomart.Attribute name='pvitticeps_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pvitticeps_homolog_subtype | <biomart.Attribute name='pvitticeps_homolog_subtype', display_name='Last common ancestor with Central bearded dragon', description=''>
pvitticeps_homolog_orthology_type | <biomart.Attribute name='pvitticeps_homolog_orthology_type', display_name='Central bearded dragon homology type', description=''>
pvitticeps_homolog_perc_id | <biomart.Attribute name='pvitticeps_homolog_perc_id', display_name='%id. target Central bearded dragon gene identical to query gene', description=''>
pvitticeps_homolog_perc_id_r1 | <biomart.Attribute name='pvitticeps_homolog_perc_id_r1', display_name='%id. query gene identical to target Central bearded dragon gene', description=''>
pvitticeps_homolog_goc_score | <biomart.Attribute name='pvitticeps_homolog_goc_score', display_name='Central bearded dragon Gene-order conservation score', description=''>
pvitticeps_homolog_wga_coverage | <biomart.Attribute name='pvitticeps_homolog_wga_coverage', display_name='Central bearded dragon Whole-genome alignment coverage', description=''>
pvitticeps_homolog_orthology_confidence | <biomart.Attribute name='pvitticeps_homolog_orthology_confidence', display_name='Central bearded dragon orthology confidence [0 low, 1 high]', description=''>
cgobio_homolog_ensembl_gene | <biomart.Attribute name='cgobio_homolog_ensembl_gene', display_name='Channel bull blenny gene stable ID', description=''>
cgobio_homolog_associated_gene_name | <biomart.Attribute name='cgobio_homolog_associated_gene_name', display_name='Channel bull blenny gene name', description=''>
cgobio_homolog_ensembl_peptide | <biomart.Attribute name='cgobio_homolog_ensembl_peptide', display_name='Channel bull blenny protein or transcript stable ID', description=''>
cgobio_homolog_chromosome | <biomart.Attribute name='cgobio_homolog_chromosome', display_name='Channel bull blenny chromosome/scaffold name', description=''>
cgobio_homolog_chrom_start | <biomart.Attribute name='cgobio_homolog_chrom_start', display_name='Channel bull blenny chromosome/scaffold start (bp)', description=''>
cgobio_homolog_chrom_end | <biomart.Attribute name='cgobio_homolog_chrom_end', display_name='Channel bull blenny chromosome/scaffold end (bp)', description=''>
cgobio_homolog_canonical_transcript_protein | <biomart.Attribute name='cgobio_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cgobio_homolog_subtype | <biomart.Attribute name='cgobio_homolog_subtype', display_name='Last common ancestor with Channel bull blenny', description=''>
cgobio_homolog_orthology_type | <biomart.Attribute name='cgobio_homolog_orthology_type', display_name='Channel bull blenny homology type', description=''>
cgobio_homolog_perc_id | <biomart.Attribute name='cgobio_homolog_perc_id', display_name='%id. target Channel bull blenny gene identical to query gene', description=''>
cgobio_homolog_perc_id_r1 | <biomart.Attribute name='cgobio_homolog_perc_id_r1', display_name='%id. query gene identical to target Channel bull blenny gene', description=''>
cgobio_homolog_goc_score | <biomart.Attribute name='cgobio_homolog_goc_score', display_name='Channel bull blenny Gene-order conservation score', description=''>
cgobio_homolog_wga_coverage | <biomart.Attribute name='cgobio_homolog_wga_coverage', display_name='Channel bull blenny Whole-genome alignment coverage', description=''>
cgobio_homolog_orthology_confidence | <biomart.Attribute name='cgobio_homolog_orthology_confidence', display_name='Channel bull blenny orthology confidence [0 low, 1 high]', description=''>
ipunctatus_homolog_ensembl_gene | <biomart.Attribute name='ipunctatus_homolog_ensembl_gene', display_name='Channel catfish gene stable ID', description=''>
ipunctatus_homolog_associated_gene_name | <biomart.Attribute name='ipunctatus_homolog_associated_gene_name', display_name='Channel catfish gene name', description=''>
ipunctatus_homolog_ensembl_peptide | <biomart.Attribute name='ipunctatus_homolog_ensembl_peptide', display_name='Channel catfish protein or transcript stable ID', description=''>
ipunctatus_homolog_chromosome | <biomart.Attribute name='ipunctatus_homolog_chromosome', display_name='Channel catfish chromosome/scaffold name', description=''>
ipunctatus_homolog_chrom_start | <biomart.Attribute name='ipunctatus_homolog_chrom_start', display_name='Channel catfish chromosome/scaffold start (bp)', description=''>
ipunctatus_homolog_chrom_end | <biomart.Attribute name='ipunctatus_homolog_chrom_end', display_name='Channel catfish chromosome/scaffold end (bp)', description=''>
ipunctatus_homolog_canonical_transcript_protein | <biomart.Attribute name='ipunctatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ipunctatus_homolog_subtype | <biomart.Attribute name='ipunctatus_homolog_subtype', display_name='Last common ancestor with Channel catfish', description=''>
ipunctatus_homolog_orthology_type | <biomart.Attribute name='ipunctatus_homolog_orthology_type', display_name='Channel catfish homology type', description=''>
ipunctatus_homolog_perc_id | <biomart.Attribute name='ipunctatus_homolog_perc_id', display_name='%id. target Channel catfish gene identical to query gene', description=''>
ipunctatus_homolog_perc_id_r1 | <biomart.Attribute name='ipunctatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Channel catfish gene', description=''>
ipunctatus_homolog_goc_score | <biomart.Attribute name='ipunctatus_homolog_goc_score', display_name='Channel catfish Gene-order conservation score', description=''>
ipunctatus_homolog_wga_coverage | <biomart.Attribute name='ipunctatus_homolog_wga_coverage', display_name='Channel catfish Whole-genome alignment coverage', description=''>
ipunctatus_homolog_orthology_confidence | <biomart.Attribute name='ipunctatus_homolog_orthology_confidence', display_name='Channel catfish orthology confidence [0 low, 1 high]', description=''>
ggallus_homolog_ensembl_gene | <biomart.Attribute name='ggallus_homolog_ensembl_gene', display_name='Chicken gene stable ID', description=''>
ggallus_homolog_associated_gene_name | <biomart.Attribute name='ggallus_homolog_associated_gene_name', display_name='Chicken gene name', description=''>
ggallus_homolog_ensembl_peptide | <biomart.Attribute name='ggallus_homolog_ensembl_peptide', display_name='Chicken protein or transcript stable ID', description=''>
ggallus_homolog_chromosome | <biomart.Attribute name='ggallus_homolog_chromosome', display_name='Chicken chromosome/scaffold name', description=''>
ggallus_homolog_chrom_start | <biomart.Attribute name='ggallus_homolog_chrom_start', display_name='Chicken chromosome/scaffold start (bp)', description=''>
ggallus_homolog_chrom_end | <biomart.Attribute name='ggallus_homolog_chrom_end', display_name='Chicken chromosome/scaffold end (bp)', description=''>
ggallus_homolog_canonical_transcript_protein | <biomart.Attribute name='ggallus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ggallus_homolog_subtype | <biomart.Attribute name='ggallus_homolog_subtype', display_name='Last common ancestor with Chicken', description=''>
ggallus_homolog_orthology_type | <biomart.Attribute name='ggallus_homolog_orthology_type', display_name='Chicken homology type', description=''>
ggallus_homolog_perc_id | <biomart.Attribute name='ggallus_homolog_perc_id', display_name='%id. target Chicken gene identical to query gene', description=''>
ggallus_homolog_perc_id_r1 | <biomart.Attribute name='ggallus_homolog_perc_id_r1', display_name='%id. query gene identical to target Chicken gene', description=''>
ggallus_homolog_goc_score | <biomart.Attribute name='ggallus_homolog_goc_score', display_name='Chicken Gene-order conservation score', description=''>
ggallus_homolog_wga_coverage | <biomart.Attribute name='ggallus_homolog_wga_coverage', display_name='Chicken Whole-genome alignment coverage', description=''>
ggallus_homolog_orthology_confidence | <biomart.Attribute name='ggallus_homolog_orthology_confidence', display_name='Chicken orthology confidence [0 low, 1 high]', description=''>
ptroglodytes_homolog_ensembl_gene | <biomart.Attribute name='ptroglodytes_homolog_ensembl_gene', display_name='Chimpanzee gene stable ID', description=''>
ptroglodytes_homolog_associated_gene_name | <biomart.Attribute name='ptroglodytes_homolog_associated_gene_name', display_name='Chimpanzee gene name', description=''>
ptroglodytes_homolog_ensembl_peptide | <biomart.Attribute name='ptroglodytes_homolog_ensembl_peptide', display_name='Chimpanzee protein or transcript stable ID', description=''>
ptroglodytes_homolog_chromosome | <biomart.Attribute name='ptroglodytes_homolog_chromosome', display_name='Chimpanzee chromosome/scaffold name', description=''>
ptroglodytes_homolog_chrom_start | <biomart.Attribute name='ptroglodytes_homolog_chrom_start', display_name='Chimpanzee chromosome/scaffold start (bp)', description=''>
ptroglodytes_homolog_chrom_end | <biomart.Attribute name='ptroglodytes_homolog_chrom_end', display_name='Chimpanzee chromosome/scaffold end (bp)', description=''>
ptroglodytes_homolog_canonical_transcript_protein | <biomart.Attribute name='ptroglodytes_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ptroglodytes_homolog_subtype | <biomart.Attribute name='ptroglodytes_homolog_subtype', display_name='Last common ancestor with Chimpanzee', description=''>
ptroglodytes_homolog_orthology_type | <biomart.Attribute name='ptroglodytes_homolog_orthology_type', display_name='Chimpanzee homology type', description=''>
ptroglodytes_homolog_perc_id | <biomart.Attribute name='ptroglodytes_homolog_perc_id', display_name='%id. target Chimpanzee gene identical to query gene', description=''>
ptroglodytes_homolog_perc_id_r1 | <biomart.Attribute name='ptroglodytes_homolog_perc_id_r1', display_name='%id. query gene identical to target Chimpanzee gene', description=''>
ptroglodytes_homolog_goc_score | <biomart.Attribute name='ptroglodytes_homolog_goc_score', display_name='Chimpanzee Gene-order conservation score', description=''>
ptroglodytes_homolog_wga_coverage | <biomart.Attribute name='ptroglodytes_homolog_wga_coverage', display_name='Chimpanzee Whole-genome alignment coverage', description=''>
ptroglodytes_homolog_orthology_confidence | <biomart.Attribute name='ptroglodytes_homolog_orthology_confidence', display_name='Chimpanzee orthology confidence [0 low, 1 high]', description=''>
cgchok1gshd_homolog_ensembl_gene | <biomart.Attribute name='cgchok1gshd_homolog_ensembl_gene', display_name='Chinese hamster CHOK1GS gene stable ID', description=''>
cgchok1gshd_homolog_associated_gene_name | <biomart.Attribute name='cgchok1gshd_homolog_associated_gene_name', display_name='Chinese hamster CHOK1GS gene name', description=''>
cgchok1gshd_homolog_ensembl_peptide | <biomart.Attribute name='cgchok1gshd_homolog_ensembl_peptide', display_name='Chinese hamster CHOK1GS protein or transcript stable ID', description=''>
cgchok1gshd_homolog_chromosome | <biomart.Attribute name='cgchok1gshd_homolog_chromosome', display_name='Chinese hamster CHOK1GS chromosome/scaffold name', description=''>
cgchok1gshd_homolog_chrom_start | <biomart.Attribute name='cgchok1gshd_homolog_chrom_start', display_name='Chinese hamster CHOK1GS chromosome/scaffold start (bp)', description=''>
cgchok1gshd_homolog_chrom_end | <biomart.Attribute name='cgchok1gshd_homolog_chrom_end', display_name='Chinese hamster CHOK1GS chromosome/scaffold end (bp)', description=''>
cgchok1gshd_homolog_canonical_transcript_protein | <biomart.Attribute name='cgchok1gshd_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cgchok1gshd_homolog_subtype | <biomart.Attribute name='cgchok1gshd_homolog_subtype', display_name='Last common ancestor with Chinese hamster CHOK1GS', description=''>
cgchok1gshd_homolog_orthology_type | <biomart.Attribute name='cgchok1gshd_homolog_orthology_type', display_name='Chinese hamster CHOK1GS homology type', description=''>
cgchok1gshd_homolog_perc_id | <biomart.Attribute name='cgchok1gshd_homolog_perc_id', display_name='%id. target Chinese hamster CHOK1GS gene identical to query gene', description=''>
cgchok1gshd_homolog_perc_id_r1 | <biomart.Attribute name='cgchok1gshd_homolog_perc_id_r1', display_name='%id. query gene identical to target Chinese hamster CHOK1GS gene', description=''>
cgchok1gshd_homolog_goc_score | <biomart.Attribute name='cgchok1gshd_homolog_goc_score', display_name='Chinese hamster CHOK1GS Gene-order conservation score', description=''>
cgchok1gshd_homolog_wga_coverage | <biomart.Attribute name='cgchok1gshd_homolog_wga_coverage', display_name='Chinese hamster CHOK1GS Whole-genome alignment coverage', description=''>
cgchok1gshd_homolog_orthology_confidence | <biomart.Attribute name='cgchok1gshd_homolog_orthology_confidence', display_name='Chinese hamster CHOK1GS orthology confidence [0 low, 1 high]', description=''>
cgcrigri_homolog_ensembl_gene | <biomart.Attribute name='cgcrigri_homolog_ensembl_gene', display_name='Chinese hamster CriGri gene stable ID', description=''>
cgcrigri_homolog_associated_gene_name | <biomart.Attribute name='cgcrigri_homolog_associated_gene_name', display_name='Chinese hamster CriGri gene name', description=''>
cgcrigri_homolog_ensembl_peptide | <biomart.Attribute name='cgcrigri_homolog_ensembl_peptide', display_name='Chinese hamster CriGri protein or transcript stable ID', description=''>
cgcrigri_homolog_chromosome | <biomart.Attribute name='cgcrigri_homolog_chromosome', display_name='Chinese hamster CriGri chromosome/scaffold name', description=''>
cgcrigri_homolog_chrom_start | <biomart.Attribute name='cgcrigri_homolog_chrom_start', display_name='Chinese hamster CriGri chromosome/scaffold start (bp)', description=''>
cgcrigri_homolog_chrom_end | <biomart.Attribute name='cgcrigri_homolog_chrom_end', display_name='Chinese hamster CriGri chromosome/scaffold end (bp)', description=''>
cgcrigri_homolog_canonical_transcript_protein | <biomart.Attribute name='cgcrigri_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cgcrigri_homolog_subtype | <biomart.Attribute name='cgcrigri_homolog_subtype', display_name='Last common ancestor with Chinese hamster CriGri', description=''>
cgcrigri_homolog_orthology_type | <biomart.Attribute name='cgcrigri_homolog_orthology_type', display_name='Chinese hamster CriGri homology type', description=''>
cgcrigri_homolog_perc_id | <biomart.Attribute name='cgcrigri_homolog_perc_id', display_name='%id. target Chinese hamster CriGri gene identical to query gene', description=''>
cgcrigri_homolog_perc_id_r1 | <biomart.Attribute name='cgcrigri_homolog_perc_id_r1', display_name='%id. query gene identical to target Chinese hamster CriGri gene', description=''>
cgcrigri_homolog_goc_score | <biomart.Attribute name='cgcrigri_homolog_goc_score', display_name='Chinese hamster CriGri Gene-order conservation score', description=''>
cgcrigri_homolog_wga_coverage | <biomart.Attribute name='cgcrigri_homolog_wga_coverage', display_name='Chinese hamster CriGri Whole-genome alignment coverage', description=''>
cgcrigri_homolog_orthology_confidence | <biomart.Attribute name='cgcrigri_homolog_orthology_confidence', display_name='Chinese hamster CriGri orthology confidence [0 low, 1 high]', description=''>
cgpicr_homolog_ensembl_gene | <biomart.Attribute name='cgpicr_homolog_ensembl_gene', display_name='Chinese hamster PICR gene stable ID', description=''>
cgpicr_homolog_associated_gene_name | <biomart.Attribute name='cgpicr_homolog_associated_gene_name', display_name='Chinese hamster PICR gene name', description=''>
cgpicr_homolog_ensembl_peptide | <biomart.Attribute name='cgpicr_homolog_ensembl_peptide', display_name='Chinese hamster PICR protein or transcript stable ID', description=''>
cgpicr_homolog_chromosome | <biomart.Attribute name='cgpicr_homolog_chromosome', display_name='Chinese hamster PICR chromosome/scaffold name', description=''>
cgpicr_homolog_chrom_start | <biomart.Attribute name='cgpicr_homolog_chrom_start', display_name='Chinese hamster PICR chromosome/scaffold start (bp)', description=''>
cgpicr_homolog_chrom_end | <biomart.Attribute name='cgpicr_homolog_chrom_end', display_name='Chinese hamster PICR chromosome/scaffold end (bp)', description=''>
cgpicr_homolog_canonical_transcript_protein | <biomart.Attribute name='cgpicr_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cgpicr_homolog_subtype | <biomart.Attribute name='cgpicr_homolog_subtype', display_name='Last common ancestor with Chinese hamster PICR', description=''>
cgpicr_homolog_orthology_type | <biomart.Attribute name='cgpicr_homolog_orthology_type', display_name='Chinese hamster PICR homology type', description=''>
cgpicr_homolog_perc_id | <biomart.Attribute name='cgpicr_homolog_perc_id', display_name='%id. target Chinese hamster PICR gene identical to query gene', description=''>
cgpicr_homolog_perc_id_r1 | <biomart.Attribute name='cgpicr_homolog_perc_id_r1', display_name='%id. query gene identical to target Chinese hamster PICR gene', description=''>
cgpicr_homolog_goc_score | <biomart.Attribute name='cgpicr_homolog_goc_score', display_name='Chinese hamster PICR Gene-order conservation score', description=''>
cgpicr_homolog_wga_coverage | <biomart.Attribute name='cgpicr_homolog_wga_coverage', display_name='Chinese hamster PICR Whole-genome alignment coverage', description=''>
cgpicr_homolog_orthology_confidence | <biomart.Attribute name='cgpicr_homolog_orthology_confidence', display_name='Chinese hamster PICR orthology confidence [0 low, 1 high]', description=''>
psinensis_homolog_ensembl_gene | <biomart.Attribute name='psinensis_homolog_ensembl_gene', display_name='Chinese softshell turtle gene stable ID', description=''>
psinensis_homolog_associated_gene_name | <biomart.Attribute name='psinensis_homolog_associated_gene_name', display_name='Chinese softshell turtle gene name', description=''>
psinensis_homolog_ensembl_peptide | <biomart.Attribute name='psinensis_homolog_ensembl_peptide', display_name='Chinese softshell turtle protein or transcript stable ID', description=''>
psinensis_homolog_chromosome | <biomart.Attribute name='psinensis_homolog_chromosome', display_name='Chinese softshell turtle chromosome/scaffold name', description=''>
psinensis_homolog_chrom_start | <biomart.Attribute name='psinensis_homolog_chrom_start', display_name='Chinese softshell turtle chromosome/scaffold start (bp)', description=''>
psinensis_homolog_chrom_end | <biomart.Attribute name='psinensis_homolog_chrom_end', display_name='Chinese softshell turtle chromosome/scaffold end (bp)', description=''>
psinensis_homolog_canonical_transcript_protein | <biomart.Attribute name='psinensis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
psinensis_homolog_subtype | <biomart.Attribute name='psinensis_homolog_subtype', display_name='Last common ancestor with Chinese softshell turtle', description=''>
psinensis_homolog_orthology_type | <biomart.Attribute name='psinensis_homolog_orthology_type', display_name='Chinese softshell turtle homology type', description=''>
psinensis_homolog_perc_id | <biomart.Attribute name='psinensis_homolog_perc_id', display_name='%id. target Chinese softshell turtle gene identical to query gene', description=''>
psinensis_homolog_perc_id_r1 | <biomart.Attribute name='psinensis_homolog_perc_id_r1', display_name='%id. query gene identical to target Chinese softshell turtle gene', description=''>
psinensis_homolog_goc_score | <biomart.Attribute name='psinensis_homolog_goc_score', display_name='Chinese softshell turtle Gene-order conservation score', description=''>
psinensis_homolog_wga_coverage | <biomart.Attribute name='psinensis_homolog_wga_coverage', display_name='Chinese softshell turtle Whole-genome alignment coverage', description=''>
psinensis_homolog_orthology_confidence | <biomart.Attribute name='psinensis_homolog_orthology_confidence', display_name='Chinese softshell turtle orthology confidence [0 low, 1 high]', description=''>
atestudineus_homolog_ensembl_gene | <biomart.Attribute name='atestudineus_homolog_ensembl_gene', display_name='Climbing perch gene stable ID', description=''>
atestudineus_homolog_associated_gene_name | <biomart.Attribute name='atestudineus_homolog_associated_gene_name', display_name='Climbing perch gene name', description=''>
atestudineus_homolog_ensembl_peptide | <biomart.Attribute name='atestudineus_homolog_ensembl_peptide', display_name='Climbing perch protein or transcript stable ID', description=''>
atestudineus_homolog_chromosome | <biomart.Attribute name='atestudineus_homolog_chromosome', display_name='Climbing perch chromosome/scaffold name', description=''>
atestudineus_homolog_chrom_start | <biomart.Attribute name='atestudineus_homolog_chrom_start', display_name='Climbing perch chromosome/scaffold start (bp)', description=''>
atestudineus_homolog_chrom_end | <biomart.Attribute name='atestudineus_homolog_chrom_end', display_name='Climbing perch chromosome/scaffold end (bp)', description=''>
atestudineus_homolog_canonical_transcript_protein | <biomart.Attribute name='atestudineus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
atestudineus_homolog_subtype | <biomart.Attribute name='atestudineus_homolog_subtype', display_name='Last common ancestor with Climbing perch', description=''>
atestudineus_homolog_orthology_type | <biomart.Attribute name='atestudineus_homolog_orthology_type', display_name='Climbing perch homology type', description=''>
atestudineus_homolog_perc_id | <biomart.Attribute name='atestudineus_homolog_perc_id', display_name='%id. target Climbing perch gene identical to query gene', description=''>
atestudineus_homolog_perc_id_r1 | <biomart.Attribute name='atestudineus_homolog_perc_id_r1', display_name='%id. query gene identical to target Climbing perch gene', description=''>
atestudineus_homolog_goc_score | <biomart.Attribute name='atestudineus_homolog_goc_score', display_name='Climbing perch Gene-order conservation score', description=''>
atestudineus_homolog_wga_coverage | <biomart.Attribute name='atestudineus_homolog_wga_coverage', display_name='Climbing perch Whole-genome alignment coverage', description=''>
atestudineus_homolog_orthology_confidence | <biomart.Attribute name='atestudineus_homolog_orthology_confidence', display_name='Climbing perch orthology confidence [0 low, 1 high]', description=''>
gmorhua_homolog_ensembl_gene | <biomart.Attribute name='gmorhua_homolog_ensembl_gene', display_name='Cod gene stable ID', description=''>
gmorhua_homolog_associated_gene_name | <biomart.Attribute name='gmorhua_homolog_associated_gene_name', display_name='Cod gene name', description=''>
gmorhua_homolog_ensembl_peptide | <biomart.Attribute name='gmorhua_homolog_ensembl_peptide', display_name='Cod protein or transcript stable ID', description=''>
gmorhua_homolog_chromosome | <biomart.Attribute name='gmorhua_homolog_chromosome', display_name='Cod chromosome/scaffold name', description=''>
gmorhua_homolog_chrom_start | <biomart.Attribute name='gmorhua_homolog_chrom_start', display_name='Cod chromosome/scaffold start (bp)', description=''>
gmorhua_homolog_chrom_end | <biomart.Attribute name='gmorhua_homolog_chrom_end', display_name='Cod chromosome/scaffold end (bp)', description=''>
gmorhua_homolog_canonical_transcript_protein | <biomart.Attribute name='gmorhua_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
gmorhua_homolog_subtype | <biomart.Attribute name='gmorhua_homolog_subtype', display_name='Last common ancestor with Cod', description=''>
gmorhua_homolog_orthology_type | <biomart.Attribute name='gmorhua_homolog_orthology_type', display_name='Cod homology type', description=''>
gmorhua_homolog_perc_id | <biomart.Attribute name='gmorhua_homolog_perc_id', display_name='%id. target Cod gene identical to query gene', description=''>
gmorhua_homolog_perc_id_r1 | <biomart.Attribute name='gmorhua_homolog_perc_id_r1', display_name='%id. query gene identical to target Cod gene', description=''>
gmorhua_homolog_goc_score | <biomart.Attribute name='gmorhua_homolog_goc_score', display_name='Cod Gene-order conservation score', description=''>
gmorhua_homolog_wga_coverage | <biomart.Attribute name='gmorhua_homolog_wga_coverage', display_name='Cod Whole-genome alignment coverage', description=''>
gmorhua_homolog_orthology_confidence | <biomart.Attribute name='gmorhua_homolog_orthology_confidence', display_name='Cod orthology confidence [0 low, 1 high]', description=''>
lchalumnae_homolog_ensembl_gene | <biomart.Attribute name='lchalumnae_homolog_ensembl_gene', display_name='Coelacanth gene stable ID', description=''>
lchalumnae_homolog_associated_gene_name | <biomart.Attribute name='lchalumnae_homolog_associated_gene_name', display_name='Coelacanth gene name', description=''>
lchalumnae_homolog_ensembl_peptide | <biomart.Attribute name='lchalumnae_homolog_ensembl_peptide', display_name='Coelacanth protein or transcript stable ID', description=''>
lchalumnae_homolog_chromosome | <biomart.Attribute name='lchalumnae_homolog_chromosome', display_name='Coelacanth chromosome/scaffold name', description=''>
lchalumnae_homolog_chrom_start | <biomart.Attribute name='lchalumnae_homolog_chrom_start', display_name='Coelacanth chromosome/scaffold start (bp)', description=''>
lchalumnae_homolog_chrom_end | <biomart.Attribute name='lchalumnae_homolog_chrom_end', display_name='Coelacanth chromosome/scaffold end (bp)', description=''>
lchalumnae_homolog_canonical_transcript_protein | <biomart.Attribute name='lchalumnae_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lchalumnae_homolog_subtype | <biomart.Attribute name='lchalumnae_homolog_subtype', display_name='Last common ancestor with Coelacanth', description=''>
lchalumnae_homolog_orthology_type | <biomart.Attribute name='lchalumnae_homolog_orthology_type', display_name='Coelacanth homology type', description=''>
lchalumnae_homolog_perc_id | <biomart.Attribute name='lchalumnae_homolog_perc_id', display_name='%id. target Coelacanth gene identical to query gene', description=''>
lchalumnae_homolog_perc_id_r1 | <biomart.Attribute name='lchalumnae_homolog_perc_id_r1', display_name='%id. query gene identical to target Coelacanth gene', description=''>
lchalumnae_homolog_goc_score | <biomart.Attribute name='lchalumnae_homolog_goc_score', display_name='Coelacanth Gene-order conservation score', description=''>
lchalumnae_homolog_wga_coverage | <biomart.Attribute name='lchalumnae_homolog_wga_coverage', display_name='Coelacanth Whole-genome alignment coverage', description=''>
lchalumnae_homolog_orthology_confidence | <biomart.Attribute name='lchalumnae_homolog_orthology_confidence', display_name='Coelacanth orthology confidence [0 low, 1 high]', description=''>
scanaria_homolog_ensembl_gene | <biomart.Attribute name='scanaria_homolog_ensembl_gene', display_name='Common canary gene stable ID', description=''>
scanaria_homolog_associated_gene_name | <biomart.Attribute name='scanaria_homolog_associated_gene_name', display_name='Common canary gene name', description=''>
scanaria_homolog_ensembl_peptide | <biomart.Attribute name='scanaria_homolog_ensembl_peptide', display_name='Common canary protein or transcript stable ID', description=''>
scanaria_homolog_chromosome | <biomart.Attribute name='scanaria_homolog_chromosome', display_name='Common canary chromosome/scaffold name', description=''>
scanaria_homolog_chrom_start | <biomart.Attribute name='scanaria_homolog_chrom_start', display_name='Common canary chromosome/scaffold start (bp)', description=''>
scanaria_homolog_chrom_end | <biomart.Attribute name='scanaria_homolog_chrom_end', display_name='Common canary chromosome/scaffold end (bp)', description=''>
scanaria_homolog_canonical_transcript_protein | <biomart.Attribute name='scanaria_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
scanaria_homolog_subtype | <biomart.Attribute name='scanaria_homolog_subtype', display_name='Last common ancestor with Common canary', description=''>
scanaria_homolog_orthology_type | <biomart.Attribute name='scanaria_homolog_orthology_type', display_name='Common canary homology type', description=''>
scanaria_homolog_perc_id | <biomart.Attribute name='scanaria_homolog_perc_id', display_name='%id. target Common canary gene identical to query gene', description=''>
scanaria_homolog_perc_id_r1 | <biomart.Attribute name='scanaria_homolog_perc_id_r1', display_name='%id. query gene identical to target Common canary gene', description=''>
scanaria_homolog_goc_score | <biomart.Attribute name='scanaria_homolog_goc_score', display_name='Common canary Gene-order conservation score', description=''>
scanaria_homolog_wga_coverage | <biomart.Attribute name='scanaria_homolog_wga_coverage', display_name='Common canary Whole-genome alignment coverage', description=''>
scanaria_homolog_orthology_confidence | <biomart.Attribute name='scanaria_homolog_orthology_confidence', display_name='Common canary orthology confidence [0 low, 1 high]', description=''>
vursinus_homolog_ensembl_gene | <biomart.Attribute name='vursinus_homolog_ensembl_gene', display_name='Common wombat gene stable ID', description=''>
vursinus_homolog_associated_gene_name | <biomart.Attribute name='vursinus_homolog_associated_gene_name', display_name='Common wombat gene name', description=''>
vursinus_homolog_ensembl_peptide | <biomart.Attribute name='vursinus_homolog_ensembl_peptide', display_name='Common wombat protein or transcript stable ID', description=''>
vursinus_homolog_chromosome | <biomart.Attribute name='vursinus_homolog_chromosome', display_name='Common wombat chromosome/scaffold name', description=''>
vursinus_homolog_chrom_start | <biomart.Attribute name='vursinus_homolog_chrom_start', display_name='Common wombat chromosome/scaffold start (bp)', description=''>
vursinus_homolog_chrom_end | <biomart.Attribute name='vursinus_homolog_chrom_end', display_name='Common wombat chromosome/scaffold end (bp)', description=''>
vursinus_homolog_canonical_transcript_protein | <biomart.Attribute name='vursinus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
vursinus_homolog_subtype | <biomart.Attribute name='vursinus_homolog_subtype', display_name='Last common ancestor with Common wombat', description=''>
vursinus_homolog_orthology_type | <biomart.Attribute name='vursinus_homolog_orthology_type', display_name='Common wombat homology type', description=''>
vursinus_homolog_perc_id | <biomart.Attribute name='vursinus_homolog_perc_id', display_name='%id. target Common wombat gene identical to query gene', description=''>
vursinus_homolog_perc_id_r1 | <biomart.Attribute name='vursinus_homolog_perc_id_r1', display_name='%id. query gene identical to target Common wombat gene', description=''>
vursinus_homolog_goc_score | <biomart.Attribute name='vursinus_homolog_goc_score', display_name='Common wombat Gene-order conservation score', description=''>
vursinus_homolog_wga_coverage | <biomart.Attribute name='vursinus_homolog_wga_coverage', display_name='Common wombat Whole-genome alignment coverage', description=''>
vursinus_homolog_orthology_confidence | <biomart.Attribute name='vursinus_homolog_orthology_confidence', display_name='Common wombat orthology confidence [0 low, 1 high]', description=''>
pcoquereli_homolog_ensembl_gene | <biomart.Attribute name='pcoquereli_homolog_ensembl_gene', display_name="Coquerel's sifaka gene stable ID", description=''>
pcoquereli_homolog_associated_gene_name | <biomart.Attribute name='pcoquereli_homolog_associated_gene_name', display_name="Coquerel's sifaka gene name", description=''>
pcoquereli_homolog_ensembl_peptide | <biomart.Attribute name='pcoquereli_homolog_ensembl_peptide', display_name="Coquerel's sifaka protein or transcript stable ID", description=''>
pcoquereli_homolog_chromosome | <biomart.Attribute name='pcoquereli_homolog_chromosome', display_name="Coquerel's sifaka chromosome/scaffold name", description=''>
pcoquereli_homolog_chrom_start | <biomart.Attribute name='pcoquereli_homolog_chrom_start', display_name="Coquerel's sifaka chromosome/scaffold start (bp)", description=''>
pcoquereli_homolog_chrom_end | <biomart.Attribute name='pcoquereli_homolog_chrom_end', display_name="Coquerel's sifaka chromosome/scaffold end (bp)", description=''>
pcoquereli_homolog_canonical_transcript_protein | <biomart.Attribute name='pcoquereli_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pcoquereli_homolog_subtype | <biomart.Attribute name='pcoquereli_homolog_subtype', display_name="Last common ancestor with Coquerel's sifaka", description=''>
pcoquereli_homolog_orthology_type | <biomart.Attribute name='pcoquereli_homolog_orthology_type', display_name="Coquerel's sifaka homology type", description=''>
pcoquereli_homolog_perc_id | <biomart.Attribute name='pcoquereli_homolog_perc_id', display_name="%id. target Coquerel's sifaka gene identical to query gene", description=''>
pcoquereli_homolog_perc_id_r1 | <biomart.Attribute name='pcoquereli_homolog_perc_id_r1', display_name="%id. query gene identical to target Coquerel's sifaka gene", description=''>
pcoquereli_homolog_goc_score | <biomart.Attribute name='pcoquereli_homolog_goc_score', display_name="Coquerel's sifaka Gene-order conservation score", description=''>
pcoquereli_homolog_wga_coverage | <biomart.Attribute name='pcoquereli_homolog_wga_coverage', display_name="Coquerel's sifaka Whole-genome alignment coverage", description=''>
pcoquereli_homolog_orthology_confidence | <biomart.Attribute name='pcoquereli_homolog_orthology_confidence', display_name="Coquerel's sifaka orthology confidence [0 low, 1 high]", description=''>
btaurus_homolog_ensembl_gene | <biomart.Attribute name='btaurus_homolog_ensembl_gene', display_name='Cow gene stable ID', description=''>
btaurus_homolog_associated_gene_name | <biomart.Attribute name='btaurus_homolog_associated_gene_name', display_name='Cow gene name', description=''>
btaurus_homolog_ensembl_peptide | <biomart.Attribute name='btaurus_homolog_ensembl_peptide', display_name='Cow protein or transcript stable ID', description=''>
btaurus_homolog_chromosome | <biomart.Attribute name='btaurus_homolog_chromosome', display_name='Cow chromosome/scaffold name', description=''>
btaurus_homolog_chrom_start | <biomart.Attribute name='btaurus_homolog_chrom_start', display_name='Cow chromosome/scaffold start (bp)', description=''>
btaurus_homolog_chrom_end | <biomart.Attribute name='btaurus_homolog_chrom_end', display_name='Cow chromosome/scaffold end (bp)', description=''>
btaurus_homolog_canonical_transcript_protein | <biomart.Attribute name='btaurus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
btaurus_homolog_subtype | <biomart.Attribute name='btaurus_homolog_subtype', display_name='Last common ancestor with Cow', description=''>
btaurus_homolog_orthology_type | <biomart.Attribute name='btaurus_homolog_orthology_type', display_name='Cow homology type', description=''>
btaurus_homolog_perc_id | <biomart.Attribute name='btaurus_homolog_perc_id', display_name='%id. target Cow gene identical to query gene', description=''>
btaurus_homolog_perc_id_r1 | <biomart.Attribute name='btaurus_homolog_perc_id_r1', display_name='%id. query gene identical to target Cow gene', description=''>
btaurus_homolog_goc_score | <biomart.Attribute name='btaurus_homolog_goc_score', display_name='Cow Gene-order conservation score', description=''>
btaurus_homolog_wga_coverage | <biomart.Attribute name='btaurus_homolog_wga_coverage', display_name='Cow Whole-genome alignment coverage', description=''>
btaurus_homolog_orthology_confidence | <biomart.Attribute name='btaurus_homolog_orthology_confidence', display_name='Cow orthology confidence [0 low, 1 high]', description=''>
mfascicularis_homolog_ensembl_gene | <biomart.Attribute name='mfascicularis_homolog_ensembl_gene', display_name='Crab-eating macaque gene stable ID', description=''>
mfascicularis_homolog_associated_gene_name | <biomart.Attribute name='mfascicularis_homolog_associated_gene_name', display_name='Crab-eating macaque gene name', description=''>
mfascicularis_homolog_ensembl_peptide | <biomart.Attribute name='mfascicularis_homolog_ensembl_peptide', display_name='Crab-eating macaque protein or transcript stable ID', description=''>
mfascicularis_homolog_chromosome | <biomart.Attribute name='mfascicularis_homolog_chromosome', display_name='Crab-eating macaque chromosome/scaffold name', description=''>
mfascicularis_homolog_chrom_start | <biomart.Attribute name='mfascicularis_homolog_chrom_start', display_name='Crab-eating macaque chromosome/scaffold start (bp)', description=''>
mfascicularis_homolog_chrom_end | <biomart.Attribute name='mfascicularis_homolog_chrom_end', display_name='Crab-eating macaque chromosome/scaffold end (bp)', description=''>
mfascicularis_homolog_canonical_transcript_protein | <biomart.Attribute name='mfascicularis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mfascicularis_homolog_subtype | <biomart.Attribute name='mfascicularis_homolog_subtype', display_name='Last common ancestor with Crab-eating macaque', description=''>
mfascicularis_homolog_orthology_type | <biomart.Attribute name='mfascicularis_homolog_orthology_type', display_name='Crab-eating macaque homology type', description=''>
mfascicularis_homolog_perc_id | <biomart.Attribute name='mfascicularis_homolog_perc_id', display_name='%id. target Crab-eating macaque gene identical to query gene', description=''>
mfascicularis_homolog_perc_id_r1 | <biomart.Attribute name='mfascicularis_homolog_perc_id_r1', display_name='%id. query gene identical to target Crab-eating macaque gene', description=''>
mfascicularis_homolog_goc_score | <biomart.Attribute name='mfascicularis_homolog_goc_score', display_name='Crab-eating macaque Gene-order conservation score', description=''>
mfascicularis_homolog_wga_coverage | <biomart.Attribute name='mfascicularis_homolog_wga_coverage', display_name='Crab-eating macaque Whole-genome alignment coverage', description=''>
mfascicularis_homolog_orthology_confidence | <biomart.Attribute name='mfascicularis_homolog_orthology_confidence', display_name='Crab-eating macaque orthology confidence [0 low, 1 high]', description=''>
fdamarensis_homolog_ensembl_gene | <biomart.Attribute name='fdamarensis_homolog_ensembl_gene', display_name='Damara mole rat gene stable ID', description=''>
fdamarensis_homolog_associated_gene_name | <biomart.Attribute name='fdamarensis_homolog_associated_gene_name', display_name='Damara mole rat gene name', description=''>
fdamarensis_homolog_ensembl_peptide | <biomart.Attribute name='fdamarensis_homolog_ensembl_peptide', display_name='Damara mole rat protein or transcript stable ID', description=''>
fdamarensis_homolog_chromosome | <biomart.Attribute name='fdamarensis_homolog_chromosome', display_name='Damara mole rat chromosome/scaffold name', description=''>
fdamarensis_homolog_chrom_start | <biomart.Attribute name='fdamarensis_homolog_chrom_start', display_name='Damara mole rat chromosome/scaffold start (bp)', description=''>
fdamarensis_homolog_chrom_end | <biomart.Attribute name='fdamarensis_homolog_chrom_end', display_name='Damara mole rat chromosome/scaffold end (bp)', description=''>
fdamarensis_homolog_canonical_transcript_protein | <biomart.Attribute name='fdamarensis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
fdamarensis_homolog_subtype | <biomart.Attribute name='fdamarensis_homolog_subtype', display_name='Last common ancestor with Damara mole rat', description=''>
fdamarensis_homolog_orthology_type | <biomart.Attribute name='fdamarensis_homolog_orthology_type', display_name='Damara mole rat homology type', description=''>
fdamarensis_homolog_perc_id | <biomart.Attribute name='fdamarensis_homolog_perc_id', display_name='%id. target Damara mole rat gene identical to query gene', description=''>
fdamarensis_homolog_perc_id_r1 | <biomart.Attribute name='fdamarensis_homolog_perc_id_r1', display_name='%id. query gene identical to target Damara mole rat gene', description=''>
fdamarensis_homolog_goc_score | <biomart.Attribute name='fdamarensis_homolog_goc_score', display_name='Damara mole rat Gene-order conservation score', description=''>
fdamarensis_homolog_wga_coverage | <biomart.Attribute name='fdamarensis_homolog_wga_coverage', display_name='Damara mole rat Whole-genome alignment coverage', description=''>
fdamarensis_homolog_orthology_confidence | <biomart.Attribute name='fdamarensis_homolog_orthology_confidence', display_name='Damara mole rat orthology confidence [0 low, 1 high]', description=''>
jhyemalis_homolog_ensembl_gene | <biomart.Attribute name='jhyemalis_homolog_ensembl_gene', display_name='Dark-eyed junco gene stable ID', description=''>
jhyemalis_homolog_associated_gene_name | <biomart.Attribute name='jhyemalis_homolog_associated_gene_name', display_name='Dark-eyed junco gene name', description=''>
jhyemalis_homolog_ensembl_peptide | <biomart.Attribute name='jhyemalis_homolog_ensembl_peptide', display_name='Dark-eyed junco protein or transcript stable ID', description=''>
jhyemalis_homolog_chromosome | <biomart.Attribute name='jhyemalis_homolog_chromosome', display_name='Dark-eyed junco chromosome/scaffold name', description=''>
jhyemalis_homolog_chrom_start | <biomart.Attribute name='jhyemalis_homolog_chrom_start', display_name='Dark-eyed junco chromosome/scaffold start (bp)', description=''>
jhyemalis_homolog_chrom_end | <biomart.Attribute name='jhyemalis_homolog_chrom_end', display_name='Dark-eyed junco chromosome/scaffold end (bp)', description=''>
jhyemalis_homolog_canonical_transcript_protein | <biomart.Attribute name='jhyemalis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
jhyemalis_homolog_subtype | <biomart.Attribute name='jhyemalis_homolog_subtype', display_name='Last common ancestor with Dark-eyed junco', description=''>
jhyemalis_homolog_orthology_type | <biomart.Attribute name='jhyemalis_homolog_orthology_type', display_name='Dark-eyed junco homology type', description=''>
jhyemalis_homolog_perc_id | <biomart.Attribute name='jhyemalis_homolog_perc_id', display_name='%id. target Dark-eyed junco gene identical to query gene', description=''>
jhyemalis_homolog_perc_id_r1 | <biomart.Attribute name='jhyemalis_homolog_perc_id_r1', display_name='%id. query gene identical to target Dark-eyed junco gene', description=''>
jhyemalis_homolog_goc_score | <biomart.Attribute name='jhyemalis_homolog_goc_score', display_name='Dark-eyed junco Gene-order conservation score', description=''>
jhyemalis_homolog_wga_coverage | <biomart.Attribute name='jhyemalis_homolog_wga_coverage', display_name='Dark-eyed junco Whole-genome alignment coverage', description=''>
jhyemalis_homolog_orthology_confidence | <biomart.Attribute name='jhyemalis_homolog_orthology_confidence', display_name='Dark-eyed junco orthology confidence [0 low, 1 high]', description=''>
odegus_homolog_ensembl_gene | <biomart.Attribute name='odegus_homolog_ensembl_gene', display_name='Degu gene stable ID', description=''>
odegus_homolog_associated_gene_name | <biomart.Attribute name='odegus_homolog_associated_gene_name', display_name='Degu gene name', description=''>
odegus_homolog_ensembl_peptide | <biomart.Attribute name='odegus_homolog_ensembl_peptide', display_name='Degu protein or transcript stable ID', description=''>
odegus_homolog_chromosome | <biomart.Attribute name='odegus_homolog_chromosome', display_name='Degu chromosome/scaffold name', description=''>
odegus_homolog_chrom_start | <biomart.Attribute name='odegus_homolog_chrom_start', display_name='Degu chromosome/scaffold start (bp)', description=''>
odegus_homolog_chrom_end | <biomart.Attribute name='odegus_homolog_chrom_end', display_name='Degu chromosome/scaffold end (bp)', description=''>
odegus_homolog_canonical_transcript_protein | <biomart.Attribute name='odegus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
odegus_homolog_subtype | <biomart.Attribute name='odegus_homolog_subtype', display_name='Last common ancestor with Degu', description=''>
odegus_homolog_orthology_type | <biomart.Attribute name='odegus_homolog_orthology_type', display_name='Degu homology type', description=''>
odegus_homolog_perc_id | <biomart.Attribute name='odegus_homolog_perc_id', display_name='%id. target Degu gene identical to query gene', description=''>
odegus_homolog_perc_id_r1 | <biomart.Attribute name='odegus_homolog_perc_id_r1', display_name='%id. query gene identical to target Degu gene', description=''>
odegus_homolog_goc_score | <biomart.Attribute name='odegus_homolog_goc_score', display_name='Degu Gene-order conservation score', description=''>
odegus_homolog_wga_coverage | <biomart.Attribute name='odegus_homolog_wga_coverage', display_name='Degu Whole-genome alignment coverage', description=''>
odegus_homolog_orthology_confidence | <biomart.Attribute name='odegus_homolog_orthology_confidence', display_name='Degu orthology confidence [0 low, 1 high]', description=''>
dclupeoides_homolog_ensembl_gene | <biomart.Attribute name='dclupeoides_homolog_ensembl_gene', display_name='Denticle herring gene stable ID', description=''>
dclupeoides_homolog_associated_gene_name | <biomart.Attribute name='dclupeoides_homolog_associated_gene_name', display_name='Denticle herring gene name', description=''>
dclupeoides_homolog_ensembl_peptide | <biomart.Attribute name='dclupeoides_homolog_ensembl_peptide', display_name='Denticle herring protein or transcript stable ID', description=''>
dclupeoides_homolog_chromosome | <biomart.Attribute name='dclupeoides_homolog_chromosome', display_name='Denticle herring chromosome/scaffold name', description=''>
dclupeoides_homolog_chrom_start | <biomart.Attribute name='dclupeoides_homolog_chrom_start', display_name='Denticle herring chromosome/scaffold start (bp)', description=''>
dclupeoides_homolog_chrom_end | <biomart.Attribute name='dclupeoides_homolog_chrom_end', display_name='Denticle herring chromosome/scaffold end (bp)', description=''>
dclupeoides_homolog_canonical_transcript_protein | <biomart.Attribute name='dclupeoides_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
dclupeoides_homolog_subtype | <biomart.Attribute name='dclupeoides_homolog_subtype', display_name='Last common ancestor with Denticle herring', description=''>
dclupeoides_homolog_orthology_type | <biomart.Attribute name='dclupeoides_homolog_orthology_type', display_name='Denticle herring homology type', description=''>
dclupeoides_homolog_perc_id | <biomart.Attribute name='dclupeoides_homolog_perc_id', display_name='%id. target Denticle herring gene identical to query gene', description=''>
dclupeoides_homolog_perc_id_r1 | <biomart.Attribute name='dclupeoides_homolog_perc_id_r1', display_name='%id. query gene identical to target Denticle herring gene', description=''>
dclupeoides_homolog_goc_score | <biomart.Attribute name='dclupeoides_homolog_goc_score', display_name='Denticle herring Gene-order conservation score', description=''>
dclupeoides_homolog_wga_coverage | <biomart.Attribute name='dclupeoides_homolog_wga_coverage', display_name='Denticle herring Whole-genome alignment coverage', description=''>
dclupeoides_homolog_orthology_confidence | <biomart.Attribute name='dclupeoides_homolog_orthology_confidence', display_name='Denticle herring orthology confidence [0 low, 1 high]', description=''>
cldingo_homolog_ensembl_gene | <biomart.Attribute name='cldingo_homolog_ensembl_gene', display_name='Dingo gene stable ID', description=''>
cldingo_homolog_associated_gene_name | <biomart.Attribute name='cldingo_homolog_associated_gene_name', display_name='Dingo gene name', description=''>
cldingo_homolog_ensembl_peptide | <biomart.Attribute name='cldingo_homolog_ensembl_peptide', display_name='Dingo protein or transcript stable ID', description=''>
cldingo_homolog_chromosome | <biomart.Attribute name='cldingo_homolog_chromosome', display_name='Dingo chromosome/scaffold name', description=''>
cldingo_homolog_chrom_start | <biomart.Attribute name='cldingo_homolog_chrom_start', display_name='Dingo chromosome/scaffold start (bp)', description=''>
cldingo_homolog_chrom_end | <biomart.Attribute name='cldingo_homolog_chrom_end', display_name='Dingo chromosome/scaffold end (bp)', description=''>
cldingo_homolog_canonical_transcript_protein | <biomart.Attribute name='cldingo_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cldingo_homolog_subtype | <biomart.Attribute name='cldingo_homolog_subtype', display_name='Last common ancestor with Dingo', description=''>
cldingo_homolog_orthology_type | <biomart.Attribute name='cldingo_homolog_orthology_type', display_name='Dingo homology type', description=''>
cldingo_homolog_perc_id | <biomart.Attribute name='cldingo_homolog_perc_id', display_name='%id. target Dingo gene identical to query gene', description=''>
cldingo_homolog_perc_id_r1 | <biomart.Attribute name='cldingo_homolog_perc_id_r1', display_name='%id. query gene identical to target Dingo gene', description=''>
cldingo_homolog_goc_score | <biomart.Attribute name='cldingo_homolog_goc_score', display_name='Dingo Gene-order conservation score', description=''>
cldingo_homolog_wga_coverage | <biomart.Attribute name='cldingo_homolog_wga_coverage', display_name='Dingo Whole-genome alignment coverage', description=''>
cldingo_homolog_orthology_confidence | <biomart.Attribute name='cldingo_homolog_orthology_confidence', display_name='Dingo orthology confidence [0 low, 1 high]', description=''>
clfamiliaris_homolog_ensembl_gene | <biomart.Attribute name='clfamiliaris_homolog_ensembl_gene', display_name='Dog gene stable ID', description=''>
clfamiliaris_homolog_associated_gene_name | <biomart.Attribute name='clfamiliaris_homolog_associated_gene_name', display_name='Dog gene name', description=''>
clfamiliaris_homolog_ensembl_peptide | <biomart.Attribute name='clfamiliaris_homolog_ensembl_peptide', display_name='Dog protein or transcript stable ID', description=''>
clfamiliaris_homolog_chromosome | <biomart.Attribute name='clfamiliaris_homolog_chromosome', display_name='Dog chromosome/scaffold name', description=''>
clfamiliaris_homolog_chrom_start | <biomart.Attribute name='clfamiliaris_homolog_chrom_start', display_name='Dog chromosome/scaffold start (bp)', description=''>
clfamiliaris_homolog_chrom_end | <biomart.Attribute name='clfamiliaris_homolog_chrom_end', display_name='Dog chromosome/scaffold end (bp)', description=''>
clfamiliaris_homolog_canonical_transcript_protein | <biomart.Attribute name='clfamiliaris_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
clfamiliaris_homolog_subtype | <biomart.Attribute name='clfamiliaris_homolog_subtype', display_name='Last common ancestor with Dog', description=''>
clfamiliaris_homolog_orthology_type | <biomart.Attribute name='clfamiliaris_homolog_orthology_type', display_name='Dog homology type', description=''>
clfamiliaris_homolog_perc_id | <biomart.Attribute name='clfamiliaris_homolog_perc_id', display_name='%id. target Dog gene identical to query gene', description=''>
clfamiliaris_homolog_perc_id_r1 | <biomart.Attribute name='clfamiliaris_homolog_perc_id_r1', display_name='%id. query gene identical to target Dog gene', description=''>
clfamiliaris_homolog_goc_score | <biomart.Attribute name='clfamiliaris_homolog_goc_score', display_name='Dog Gene-order conservation score', description=''>
clfamiliaris_homolog_wga_coverage | <biomart.Attribute name='clfamiliaris_homolog_wga_coverage', display_name='Dog Whole-genome alignment coverage', description=''>
clfamiliaris_homolog_orthology_confidence | <biomart.Attribute name='clfamiliaris_homolog_orthology_confidence', display_name='Dog orthology confidence [0 low, 1 high]', description=''>
ttruncatus_homolog_ensembl_gene | <biomart.Attribute name='ttruncatus_homolog_ensembl_gene', display_name='Dolphin gene stable ID', description=''>
ttruncatus_homolog_associated_gene_name | <biomart.Attribute name='ttruncatus_homolog_associated_gene_name', display_name='Dolphin gene name', description=''>
ttruncatus_homolog_ensembl_peptide | <biomart.Attribute name='ttruncatus_homolog_ensembl_peptide', display_name='Dolphin protein or transcript stable ID', description=''>
ttruncatus_homolog_chromosome | <biomart.Attribute name='ttruncatus_homolog_chromosome', display_name='Dolphin chromosome/scaffold name', description=''>
ttruncatus_homolog_chrom_start | <biomart.Attribute name='ttruncatus_homolog_chrom_start', display_name='Dolphin chromosome/scaffold start (bp)', description=''>
ttruncatus_homolog_chrom_end | <biomart.Attribute name='ttruncatus_homolog_chrom_end', display_name='Dolphin chromosome/scaffold end (bp)', description=''>
ttruncatus_homolog_canonical_transcript_protein | <biomart.Attribute name='ttruncatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ttruncatus_homolog_subtype | <biomart.Attribute name='ttruncatus_homolog_subtype', display_name='Last common ancestor with Dolphin', description=''>
ttruncatus_homolog_orthology_type | <biomart.Attribute name='ttruncatus_homolog_orthology_type', display_name='Dolphin homology type', description=''>
ttruncatus_homolog_perc_id | <biomart.Attribute name='ttruncatus_homolog_perc_id', display_name='%id. target Dolphin gene identical to query gene', description=''>
ttruncatus_homolog_perc_id_r1 | <biomart.Attribute name='ttruncatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Dolphin gene', description=''>
ttruncatus_homolog_goc_score | <biomart.Attribute name='ttruncatus_homolog_goc_score', display_name='Dolphin Gene-order conservation score', description=''>
ttruncatus_homolog_wga_coverage | <biomart.Attribute name='ttruncatus_homolog_wga_coverage', display_name='Dolphin Whole-genome alignment coverage', description=''>
ttruncatus_homolog_orthology_confidence | <biomart.Attribute name='ttruncatus_homolog_orthology_confidence', display_name='Dolphin orthology confidence [0 low, 1 high]', description=''>
bgrunniens_homolog_ensembl_gene | <biomart.Attribute name='bgrunniens_homolog_ensembl_gene', display_name='Domestic yak gene stable ID', description=''>
bgrunniens_homolog_associated_gene_name | <biomart.Attribute name='bgrunniens_homolog_associated_gene_name', display_name='Domestic yak gene name', description=''>
bgrunniens_homolog_ensembl_peptide | <biomart.Attribute name='bgrunniens_homolog_ensembl_peptide', display_name='Domestic yak protein or transcript stable ID', description=''>
bgrunniens_homolog_chromosome | <biomart.Attribute name='bgrunniens_homolog_chromosome', display_name='Domestic yak chromosome/scaffold name', description=''>
bgrunniens_homolog_chrom_start | <biomart.Attribute name='bgrunniens_homolog_chrom_start', display_name='Domestic yak chromosome/scaffold start (bp)', description=''>
bgrunniens_homolog_chrom_end | <biomart.Attribute name='bgrunniens_homolog_chrom_end', display_name='Domestic yak chromosome/scaffold end (bp)', description=''>
bgrunniens_homolog_canonical_transcript_protein | <biomart.Attribute name='bgrunniens_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
bgrunniens_homolog_subtype | <biomart.Attribute name='bgrunniens_homolog_subtype', display_name='Last common ancestor with Domestic yak', description=''>
bgrunniens_homolog_orthology_type | <biomart.Attribute name='bgrunniens_homolog_orthology_type', display_name='Domestic yak homology type', description=''>
bgrunniens_homolog_perc_id | <biomart.Attribute name='bgrunniens_homolog_perc_id', display_name='%id. target Domestic yak gene identical to query gene', description=''>
bgrunniens_homolog_perc_id_r1 | <biomart.Attribute name='bgrunniens_homolog_perc_id_r1', display_name='%id. query gene identical to target Domestic yak gene', description=''>
bgrunniens_homolog_goc_score | <biomart.Attribute name='bgrunniens_homolog_goc_score', display_name='Domestic yak Gene-order conservation score', description=''>
bgrunniens_homolog_wga_coverage | <biomart.Attribute name='bgrunniens_homolog_wga_coverage', display_name='Domestic yak Whole-genome alignment coverage', description=''>
bgrunniens_homolog_orthology_confidence | <biomart.Attribute name='bgrunniens_homolog_orthology_confidence', display_name='Domestic yak orthology confidence [0 low, 1 high]', description=''>
eaasinus_homolog_ensembl_gene | <biomart.Attribute name='eaasinus_homolog_ensembl_gene', display_name='Donkey gene stable ID', description=''>
eaasinus_homolog_associated_gene_name | <biomart.Attribute name='eaasinus_homolog_associated_gene_name', display_name='Donkey gene name', description=''>
eaasinus_homolog_ensembl_peptide | <biomart.Attribute name='eaasinus_homolog_ensembl_peptide', display_name='Donkey protein or transcript stable ID', description=''>
eaasinus_homolog_chromosome | <biomart.Attribute name='eaasinus_homolog_chromosome', display_name='Donkey chromosome/scaffold name', description=''>
eaasinus_homolog_chrom_start | <biomart.Attribute name='eaasinus_homolog_chrom_start', display_name='Donkey chromosome/scaffold start (bp)', description=''>
eaasinus_homolog_chrom_end | <biomart.Attribute name='eaasinus_homolog_chrom_end', display_name='Donkey chromosome/scaffold end (bp)', description=''>
eaasinus_homolog_canonical_transcript_protein | <biomart.Attribute name='eaasinus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
eaasinus_homolog_subtype | <biomart.Attribute name='eaasinus_homolog_subtype', display_name='Last common ancestor with Donkey', description=''>
eaasinus_homolog_orthology_type | <biomart.Attribute name='eaasinus_homolog_orthology_type', display_name='Donkey homology type', description=''>
eaasinus_homolog_perc_id | <biomart.Attribute name='eaasinus_homolog_perc_id', display_name='%id. target Donkey gene identical to query gene', description=''>
eaasinus_homolog_perc_id_r1 | <biomart.Attribute name='eaasinus_homolog_perc_id_r1', display_name='%id. query gene identical to target Donkey gene', description=''>
eaasinus_homolog_goc_score | <biomart.Attribute name='eaasinus_homolog_goc_score', display_name='Donkey Gene-order conservation score', description=''>
eaasinus_homolog_wga_coverage | <biomart.Attribute name='eaasinus_homolog_wga_coverage', display_name='Donkey Whole-genome alignment coverage', description=''>
eaasinus_homolog_orthology_confidence | <biomart.Attribute name='eaasinus_homolog_orthology_confidence', display_name='Donkey orthology confidence [0 low, 1 high]', description=''>
mleucophaeus_homolog_ensembl_gene | <biomart.Attribute name='mleucophaeus_homolog_ensembl_gene', display_name='Drill gene stable ID', description=''>
mleucophaeus_homolog_associated_gene_name | <biomart.Attribute name='mleucophaeus_homolog_associated_gene_name', display_name='Drill gene name', description=''>
mleucophaeus_homolog_ensembl_peptide | <biomart.Attribute name='mleucophaeus_homolog_ensembl_peptide', display_name='Drill protein or transcript stable ID', description=''>
mleucophaeus_homolog_chromosome | <biomart.Attribute name='mleucophaeus_homolog_chromosome', display_name='Drill chromosome/scaffold name', description=''>
mleucophaeus_homolog_chrom_start | <biomart.Attribute name='mleucophaeus_homolog_chrom_start', display_name='Drill chromosome/scaffold start (bp)', description=''>
mleucophaeus_homolog_chrom_end | <biomart.Attribute name='mleucophaeus_homolog_chrom_end', display_name='Drill chromosome/scaffold end (bp)', description=''>
mleucophaeus_homolog_canonical_transcript_protein | <biomart.Attribute name='mleucophaeus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mleucophaeus_homolog_subtype | <biomart.Attribute name='mleucophaeus_homolog_subtype', display_name='Last common ancestor with Drill', description=''>
mleucophaeus_homolog_orthology_type | <biomart.Attribute name='mleucophaeus_homolog_orthology_type', display_name='Drill homology type', description=''>
mleucophaeus_homolog_perc_id | <biomart.Attribute name='mleucophaeus_homolog_perc_id', display_name='%id. target Drill gene identical to query gene', description=''>
mleucophaeus_homolog_perc_id_r1 | <biomart.Attribute name='mleucophaeus_homolog_perc_id_r1', display_name='%id. query gene identical to target Drill gene', description=''>
mleucophaeus_homolog_goc_score | <biomart.Attribute name='mleucophaeus_homolog_goc_score', display_name='Drill Gene-order conservation score', description=''>
mleucophaeus_homolog_wga_coverage | <biomart.Attribute name='mleucophaeus_homolog_wga_coverage', display_name='Drill Whole-genome alignment coverage', description=''>
mleucophaeus_homolog_orthology_confidence | <biomart.Attribute name='mleucophaeus_homolog_orthology_confidence', display_name='Drill orthology confidence [0 low, 1 high]', description=''>
dmelanogaster_homolog_ensembl_gene | <biomart.Attribute name='dmelanogaster_homolog_ensembl_gene', display_name='Drosophila melanogaster gene stable ID', description=''>
dmelanogaster_homolog_associated_gene_name | <biomart.Attribute name='dmelanogaster_homolog_associated_gene_name', display_name='Drosophila melanogaster gene name', description=''>
dmelanogaster_homolog_ensembl_peptide | <biomart.Attribute name='dmelanogaster_homolog_ensembl_peptide', display_name='Drosophila melanogaster protein or transcript stable ID', description=''>
dmelanogaster_homolog_chromosome | <biomart.Attribute name='dmelanogaster_homolog_chromosome', display_name='Drosophila melanogaster chromosome/scaffold name', description=''>
dmelanogaster_homolog_chrom_start | <biomart.Attribute name='dmelanogaster_homolog_chrom_start', display_name='Drosophila melanogaster chromosome/scaffold start (bp)', description=''>
dmelanogaster_homolog_chrom_end | <biomart.Attribute name='dmelanogaster_homolog_chrom_end', display_name='Drosophila melanogaster chromosome/scaffold end (bp)', description=''>
dmelanogaster_homolog_canonical_transcript_protein | <biomart.Attribute name='dmelanogaster_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
dmelanogaster_homolog_subtype | <biomart.Attribute name='dmelanogaster_homolog_subtype', display_name='Last common ancestor with Drosophila melanogaster', description=''>
dmelanogaster_homolog_orthology_type | <biomart.Attribute name='dmelanogaster_homolog_orthology_type', display_name='Drosophila melanogaster homology type', description=''>
dmelanogaster_homolog_perc_id | <biomart.Attribute name='dmelanogaster_homolog_perc_id', display_name='%id. target Drosophila melanogaster gene identical to query gene', description=''>
dmelanogaster_homolog_perc_id_r1 | <biomart.Attribute name='dmelanogaster_homolog_perc_id_r1', display_name='%id. query gene identical to target Drosophila melanogaster gene', description=''>
dmelanogaster_homolog_orthology_confidence | <biomart.Attribute name='dmelanogaster_homolog_orthology_confidence', display_name='Drosophila melanogaster orthology confidence [0 low, 1 high]', description=''>
applatyrhynchos_homolog_ensembl_gene | <biomart.Attribute name='applatyrhynchos_homolog_ensembl_gene', display_name='Duck gene stable ID', description=''>
applatyrhynchos_homolog_associated_gene_name | <biomart.Attribute name='applatyrhynchos_homolog_associated_gene_name', display_name='Duck gene name', description=''>
applatyrhynchos_homolog_ensembl_peptide | <biomart.Attribute name='applatyrhynchos_homolog_ensembl_peptide', display_name='Duck protein or transcript stable ID', description=''>
applatyrhynchos_homolog_chromosome | <biomart.Attribute name='applatyrhynchos_homolog_chromosome', display_name='Duck chromosome/scaffold name', description=''>
applatyrhynchos_homolog_chrom_start | <biomart.Attribute name='applatyrhynchos_homolog_chrom_start', display_name='Duck chromosome/scaffold start (bp)', description=''>
applatyrhynchos_homolog_chrom_end | <biomart.Attribute name='applatyrhynchos_homolog_chrom_end', display_name='Duck chromosome/scaffold end (bp)', description=''>
applatyrhynchos_homolog_canonical_transcript_protein | <biomart.Attribute name='applatyrhynchos_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
applatyrhynchos_homolog_subtype | <biomart.Attribute name='applatyrhynchos_homolog_subtype', display_name='Last common ancestor with Duck', description=''>
applatyrhynchos_homolog_orthology_type | <biomart.Attribute name='applatyrhynchos_homolog_orthology_type', display_name='Duck homology type', description=''>
applatyrhynchos_homolog_perc_id | <biomart.Attribute name='applatyrhynchos_homolog_perc_id', display_name='%id. target Duck gene identical to query gene', description=''>
applatyrhynchos_homolog_perc_id_r1 | <biomart.Attribute name='applatyrhynchos_homolog_perc_id_r1', display_name='%id. query gene identical to target Duck gene', description=''>
applatyrhynchos_homolog_goc_score | <biomart.Attribute name='applatyrhynchos_homolog_goc_score', display_name='Duck Gene-order conservation score', description=''>
applatyrhynchos_homolog_wga_coverage | <biomart.Attribute name='applatyrhynchos_homolog_wga_coverage', display_name='Duck Whole-genome alignment coverage', description=''>
applatyrhynchos_homolog_orthology_confidence | <biomart.Attribute name='applatyrhynchos_homolog_orthology_confidence', display_name='Duck orthology confidence [0 low, 1 high]', description=''>
acalliptera_homolog_ensembl_gene | <biomart.Attribute name='acalliptera_homolog_ensembl_gene', display_name='Eastern happy gene stable ID', description=''>
acalliptera_homolog_associated_gene_name | <biomart.Attribute name='acalliptera_homolog_associated_gene_name', display_name='Eastern happy gene name', description=''>
acalliptera_homolog_ensembl_peptide | <biomart.Attribute name='acalliptera_homolog_ensembl_peptide', display_name='Eastern happy protein or transcript stable ID', description=''>
acalliptera_homolog_chromosome | <biomart.Attribute name='acalliptera_homolog_chromosome', display_name='Eastern happy chromosome/scaffold name', description=''>
acalliptera_homolog_chrom_start | <biomart.Attribute name='acalliptera_homolog_chrom_start', display_name='Eastern happy chromosome/scaffold start (bp)', description=''>
acalliptera_homolog_chrom_end | <biomart.Attribute name='acalliptera_homolog_chrom_end', display_name='Eastern happy chromosome/scaffold end (bp)', description=''>
acalliptera_homolog_canonical_transcript_protein | <biomart.Attribute name='acalliptera_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
acalliptera_homolog_subtype | <biomart.Attribute name='acalliptera_homolog_subtype', display_name='Last common ancestor with Eastern happy', description=''>
acalliptera_homolog_orthology_type | <biomart.Attribute name='acalliptera_homolog_orthology_type', display_name='Eastern happy homology type', description=''>
acalliptera_homolog_perc_id | <biomart.Attribute name='acalliptera_homolog_perc_id', display_name='%id. target Eastern happy gene identical to query gene', description=''>
acalliptera_homolog_perc_id_r1 | <biomart.Attribute name='acalliptera_homolog_perc_id_r1', display_name='%id. query gene identical to target Eastern happy gene', description=''>
acalliptera_homolog_goc_score | <biomart.Attribute name='acalliptera_homolog_goc_score', display_name='Eastern happy Gene-order conservation score', description=''>
acalliptera_homolog_wga_coverage | <biomart.Attribute name='acalliptera_homolog_wga_coverage', display_name='Eastern happy Whole-genome alignment coverage', description=''>
acalliptera_homolog_orthology_confidence | <biomart.Attribute name='acalliptera_homolog_orthology_confidence', display_name='Eastern happy orthology confidence [0 low, 1 high]', description=''>
eelectricus_homolog_ensembl_gene | <biomart.Attribute name='eelectricus_homolog_ensembl_gene', display_name='Electric eel gene stable ID', description=''>
eelectricus_homolog_associated_gene_name | <biomart.Attribute name='eelectricus_homolog_associated_gene_name', display_name='Electric eel gene name', description=''>
eelectricus_homolog_ensembl_peptide | <biomart.Attribute name='eelectricus_homolog_ensembl_peptide', display_name='Electric eel protein or transcript stable ID', description=''>
eelectricus_homolog_chromosome | <biomart.Attribute name='eelectricus_homolog_chromosome', display_name='Electric eel chromosome/scaffold name', description=''>
eelectricus_homolog_chrom_start | <biomart.Attribute name='eelectricus_homolog_chrom_start', display_name='Electric eel chromosome/scaffold start (bp)', description=''>
eelectricus_homolog_chrom_end | <biomart.Attribute name='eelectricus_homolog_chrom_end', display_name='Electric eel chromosome/scaffold end (bp)', description=''>
eelectricus_homolog_canonical_transcript_protein | <biomart.Attribute name='eelectricus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
eelectricus_homolog_subtype | <biomart.Attribute name='eelectricus_homolog_subtype', display_name='Last common ancestor with Electric eel', description=''>
eelectricus_homolog_orthology_type | <biomart.Attribute name='eelectricus_homolog_orthology_type', display_name='Electric eel homology type', description=''>
eelectricus_homolog_perc_id | <biomart.Attribute name='eelectricus_homolog_perc_id', display_name='%id. target Electric eel gene identical to query gene', description=''>
eelectricus_homolog_perc_id_r1 | <biomart.Attribute name='eelectricus_homolog_perc_id_r1', display_name='%id. query gene identical to target Electric eel gene', description=''>
eelectricus_homolog_goc_score | <biomart.Attribute name='eelectricus_homolog_goc_score', display_name='Electric eel Gene-order conservation score', description=''>
eelectricus_homolog_wga_coverage | <biomart.Attribute name='eelectricus_homolog_wga_coverage', display_name='Electric eel Whole-genome alignment coverage', description=''>
eelectricus_homolog_orthology_confidence | <biomart.Attribute name='eelectricus_homolog_orthology_confidence', display_name='Electric eel orthology confidence [0 low, 1 high]', description=''>
lafricana_homolog_ensembl_gene | <biomart.Attribute name='lafricana_homolog_ensembl_gene', display_name='Elephant gene stable ID', description=''>
lafricana_homolog_associated_gene_name | <biomart.Attribute name='lafricana_homolog_associated_gene_name', display_name='Elephant gene name', description=''>
lafricana_homolog_ensembl_peptide | <biomart.Attribute name='lafricana_homolog_ensembl_peptide', display_name='Elephant protein or transcript stable ID', description=''>
lafricana_homolog_chromosome | <biomart.Attribute name='lafricana_homolog_chromosome', display_name='Elephant chromosome/scaffold name', description=''>
lafricana_homolog_chrom_start | <biomart.Attribute name='lafricana_homolog_chrom_start', display_name='Elephant chromosome/scaffold start (bp)', description=''>
lafricana_homolog_chrom_end | <biomart.Attribute name='lafricana_homolog_chrom_end', display_name='Elephant chromosome/scaffold end (bp)', description=''>
lafricana_homolog_canonical_transcript_protein | <biomart.Attribute name='lafricana_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lafricana_homolog_subtype | <biomart.Attribute name='lafricana_homolog_subtype', display_name='Last common ancestor with Elephant', description=''>
lafricana_homolog_orthology_type | <biomart.Attribute name='lafricana_homolog_orthology_type', display_name='Elephant homology type', description=''>
lafricana_homolog_perc_id | <biomart.Attribute name='lafricana_homolog_perc_id', display_name='%id. target Elephant gene identical to query gene', description=''>
lafricana_homolog_perc_id_r1 | <biomart.Attribute name='lafricana_homolog_perc_id_r1', display_name='%id. query gene identical to target Elephant gene', description=''>
lafricana_homolog_goc_score | <biomart.Attribute name='lafricana_homolog_goc_score', display_name='Elephant Gene-order conservation score', description=''>
lafricana_homolog_wga_coverage | <biomart.Attribute name='lafricana_homolog_wga_coverage', display_name='Elephant Whole-genome alignment coverage', description=''>
lafricana_homolog_orthology_confidence | <biomart.Attribute name='lafricana_homolog_orthology_confidence', display_name='Elephant orthology confidence [0 low, 1 high]', description=''>
cmilii_homolog_ensembl_gene | <biomart.Attribute name='cmilii_homolog_ensembl_gene', display_name='Elephant shark gene stable ID', description=''>
cmilii_homolog_associated_gene_name | <biomart.Attribute name='cmilii_homolog_associated_gene_name', display_name='Elephant shark gene name', description=''>
cmilii_homolog_ensembl_peptide | <biomart.Attribute name='cmilii_homolog_ensembl_peptide', display_name='Elephant shark protein or transcript stable ID', description=''>
cmilii_homolog_chromosome | <biomart.Attribute name='cmilii_homolog_chromosome', display_name='Elephant shark chromosome/scaffold name', description=''>
cmilii_homolog_chrom_start | <biomart.Attribute name='cmilii_homolog_chrom_start', display_name='Elephant shark chromosome/scaffold start (bp)', description=''>
cmilii_homolog_chrom_end | <biomart.Attribute name='cmilii_homolog_chrom_end', display_name='Elephant shark chromosome/scaffold end (bp)', description=''>
cmilii_homolog_canonical_transcript_protein | <biomart.Attribute name='cmilii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cmilii_homolog_subtype | <biomart.Attribute name='cmilii_homolog_subtype', display_name='Last common ancestor with Elephant shark', description=''>
cmilii_homolog_orthology_type | <biomart.Attribute name='cmilii_homolog_orthology_type', display_name='Elephant shark homology type', description=''>
cmilii_homolog_perc_id | <biomart.Attribute name='cmilii_homolog_perc_id', display_name='%id. target Elephant shark gene identical to query gene', description=''>
cmilii_homolog_perc_id_r1 | <biomart.Attribute name='cmilii_homolog_perc_id_r1', display_name='%id. query gene identical to target Elephant shark gene', description=''>
cmilii_homolog_wga_coverage | <biomart.Attribute name='cmilii_homolog_wga_coverage', display_name='Elephant shark Whole-genome alignment coverage', description=''>
cmilii_homolog_orthology_confidence | <biomart.Attribute name='cmilii_homolog_orthology_confidence', display_name='Elephant shark orthology confidence [0 low, 1 high]', description=''>
dnovaehollandiae_homolog_ensembl_gene | <biomart.Attribute name='dnovaehollandiae_homolog_ensembl_gene', display_name='Emu gene stable ID', description=''>
dnovaehollandiae_homolog_associated_gene_name | <biomart.Attribute name='dnovaehollandiae_homolog_associated_gene_name', display_name='Emu gene name', description=''>
dnovaehollandiae_homolog_ensembl_peptide | <biomart.Attribute name='dnovaehollandiae_homolog_ensembl_peptide', display_name='Emu protein or transcript stable ID', description=''>
dnovaehollandiae_homolog_chromosome | <biomart.Attribute name='dnovaehollandiae_homolog_chromosome', display_name='Emu chromosome/scaffold name', description=''>
dnovaehollandiae_homolog_chrom_start | <biomart.Attribute name='dnovaehollandiae_homolog_chrom_start', display_name='Emu chromosome/scaffold start (bp)', description=''>
dnovaehollandiae_homolog_chrom_end | <biomart.Attribute name='dnovaehollandiae_homolog_chrom_end', display_name='Emu chromosome/scaffold end (bp)', description=''>
dnovaehollandiae_homolog_canonical_transcript_protein | <biomart.Attribute name='dnovaehollandiae_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
dnovaehollandiae_homolog_subtype | <biomart.Attribute name='dnovaehollandiae_homolog_subtype', display_name='Last common ancestor with Emu', description=''>
dnovaehollandiae_homolog_orthology_type | <biomart.Attribute name='dnovaehollandiae_homolog_orthology_type', display_name='Emu homology type', description=''>
dnovaehollandiae_homolog_perc_id | <biomart.Attribute name='dnovaehollandiae_homolog_perc_id', display_name='%id. target Emu gene identical to query gene', description=''>
dnovaehollandiae_homolog_perc_id_r1 | <biomart.Attribute name='dnovaehollandiae_homolog_perc_id_r1', display_name='%id. query gene identical to target Emu gene', description=''>
dnovaehollandiae_homolog_goc_score | <biomart.Attribute name='dnovaehollandiae_homolog_goc_score', display_name='Emu Gene-order conservation score', description=''>
dnovaehollandiae_homolog_wga_coverage | <biomart.Attribute name='dnovaehollandiae_homolog_wga_coverage', display_name='Emu Whole-genome alignment coverage', description=''>
dnovaehollandiae_homolog_orthology_confidence | <biomart.Attribute name='dnovaehollandiae_homolog_orthology_confidence', display_name='Emu orthology confidence [0 low, 1 high]', description=''>
mpfuro_homolog_ensembl_gene | <biomart.Attribute name='mpfuro_homolog_ensembl_gene', display_name='Ferret gene stable ID', description=''>
mpfuro_homolog_associated_gene_name | <biomart.Attribute name='mpfuro_homolog_associated_gene_name', display_name='Ferret gene name', description=''>
mpfuro_homolog_ensembl_peptide | <biomart.Attribute name='mpfuro_homolog_ensembl_peptide', display_name='Ferret protein or transcript stable ID', description=''>
mpfuro_homolog_chromosome | <biomart.Attribute name='mpfuro_homolog_chromosome', display_name='Ferret chromosome/scaffold name', description=''>
mpfuro_homolog_chrom_start | <biomart.Attribute name='mpfuro_homolog_chrom_start', display_name='Ferret chromosome/scaffold start (bp)', description=''>
mpfuro_homolog_chrom_end | <biomart.Attribute name='mpfuro_homolog_chrom_end', display_name='Ferret chromosome/scaffold end (bp)', description=''>
mpfuro_homolog_canonical_transcript_protein | <biomart.Attribute name='mpfuro_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mpfuro_homolog_subtype | <biomart.Attribute name='mpfuro_homolog_subtype', display_name='Last common ancestor with Ferret', description=''>
mpfuro_homolog_orthology_type | <biomart.Attribute name='mpfuro_homolog_orthology_type', display_name='Ferret homology type', description=''>
mpfuro_homolog_perc_id | <biomart.Attribute name='mpfuro_homolog_perc_id', display_name='%id. target Ferret gene identical to query gene', description=''>
mpfuro_homolog_perc_id_r1 | <biomart.Attribute name='mpfuro_homolog_perc_id_r1', display_name='%id. query gene identical to target Ferret gene', description=''>
mpfuro_homolog_goc_score | <biomart.Attribute name='mpfuro_homolog_goc_score', display_name='Ferret Gene-order conservation score', description=''>
mpfuro_homolog_wga_coverage | <biomart.Attribute name='mpfuro_homolog_wga_coverage', display_name='Ferret Whole-genome alignment coverage', description=''>
mpfuro_homolog_orthology_confidence | <biomart.Attribute name='mpfuro_homolog_orthology_confidence', display_name='Ferret orthology confidence [0 low, 1 high]', description=''>
falbicollis_homolog_ensembl_gene | <biomart.Attribute name='falbicollis_homolog_ensembl_gene', display_name='Flycatcher gene stable ID', description=''>
falbicollis_homolog_associated_gene_name | <biomart.Attribute name='falbicollis_homolog_associated_gene_name', display_name='Flycatcher gene name', description=''>
falbicollis_homolog_ensembl_peptide | <biomart.Attribute name='falbicollis_homolog_ensembl_peptide', display_name='Flycatcher protein or transcript stable ID', description=''>
falbicollis_homolog_chromosome | <biomart.Attribute name='falbicollis_homolog_chromosome', display_name='Flycatcher chromosome/scaffold name', description=''>
falbicollis_homolog_chrom_start | <biomart.Attribute name='falbicollis_homolog_chrom_start', display_name='Flycatcher chromosome/scaffold start (bp)', description=''>
falbicollis_homolog_chrom_end | <biomart.Attribute name='falbicollis_homolog_chrom_end', display_name='Flycatcher chromosome/scaffold end (bp)', description=''>
falbicollis_homolog_canonical_transcript_protein | <biomart.Attribute name='falbicollis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
falbicollis_homolog_subtype | <biomart.Attribute name='falbicollis_homolog_subtype', display_name='Last common ancestor with Flycatcher', description=''>
falbicollis_homolog_orthology_type | <biomart.Attribute name='falbicollis_homolog_orthology_type', display_name='Flycatcher homology type', description=''>
falbicollis_homolog_perc_id | <biomart.Attribute name='falbicollis_homolog_perc_id', display_name='%id. target Flycatcher gene identical to query gene', description=''>
falbicollis_homolog_perc_id_r1 | <biomart.Attribute name='falbicollis_homolog_perc_id_r1', display_name='%id. query gene identical to target Flycatcher gene', description=''>
falbicollis_homolog_goc_score | <biomart.Attribute name='falbicollis_homolog_goc_score', display_name='Flycatcher Gene-order conservation score', description=''>
falbicollis_homolog_wga_coverage | <biomart.Attribute name='falbicollis_homolog_wga_coverage', display_name='Flycatcher Whole-genome alignment coverage', description=''>
falbicollis_homolog_orthology_confidence | <biomart.Attribute name='falbicollis_homolog_orthology_confidence', display_name='Flycatcher orthology confidence [0 low, 1 high]', description=''>
trubripes_homolog_ensembl_gene | <biomart.Attribute name='trubripes_homolog_ensembl_gene', display_name='Fugu gene stable ID', description=''>
trubripes_homolog_associated_gene_name | <biomart.Attribute name='trubripes_homolog_associated_gene_name', display_name='Fugu gene name', description=''>
trubripes_homolog_ensembl_peptide | <biomart.Attribute name='trubripes_homolog_ensembl_peptide', display_name='Fugu protein or transcript stable ID', description=''>
trubripes_homolog_chromosome | <biomart.Attribute name='trubripes_homolog_chromosome', display_name='Fugu chromosome/scaffold name', description=''>
trubripes_homolog_chrom_start | <biomart.Attribute name='trubripes_homolog_chrom_start', display_name='Fugu chromosome/scaffold start (bp)', description=''>
trubripes_homolog_chrom_end | <biomart.Attribute name='trubripes_homolog_chrom_end', display_name='Fugu chromosome/scaffold end (bp)', description=''>
trubripes_homolog_canonical_transcript_protein | <biomart.Attribute name='trubripes_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
trubripes_homolog_subtype | <biomart.Attribute name='trubripes_homolog_subtype', display_name='Last common ancestor with Fugu', description=''>
trubripes_homolog_orthology_type | <biomart.Attribute name='trubripes_homolog_orthology_type', display_name='Fugu homology type', description=''>
trubripes_homolog_perc_id | <biomart.Attribute name='trubripes_homolog_perc_id', display_name='%id. target Fugu gene identical to query gene', description=''>
trubripes_homolog_perc_id_r1 | <biomart.Attribute name='trubripes_homolog_perc_id_r1', display_name='%id. query gene identical to target Fugu gene', description=''>
trubripes_homolog_goc_score | <biomart.Attribute name='trubripes_homolog_goc_score', display_name='Fugu Gene-order conservation score', description=''>
trubripes_homolog_wga_coverage | <biomart.Attribute name='trubripes_homolog_wga_coverage', display_name='Fugu Whole-genome alignment coverage', description=''>
trubripes_homolog_orthology_confidence | <biomart.Attribute name='trubripes_homolog_orthology_confidence', display_name='Fugu orthology confidence [0 low, 1 high]', description=''>
tgelada_homolog_ensembl_gene | <biomart.Attribute name='tgelada_homolog_ensembl_gene', display_name='Gelada gene stable ID', description=''>
tgelada_homolog_associated_gene_name | <biomart.Attribute name='tgelada_homolog_associated_gene_name', display_name='Gelada gene name', description=''>
tgelada_homolog_ensembl_peptide | <biomart.Attribute name='tgelada_homolog_ensembl_peptide', display_name='Gelada protein or transcript stable ID', description=''>
tgelada_homolog_chromosome | <biomart.Attribute name='tgelada_homolog_chromosome', display_name='Gelada chromosome/scaffold name', description=''>
tgelada_homolog_chrom_start | <biomart.Attribute name='tgelada_homolog_chrom_start', display_name='Gelada chromosome/scaffold start (bp)', description=''>
tgelada_homolog_chrom_end | <biomart.Attribute name='tgelada_homolog_chrom_end', display_name='Gelada chromosome/scaffold end (bp)', description=''>
tgelada_homolog_canonical_transcript_protein | <biomart.Attribute name='tgelada_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
tgelada_homolog_subtype | <biomart.Attribute name='tgelada_homolog_subtype', display_name='Last common ancestor with Gelada', description=''>
tgelada_homolog_orthology_type | <biomart.Attribute name='tgelada_homolog_orthology_type', display_name='Gelada homology type', description=''>
tgelada_homolog_perc_id | <biomart.Attribute name='tgelada_homolog_perc_id', display_name='%id. target Gelada gene identical to query gene', description=''>
tgelada_homolog_perc_id_r1 | <biomart.Attribute name='tgelada_homolog_perc_id_r1', display_name='%id. query gene identical to target Gelada gene', description=''>
tgelada_homolog_goc_score | <biomart.Attribute name='tgelada_homolog_goc_score', display_name='Gelada Gene-order conservation score', description=''>
tgelada_homolog_wga_coverage | <biomart.Attribute name='tgelada_homolog_wga_coverage', display_name='Gelada Whole-genome alignment coverage', description=''>
tgelada_homolog_orthology_confidence | <biomart.Attribute name='tgelada_homolog_orthology_confidence', display_name='Gelada orthology confidence [0 low, 1 high]', description=''>
nleucogenys_homolog_ensembl_gene | <biomart.Attribute name='nleucogenys_homolog_ensembl_gene', display_name='Gibbon gene stable ID', description=''>
nleucogenys_homolog_associated_gene_name | <biomart.Attribute name='nleucogenys_homolog_associated_gene_name', display_name='Gibbon gene name', description=''>
nleucogenys_homolog_ensembl_peptide | <biomart.Attribute name='nleucogenys_homolog_ensembl_peptide', display_name='Gibbon protein or transcript stable ID', description=''>
nleucogenys_homolog_chromosome | <biomart.Attribute name='nleucogenys_homolog_chromosome', display_name='Gibbon chromosome/scaffold name', description=''>
nleucogenys_homolog_chrom_start | <biomart.Attribute name='nleucogenys_homolog_chrom_start', display_name='Gibbon chromosome/scaffold start (bp)', description=''>
nleucogenys_homolog_chrom_end | <biomart.Attribute name='nleucogenys_homolog_chrom_end', display_name='Gibbon chromosome/scaffold end (bp)', description=''>
nleucogenys_homolog_canonical_transcript_protein | <biomart.Attribute name='nleucogenys_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
nleucogenys_homolog_subtype | <biomart.Attribute name='nleucogenys_homolog_subtype', display_name='Last common ancestor with Gibbon', description=''>
nleucogenys_homolog_orthology_type | <biomart.Attribute name='nleucogenys_homolog_orthology_type', display_name='Gibbon homology type', description=''>
nleucogenys_homolog_perc_id | <biomart.Attribute name='nleucogenys_homolog_perc_id', display_name='%id. target Gibbon gene identical to query gene', description=''>
nleucogenys_homolog_perc_id_r1 | <biomart.Attribute name='nleucogenys_homolog_perc_id_r1', display_name='%id. query gene identical to target Gibbon gene', description=''>
nleucogenys_homolog_goc_score | <biomart.Attribute name='nleucogenys_homolog_goc_score', display_name='Gibbon Gene-order conservation score', description=''>
nleucogenys_homolog_wga_coverage | <biomart.Attribute name='nleucogenys_homolog_wga_coverage', display_name='Gibbon Whole-genome alignment coverage', description=''>
nleucogenys_homolog_orthology_confidence | <biomart.Attribute name='nleucogenys_homolog_orthology_confidence', display_name='Gibbon orthology confidence [0 low, 1 high]', description=''>
saurata_homolog_ensembl_gene | <biomart.Attribute name='saurata_homolog_ensembl_gene', display_name='Gilthead seabream gene stable ID', description=''>
saurata_homolog_associated_gene_name | <biomart.Attribute name='saurata_homolog_associated_gene_name', display_name='Gilthead seabream gene name', description=''>
saurata_homolog_ensembl_peptide | <biomart.Attribute name='saurata_homolog_ensembl_peptide', display_name='Gilthead seabream protein or transcript stable ID', description=''>
saurata_homolog_chromosome | <biomart.Attribute name='saurata_homolog_chromosome', display_name='Gilthead seabream chromosome/scaffold name', description=''>
saurata_homolog_chrom_start | <biomart.Attribute name='saurata_homolog_chrom_start', display_name='Gilthead seabream chromosome/scaffold start (bp)', description=''>
saurata_homolog_chrom_end | <biomart.Attribute name='saurata_homolog_chrom_end', display_name='Gilthead seabream chromosome/scaffold end (bp)', description=''>
saurata_homolog_canonical_transcript_protein | <biomart.Attribute name='saurata_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
saurata_homolog_subtype | <biomart.Attribute name='saurata_homolog_subtype', display_name='Last common ancestor with Gilthead seabream', description=''>
saurata_homolog_orthology_type | <biomart.Attribute name='saurata_homolog_orthology_type', display_name='Gilthead seabream homology type', description=''>
saurata_homolog_perc_id | <biomart.Attribute name='saurata_homolog_perc_id', display_name='%id. target Gilthead seabream gene identical to query gene', description=''>
saurata_homolog_perc_id_r1 | <biomart.Attribute name='saurata_homolog_perc_id_r1', display_name='%id. query gene identical to target Gilthead seabream gene', description=''>
saurata_homolog_goc_score | <biomart.Attribute name='saurata_homolog_goc_score', display_name='Gilthead seabream Gene-order conservation score', description=''>
saurata_homolog_wga_coverage | <biomart.Attribute name='saurata_homolog_wga_coverage', display_name='Gilthead seabream Whole-genome alignment coverage', description=''>
saurata_homolog_orthology_confidence | <biomart.Attribute name='saurata_homolog_orthology_confidence', display_name='Gilthead seabream orthology confidence [0 low, 1 high]', description=''>
chircus_homolog_ensembl_gene | <biomart.Attribute name='chircus_homolog_ensembl_gene', display_name='Goat gene stable ID', description=''>
chircus_homolog_associated_gene_name | <biomart.Attribute name='chircus_homolog_associated_gene_name', display_name='Goat gene name', description=''>
chircus_homolog_ensembl_peptide | <biomart.Attribute name='chircus_homolog_ensembl_peptide', display_name='Goat protein or transcript stable ID', description=''>
chircus_homolog_chromosome | <biomart.Attribute name='chircus_homolog_chromosome', display_name='Goat chromosome/scaffold name', description=''>
chircus_homolog_chrom_start | <biomart.Attribute name='chircus_homolog_chrom_start', display_name='Goat chromosome/scaffold start (bp)', description=''>
chircus_homolog_chrom_end | <biomart.Attribute name='chircus_homolog_chrom_end', display_name='Goat chromosome/scaffold end (bp)', description=''>
chircus_homolog_canonical_transcript_protein | <biomart.Attribute name='chircus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
chircus_homolog_subtype | <biomart.Attribute name='chircus_homolog_subtype', display_name='Last common ancestor with Goat', description=''>
chircus_homolog_orthology_type | <biomart.Attribute name='chircus_homolog_orthology_type', display_name='Goat homology type', description=''>
chircus_homolog_perc_id | <biomart.Attribute name='chircus_homolog_perc_id', display_name='%id. target Goat gene identical to query gene', description=''>
chircus_homolog_perc_id_r1 | <biomart.Attribute name='chircus_homolog_perc_id_r1', display_name='%id. query gene identical to target Goat gene', description=''>
chircus_homolog_goc_score | <biomart.Attribute name='chircus_homolog_goc_score', display_name='Goat Gene-order conservation score', description=''>
chircus_homolog_wga_coverage | <biomart.Attribute name='chircus_homolog_wga_coverage', display_name='Goat Whole-genome alignment coverage', description=''>
chircus_homolog_orthology_confidence | <biomart.Attribute name='chircus_homolog_orthology_confidence', display_name='Goat orthology confidence [0 low, 1 high]', description=''>
mauratus_homolog_ensembl_gene | <biomart.Attribute name='mauratus_homolog_ensembl_gene', display_name='Golden Hamster gene stable ID', description=''>
mauratus_homolog_associated_gene_name | <biomart.Attribute name='mauratus_homolog_associated_gene_name', display_name='Golden Hamster gene name', description=''>
mauratus_homolog_ensembl_peptide | <biomart.Attribute name='mauratus_homolog_ensembl_peptide', display_name='Golden Hamster protein or transcript stable ID', description=''>
mauratus_homolog_chromosome | <biomart.Attribute name='mauratus_homolog_chromosome', display_name='Golden Hamster chromosome/scaffold name', description=''>
mauratus_homolog_chrom_start | <biomart.Attribute name='mauratus_homolog_chrom_start', display_name='Golden Hamster chromosome/scaffold start (bp)', description=''>
mauratus_homolog_chrom_end | <biomart.Attribute name='mauratus_homolog_chrom_end', display_name='Golden Hamster chromosome/scaffold end (bp)', description=''>
mauratus_homolog_canonical_transcript_protein | <biomart.Attribute name='mauratus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mauratus_homolog_subtype | <biomart.Attribute name='mauratus_homolog_subtype', display_name='Last common ancestor with Golden Hamster', description=''>
mauratus_homolog_orthology_type | <biomart.Attribute name='mauratus_homolog_orthology_type', display_name='Golden Hamster homology type', description=''>
mauratus_homolog_perc_id | <biomart.Attribute name='mauratus_homolog_perc_id', display_name='%id. target Golden Hamster gene identical to query gene', description=''>
mauratus_homolog_perc_id_r1 | <biomart.Attribute name='mauratus_homolog_perc_id_r1', display_name='%id. query gene identical to target Golden Hamster gene', description=''>
mauratus_homolog_goc_score | <biomart.Attribute name='mauratus_homolog_goc_score', display_name='Golden Hamster Gene-order conservation score', description=''>
mauratus_homolog_wga_coverage | <biomart.Attribute name='mauratus_homolog_wga_coverage', display_name='Golden Hamster Whole-genome alignment coverage', description=''>
mauratus_homolog_orthology_confidence | <biomart.Attribute name='mauratus_homolog_orthology_confidence', display_name='Golden Hamster orthology confidence [0 low, 1 high]', description=''>
acchrysaetos_homolog_ensembl_gene | <biomart.Attribute name='acchrysaetos_homolog_ensembl_gene', display_name='Golden eagle gene stable ID', description=''>
acchrysaetos_homolog_associated_gene_name | <biomart.Attribute name='acchrysaetos_homolog_associated_gene_name', display_name='Golden eagle gene name', description=''>
acchrysaetos_homolog_ensembl_peptide | <biomart.Attribute name='acchrysaetos_homolog_ensembl_peptide', display_name='Golden eagle protein or transcript stable ID', description=''>
acchrysaetos_homolog_chromosome | <biomart.Attribute name='acchrysaetos_homolog_chromosome', display_name='Golden eagle chromosome/scaffold name', description=''>
acchrysaetos_homolog_chrom_start | <biomart.Attribute name='acchrysaetos_homolog_chrom_start', display_name='Golden eagle chromosome/scaffold start (bp)', description=''>
acchrysaetos_homolog_chrom_end | <biomart.Attribute name='acchrysaetos_homolog_chrom_end', display_name='Golden eagle chromosome/scaffold end (bp)', description=''>
acchrysaetos_homolog_canonical_transcript_protein | <biomart.Attribute name='acchrysaetos_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
acchrysaetos_homolog_subtype | <biomart.Attribute name='acchrysaetos_homolog_subtype', display_name='Last common ancestor with Golden eagle', description=''>
acchrysaetos_homolog_orthology_type | <biomart.Attribute name='acchrysaetos_homolog_orthology_type', display_name='Golden eagle homology type', description=''>
acchrysaetos_homolog_perc_id | <biomart.Attribute name='acchrysaetos_homolog_perc_id', display_name='%id. target Golden eagle gene identical to query gene', description=''>
acchrysaetos_homolog_perc_id_r1 | <biomart.Attribute name='acchrysaetos_homolog_perc_id_r1', display_name='%id. query gene identical to target Golden eagle gene', description=''>
acchrysaetos_homolog_goc_score | <biomart.Attribute name='acchrysaetos_homolog_goc_score', display_name='Golden eagle Gene-order conservation score', description=''>
acchrysaetos_homolog_wga_coverage | <biomart.Attribute name='acchrysaetos_homolog_wga_coverage', display_name='Golden eagle Whole-genome alignment coverage', description=''>
acchrysaetos_homolog_orthology_confidence | <biomart.Attribute name='acchrysaetos_homolog_orthology_confidence', display_name='Golden eagle orthology confidence [0 low, 1 high]', description=''>
rroxellana_homolog_ensembl_gene | <biomart.Attribute name='rroxellana_homolog_ensembl_gene', display_name='Golden snub-nosed monkey gene stable ID', description=''>
rroxellana_homolog_associated_gene_name | <biomart.Attribute name='rroxellana_homolog_associated_gene_name', display_name='Golden snub-nosed monkey gene name', description=''>
rroxellana_homolog_ensembl_peptide | <biomart.Attribute name='rroxellana_homolog_ensembl_peptide', display_name='Golden snub-nosed monkey protein or transcript stable ID', description=''>
rroxellana_homolog_chromosome | <biomart.Attribute name='rroxellana_homolog_chromosome', display_name='Golden snub-nosed monkey chromosome/scaffold name', description=''>
rroxellana_homolog_chrom_start | <biomart.Attribute name='rroxellana_homolog_chrom_start', display_name='Golden snub-nosed monkey chromosome/scaffold start (bp)', description=''>
rroxellana_homolog_chrom_end | <biomart.Attribute name='rroxellana_homolog_chrom_end', display_name='Golden snub-nosed monkey chromosome/scaffold end (bp)', description=''>
rroxellana_homolog_canonical_transcript_protein | <biomart.Attribute name='rroxellana_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
rroxellana_homolog_subtype | <biomart.Attribute name='rroxellana_homolog_subtype', display_name='Last common ancestor with Golden snub-nosed monkey', description=''>
rroxellana_homolog_orthology_type | <biomart.Attribute name='rroxellana_homolog_orthology_type', display_name='Golden snub-nosed monkey homology type', description=''>
rroxellana_homolog_perc_id | <biomart.Attribute name='rroxellana_homolog_perc_id', display_name='%id. target Golden snub-nosed monkey gene identical to query gene', description=''>
rroxellana_homolog_perc_id_r1 | <biomart.Attribute name='rroxellana_homolog_perc_id_r1', display_name='%id. query gene identical to target Golden snub-nosed monkey gene', description=''>
rroxellana_homolog_goc_score | <biomart.Attribute name='rroxellana_homolog_goc_score', display_name='Golden snub-nosed monkey Gene-order conservation score', description=''>
rroxellana_homolog_wga_coverage | <biomart.Attribute name='rroxellana_homolog_wga_coverage', display_name='Golden snub-nosed monkey Whole-genome alignment coverage', description=''>
rroxellana_homolog_orthology_confidence | <biomart.Attribute name='rroxellana_homolog_orthology_confidence', display_name='Golden snub-nosed monkey orthology confidence [0 low, 1 high]', description=''>
mvitellinus_homolog_ensembl_gene | <biomart.Attribute name='mvitellinus_homolog_ensembl_gene', display_name='Golden-collared manakin gene stable ID', description=''>
mvitellinus_homolog_associated_gene_name | <biomart.Attribute name='mvitellinus_homolog_associated_gene_name', display_name='Golden-collared manakin gene name', description=''>
mvitellinus_homolog_ensembl_peptide | <biomart.Attribute name='mvitellinus_homolog_ensembl_peptide', display_name='Golden-collared manakin protein or transcript stable ID', description=''>
mvitellinus_homolog_chromosome | <biomart.Attribute name='mvitellinus_homolog_chromosome', display_name='Golden-collared manakin chromosome/scaffold name', description=''>
mvitellinus_homolog_chrom_start | <biomart.Attribute name='mvitellinus_homolog_chrom_start', display_name='Golden-collared manakin chromosome/scaffold start (bp)', description=''>
mvitellinus_homolog_chrom_end | <biomart.Attribute name='mvitellinus_homolog_chrom_end', display_name='Golden-collared manakin chromosome/scaffold end (bp)', description=''>
mvitellinus_homolog_canonical_transcript_protein | <biomart.Attribute name='mvitellinus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mvitellinus_homolog_subtype | <biomart.Attribute name='mvitellinus_homolog_subtype', display_name='Last common ancestor with Golden-collared manakin', description=''>
mvitellinus_homolog_orthology_type | <biomart.Attribute name='mvitellinus_homolog_orthology_type', display_name='Golden-collared manakin homology type', description=''>
mvitellinus_homolog_perc_id | <biomart.Attribute name='mvitellinus_homolog_perc_id', display_name='%id. target Golden-collared manakin gene identical to query gene', description=''>
mvitellinus_homolog_perc_id_r1 | <biomart.Attribute name='mvitellinus_homolog_perc_id_r1', display_name='%id. query gene identical to target Golden-collared manakin gene', description=''>
mvitellinus_homolog_goc_score | <biomart.Attribute name='mvitellinus_homolog_goc_score', display_name='Golden-collared manakin Gene-order conservation score', description=''>
mvitellinus_homolog_wga_coverage | <biomart.Attribute name='mvitellinus_homolog_wga_coverage', display_name='Golden-collared manakin Whole-genome alignment coverage', description=''>
mvitellinus_homolog_orthology_confidence | <biomart.Attribute name='mvitellinus_homolog_orthology_confidence', display_name='Golden-collared manakin orthology confidence [0 low, 1 high]', description=''>
ggorilla_homolog_ensembl_gene | <biomart.Attribute name='ggorilla_homolog_ensembl_gene', display_name='Gorilla gene stable ID', description=''>
ggorilla_homolog_associated_gene_name | <biomart.Attribute name='ggorilla_homolog_associated_gene_name', display_name='Gorilla gene name', description=''>
ggorilla_homolog_ensembl_peptide | <biomart.Attribute name='ggorilla_homolog_ensembl_peptide', display_name='Gorilla protein or transcript stable ID', description=''>
ggorilla_homolog_chromosome | <biomart.Attribute name='ggorilla_homolog_chromosome', display_name='Gorilla chromosome/scaffold name', description=''>
ggorilla_homolog_chrom_start | <biomart.Attribute name='ggorilla_homolog_chrom_start', display_name='Gorilla chromosome/scaffold start (bp)', description=''>
ggorilla_homolog_chrom_end | <biomart.Attribute name='ggorilla_homolog_chrom_end', display_name='Gorilla chromosome/scaffold end (bp)', description=''>
ggorilla_homolog_canonical_transcript_protein | <biomart.Attribute name='ggorilla_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ggorilla_homolog_subtype | <biomart.Attribute name='ggorilla_homolog_subtype', display_name='Last common ancestor with Gorilla', description=''>
ggorilla_homolog_orthology_type | <biomart.Attribute name='ggorilla_homolog_orthology_type', display_name='Gorilla homology type', description=''>
ggorilla_homolog_perc_id | <biomart.Attribute name='ggorilla_homolog_perc_id', display_name='%id. target Gorilla gene identical to query gene', description=''>
ggorilla_homolog_perc_id_r1 | <biomart.Attribute name='ggorilla_homolog_perc_id_r1', display_name='%id. query gene identical to target Gorilla gene', description=''>
ggorilla_homolog_goc_score | <biomart.Attribute name='ggorilla_homolog_goc_score', display_name='Gorilla Gene-order conservation score', description=''>
ggorilla_homolog_wga_coverage | <biomart.Attribute name='ggorilla_homolog_wga_coverage', display_name='Gorilla Whole-genome alignment coverage', description=''>
ggorilla_homolog_orthology_confidence | <biomart.Attribute name='ggorilla_homolog_orthology_confidence', display_name='Gorilla orthology confidence [0 low, 1 high]', description=''>
sdumerili_homolog_ensembl_gene | <biomart.Attribute name='sdumerili_homolog_ensembl_gene', display_name='Greater amberjack gene stable ID', description=''>
sdumerili_homolog_associated_gene_name | <biomart.Attribute name='sdumerili_homolog_associated_gene_name', display_name='Greater amberjack gene name', description=''>
sdumerili_homolog_ensembl_peptide | <biomart.Attribute name='sdumerili_homolog_ensembl_peptide', display_name='Greater amberjack protein or transcript stable ID', description=''>
sdumerili_homolog_chromosome | <biomart.Attribute name='sdumerili_homolog_chromosome', display_name='Greater amberjack chromosome/scaffold name', description=''>
sdumerili_homolog_chrom_start | <biomart.Attribute name='sdumerili_homolog_chrom_start', display_name='Greater amberjack chromosome/scaffold start (bp)', description=''>
sdumerili_homolog_chrom_end | <biomart.Attribute name='sdumerili_homolog_chrom_end', display_name='Greater amberjack chromosome/scaffold end (bp)', description=''>
sdumerili_homolog_canonical_transcript_protein | <biomart.Attribute name='sdumerili_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sdumerili_homolog_subtype | <biomart.Attribute name='sdumerili_homolog_subtype', display_name='Last common ancestor with Greater amberjack', description=''>
sdumerili_homolog_orthology_type | <biomart.Attribute name='sdumerili_homolog_orthology_type', display_name='Greater amberjack homology type', description=''>
sdumerili_homolog_perc_id | <biomart.Attribute name='sdumerili_homolog_perc_id', display_name='%id. target Greater amberjack gene identical to query gene', description=''>
sdumerili_homolog_perc_id_r1 | <biomart.Attribute name='sdumerili_homolog_perc_id_r1', display_name='%id. query gene identical to target Greater amberjack gene', description=''>
sdumerili_homolog_goc_score | <biomart.Attribute name='sdumerili_homolog_goc_score', display_name='Greater amberjack Gene-order conservation score', description=''>
sdumerili_homolog_wga_coverage | <biomart.Attribute name='sdumerili_homolog_wga_coverage', display_name='Greater amberjack Whole-genome alignment coverage', description=''>
sdumerili_homolog_orthology_confidence | <biomart.Attribute name='sdumerili_homolog_orthology_confidence', display_name='Greater amberjack orthology confidence [0 low, 1 high]', description=''>
psimus_homolog_ensembl_gene | <biomart.Attribute name='psimus_homolog_ensembl_gene', display_name='Greater bamboo lemur gene stable ID', description=''>
psimus_homolog_associated_gene_name | <biomart.Attribute name='psimus_homolog_associated_gene_name', display_name='Greater bamboo lemur gene name', description=''>
psimus_homolog_ensembl_peptide | <biomart.Attribute name='psimus_homolog_ensembl_peptide', display_name='Greater bamboo lemur protein or transcript stable ID', description=''>
psimus_homolog_chromosome | <biomart.Attribute name='psimus_homolog_chromosome', display_name='Greater bamboo lemur chromosome/scaffold name', description=''>
psimus_homolog_chrom_start | <biomart.Attribute name='psimus_homolog_chrom_start', display_name='Greater bamboo lemur chromosome/scaffold start (bp)', description=''>
psimus_homolog_chrom_end | <biomart.Attribute name='psimus_homolog_chrom_end', display_name='Greater bamboo lemur chromosome/scaffold end (bp)', description=''>
psimus_homolog_canonical_transcript_protein | <biomart.Attribute name='psimus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
psimus_homolog_subtype | <biomart.Attribute name='psimus_homolog_subtype', display_name='Last common ancestor with Greater bamboo lemur', description=''>
psimus_homolog_orthology_type | <biomart.Attribute name='psimus_homolog_orthology_type', display_name='Greater bamboo lemur homology type', description=''>
psimus_homolog_perc_id | <biomart.Attribute name='psimus_homolog_perc_id', display_name='%id. target Greater bamboo lemur gene identical to query gene', description=''>
psimus_homolog_perc_id_r1 | <biomart.Attribute name='psimus_homolog_perc_id_r1', display_name='%id. query gene identical to target Greater bamboo lemur gene', description=''>
psimus_homolog_goc_score | <biomart.Attribute name='psimus_homolog_goc_score', display_name='Greater bamboo lemur Gene-order conservation score', description=''>
psimus_homolog_wga_coverage | <biomart.Attribute name='psimus_homolog_wga_coverage', display_name='Greater bamboo lemur Whole-genome alignment coverage', description=''>
psimus_homolog_orthology_confidence | <biomart.Attribute name='psimus_homolog_orthology_confidence', display_name='Greater bamboo lemur orthology confidence [0 low, 1 high]', description=''>
rferrumequinum_homolog_ensembl_gene | <biomart.Attribute name='rferrumequinum_homolog_ensembl_gene', display_name='Greater horseshoe bat gene stable ID', description=''>
rferrumequinum_homolog_associated_gene_name | <biomart.Attribute name='rferrumequinum_homolog_associated_gene_name', display_name='Greater horseshoe bat gene name', description=''>
rferrumequinum_homolog_ensembl_peptide | <biomart.Attribute name='rferrumequinum_homolog_ensembl_peptide', display_name='Greater horseshoe bat protein or transcript stable ID', description=''>
rferrumequinum_homolog_chromosome | <biomart.Attribute name='rferrumequinum_homolog_chromosome', display_name='Greater horseshoe bat chromosome/scaffold name', description=''>
rferrumequinum_homolog_chrom_start | <biomart.Attribute name='rferrumequinum_homolog_chrom_start', display_name='Greater horseshoe bat chromosome/scaffold start (bp)', description=''>
rferrumequinum_homolog_chrom_end | <biomart.Attribute name='rferrumequinum_homolog_chrom_end', display_name='Greater horseshoe bat chromosome/scaffold end (bp)', description=''>
rferrumequinum_homolog_canonical_transcript_protein | <biomart.Attribute name='rferrumequinum_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
rferrumequinum_homolog_subtype | <biomart.Attribute name='rferrumequinum_homolog_subtype', display_name='Last common ancestor with Greater horseshoe bat', description=''>
rferrumequinum_homolog_orthology_type | <biomart.Attribute name='rferrumequinum_homolog_orthology_type', display_name='Greater horseshoe bat homology type', description=''>
rferrumequinum_homolog_perc_id | <biomart.Attribute name='rferrumequinum_homolog_perc_id', display_name='%id. target Greater horseshoe bat gene identical to query gene', description=''>
rferrumequinum_homolog_perc_id_r1 | <biomart.Attribute name='rferrumequinum_homolog_perc_id_r1', display_name='%id. query gene identical to target Greater horseshoe bat gene', description=''>
rferrumequinum_homolog_goc_score | <biomart.Attribute name='rferrumequinum_homolog_goc_score', display_name='Greater horseshoe bat Gene-order conservation score', description=''>
rferrumequinum_homolog_wga_coverage | <biomart.Attribute name='rferrumequinum_homolog_wga_coverage', display_name='Greater horseshoe bat Whole-genome alignment coverage', description=''>
rferrumequinum_homolog_orthology_confidence | <biomart.Attribute name='rferrumequinum_homolog_orthology_confidence', display_name='Greater horseshoe bat orthology confidence [0 low, 1 high]', description=''>
cporcellus_homolog_ensembl_gene | <biomart.Attribute name='cporcellus_homolog_ensembl_gene', display_name='Guinea Pig gene stable ID', description=''>
cporcellus_homolog_associated_gene_name | <biomart.Attribute name='cporcellus_homolog_associated_gene_name', display_name='Guinea Pig gene name', description=''>
cporcellus_homolog_ensembl_peptide | <biomart.Attribute name='cporcellus_homolog_ensembl_peptide', display_name='Guinea Pig protein or transcript stable ID', description=''>
cporcellus_homolog_chromosome | <biomart.Attribute name='cporcellus_homolog_chromosome', display_name='Guinea Pig chromosome/scaffold name', description=''>
cporcellus_homolog_chrom_start | <biomart.Attribute name='cporcellus_homolog_chrom_start', display_name='Guinea Pig chromosome/scaffold start (bp)', description=''>
cporcellus_homolog_chrom_end | <biomart.Attribute name='cporcellus_homolog_chrom_end', display_name='Guinea Pig chromosome/scaffold end (bp)', description=''>
cporcellus_homolog_canonical_transcript_protein | <biomart.Attribute name='cporcellus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cporcellus_homolog_subtype | <biomart.Attribute name='cporcellus_homolog_subtype', display_name='Last common ancestor with Guinea Pig', description=''>
cporcellus_homolog_orthology_type | <biomart.Attribute name='cporcellus_homolog_orthology_type', display_name='Guinea Pig homology type', description=''>
cporcellus_homolog_perc_id | <biomart.Attribute name='cporcellus_homolog_perc_id', display_name='%id. target Guinea Pig gene identical to query gene', description=''>
cporcellus_homolog_perc_id_r1 | <biomart.Attribute name='cporcellus_homolog_perc_id_r1', display_name='%id. query gene identical to target Guinea Pig gene', description=''>
cporcellus_homolog_goc_score | <biomart.Attribute name='cporcellus_homolog_goc_score', display_name='Guinea Pig Gene-order conservation score', description=''>
cporcellus_homolog_wga_coverage | <biomart.Attribute name='cporcellus_homolog_wga_coverage', display_name='Guinea Pig Whole-genome alignment coverage', description=''>
cporcellus_homolog_orthology_confidence | <biomart.Attribute name='cporcellus_homolog_orthology_confidence', display_name='Guinea Pig orthology confidence [0 low, 1 high]', description=''>
preticulata_homolog_ensembl_gene | <biomart.Attribute name='preticulata_homolog_ensembl_gene', display_name='Guppy gene stable ID', description=''>
preticulata_homolog_associated_gene_name | <biomart.Attribute name='preticulata_homolog_associated_gene_name', display_name='Guppy gene name', description=''>
preticulata_homolog_ensembl_peptide | <biomart.Attribute name='preticulata_homolog_ensembl_peptide', display_name='Guppy protein or transcript stable ID', description=''>
preticulata_homolog_chromosome | <biomart.Attribute name='preticulata_homolog_chromosome', display_name='Guppy chromosome/scaffold name', description=''>
preticulata_homolog_chrom_start | <biomart.Attribute name='preticulata_homolog_chrom_start', display_name='Guppy chromosome/scaffold start (bp)', description=''>
preticulata_homolog_chrom_end | <biomart.Attribute name='preticulata_homolog_chrom_end', display_name='Guppy chromosome/scaffold end (bp)', description=''>
preticulata_homolog_canonical_transcript_protein | <biomart.Attribute name='preticulata_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
preticulata_homolog_subtype | <biomart.Attribute name='preticulata_homolog_subtype', display_name='Last common ancestor with Guppy', description=''>
preticulata_homolog_orthology_type | <biomart.Attribute name='preticulata_homolog_orthology_type', display_name='Guppy homology type', description=''>
preticulata_homolog_perc_id | <biomart.Attribute name='preticulata_homolog_perc_id', display_name='%id. target Guppy gene identical to query gene', description=''>
preticulata_homolog_perc_id_r1 | <biomart.Attribute name='preticulata_homolog_perc_id_r1', display_name='%id. query gene identical to target Guppy gene', description=''>
preticulata_homolog_goc_score | <biomart.Attribute name='preticulata_homolog_goc_score', display_name='Guppy Gene-order conservation score', description=''>
preticulata_homolog_wga_coverage | <biomart.Attribute name='preticulata_homolog_wga_coverage', display_name='Guppy Whole-genome alignment coverage', description=''>
preticulata_homolog_orthology_confidence | <biomart.Attribute name='preticulata_homolog_orthology_confidence', display_name='Guppy orthology confidence [0 low, 1 high]', description=''>
eburgeri_homolog_ensembl_gene | <biomart.Attribute name='eburgeri_homolog_ensembl_gene', display_name='Hagfish gene stable ID', description=''>
eburgeri_homolog_associated_gene_name | <biomart.Attribute name='eburgeri_homolog_associated_gene_name', display_name='Hagfish gene name', description=''>
eburgeri_homolog_ensembl_peptide | <biomart.Attribute name='eburgeri_homolog_ensembl_peptide', display_name='Hagfish protein or transcript stable ID', description=''>
eburgeri_homolog_chromosome | <biomart.Attribute name='eburgeri_homolog_chromosome', display_name='Hagfish chromosome/scaffold name', description=''>
eburgeri_homolog_chrom_start | <biomart.Attribute name='eburgeri_homolog_chrom_start', display_name='Hagfish chromosome/scaffold start (bp)', description=''>
eburgeri_homolog_chrom_end | <biomart.Attribute name='eburgeri_homolog_chrom_end', display_name='Hagfish chromosome/scaffold end (bp)', description=''>
eburgeri_homolog_canonical_transcript_protein | <biomart.Attribute name='eburgeri_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
eburgeri_homolog_subtype | <biomart.Attribute name='eburgeri_homolog_subtype', display_name='Last common ancestor with Hagfish', description=''>
eburgeri_homolog_orthology_type | <biomart.Attribute name='eburgeri_homolog_orthology_type', display_name='Hagfish homology type', description=''>
eburgeri_homolog_perc_id | <biomart.Attribute name='eburgeri_homolog_perc_id', display_name='%id. target Hagfish gene identical to query gene', description=''>
eburgeri_homolog_perc_id_r1 | <biomart.Attribute name='eburgeri_homolog_perc_id_r1', display_name='%id. query gene identical to target Hagfish gene', description=''>
eburgeri_homolog_wga_coverage | <biomart.Attribute name='eburgeri_homolog_wga_coverage', display_name='Hagfish Whole-genome alignment coverage', description=''>
eburgeri_homolog_orthology_confidence | <biomart.Attribute name='eburgeri_homolog_orthology_confidence', display_name='Hagfish orthology confidence [0 low, 1 high]', description=''>
eeuropaeus_homolog_ensembl_gene | <biomart.Attribute name='eeuropaeus_homolog_ensembl_gene', display_name='Hedgehog gene stable ID', description=''>
eeuropaeus_homolog_associated_gene_name | <biomart.Attribute name='eeuropaeus_homolog_associated_gene_name', display_name='Hedgehog gene name', description=''>
eeuropaeus_homolog_ensembl_peptide | <biomart.Attribute name='eeuropaeus_homolog_ensembl_peptide', display_name='Hedgehog protein or transcript stable ID', description=''>
eeuropaeus_homolog_chromosome | <biomart.Attribute name='eeuropaeus_homolog_chromosome', display_name='Hedgehog chromosome/scaffold name', description=''>
eeuropaeus_homolog_chrom_start | <biomart.Attribute name='eeuropaeus_homolog_chrom_start', display_name='Hedgehog chromosome/scaffold start (bp)', description=''>
eeuropaeus_homolog_chrom_end | <biomart.Attribute name='eeuropaeus_homolog_chrom_end', display_name='Hedgehog chromosome/scaffold end (bp)', description=''>
eeuropaeus_homolog_canonical_transcript_protein | <biomart.Attribute name='eeuropaeus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
eeuropaeus_homolog_subtype | <biomart.Attribute name='eeuropaeus_homolog_subtype', display_name='Last common ancestor with Hedgehog', description=''>
eeuropaeus_homolog_orthology_type | <biomart.Attribute name='eeuropaeus_homolog_orthology_type', display_name='Hedgehog homology type', description=''>
eeuropaeus_homolog_perc_id | <biomart.Attribute name='eeuropaeus_homolog_perc_id', display_name='%id. target Hedgehog gene identical to query gene', description=''>
eeuropaeus_homolog_perc_id_r1 | <biomart.Attribute name='eeuropaeus_homolog_perc_id_r1', display_name='%id. query gene identical to target Hedgehog gene', description=''>
eeuropaeus_homolog_goc_score | <biomart.Attribute name='eeuropaeus_homolog_goc_score', display_name='Hedgehog Gene-order conservation score', description=''>
eeuropaeus_homolog_wga_coverage | <biomart.Attribute name='eeuropaeus_homolog_wga_coverage', display_name='Hedgehog Whole-genome alignment coverage', description=''>
eeuropaeus_homolog_orthology_confidence | <biomart.Attribute name='eeuropaeus_homolog_orthology_confidence', display_name='Hedgehog orthology confidence [0 low, 1 high]', description=''>
ecaballus_homolog_ensembl_gene | <biomart.Attribute name='ecaballus_homolog_ensembl_gene', display_name='Horse gene stable ID', description=''>
ecaballus_homolog_associated_gene_name | <biomart.Attribute name='ecaballus_homolog_associated_gene_name', display_name='Horse gene name', description=''>
ecaballus_homolog_ensembl_peptide | <biomart.Attribute name='ecaballus_homolog_ensembl_peptide', display_name='Horse protein or transcript stable ID', description=''>
ecaballus_homolog_chromosome | <biomart.Attribute name='ecaballus_homolog_chromosome', display_name='Horse chromosome/scaffold name', description=''>
ecaballus_homolog_chrom_start | <biomart.Attribute name='ecaballus_homolog_chrom_start', display_name='Horse chromosome/scaffold start (bp)', description=''>
ecaballus_homolog_chrom_end | <biomart.Attribute name='ecaballus_homolog_chrom_end', display_name='Horse chromosome/scaffold end (bp)', description=''>
ecaballus_homolog_canonical_transcript_protein | <biomart.Attribute name='ecaballus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ecaballus_homolog_subtype | <biomart.Attribute name='ecaballus_homolog_subtype', display_name='Last common ancestor with Horse', description=''>
ecaballus_homolog_orthology_type | <biomart.Attribute name='ecaballus_homolog_orthology_type', display_name='Horse homology type', description=''>
ecaballus_homolog_perc_id | <biomart.Attribute name='ecaballus_homolog_perc_id', display_name='%id. target Horse gene identical to query gene', description=''>
ecaballus_homolog_perc_id_r1 | <biomart.Attribute name='ecaballus_homolog_perc_id_r1', display_name='%id. query gene identical to target Horse gene', description=''>
ecaballus_homolog_goc_score | <biomart.Attribute name='ecaballus_homolog_goc_score', display_name='Horse Gene-order conservation score', description=''>
ecaballus_homolog_wga_coverage | <biomart.Attribute name='ecaballus_homolog_wga_coverage', display_name='Horse Whole-genome alignment coverage', description=''>
ecaballus_homolog_orthology_confidence | <biomart.Attribute name='ecaballus_homolog_orthology_confidence', display_name='Horse orthology confidence [0 low, 1 high]', description=''>
hhucho_homolog_ensembl_gene | <biomart.Attribute name='hhucho_homolog_ensembl_gene', display_name='Huchen gene stable ID', description=''>
hhucho_homolog_associated_gene_name | <biomart.Attribute name='hhucho_homolog_associated_gene_name', display_name='Huchen gene name', description=''>
hhucho_homolog_ensembl_peptide | <biomart.Attribute name='hhucho_homolog_ensembl_peptide', display_name='Huchen protein or transcript stable ID', description=''>
hhucho_homolog_chromosome | <biomart.Attribute name='hhucho_homolog_chromosome', display_name='Huchen chromosome/scaffold name', description=''>
hhucho_homolog_chrom_start | <biomart.Attribute name='hhucho_homolog_chrom_start', display_name='Huchen chromosome/scaffold start (bp)', description=''>
hhucho_homolog_chrom_end | <biomart.Attribute name='hhucho_homolog_chrom_end', display_name='Huchen chromosome/scaffold end (bp)', description=''>
hhucho_homolog_canonical_transcript_protein | <biomart.Attribute name='hhucho_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
hhucho_homolog_subtype | <biomart.Attribute name='hhucho_homolog_subtype', display_name='Last common ancestor with Huchen', description=''>
hhucho_homolog_orthology_type | <biomart.Attribute name='hhucho_homolog_orthology_type', display_name='Huchen homology type', description=''>
hhucho_homolog_perc_id | <biomart.Attribute name='hhucho_homolog_perc_id', display_name='%id. target Huchen gene identical to query gene', description=''>
hhucho_homolog_perc_id_r1 | <biomart.Attribute name='hhucho_homolog_perc_id_r1', display_name='%id. query gene identical to target Huchen gene', description=''>
hhucho_homolog_goc_score | <biomart.Attribute name='hhucho_homolog_goc_score', display_name='Huchen Gene-order conservation score', description=''>
hhucho_homolog_wga_coverage | <biomart.Attribute name='hhucho_homolog_wga_coverage', display_name='Huchen Whole-genome alignment coverage', description=''>
hhucho_homolog_orthology_confidence | <biomart.Attribute name='hhucho_homolog_orthology_confidence', display_name='Huchen orthology confidence [0 low, 1 high]', description=''>
bihybrid_homolog_ensembl_gene | <biomart.Attribute name='bihybrid_homolog_ensembl_gene', display_name='Hybrid - Bos Indicus gene stable ID', description=''>
bihybrid_homolog_associated_gene_name | <biomart.Attribute name='bihybrid_homolog_associated_gene_name', display_name='Hybrid - Bos Indicus gene name', description=''>
bihybrid_homolog_ensembl_peptide | <biomart.Attribute name='bihybrid_homolog_ensembl_peptide', display_name='Hybrid - Bos Indicus protein or transcript stable ID', description=''>
bihybrid_homolog_chromosome | <biomart.Attribute name='bihybrid_homolog_chromosome', display_name='Hybrid - Bos Indicus chromosome/scaffold name', description=''>
bihybrid_homolog_chrom_start | <biomart.Attribute name='bihybrid_homolog_chrom_start', display_name='Hybrid - Bos Indicus chromosome/scaffold start (bp)', description=''>
bihybrid_homolog_chrom_end | <biomart.Attribute name='bihybrid_homolog_chrom_end', display_name='Hybrid - Bos Indicus chromosome/scaffold end (bp)', description=''>
bihybrid_homolog_canonical_transcript_protein | <biomart.Attribute name='bihybrid_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
bihybrid_homolog_subtype | <biomart.Attribute name='bihybrid_homolog_subtype', display_name='Last common ancestor with Hybrid - Bos Indicus', description=''>
bihybrid_homolog_orthology_type | <biomart.Attribute name='bihybrid_homolog_orthology_type', display_name='Hybrid - Bos Indicus homology type', description=''>
bihybrid_homolog_perc_id | <biomart.Attribute name='bihybrid_homolog_perc_id', display_name='%id. target Hybrid - Bos Indicus gene identical to query gene', description=''>
bihybrid_homolog_perc_id_r1 | <biomart.Attribute name='bihybrid_homolog_perc_id_r1', display_name='%id. query gene identical to target Hybrid - Bos Indicus gene', description=''>
bihybrid_homolog_goc_score | <biomart.Attribute name='bihybrid_homolog_goc_score', display_name='Hybrid - Bos Indicus Gene-order conservation score', description=''>
bihybrid_homolog_wga_coverage | <biomart.Attribute name='bihybrid_homolog_wga_coverage', display_name='Hybrid - Bos Indicus Whole-genome alignment coverage', description=''>
bihybrid_homolog_orthology_confidence | <biomart.Attribute name='bihybrid_homolog_orthology_confidence', display_name='Hybrid - Bos Indicus orthology confidence [0 low, 1 high]', description=''>
bthybrid_homolog_ensembl_gene | <biomart.Attribute name='bthybrid_homolog_ensembl_gene', display_name='Hybrid - Bos Taurus gene stable ID', description=''>
bthybrid_homolog_associated_gene_name | <biomart.Attribute name='bthybrid_homolog_associated_gene_name', display_name='Hybrid - Bos Taurus gene name', description=''>
bthybrid_homolog_ensembl_peptide | <biomart.Attribute name='bthybrid_homolog_ensembl_peptide', display_name='Hybrid - Bos Taurus protein or transcript stable ID', description=''>
bthybrid_homolog_chromosome | <biomart.Attribute name='bthybrid_homolog_chromosome', display_name='Hybrid - Bos Taurus chromosome/scaffold name', description=''>
bthybrid_homolog_chrom_start | <biomart.Attribute name='bthybrid_homolog_chrom_start', display_name='Hybrid - Bos Taurus chromosome/scaffold start (bp)', description=''>
bthybrid_homolog_chrom_end | <biomart.Attribute name='bthybrid_homolog_chrom_end', display_name='Hybrid - Bos Taurus chromosome/scaffold end (bp)', description=''>
bthybrid_homolog_canonical_transcript_protein | <biomart.Attribute name='bthybrid_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
bthybrid_homolog_subtype | <biomart.Attribute name='bthybrid_homolog_subtype', display_name='Last common ancestor with Hybrid - Bos Taurus', description=''>
bthybrid_homolog_orthology_type | <biomart.Attribute name='bthybrid_homolog_orthology_type', display_name='Hybrid - Bos Taurus homology type', description=''>
bthybrid_homolog_perc_id | <biomart.Attribute name='bthybrid_homolog_perc_id', display_name='%id. target Hybrid - Bos Taurus gene identical to query gene', description=''>
bthybrid_homolog_perc_id_r1 | <biomart.Attribute name='bthybrid_homolog_perc_id_r1', display_name='%id. query gene identical to target Hybrid - Bos Taurus gene', description=''>
bthybrid_homolog_goc_score | <biomart.Attribute name='bthybrid_homolog_goc_score', display_name='Hybrid - Bos Taurus Gene-order conservation score', description=''>
bthybrid_homolog_wga_coverage | <biomart.Attribute name='bthybrid_homolog_wga_coverage', display_name='Hybrid - Bos Taurus Whole-genome alignment coverage', description=''>
bthybrid_homolog_orthology_confidence | <biomart.Attribute name='bthybrid_homolog_orthology_confidence', display_name='Hybrid - Bos Taurus orthology confidence [0 low, 1 high]', description=''>
pcapensis_homolog_ensembl_gene | <biomart.Attribute name='pcapensis_homolog_ensembl_gene', display_name='Hyrax gene stable ID', description=''>
pcapensis_homolog_associated_gene_name | <biomart.Attribute name='pcapensis_homolog_associated_gene_name', display_name='Hyrax gene name', description=''>
pcapensis_homolog_ensembl_peptide | <biomart.Attribute name='pcapensis_homolog_ensembl_peptide', display_name='Hyrax protein or transcript stable ID', description=''>
pcapensis_homolog_chromosome | <biomart.Attribute name='pcapensis_homolog_chromosome', display_name='Hyrax chromosome/scaffold name', description=''>
pcapensis_homolog_chrom_start | <biomart.Attribute name='pcapensis_homolog_chrom_start', display_name='Hyrax chromosome/scaffold start (bp)', description=''>
pcapensis_homolog_chrom_end | <biomart.Attribute name='pcapensis_homolog_chrom_end', display_name='Hyrax chromosome/scaffold end (bp)', description=''>
pcapensis_homolog_canonical_transcript_protein | <biomart.Attribute name='pcapensis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pcapensis_homolog_subtype | <biomart.Attribute name='pcapensis_homolog_subtype', display_name='Last common ancestor with Hyrax', description=''>
pcapensis_homolog_orthology_type | <biomart.Attribute name='pcapensis_homolog_orthology_type', display_name='Hyrax homology type', description=''>
pcapensis_homolog_perc_id | <biomart.Attribute name='pcapensis_homolog_perc_id', display_name='%id. target Hyrax gene identical to query gene', description=''>
pcapensis_homolog_perc_id_r1 | <biomart.Attribute name='pcapensis_homolog_perc_id_r1', display_name='%id. query gene identical to target Hyrax gene', description=''>
pcapensis_homolog_goc_score | <biomart.Attribute name='pcapensis_homolog_goc_score', display_name='Hyrax Gene-order conservation score', description=''>
pcapensis_homolog_wga_coverage | <biomart.Attribute name='pcapensis_homolog_wga_coverage', display_name='Hyrax Whole-genome alignment coverage', description=''>
pcapensis_homolog_orthology_confidence | <biomart.Attribute name='pcapensis_homolog_orthology_confidence', display_name='Hyrax orthology confidence [0 low, 1 high]', description=''>
pranga_homolog_ensembl_gene | <biomart.Attribute name='pranga_homolog_ensembl_gene', display_name='Indian glassy fish gene stable ID', description=''>
pranga_homolog_associated_gene_name | <biomart.Attribute name='pranga_homolog_associated_gene_name', display_name='Indian glassy fish gene name', description=''>
pranga_homolog_ensembl_peptide | <biomart.Attribute name='pranga_homolog_ensembl_peptide', display_name='Indian glassy fish protein or transcript stable ID', description=''>
pranga_homolog_chromosome | <biomart.Attribute name='pranga_homolog_chromosome', display_name='Indian glassy fish chromosome/scaffold name', description=''>
pranga_homolog_chrom_start | <biomart.Attribute name='pranga_homolog_chrom_start', display_name='Indian glassy fish chromosome/scaffold start (bp)', description=''>
pranga_homolog_chrom_end | <biomart.Attribute name='pranga_homolog_chrom_end', display_name='Indian glassy fish chromosome/scaffold end (bp)', description=''>
pranga_homolog_canonical_transcript_protein | <biomart.Attribute name='pranga_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pranga_homolog_subtype | <biomart.Attribute name='pranga_homolog_subtype', display_name='Last common ancestor with Indian glassy fish', description=''>
pranga_homolog_orthology_type | <biomart.Attribute name='pranga_homolog_orthology_type', display_name='Indian glassy fish homology type', description=''>
pranga_homolog_perc_id | <biomart.Attribute name='pranga_homolog_perc_id', display_name='%id. target Indian glassy fish gene identical to query gene', description=''>
pranga_homolog_perc_id_r1 | <biomart.Attribute name='pranga_homolog_perc_id_r1', display_name='%id. query gene identical to target Indian glassy fish gene', description=''>
pranga_homolog_goc_score | <biomart.Attribute name='pranga_homolog_goc_score', display_name='Indian glassy fish Gene-order conservation score', description=''>
pranga_homolog_wga_coverage | <biomart.Attribute name='pranga_homolog_wga_coverage', display_name='Indian glassy fish Whole-genome alignment coverage', description=''>
pranga_homolog_orthology_confidence | <biomart.Attribute name='pranga_homolog_orthology_confidence', display_name='Indian glassy fish orthology confidence [0 low, 1 high]', description=''>
omelastigma_homolog_ensembl_gene | <biomart.Attribute name='omelastigma_homolog_ensembl_gene', display_name='Indian medaka gene stable ID', description=''>
omelastigma_homolog_associated_gene_name | <biomart.Attribute name='omelastigma_homolog_associated_gene_name', display_name='Indian medaka gene name', description=''>
omelastigma_homolog_ensembl_peptide | <biomart.Attribute name='omelastigma_homolog_ensembl_peptide', display_name='Indian medaka protein or transcript stable ID', description=''>
omelastigma_homolog_chromosome | <biomart.Attribute name='omelastigma_homolog_chromosome', display_name='Indian medaka chromosome/scaffold name', description=''>
omelastigma_homolog_chrom_start | <biomart.Attribute name='omelastigma_homolog_chrom_start', display_name='Indian medaka chromosome/scaffold start (bp)', description=''>
omelastigma_homolog_chrom_end | <biomart.Attribute name='omelastigma_homolog_chrom_end', display_name='Indian medaka chromosome/scaffold end (bp)', description=''>
omelastigma_homolog_canonical_transcript_protein | <biomart.Attribute name='omelastigma_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
omelastigma_homolog_subtype | <biomart.Attribute name='omelastigma_homolog_subtype', display_name='Last common ancestor with Indian medaka', description=''>
omelastigma_homolog_orthology_type | <biomart.Attribute name='omelastigma_homolog_orthology_type', display_name='Indian medaka homology type', description=''>
omelastigma_homolog_perc_id | <biomart.Attribute name='omelastigma_homolog_perc_id', display_name='%id. target Indian medaka gene identical to query gene', description=''>
omelastigma_homolog_perc_id_r1 | <biomart.Attribute name='omelastigma_homolog_perc_id_r1', display_name='%id. query gene identical to target Indian medaka gene', description=''>
omelastigma_homolog_goc_score | <biomart.Attribute name='omelastigma_homolog_goc_score', display_name='Indian medaka Gene-order conservation score', description=''>
omelastigma_homolog_wga_coverage | <biomart.Attribute name='omelastigma_homolog_wga_coverage', display_name='Indian medaka Whole-genome alignment coverage', description=''>
omelastigma_homolog_orthology_confidence | <biomart.Attribute name='omelastigma_homolog_orthology_confidence', display_name='Indian medaka orthology confidence [0 low, 1 high]', description=''>
olhni_homolog_ensembl_gene | <biomart.Attribute name='olhni_homolog_ensembl_gene', display_name='Japanese medaka HNI gene stable ID', description=''>
olhni_homolog_associated_gene_name | <biomart.Attribute name='olhni_homolog_associated_gene_name', display_name='Japanese medaka HNI gene name', description=''>
olhni_homolog_ensembl_peptide | <biomart.Attribute name='olhni_homolog_ensembl_peptide', display_name='Japanese medaka HNI protein or transcript stable ID', description=''>
olhni_homolog_chromosome | <biomart.Attribute name='olhni_homolog_chromosome', display_name='Japanese medaka HNI chromosome/scaffold name', description=''>
olhni_homolog_chrom_start | <biomart.Attribute name='olhni_homolog_chrom_start', display_name='Japanese medaka HNI chromosome/scaffold start (bp)', description=''>
olhni_homolog_chrom_end | <biomart.Attribute name='olhni_homolog_chrom_end', display_name='Japanese medaka HNI chromosome/scaffold end (bp)', description=''>
olhni_homolog_canonical_transcript_protein | <biomart.Attribute name='olhni_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
olhni_homolog_subtype | <biomart.Attribute name='olhni_homolog_subtype', display_name='Last common ancestor with Japanese medaka HNI', description=''>
olhni_homolog_orthology_type | <biomart.Attribute name='olhni_homolog_orthology_type', display_name='Japanese medaka HNI homology type', description=''>
olhni_homolog_perc_id | <biomart.Attribute name='olhni_homolog_perc_id', display_name='%id. target Japanese medaka HNI gene identical to query gene', description=''>
olhni_homolog_perc_id_r1 | <biomart.Attribute name='olhni_homolog_perc_id_r1', display_name='%id. query gene identical to target Japanese medaka HNI gene', description=''>
olhni_homolog_goc_score | <biomart.Attribute name='olhni_homolog_goc_score', display_name='Japanese medaka HNI Gene-order conservation score', description=''>
olhni_homolog_wga_coverage | <biomart.Attribute name='olhni_homolog_wga_coverage', display_name='Japanese medaka HNI Whole-genome alignment coverage', description=''>
olhni_homolog_orthology_confidence | <biomart.Attribute name='olhni_homolog_orthology_confidence', display_name='Japanese medaka HNI orthology confidence [0 low, 1 high]', description=''>
olhsok_homolog_ensembl_gene | <biomart.Attribute name='olhsok_homolog_ensembl_gene', display_name='Japanese medaka HSOK gene stable ID', description=''>
olhsok_homolog_associated_gene_name | <biomart.Attribute name='olhsok_homolog_associated_gene_name', display_name='Japanese medaka HSOK gene name', description=''>
olhsok_homolog_ensembl_peptide | <biomart.Attribute name='olhsok_homolog_ensembl_peptide', display_name='Japanese medaka HSOK protein or transcript stable ID', description=''>
olhsok_homolog_chromosome | <biomart.Attribute name='olhsok_homolog_chromosome', display_name='Japanese medaka HSOK chromosome/scaffold name', description=''>
olhsok_homolog_chrom_start | <biomart.Attribute name='olhsok_homolog_chrom_start', display_name='Japanese medaka HSOK chromosome/scaffold start (bp)', description=''>
olhsok_homolog_chrom_end | <biomart.Attribute name='olhsok_homolog_chrom_end', display_name='Japanese medaka HSOK chromosome/scaffold end (bp)', description=''>
olhsok_homolog_canonical_transcript_protein | <biomart.Attribute name='olhsok_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
olhsok_homolog_subtype | <biomart.Attribute name='olhsok_homolog_subtype', display_name='Last common ancestor with Japanese medaka HSOK', description=''>
olhsok_homolog_orthology_type | <biomart.Attribute name='olhsok_homolog_orthology_type', display_name='Japanese medaka HSOK homology type', description=''>
olhsok_homolog_perc_id | <biomart.Attribute name='olhsok_homolog_perc_id', display_name='%id. target Japanese medaka HSOK gene identical to query gene', description=''>
olhsok_homolog_perc_id_r1 | <biomart.Attribute name='olhsok_homolog_perc_id_r1', display_name='%id. query gene identical to target Japanese medaka HSOK gene', description=''>
olhsok_homolog_goc_score | <biomart.Attribute name='olhsok_homolog_goc_score', display_name='Japanese medaka HSOK Gene-order conservation score', description=''>
olhsok_homolog_wga_coverage | <biomart.Attribute name='olhsok_homolog_wga_coverage', display_name='Japanese medaka HSOK Whole-genome alignment coverage', description=''>
olhsok_homolog_orthology_confidence | <biomart.Attribute name='olhsok_homolog_orthology_confidence', display_name='Japanese medaka HSOK orthology confidence [0 low, 1 high]', description=''>
olatipes_homolog_ensembl_gene | <biomart.Attribute name='olatipes_homolog_ensembl_gene', display_name='Japanese medaka HdrR gene stable ID', description=''>
olatipes_homolog_associated_gene_name | <biomart.Attribute name='olatipes_homolog_associated_gene_name', display_name='Japanese medaka HdrR gene name', description=''>
olatipes_homolog_ensembl_peptide | <biomart.Attribute name='olatipes_homolog_ensembl_peptide', display_name='Japanese medaka HdrR protein or transcript stable ID', description=''>
olatipes_homolog_chromosome | <biomart.Attribute name='olatipes_homolog_chromosome', display_name='Japanese medaka HdrR chromosome/scaffold name', description=''>
olatipes_homolog_chrom_start | <biomart.Attribute name='olatipes_homolog_chrom_start', display_name='Japanese medaka HdrR chromosome/scaffold start (bp)', description=''>
olatipes_homolog_chrom_end | <biomart.Attribute name='olatipes_homolog_chrom_end', display_name='Japanese medaka HdrR chromosome/scaffold end (bp)', description=''>
olatipes_homolog_canonical_transcript_protein | <biomart.Attribute name='olatipes_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
olatipes_homolog_subtype | <biomart.Attribute name='olatipes_homolog_subtype', display_name='Last common ancestor with Japanese medaka HdrR', description=''>
olatipes_homolog_orthology_type | <biomart.Attribute name='olatipes_homolog_orthology_type', display_name='Japanese medaka HdrR homology type', description=''>
olatipes_homolog_perc_id | <biomart.Attribute name='olatipes_homolog_perc_id', display_name='%id. target Japanese medaka HdrR gene identical to query gene', description=''>
olatipes_homolog_perc_id_r1 | <biomart.Attribute name='olatipes_homolog_perc_id_r1', display_name='%id. query gene identical to target Japanese medaka HdrR gene', description=''>
olatipes_homolog_goc_score | <biomart.Attribute name='olatipes_homolog_goc_score', display_name='Japanese medaka HdrR Gene-order conservation score', description=''>
olatipes_homolog_wga_coverage | <biomart.Attribute name='olatipes_homolog_wga_coverage', display_name='Japanese medaka HdrR Whole-genome alignment coverage', description=''>
olatipes_homolog_orthology_confidence | <biomart.Attribute name='olatipes_homolog_orthology_confidence', display_name='Japanese medaka HdrR orthology confidence [0 low, 1 high]', description=''>
cjaponica_homolog_ensembl_gene | <biomart.Attribute name='cjaponica_homolog_ensembl_gene', display_name='Japanese quail gene stable ID', description=''>
cjaponica_homolog_associated_gene_name | <biomart.Attribute name='cjaponica_homolog_associated_gene_name', display_name='Japanese quail gene name', description=''>
cjaponica_homolog_ensembl_peptide | <biomart.Attribute name='cjaponica_homolog_ensembl_peptide', display_name='Japanese quail protein or transcript stable ID', description=''>
cjaponica_homolog_chromosome | <biomart.Attribute name='cjaponica_homolog_chromosome', display_name='Japanese quail chromosome/scaffold name', description=''>
cjaponica_homolog_chrom_start | <biomart.Attribute name='cjaponica_homolog_chrom_start', display_name='Japanese quail chromosome/scaffold start (bp)', description=''>
cjaponica_homolog_chrom_end | <biomart.Attribute name='cjaponica_homolog_chrom_end', display_name='Japanese quail chromosome/scaffold end (bp)', description=''>
cjaponica_homolog_canonical_transcript_protein | <biomart.Attribute name='cjaponica_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cjaponica_homolog_subtype | <biomart.Attribute name='cjaponica_homolog_subtype', display_name='Last common ancestor with Japanese quail', description=''>
cjaponica_homolog_orthology_type | <biomart.Attribute name='cjaponica_homolog_orthology_type', display_name='Japanese quail homology type', description=''>
cjaponica_homolog_perc_id | <biomart.Attribute name='cjaponica_homolog_perc_id', display_name='%id. target Japanese quail gene identical to query gene', description=''>
cjaponica_homolog_perc_id_r1 | <biomart.Attribute name='cjaponica_homolog_perc_id_r1', display_name='%id. query gene identical to target Japanese quail gene', description=''>
cjaponica_homolog_goc_score | <biomart.Attribute name='cjaponica_homolog_goc_score', display_name='Japanese quail Gene-order conservation score', description=''>
cjaponica_homolog_wga_coverage | <biomart.Attribute name='cjaponica_homolog_wga_coverage', display_name='Japanese quail Whole-genome alignment coverage', description=''>
cjaponica_homolog_orthology_confidence | <biomart.Attribute name='cjaponica_homolog_orthology_confidence', display_name='Japanese quail orthology confidence [0 low, 1 high]', description=''>
sfasciatus_homolog_ensembl_gene | <biomart.Attribute name='sfasciatus_homolog_ensembl_gene', display_name='Jewelled blenny gene stable ID', description=''>
sfasciatus_homolog_associated_gene_name | <biomart.Attribute name='sfasciatus_homolog_associated_gene_name', display_name='Jewelled blenny gene name', description=''>
sfasciatus_homolog_ensembl_peptide | <biomart.Attribute name='sfasciatus_homolog_ensembl_peptide', display_name='Jewelled blenny protein or transcript stable ID', description=''>
sfasciatus_homolog_chromosome | <biomart.Attribute name='sfasciatus_homolog_chromosome', display_name='Jewelled blenny chromosome/scaffold name', description=''>
sfasciatus_homolog_chrom_start | <biomart.Attribute name='sfasciatus_homolog_chrom_start', display_name='Jewelled blenny chromosome/scaffold start (bp)', description=''>
sfasciatus_homolog_chrom_end | <biomart.Attribute name='sfasciatus_homolog_chrom_end', display_name='Jewelled blenny chromosome/scaffold end (bp)', description=''>
sfasciatus_homolog_canonical_transcript_protein | <biomart.Attribute name='sfasciatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sfasciatus_homolog_subtype | <biomart.Attribute name='sfasciatus_homolog_subtype', display_name='Last common ancestor with Jewelled blenny', description=''>
sfasciatus_homolog_orthology_type | <biomart.Attribute name='sfasciatus_homolog_orthology_type', display_name='Jewelled blenny homology type', description=''>
sfasciatus_homolog_perc_id | <biomart.Attribute name='sfasciatus_homolog_perc_id', display_name='%id. target Jewelled blenny gene identical to query gene', description=''>
sfasciatus_homolog_perc_id_r1 | <biomart.Attribute name='sfasciatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Jewelled blenny gene', description=''>
sfasciatus_homolog_goc_score | <biomart.Attribute name='sfasciatus_homolog_goc_score', display_name='Jewelled blenny Gene-order conservation score', description=''>
sfasciatus_homolog_wga_coverage | <biomart.Attribute name='sfasciatus_homolog_wga_coverage', display_name='Jewelled blenny Whole-genome alignment coverage', description=''>
sfasciatus_homolog_orthology_confidence | <biomart.Attribute name='sfasciatus_homolog_orthology_confidence', display_name='Jewelled blenny orthology confidence [0 low, 1 high]', description=''>
shabroptila_homolog_ensembl_gene | <biomart.Attribute name='shabroptila_homolog_ensembl_gene', display_name='Kakapo gene stable ID', description=''>
shabroptila_homolog_associated_gene_name | <biomart.Attribute name='shabroptila_homolog_associated_gene_name', display_name='Kakapo gene name', description=''>
shabroptila_homolog_ensembl_peptide | <biomart.Attribute name='shabroptila_homolog_ensembl_peptide', display_name='Kakapo protein or transcript stable ID', description=''>
shabroptila_homolog_chromosome | <biomart.Attribute name='shabroptila_homolog_chromosome', display_name='Kakapo chromosome/scaffold name', description=''>
shabroptila_homolog_chrom_start | <biomart.Attribute name='shabroptila_homolog_chrom_start', display_name='Kakapo chromosome/scaffold start (bp)', description=''>
shabroptila_homolog_chrom_end | <biomart.Attribute name='shabroptila_homolog_chrom_end', display_name='Kakapo chromosome/scaffold end (bp)', description=''>
shabroptila_homolog_canonical_transcript_protein | <biomart.Attribute name='shabroptila_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
shabroptila_homolog_subtype | <biomart.Attribute name='shabroptila_homolog_subtype', display_name='Last common ancestor with Kakapo', description=''>
shabroptila_homolog_orthology_type | <biomart.Attribute name='shabroptila_homolog_orthology_type', display_name='Kakapo homology type', description=''>
shabroptila_homolog_perc_id | <biomart.Attribute name='shabroptila_homolog_perc_id', display_name='%id. target Kakapo gene identical to query gene', description=''>
shabroptila_homolog_perc_id_r1 | <biomart.Attribute name='shabroptila_homolog_perc_id_r1', display_name='%id. query gene identical to target Kakapo gene', description=''>
shabroptila_homolog_goc_score | <biomart.Attribute name='shabroptila_homolog_goc_score', display_name='Kakapo Gene-order conservation score', description=''>
shabroptila_homolog_wga_coverage | <biomart.Attribute name='shabroptila_homolog_wga_coverage', display_name='Kakapo Whole-genome alignment coverage', description=''>
shabroptila_homolog_orthology_confidence | <biomart.Attribute name='shabroptila_homolog_orthology_confidence', display_name='Kakapo orthology confidence [0 low, 1 high]', description=''>
dordii_homolog_ensembl_gene | <biomart.Attribute name='dordii_homolog_ensembl_gene', display_name='Kangaroo rat gene stable ID', description=''>
dordii_homolog_associated_gene_name | <biomart.Attribute name='dordii_homolog_associated_gene_name', display_name='Kangaroo rat gene name', description=''>
dordii_homolog_ensembl_peptide | <biomart.Attribute name='dordii_homolog_ensembl_peptide', display_name='Kangaroo rat protein or transcript stable ID', description=''>
dordii_homolog_chromosome | <biomart.Attribute name='dordii_homolog_chromosome', display_name='Kangaroo rat chromosome/scaffold name', description=''>
dordii_homolog_chrom_start | <biomart.Attribute name='dordii_homolog_chrom_start', display_name='Kangaroo rat chromosome/scaffold start (bp)', description=''>
dordii_homolog_chrom_end | <biomart.Attribute name='dordii_homolog_chrom_end', display_name='Kangaroo rat chromosome/scaffold end (bp)', description=''>
dordii_homolog_canonical_transcript_protein | <biomart.Attribute name='dordii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
dordii_homolog_subtype | <biomart.Attribute name='dordii_homolog_subtype', display_name='Last common ancestor with Kangaroo rat', description=''>
dordii_homolog_orthology_type | <biomart.Attribute name='dordii_homolog_orthology_type', display_name='Kangaroo rat homology type', description=''>
dordii_homolog_perc_id | <biomart.Attribute name='dordii_homolog_perc_id', display_name='%id. target Kangaroo rat gene identical to query gene', description=''>
dordii_homolog_perc_id_r1 | <biomart.Attribute name='dordii_homolog_perc_id_r1', display_name='%id. query gene identical to target Kangaroo rat gene', description=''>
dordii_homolog_goc_score | <biomart.Attribute name='dordii_homolog_goc_score', display_name='Kangaroo rat Gene-order conservation score', description=''>
dordii_homolog_wga_coverage | <biomart.Attribute name='dordii_homolog_wga_coverage', display_name='Kangaroo rat Whole-genome alignment coverage', description=''>
dordii_homolog_orthology_confidence | <biomart.Attribute name='dordii_homolog_orthology_confidence', display_name='Kangaroo rat orthology confidence [0 low, 1 high]', description=''>
pcinereus_homolog_ensembl_gene | <biomart.Attribute name='pcinereus_homolog_ensembl_gene', display_name='Koala gene stable ID', description=''>
pcinereus_homolog_associated_gene_name | <biomart.Attribute name='pcinereus_homolog_associated_gene_name', display_name='Koala gene name', description=''>
pcinereus_homolog_ensembl_peptide | <biomart.Attribute name='pcinereus_homolog_ensembl_peptide', display_name='Koala protein or transcript stable ID', description=''>
pcinereus_homolog_chromosome | <biomart.Attribute name='pcinereus_homolog_chromosome', display_name='Koala chromosome/scaffold name', description=''>
pcinereus_homolog_chrom_start | <biomart.Attribute name='pcinereus_homolog_chrom_start', display_name='Koala chromosome/scaffold start (bp)', description=''>
pcinereus_homolog_chrom_end | <biomart.Attribute name='pcinereus_homolog_chrom_end', display_name='Koala chromosome/scaffold end (bp)', description=''>
pcinereus_homolog_canonical_transcript_protein | <biomart.Attribute name='pcinereus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pcinereus_homolog_subtype | <biomart.Attribute name='pcinereus_homolog_subtype', display_name='Last common ancestor with Koala', description=''>
pcinereus_homolog_orthology_type | <biomart.Attribute name='pcinereus_homolog_orthology_type', display_name='Koala homology type', description=''>
pcinereus_homolog_perc_id | <biomart.Attribute name='pcinereus_homolog_perc_id', display_name='%id. target Koala gene identical to query gene', description=''>
pcinereus_homolog_perc_id_r1 | <biomart.Attribute name='pcinereus_homolog_perc_id_r1', display_name='%id. query gene identical to target Koala gene', description=''>
pcinereus_homolog_goc_score | <biomart.Attribute name='pcinereus_homolog_goc_score', display_name='Koala Gene-order conservation score', description=''>
pcinereus_homolog_wga_coverage | <biomart.Attribute name='pcinereus_homolog_wga_coverage', display_name='Koala Whole-genome alignment coverage', description=''>
pcinereus_homolog_orthology_confidence | <biomart.Attribute name='pcinereus_homolog_orthology_confidence', display_name='Koala orthology confidence [0 low, 1 high]', description=''>
pmarinus_homolog_ensembl_gene | <biomart.Attribute name='pmarinus_homolog_ensembl_gene', display_name='Lamprey gene stable ID', description=''>
pmarinus_homolog_associated_gene_name | <biomart.Attribute name='pmarinus_homolog_associated_gene_name', display_name='Lamprey gene name', description=''>
pmarinus_homolog_ensembl_peptide | <biomart.Attribute name='pmarinus_homolog_ensembl_peptide', display_name='Lamprey protein or transcript stable ID', description=''>
pmarinus_homolog_chromosome | <biomart.Attribute name='pmarinus_homolog_chromosome', display_name='Lamprey chromosome/scaffold name', description=''>
pmarinus_homolog_chrom_start | <biomart.Attribute name='pmarinus_homolog_chrom_start', display_name='Lamprey chromosome/scaffold start (bp)', description=''>
pmarinus_homolog_chrom_end | <biomart.Attribute name='pmarinus_homolog_chrom_end', display_name='Lamprey chromosome/scaffold end (bp)', description=''>
pmarinus_homolog_canonical_transcript_protein | <biomart.Attribute name='pmarinus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pmarinus_homolog_subtype | <biomart.Attribute name='pmarinus_homolog_subtype', display_name='Last common ancestor with Lamprey', description=''>
pmarinus_homolog_orthology_type | <biomart.Attribute name='pmarinus_homolog_orthology_type', display_name='Lamprey homology type', description=''>
pmarinus_homolog_perc_id | <biomart.Attribute name='pmarinus_homolog_perc_id', display_name='%id. target Lamprey gene identical to query gene', description=''>
pmarinus_homolog_perc_id_r1 | <biomart.Attribute name='pmarinus_homolog_perc_id_r1', display_name='%id. query gene identical to target Lamprey gene', description=''>
pmarinus_homolog_wga_coverage | <biomart.Attribute name='pmarinus_homolog_wga_coverage', display_name='Lamprey Whole-genome alignment coverage', description=''>
pmarinus_homolog_orthology_confidence | <biomart.Attribute name='pmarinus_homolog_orthology_confidence', display_name='Lamprey orthology confidence [0 low, 1 high]', description=''>
lcrocea_homolog_ensembl_gene | <biomart.Attribute name='lcrocea_homolog_ensembl_gene', display_name='Large yellow croaker gene stable ID', description=''>
lcrocea_homolog_associated_gene_name | <biomart.Attribute name='lcrocea_homolog_associated_gene_name', display_name='Large yellow croaker gene name', description=''>
lcrocea_homolog_ensembl_peptide | <biomart.Attribute name='lcrocea_homolog_ensembl_peptide', display_name='Large yellow croaker protein or transcript stable ID', description=''>
lcrocea_homolog_chromosome | <biomart.Attribute name='lcrocea_homolog_chromosome', display_name='Large yellow croaker chromosome/scaffold name', description=''>
lcrocea_homolog_chrom_start | <biomart.Attribute name='lcrocea_homolog_chrom_start', display_name='Large yellow croaker chromosome/scaffold start (bp)', description=''>
lcrocea_homolog_chrom_end | <biomart.Attribute name='lcrocea_homolog_chrom_end', display_name='Large yellow croaker chromosome/scaffold end (bp)', description=''>
lcrocea_homolog_canonical_transcript_protein | <biomart.Attribute name='lcrocea_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
lcrocea_homolog_subtype | <biomart.Attribute name='lcrocea_homolog_subtype', display_name='Last common ancestor with Large yellow croaker', description=''>
lcrocea_homolog_orthology_type | <biomart.Attribute name='lcrocea_homolog_orthology_type', display_name='Large yellow croaker homology type', description=''>
lcrocea_homolog_perc_id | <biomart.Attribute name='lcrocea_homolog_perc_id', display_name='%id. target Large yellow croaker gene identical to query gene', description=''>
lcrocea_homolog_perc_id_r1 | <biomart.Attribute name='lcrocea_homolog_perc_id_r1', display_name='%id. query gene identical to target Large yellow croaker gene', description=''>
lcrocea_homolog_goc_score | <biomart.Attribute name='lcrocea_homolog_goc_score', display_name='Large yellow croaker Gene-order conservation score', description=''>
lcrocea_homolog_wga_coverage | <biomart.Attribute name='lcrocea_homolog_wga_coverage', display_name='Large yellow croaker Whole-genome alignment coverage', description=''>
lcrocea_homolog_orthology_confidence | <biomart.Attribute name='lcrocea_homolog_orthology_confidence', display_name='Large yellow croaker orthology confidence [0 low, 1 high]', description=''>
ppardus_homolog_ensembl_gene | <biomart.Attribute name='ppardus_homolog_ensembl_gene', display_name='Leopard gene stable ID', description=''>
ppardus_homolog_associated_gene_name | <biomart.Attribute name='ppardus_homolog_associated_gene_name', display_name='Leopard gene name', description=''>
ppardus_homolog_ensembl_peptide | <biomart.Attribute name='ppardus_homolog_ensembl_peptide', display_name='Leopard protein or transcript stable ID', description=''>
ppardus_homolog_chromosome | <biomart.Attribute name='ppardus_homolog_chromosome', display_name='Leopard chromosome/scaffold name', description=''>
ppardus_homolog_chrom_start | <biomart.Attribute name='ppardus_homolog_chrom_start', display_name='Leopard chromosome/scaffold start (bp)', description=''>
ppardus_homolog_chrom_end | <biomart.Attribute name='ppardus_homolog_chrom_end', display_name='Leopard chromosome/scaffold end (bp)', description=''>
ppardus_homolog_canonical_transcript_protein | <biomart.Attribute name='ppardus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ppardus_homolog_subtype | <biomart.Attribute name='ppardus_homolog_subtype', display_name='Last common ancestor with Leopard', description=''>
ppardus_homolog_orthology_type | <biomart.Attribute name='ppardus_homolog_orthology_type', display_name='Leopard homology type', description=''>
ppardus_homolog_perc_id | <biomart.Attribute name='ppardus_homolog_perc_id', display_name='%id. target Leopard gene identical to query gene', description=''>
ppardus_homolog_perc_id_r1 | <biomart.Attribute name='ppardus_homolog_perc_id_r1', display_name='%id. query gene identical to target Leopard gene', description=''>
ppardus_homolog_goc_score | <biomart.Attribute name='ppardus_homolog_goc_score', display_name='Leopard Gene-order conservation score', description=''>
ppardus_homolog_wga_coverage | <biomart.Attribute name='ppardus_homolog_wga_coverage', display_name='Leopard Whole-genome alignment coverage', description=''>
ppardus_homolog_orthology_confidence | <biomart.Attribute name='ppardus_homolog_orthology_confidence', display_name='Leopard orthology confidence [0 low, 1 high]', description=''>
jjaculus_homolog_ensembl_gene | <biomart.Attribute name='jjaculus_homolog_ensembl_gene', display_name='Lesser Egyptian jerboa gene stable ID', description=''>
jjaculus_homolog_associated_gene_name | <biomart.Attribute name='jjaculus_homolog_associated_gene_name', display_name='Lesser Egyptian jerboa gene name', description=''>
jjaculus_homolog_ensembl_peptide | <biomart.Attribute name='jjaculus_homolog_ensembl_peptide', display_name='Lesser Egyptian jerboa protein or transcript stable ID', description=''>
jjaculus_homolog_chromosome | <biomart.Attribute name='jjaculus_homolog_chromosome', display_name='Lesser Egyptian jerboa chromosome/scaffold name', description=''>
jjaculus_homolog_chrom_start | <biomart.Attribute name='jjaculus_homolog_chrom_start', display_name='Lesser Egyptian jerboa chromosome/scaffold start (bp)', description=''>
jjaculus_homolog_chrom_end | <biomart.Attribute name='jjaculus_homolog_chrom_end', display_name='Lesser Egyptian jerboa chromosome/scaffold end (bp)', description=''>
jjaculus_homolog_canonical_transcript_protein | <biomart.Attribute name='jjaculus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
jjaculus_homolog_subtype | <biomart.Attribute name='jjaculus_homolog_subtype', display_name='Last common ancestor with Lesser Egyptian jerboa', description=''>
jjaculus_homolog_orthology_type | <biomart.Attribute name='jjaculus_homolog_orthology_type', display_name='Lesser Egyptian jerboa homology type', description=''>
jjaculus_homolog_perc_id | <biomart.Attribute name='jjaculus_homolog_perc_id', display_name='%id. target Lesser Egyptian jerboa gene identical to query gene', description=''>
jjaculus_homolog_perc_id_r1 | <biomart.Attribute name='jjaculus_homolog_perc_id_r1', display_name='%id. query gene identical to target Lesser Egyptian jerboa gene', description=''>
jjaculus_homolog_goc_score | <biomart.Attribute name='jjaculus_homolog_goc_score', display_name='Lesser Egyptian jerboa Gene-order conservation score', description=''>
jjaculus_homolog_wga_coverage | <biomart.Attribute name='jjaculus_homolog_wga_coverage', display_name='Lesser Egyptian jerboa Whole-genome alignment coverage', description=''>
jjaculus_homolog_orthology_confidence | <biomart.Attribute name='jjaculus_homolog_orthology_confidence', display_name='Lesser Egyptian jerboa orthology confidence [0 low, 1 high]', description=''>
etelfairi_homolog_ensembl_gene | <biomart.Attribute name='etelfairi_homolog_ensembl_gene', display_name='Lesser hedgehog tenrec gene stable ID', description=''>
etelfairi_homolog_associated_gene_name | <biomart.Attribute name='etelfairi_homolog_associated_gene_name', display_name='Lesser hedgehog tenrec gene name', description=''>
etelfairi_homolog_ensembl_peptide | <biomart.Attribute name='etelfairi_homolog_ensembl_peptide', display_name='Lesser hedgehog tenrec protein or transcript stable ID', description=''>
etelfairi_homolog_chromosome | <biomart.Attribute name='etelfairi_homolog_chromosome', display_name='Lesser hedgehog tenrec chromosome/scaffold name', description=''>
etelfairi_homolog_chrom_start | <biomart.Attribute name='etelfairi_homolog_chrom_start', display_name='Lesser hedgehog tenrec chromosome/scaffold start (bp)', description=''>
etelfairi_homolog_chrom_end | <biomart.Attribute name='etelfairi_homolog_chrom_end', display_name='Lesser hedgehog tenrec chromosome/scaffold end (bp)', description=''>
etelfairi_homolog_canonical_transcript_protein | <biomart.Attribute name='etelfairi_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
etelfairi_homolog_subtype | <biomart.Attribute name='etelfairi_homolog_subtype', display_name='Last common ancestor with Lesser hedgehog tenrec', description=''>
etelfairi_homolog_orthology_type | <biomart.Attribute name='etelfairi_homolog_orthology_type', display_name='Lesser hedgehog tenrec homology type', description=''>
etelfairi_homolog_perc_id | <biomart.Attribute name='etelfairi_homolog_perc_id', display_name='%id. target Lesser hedgehog tenrec gene identical to query gene', description=''>
etelfairi_homolog_perc_id_r1 | <biomart.Attribute name='etelfairi_homolog_perc_id_r1', display_name='%id. query gene identical to target Lesser hedgehog tenrec gene', description=''>
etelfairi_homolog_goc_score | <biomart.Attribute name='etelfairi_homolog_goc_score', display_name='Lesser hedgehog tenrec Gene-order conservation score', description=''>
etelfairi_homolog_wga_coverage | <biomart.Attribute name='etelfairi_homolog_wga_coverage', display_name='Lesser hedgehog tenrec Whole-genome alignment coverage', description=''>
etelfairi_homolog_orthology_confidence | <biomart.Attribute name='etelfairi_homolog_orthology_confidence', display_name='Lesser hedgehog tenrec orthology confidence [0 low, 1 high]', description=''>
enaucrates_homolog_ensembl_gene | <biomart.Attribute name='enaucrates_homolog_ensembl_gene', display_name='Live sharksucker gene stable ID', description=''>
enaucrates_homolog_associated_gene_name | <biomart.Attribute name='enaucrates_homolog_associated_gene_name', display_name='Live sharksucker gene name', description=''>
enaucrates_homolog_ensembl_peptide | <biomart.Attribute name='enaucrates_homolog_ensembl_peptide', display_name='Live sharksucker protein or transcript stable ID', description=''>
enaucrates_homolog_chromosome | <biomart.Attribute name='enaucrates_homolog_chromosome', display_name='Live sharksucker chromosome/scaffold name', description=''>
enaucrates_homolog_chrom_start | <biomart.Attribute name='enaucrates_homolog_chrom_start', display_name='Live sharksucker chromosome/scaffold start (bp)', description=''>
enaucrates_homolog_chrom_end | <biomart.Attribute name='enaucrates_homolog_chrom_end', display_name='Live sharksucker chromosome/scaffold end (bp)', description=''>
enaucrates_homolog_canonical_transcript_protein | <biomart.Attribute name='enaucrates_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
enaucrates_homolog_subtype | <biomart.Attribute name='enaucrates_homolog_subtype', display_name='Last common ancestor with Live sharksucker', description=''>
enaucrates_homolog_orthology_type | <biomart.Attribute name='enaucrates_homolog_orthology_type', display_name='Live sharksucker homology type', description=''>
enaucrates_homolog_perc_id | <biomart.Attribute name='enaucrates_homolog_perc_id', display_name='%id. target Live sharksucker gene identical to query gene', description=''>
enaucrates_homolog_perc_id_r1 | <biomart.Attribute name='enaucrates_homolog_perc_id_r1', display_name='%id. query gene identical to target Live sharksucker gene', description=''>
enaucrates_homolog_goc_score | <biomart.Attribute name='enaucrates_homolog_goc_score', display_name='Live sharksucker Gene-order conservation score', description=''>
enaucrates_homolog_wga_coverage | <biomart.Attribute name='enaucrates_homolog_wga_coverage', display_name='Live sharksucker Whole-genome alignment coverage', description=''>
enaucrates_homolog_orthology_confidence | <biomart.Attribute name='enaucrates_homolog_orthology_confidence', display_name='Live sharksucker orthology confidence [0 low, 1 high]', description=''>
clanigera_homolog_ensembl_gene | <biomart.Attribute name='clanigera_homolog_ensembl_gene', display_name='Long-tailed chinchilla gene stable ID', description=''>
clanigera_homolog_associated_gene_name | <biomart.Attribute name='clanigera_homolog_associated_gene_name', display_name='Long-tailed chinchilla gene name', description=''>
clanigera_homolog_ensembl_peptide | <biomart.Attribute name='clanigera_homolog_ensembl_peptide', display_name='Long-tailed chinchilla protein or transcript stable ID', description=''>
clanigera_homolog_chromosome | <biomart.Attribute name='clanigera_homolog_chromosome', display_name='Long-tailed chinchilla chromosome/scaffold name', description=''>
clanigera_homolog_chrom_start | <biomart.Attribute name='clanigera_homolog_chrom_start', display_name='Long-tailed chinchilla chromosome/scaffold start (bp)', description=''>
clanigera_homolog_chrom_end | <biomart.Attribute name='clanigera_homolog_chrom_end', display_name='Long-tailed chinchilla chromosome/scaffold end (bp)', description=''>
clanigera_homolog_canonical_transcript_protein | <biomart.Attribute name='clanigera_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
clanigera_homolog_subtype | <biomart.Attribute name='clanigera_homolog_subtype', display_name='Last common ancestor with Long-tailed chinchilla', description=''>
clanigera_homolog_orthology_type | <biomart.Attribute name='clanigera_homolog_orthology_type', display_name='Long-tailed chinchilla homology type', description=''>
clanigera_homolog_perc_id | <biomart.Attribute name='clanigera_homolog_perc_id', display_name='%id. target Long-tailed chinchilla gene identical to query gene', description=''>
clanigera_homolog_perc_id_r1 | <biomart.Attribute name='clanigera_homolog_perc_id_r1', display_name='%id. query gene identical to target Long-tailed chinchilla gene', description=''>
clanigera_homolog_goc_score | <biomart.Attribute name='clanigera_homolog_goc_score', display_name='Long-tailed chinchilla Gene-order conservation score', description=''>
clanigera_homolog_wga_coverage | <biomart.Attribute name='clanigera_homolog_wga_coverage', display_name='Long-tailed chinchilla Whole-genome alignment coverage', description=''>
clanigera_homolog_orthology_confidence | <biomart.Attribute name='clanigera_homolog_orthology_confidence', display_name='Long-tailed chinchilla orthology confidence [0 low, 1 high]', description=''>
anancymaae_homolog_ensembl_gene | <biomart.Attribute name='anancymaae_homolog_ensembl_gene', display_name="Ma's night monkey gene stable ID", description=''>
anancymaae_homolog_associated_gene_name | <biomart.Attribute name='anancymaae_homolog_associated_gene_name', display_name="Ma's night monkey gene name", description=''>
anancymaae_homolog_ensembl_peptide | <biomart.Attribute name='anancymaae_homolog_ensembl_peptide', display_name="Ma's night monkey protein or transcript stable ID", description=''>
anancymaae_homolog_chromosome | <biomart.Attribute name='anancymaae_homolog_chromosome', display_name="Ma's night monkey chromosome/scaffold name", description=''>
anancymaae_homolog_chrom_start | <biomart.Attribute name='anancymaae_homolog_chrom_start', display_name="Ma's night monkey chromosome/scaffold start (bp)", description=''>
anancymaae_homolog_chrom_end | <biomart.Attribute name='anancymaae_homolog_chrom_end', display_name="Ma's night monkey chromosome/scaffold end (bp)", description=''>
anancymaae_homolog_canonical_transcript_protein | <biomart.Attribute name='anancymaae_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
anancymaae_homolog_subtype | <biomart.Attribute name='anancymaae_homolog_subtype', display_name="Last common ancestor with Ma's night monkey", description=''>
anancymaae_homolog_orthology_type | <biomart.Attribute name='anancymaae_homolog_orthology_type', display_name="Ma's night monkey homology type", description=''>
anancymaae_homolog_perc_id | <biomart.Attribute name='anancymaae_homolog_perc_id', display_name="%id. target Ma's night monkey gene identical to query gene", description=''>
anancymaae_homolog_perc_id_r1 | <biomart.Attribute name='anancymaae_homolog_perc_id_r1', display_name="%id. query gene identical to target Ma's night monkey gene", description=''>
anancymaae_homolog_goc_score | <biomart.Attribute name='anancymaae_homolog_goc_score', display_name="Ma's night monkey Gene-order conservation score", description=''>
anancymaae_homolog_wga_coverage | <biomart.Attribute name='anancymaae_homolog_wga_coverage', display_name="Ma's night monkey Whole-genome alignment coverage", description=''>
anancymaae_homolog_orthology_confidence | <biomart.Attribute name='anancymaae_homolog_orthology_confidence', display_name="Ma's night monkey orthology confidence [0 low, 1 high]", description=''>
mmulatta_homolog_ensembl_gene | <biomart.Attribute name='mmulatta_homolog_ensembl_gene', display_name='Macaque gene stable ID', description=''>
mmulatta_homolog_associated_gene_name | <biomart.Attribute name='mmulatta_homolog_associated_gene_name', display_name='Macaque gene name', description=''>
mmulatta_homolog_ensembl_peptide | <biomart.Attribute name='mmulatta_homolog_ensembl_peptide', display_name='Macaque protein or transcript stable ID', description=''>
mmulatta_homolog_chromosome | <biomart.Attribute name='mmulatta_homolog_chromosome', display_name='Macaque chromosome/scaffold name', description=''>
mmulatta_homolog_chrom_start | <biomart.Attribute name='mmulatta_homolog_chrom_start', display_name='Macaque chromosome/scaffold start (bp)', description=''>
mmulatta_homolog_chrom_end | <biomart.Attribute name='mmulatta_homolog_chrom_end', display_name='Macaque chromosome/scaffold end (bp)', description=''>
mmulatta_homolog_canonical_transcript_protein | <biomart.Attribute name='mmulatta_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mmulatta_homolog_subtype | <biomart.Attribute name='mmulatta_homolog_subtype', display_name='Last common ancestor with Macaque', description=''>
mmulatta_homolog_orthology_type | <biomart.Attribute name='mmulatta_homolog_orthology_type', display_name='Macaque homology type', description=''>
mmulatta_homolog_perc_id | <biomart.Attribute name='mmulatta_homolog_perc_id', display_name='%id. target Macaque gene identical to query gene', description=''>
mmulatta_homolog_perc_id_r1 | <biomart.Attribute name='mmulatta_homolog_perc_id_r1', display_name='%id. query gene identical to target Macaque gene', description=''>
mmulatta_homolog_goc_score | <biomart.Attribute name='mmulatta_homolog_goc_score', display_name='Macaque Gene-order conservation score', description=''>
mmulatta_homolog_wga_coverage | <biomart.Attribute name='mmulatta_homolog_wga_coverage', display_name='Macaque Whole-genome alignment coverage', description=''>
mmulatta_homolog_orthology_confidence | <biomart.Attribute name='mmulatta_homolog_orthology_confidence', display_name='Macaque orthology confidence [0 low, 1 high]', description=''>
pnyererei_homolog_ensembl_gene | <biomart.Attribute name='pnyererei_homolog_ensembl_gene', display_name='Makobe Island cichlid gene stable ID', description=''>
pnyererei_homolog_associated_gene_name | <biomart.Attribute name='pnyererei_homolog_associated_gene_name', display_name='Makobe Island cichlid gene name', description=''>
pnyererei_homolog_ensembl_peptide | <biomart.Attribute name='pnyererei_homolog_ensembl_peptide', display_name='Makobe Island cichlid protein or transcript stable ID', description=''>
pnyererei_homolog_chromosome | <biomart.Attribute name='pnyererei_homolog_chromosome', display_name='Makobe Island cichlid chromosome/scaffold name', description=''>
pnyererei_homolog_chrom_start | <biomart.Attribute name='pnyererei_homolog_chrom_start', display_name='Makobe Island cichlid chromosome/scaffold start (bp)', description=''>
pnyererei_homolog_chrom_end | <biomart.Attribute name='pnyererei_homolog_chrom_end', display_name='Makobe Island cichlid chromosome/scaffold end (bp)', description=''>
pnyererei_homolog_canonical_transcript_protein | <biomart.Attribute name='pnyererei_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pnyererei_homolog_subtype | <biomart.Attribute name='pnyererei_homolog_subtype', display_name='Last common ancestor with Makobe Island cichlid', description=''>
pnyererei_homolog_orthology_type | <biomart.Attribute name='pnyererei_homolog_orthology_type', display_name='Makobe Island cichlid homology type', description=''>
pnyererei_homolog_perc_id | <biomart.Attribute name='pnyererei_homolog_perc_id', display_name='%id. target Makobe Island cichlid gene identical to query gene', description=''>
pnyererei_homolog_perc_id_r1 | <biomart.Attribute name='pnyererei_homolog_perc_id_r1', display_name='%id. query gene identical to target Makobe Island cichlid gene', description=''>
pnyererei_homolog_goc_score | <biomart.Attribute name='pnyererei_homolog_goc_score', display_name='Makobe Island cichlid Gene-order conservation score', description=''>
pnyererei_homolog_wga_coverage | <biomart.Attribute name='pnyererei_homolog_wga_coverage', display_name='Makobe Island cichlid Whole-genome alignment coverage', description=''>
pnyererei_homolog_orthology_confidence | <biomart.Attribute name='pnyererei_homolog_orthology_confidence', display_name='Makobe Island cichlid orthology confidence [0 low, 1 high]', description=''>
aplatyrhynchos_homolog_ensembl_gene | <biomart.Attribute name='aplatyrhynchos_homolog_ensembl_gene', display_name='Mallard gene stable ID', description=''>
aplatyrhynchos_homolog_associated_gene_name | <biomart.Attribute name='aplatyrhynchos_homolog_associated_gene_name', display_name='Mallard gene name', description=''>
aplatyrhynchos_homolog_ensembl_peptide | <biomart.Attribute name='aplatyrhynchos_homolog_ensembl_peptide', display_name='Mallard protein or transcript stable ID', description=''>
aplatyrhynchos_homolog_chromosome | <biomart.Attribute name='aplatyrhynchos_homolog_chromosome', display_name='Mallard chromosome/scaffold name', description=''>
aplatyrhynchos_homolog_chrom_start | <biomart.Attribute name='aplatyrhynchos_homolog_chrom_start', display_name='Mallard chromosome/scaffold start (bp)', description=''>
aplatyrhynchos_homolog_chrom_end | <biomart.Attribute name='aplatyrhynchos_homolog_chrom_end', display_name='Mallard chromosome/scaffold end (bp)', description=''>
aplatyrhynchos_homolog_canonical_transcript_protein | <biomart.Attribute name='aplatyrhynchos_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
aplatyrhynchos_homolog_subtype | <biomart.Attribute name='aplatyrhynchos_homolog_subtype', display_name='Last common ancestor with Mallard', description=''>
aplatyrhynchos_homolog_orthology_type | <biomart.Attribute name='aplatyrhynchos_homolog_orthology_type', display_name='Mallard homology type', description=''>
aplatyrhynchos_homolog_perc_id | <biomart.Attribute name='aplatyrhynchos_homolog_perc_id', display_name='%id. target Mallard gene identical to query gene', description=''>
aplatyrhynchos_homolog_perc_id_r1 | <biomart.Attribute name='aplatyrhynchos_homolog_perc_id_r1', display_name='%id. query gene identical to target Mallard gene', description=''>
aplatyrhynchos_homolog_goc_score | <biomart.Attribute name='aplatyrhynchos_homolog_goc_score', display_name='Mallard Gene-order conservation score', description=''>
aplatyrhynchos_homolog_wga_coverage | <biomart.Attribute name='aplatyrhynchos_homolog_wga_coverage', display_name='Mallard Whole-genome alignment coverage', description=''>
aplatyrhynchos_homolog_orthology_confidence | <biomart.Attribute name='aplatyrhynchos_homolog_orthology_confidence', display_name='Mallard orthology confidence [0 low, 1 high]', description=''>
cjacchus_homolog_ensembl_gene | <biomart.Attribute name='cjacchus_homolog_ensembl_gene', display_name='Marmoset gene stable ID', description=''>
cjacchus_homolog_associated_gene_name | <biomart.Attribute name='cjacchus_homolog_associated_gene_name', display_name='Marmoset gene name', description=''>
cjacchus_homolog_ensembl_peptide | <biomart.Attribute name='cjacchus_homolog_ensembl_peptide', display_name='Marmoset protein or transcript stable ID', description=''>
cjacchus_homolog_chromosome | <biomart.Attribute name='cjacchus_homolog_chromosome', display_name='Marmoset chromosome/scaffold name', description=''>
cjacchus_homolog_chrom_start | <biomart.Attribute name='cjacchus_homolog_chrom_start', display_name='Marmoset chromosome/scaffold start (bp)', description=''>
cjacchus_homolog_chrom_end | <biomart.Attribute name='cjacchus_homolog_chrom_end', display_name='Marmoset chromosome/scaffold end (bp)', description=''>
cjacchus_homolog_canonical_transcript_protein | <biomart.Attribute name='cjacchus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cjacchus_homolog_subtype | <biomart.Attribute name='cjacchus_homolog_subtype', display_name='Last common ancestor with Marmoset', description=''>
cjacchus_homolog_orthology_type | <biomart.Attribute name='cjacchus_homolog_orthology_type', display_name='Marmoset homology type', description=''>
cjacchus_homolog_perc_id | <biomart.Attribute name='cjacchus_homolog_perc_id', display_name='%id. target Marmoset gene identical to query gene', description=''>
cjacchus_homolog_perc_id_r1 | <biomart.Attribute name='cjacchus_homolog_perc_id_r1', display_name='%id. query gene identical to target Marmoset gene', description=''>
cjacchus_homolog_goc_score | <biomart.Attribute name='cjacchus_homolog_goc_score', display_name='Marmoset Gene-order conservation score', description=''>
cjacchus_homolog_wga_coverage | <biomart.Attribute name='cjacchus_homolog_wga_coverage', display_name='Marmoset Whole-genome alignment coverage', description=''>
cjacchus_homolog_orthology_confidence | <biomart.Attribute name='cjacchus_homolog_orthology_confidence', display_name='Marmoset orthology confidence [0 low, 1 high]', description=''>
pvampyrus_homolog_ensembl_gene | <biomart.Attribute name='pvampyrus_homolog_ensembl_gene', display_name='Megabat gene stable ID', description=''>
pvampyrus_homolog_associated_gene_name | <biomart.Attribute name='pvampyrus_homolog_associated_gene_name', display_name='Megabat gene name', description=''>
pvampyrus_homolog_ensembl_peptide | <biomart.Attribute name='pvampyrus_homolog_ensembl_peptide', display_name='Megabat protein or transcript stable ID', description=''>
pvampyrus_homolog_chromosome | <biomart.Attribute name='pvampyrus_homolog_chromosome', display_name='Megabat chromosome/scaffold name', description=''>
pvampyrus_homolog_chrom_start | <biomart.Attribute name='pvampyrus_homolog_chrom_start', display_name='Megabat chromosome/scaffold start (bp)', description=''>
pvampyrus_homolog_chrom_end | <biomart.Attribute name='pvampyrus_homolog_chrom_end', display_name='Megabat chromosome/scaffold end (bp)', description=''>
pvampyrus_homolog_canonical_transcript_protein | <biomart.Attribute name='pvampyrus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pvampyrus_homolog_subtype | <biomart.Attribute name='pvampyrus_homolog_subtype', display_name='Last common ancestor with Megabat', description=''>
pvampyrus_homolog_orthology_type | <biomart.Attribute name='pvampyrus_homolog_orthology_type', display_name='Megabat homology type', description=''>
pvampyrus_homolog_perc_id | <biomart.Attribute name='pvampyrus_homolog_perc_id', display_name='%id. target Megabat gene identical to query gene', description=''>
pvampyrus_homolog_perc_id_r1 | <biomart.Attribute name='pvampyrus_homolog_perc_id_r1', display_name='%id. query gene identical to target Megabat gene', description=''>
pvampyrus_homolog_goc_score | <biomart.Attribute name='pvampyrus_homolog_goc_score', display_name='Megabat Gene-order conservation score', description=''>
pvampyrus_homolog_wga_coverage | <biomart.Attribute name='pvampyrus_homolog_wga_coverage', display_name='Megabat Whole-genome alignment coverage', description=''>
pvampyrus_homolog_orthology_confidence | <biomart.Attribute name='pvampyrus_homolog_orthology_confidence', display_name='Megabat orthology confidence [0 low, 1 high]', description=''>
amexicanus_homolog_ensembl_gene | <biomart.Attribute name='amexicanus_homolog_ensembl_gene', display_name='Mexican tetra gene stable ID', description=''>
amexicanus_homolog_associated_gene_name | <biomart.Attribute name='amexicanus_homolog_associated_gene_name', display_name='Mexican tetra gene name', description=''>
amexicanus_homolog_ensembl_peptide | <biomart.Attribute name='amexicanus_homolog_ensembl_peptide', display_name='Mexican tetra protein or transcript stable ID', description=''>
amexicanus_homolog_chromosome | <biomart.Attribute name='amexicanus_homolog_chromosome', display_name='Mexican tetra chromosome/scaffold name', description=''>
amexicanus_homolog_chrom_start | <biomart.Attribute name='amexicanus_homolog_chrom_start', display_name='Mexican tetra chromosome/scaffold start (bp)', description=''>
amexicanus_homolog_chrom_end | <biomart.Attribute name='amexicanus_homolog_chrom_end', display_name='Mexican tetra chromosome/scaffold end (bp)', description=''>
amexicanus_homolog_canonical_transcript_protein | <biomart.Attribute name='amexicanus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
amexicanus_homolog_subtype | <biomart.Attribute name='amexicanus_homolog_subtype', display_name='Last common ancestor with Mexican tetra', description=''>
amexicanus_homolog_orthology_type | <biomart.Attribute name='amexicanus_homolog_orthology_type', display_name='Mexican tetra homology type', description=''>
amexicanus_homolog_perc_id | <biomart.Attribute name='amexicanus_homolog_perc_id', display_name='%id. target Mexican tetra gene identical to query gene', description=''>
amexicanus_homolog_perc_id_r1 | <biomart.Attribute name='amexicanus_homolog_perc_id_r1', display_name='%id. query gene identical to target Mexican tetra gene', description=''>
amexicanus_homolog_goc_score | <biomart.Attribute name='amexicanus_homolog_goc_score', display_name='Mexican tetra Gene-order conservation score', description=''>
amexicanus_homolog_wga_coverage | <biomart.Attribute name='amexicanus_homolog_wga_coverage', display_name='Mexican tetra Whole-genome alignment coverage', description=''>
amexicanus_homolog_orthology_confidence | <biomart.Attribute name='amexicanus_homolog_orthology_confidence', display_name='Mexican tetra orthology confidence [0 low, 1 high]', description=''>
mlucifugus_homolog_ensembl_gene | <biomart.Attribute name='mlucifugus_homolog_ensembl_gene', display_name='Microbat gene stable ID', description=''>
mlucifugus_homolog_associated_gene_name | <biomart.Attribute name='mlucifugus_homolog_associated_gene_name', display_name='Microbat gene name', description=''>
mlucifugus_homolog_ensembl_peptide | <biomart.Attribute name='mlucifugus_homolog_ensembl_peptide', display_name='Microbat protein or transcript stable ID', description=''>
mlucifugus_homolog_chromosome | <biomart.Attribute name='mlucifugus_homolog_chromosome', display_name='Microbat chromosome/scaffold name', description=''>
mlucifugus_homolog_chrom_start | <biomart.Attribute name='mlucifugus_homolog_chrom_start', display_name='Microbat chromosome/scaffold start (bp)', description=''>
mlucifugus_homolog_chrom_end | <biomart.Attribute name='mlucifugus_homolog_chrom_end', display_name='Microbat chromosome/scaffold end (bp)', description=''>
mlucifugus_homolog_canonical_transcript_protein | <biomart.Attribute name='mlucifugus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mlucifugus_homolog_subtype | <biomart.Attribute name='mlucifugus_homolog_subtype', display_name='Last common ancestor with Microbat', description=''>
mlucifugus_homolog_orthology_type | <biomart.Attribute name='mlucifugus_homolog_orthology_type', display_name='Microbat homology type', description=''>
mlucifugus_homolog_perc_id | <biomart.Attribute name='mlucifugus_homolog_perc_id', display_name='%id. target Microbat gene identical to query gene', description=''>
mlucifugus_homolog_perc_id_r1 | <biomart.Attribute name='mlucifugus_homolog_perc_id_r1', display_name='%id. query gene identical to target Microbat gene', description=''>
mlucifugus_homolog_goc_score | <biomart.Attribute name='mlucifugus_homolog_goc_score', display_name='Microbat Gene-order conservation score', description=''>
mlucifugus_homolog_wga_coverage | <biomart.Attribute name='mlucifugus_homolog_wga_coverage', display_name='Microbat Whole-genome alignment coverage', description=''>
mlucifugus_homolog_orthology_confidence | <biomart.Attribute name='mlucifugus_homolog_orthology_confidence', display_name='Microbat orthology confidence [0 low, 1 high]', description=''>
acitrinellus_homolog_ensembl_gene | <biomart.Attribute name='acitrinellus_homolog_ensembl_gene', display_name='Midas cichlid gene stable ID', description=''>
acitrinellus_homolog_associated_gene_name | <biomart.Attribute name='acitrinellus_homolog_associated_gene_name', display_name='Midas cichlid gene name', description=''>
acitrinellus_homolog_ensembl_peptide | <biomart.Attribute name='acitrinellus_homolog_ensembl_peptide', display_name='Midas cichlid protein or transcript stable ID', description=''>
acitrinellus_homolog_chromosome | <biomart.Attribute name='acitrinellus_homolog_chromosome', display_name='Midas cichlid chromosome/scaffold name', description=''>
acitrinellus_homolog_chrom_start | <biomart.Attribute name='acitrinellus_homolog_chrom_start', display_name='Midas cichlid chromosome/scaffold start (bp)', description=''>
acitrinellus_homolog_chrom_end | <biomart.Attribute name='acitrinellus_homolog_chrom_end', display_name='Midas cichlid chromosome/scaffold end (bp)', description=''>
acitrinellus_homolog_canonical_transcript_protein | <biomart.Attribute name='acitrinellus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
acitrinellus_homolog_subtype | <biomart.Attribute name='acitrinellus_homolog_subtype', display_name='Last common ancestor with Midas cichlid', description=''>
acitrinellus_homolog_orthology_type | <biomart.Attribute name='acitrinellus_homolog_orthology_type', display_name='Midas cichlid homology type', description=''>
acitrinellus_homolog_perc_id | <biomart.Attribute name='acitrinellus_homolog_perc_id', display_name='%id. target Midas cichlid gene identical to query gene', description=''>
acitrinellus_homolog_perc_id_r1 | <biomart.Attribute name='acitrinellus_homolog_perc_id_r1', display_name='%id. query gene identical to target Midas cichlid gene', description=''>
acitrinellus_homolog_goc_score | <biomart.Attribute name='acitrinellus_homolog_goc_score', display_name='Midas cichlid Gene-order conservation score', description=''>
acitrinellus_homolog_wga_coverage | <biomart.Attribute name='acitrinellus_homolog_wga_coverage', display_name='Midas cichlid Whole-genome alignment coverage', description=''>
acitrinellus_homolog_orthology_confidence | <biomart.Attribute name='acitrinellus_homolog_orthology_confidence', display_name='Midas cichlid orthology confidence [0 low, 1 high]', description=''>
munguiculatus_homolog_ensembl_gene | <biomart.Attribute name='munguiculatus_homolog_ensembl_gene', display_name='Mongolian gerbil gene stable ID', description=''>
munguiculatus_homolog_associated_gene_name | <biomart.Attribute name='munguiculatus_homolog_associated_gene_name', display_name='Mongolian gerbil gene name', description=''>
munguiculatus_homolog_ensembl_peptide | <biomart.Attribute name='munguiculatus_homolog_ensembl_peptide', display_name='Mongolian gerbil protein or transcript stable ID', description=''>
munguiculatus_homolog_chromosome | <biomart.Attribute name='munguiculatus_homolog_chromosome', display_name='Mongolian gerbil chromosome/scaffold name', description=''>
munguiculatus_homolog_chrom_start | <biomart.Attribute name='munguiculatus_homolog_chrom_start', display_name='Mongolian gerbil chromosome/scaffold start (bp)', description=''>
munguiculatus_homolog_chrom_end | <biomart.Attribute name='munguiculatus_homolog_chrom_end', display_name='Mongolian gerbil chromosome/scaffold end (bp)', description=''>
munguiculatus_homolog_canonical_transcript_protein | <biomart.Attribute name='munguiculatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
munguiculatus_homolog_subtype | <biomart.Attribute name='munguiculatus_homolog_subtype', display_name='Last common ancestor with Mongolian gerbil', description=''>
munguiculatus_homolog_orthology_type | <biomart.Attribute name='munguiculatus_homolog_orthology_type', display_name='Mongolian gerbil homology type', description=''>
munguiculatus_homolog_perc_id | <biomart.Attribute name='munguiculatus_homolog_perc_id', display_name='%id. target Mongolian gerbil gene identical to query gene', description=''>
munguiculatus_homolog_perc_id_r1 | <biomart.Attribute name='munguiculatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Mongolian gerbil gene', description=''>
munguiculatus_homolog_goc_score | <biomart.Attribute name='munguiculatus_homolog_goc_score', display_name='Mongolian gerbil Gene-order conservation score', description=''>
munguiculatus_homolog_wga_coverage | <biomart.Attribute name='munguiculatus_homolog_wga_coverage', display_name='Mongolian gerbil Whole-genome alignment coverage', description=''>
munguiculatus_homolog_orthology_confidence | <biomart.Attribute name='munguiculatus_homolog_orthology_confidence', display_name='Mongolian gerbil orthology confidence [0 low, 1 high]', description=''>
mmusculus_homolog_ensembl_gene | <biomart.Attribute name='mmusculus_homolog_ensembl_gene', display_name='Mouse gene stable ID', description=''>
mmusculus_homolog_associated_gene_name | <biomart.Attribute name='mmusculus_homolog_associated_gene_name', display_name='Mouse gene name', description=''>
mmusculus_homolog_ensembl_peptide | <biomart.Attribute name='mmusculus_homolog_ensembl_peptide', display_name='Mouse protein or transcript stable ID', description=''>
mmusculus_homolog_chromosome | <biomart.Attribute name='mmusculus_homolog_chromosome', display_name='Mouse chromosome/scaffold name', description=''>
mmusculus_homolog_chrom_start | <biomart.Attribute name='mmusculus_homolog_chrom_start', display_name='Mouse chromosome/scaffold start (bp)', description=''>
mmusculus_homolog_chrom_end | <biomart.Attribute name='mmusculus_homolog_chrom_end', display_name='Mouse chromosome/scaffold end (bp)', description=''>
mmusculus_homolog_canonical_transcript_protein | <biomart.Attribute name='mmusculus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mmusculus_homolog_subtype | <biomart.Attribute name='mmusculus_homolog_subtype', display_name='Last common ancestor with Mouse', description=''>
mmusculus_homolog_orthology_type | <biomart.Attribute name='mmusculus_homolog_orthology_type', display_name='Mouse homology type', description=''>
mmusculus_homolog_perc_id | <biomart.Attribute name='mmusculus_homolog_perc_id', display_name='%id. target Mouse gene identical to query gene', description=''>
mmusculus_homolog_perc_id_r1 | <biomart.Attribute name='mmusculus_homolog_perc_id_r1', display_name='%id. query gene identical to target Mouse gene', description=''>
mmusculus_homolog_goc_score | <biomart.Attribute name='mmusculus_homolog_goc_score', display_name='Mouse Gene-order conservation score', description=''>
mmusculus_homolog_wga_coverage | <biomart.Attribute name='mmusculus_homolog_wga_coverage', display_name='Mouse Whole-genome alignment coverage', description=''>
mmusculus_homolog_orthology_confidence | <biomart.Attribute name='mmusculus_homolog_orthology_confidence', display_name='Mouse orthology confidence [0 low, 1 high]', description=''>
mmurinus_homolog_ensembl_gene | <biomart.Attribute name='mmurinus_homolog_ensembl_gene', display_name='Mouse Lemur gene stable ID', description=''>
mmurinus_homolog_associated_gene_name | <biomart.Attribute name='mmurinus_homolog_associated_gene_name', display_name='Mouse Lemur gene name', description=''>
mmurinus_homolog_ensembl_peptide | <biomart.Attribute name='mmurinus_homolog_ensembl_peptide', display_name='Mouse Lemur protein or transcript stable ID', description=''>
mmurinus_homolog_chromosome | <biomart.Attribute name='mmurinus_homolog_chromosome', display_name='Mouse Lemur chromosome/scaffold name', description=''>
mmurinus_homolog_chrom_start | <biomart.Attribute name='mmurinus_homolog_chrom_start', display_name='Mouse Lemur chromosome/scaffold start (bp)', description=''>
mmurinus_homolog_chrom_end | <biomart.Attribute name='mmurinus_homolog_chrom_end', display_name='Mouse Lemur chromosome/scaffold end (bp)', description=''>
mmurinus_homolog_canonical_transcript_protein | <biomart.Attribute name='mmurinus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mmurinus_homolog_subtype | <biomart.Attribute name='mmurinus_homolog_subtype', display_name='Last common ancestor with Mouse Lemur', description=''>
mmurinus_homolog_orthology_type | <biomart.Attribute name='mmurinus_homolog_orthology_type', display_name='Mouse Lemur homology type', description=''>
mmurinus_homolog_perc_id | <biomart.Attribute name='mmurinus_homolog_perc_id', display_name='%id. target Mouse Lemur gene identical to query gene', description=''>
mmurinus_homolog_perc_id_r1 | <biomart.Attribute name='mmurinus_homolog_perc_id_r1', display_name='%id. query gene identical to target Mouse Lemur gene', description=''>
mmurinus_homolog_goc_score | <biomart.Attribute name='mmurinus_homolog_goc_score', display_name='Mouse Lemur Gene-order conservation score', description=''>
mmurinus_homolog_wga_coverage | <biomart.Attribute name='mmurinus_homolog_wga_coverage', display_name='Mouse Lemur Whole-genome alignment coverage', description=''>
mmurinus_homolog_orthology_confidence | <biomart.Attribute name='mmurinus_homolog_orthology_confidence', display_name='Mouse Lemur orthology confidence [0 low, 1 high]', description=''>
fheteroclitus_homolog_ensembl_gene | <biomart.Attribute name='fheteroclitus_homolog_ensembl_gene', display_name='Mummichog gene stable ID', description=''>
fheteroclitus_homolog_associated_gene_name | <biomart.Attribute name='fheteroclitus_homolog_associated_gene_name', display_name='Mummichog gene name', description=''>
fheteroclitus_homolog_ensembl_peptide | <biomart.Attribute name='fheteroclitus_homolog_ensembl_peptide', display_name='Mummichog protein or transcript stable ID', description=''>
fheteroclitus_homolog_chromosome | <biomart.Attribute name='fheteroclitus_homolog_chromosome', display_name='Mummichog chromosome/scaffold name', description=''>
fheteroclitus_homolog_chrom_start | <biomart.Attribute name='fheteroclitus_homolog_chrom_start', display_name='Mummichog chromosome/scaffold start (bp)', description=''>
fheteroclitus_homolog_chrom_end | <biomart.Attribute name='fheteroclitus_homolog_chrom_end', display_name='Mummichog chromosome/scaffold end (bp)', description=''>
fheteroclitus_homolog_canonical_transcript_protein | <biomart.Attribute name='fheteroclitus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
fheteroclitus_homolog_subtype | <biomart.Attribute name='fheteroclitus_homolog_subtype', display_name='Last common ancestor with Mummichog', description=''>
fheteroclitus_homolog_orthology_type | <biomart.Attribute name='fheteroclitus_homolog_orthology_type', display_name='Mummichog homology type', description=''>
fheteroclitus_homolog_perc_id | <biomart.Attribute name='fheteroclitus_homolog_perc_id', display_name='%id. target Mummichog gene identical to query gene', description=''>
fheteroclitus_homolog_perc_id_r1 | <biomart.Attribute name='fheteroclitus_homolog_perc_id_r1', display_name='%id. query gene identical to target Mummichog gene', description=''>
fheteroclitus_homolog_goc_score | <biomart.Attribute name='fheteroclitus_homolog_goc_score', display_name='Mummichog Gene-order conservation score', description=''>
fheteroclitus_homolog_wga_coverage | <biomart.Attribute name='fheteroclitus_homolog_wga_coverage', display_name='Mummichog Whole-genome alignment coverage', description=''>
fheteroclitus_homolog_orthology_confidence | <biomart.Attribute name='fheteroclitus_homolog_orthology_confidence', display_name='Mummichog orthology confidence [0 low, 1 high]', description=''>
hgfemale_homolog_ensembl_gene | <biomart.Attribute name='hgfemale_homolog_ensembl_gene', display_name='Naked mole-rat female gene stable ID', description=''>
hgfemale_homolog_associated_gene_name | <biomart.Attribute name='hgfemale_homolog_associated_gene_name', display_name='Naked mole-rat female gene name', description=''>
hgfemale_homolog_ensembl_peptide | <biomart.Attribute name='hgfemale_homolog_ensembl_peptide', display_name='Naked mole-rat female protein or transcript stable ID', description=''>
hgfemale_homolog_chromosome | <biomart.Attribute name='hgfemale_homolog_chromosome', display_name='Naked mole-rat female chromosome/scaffold name', description=''>
hgfemale_homolog_chrom_start | <biomart.Attribute name='hgfemale_homolog_chrom_start', display_name='Naked mole-rat female chromosome/scaffold start (bp)', description=''>
hgfemale_homolog_chrom_end | <biomart.Attribute name='hgfemale_homolog_chrom_end', display_name='Naked mole-rat female chromosome/scaffold end (bp)', description=''>
hgfemale_homolog_canonical_transcript_protein | <biomart.Attribute name='hgfemale_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
hgfemale_homolog_subtype | <biomart.Attribute name='hgfemale_homolog_subtype', display_name='Last common ancestor with Naked mole-rat female', description=''>
hgfemale_homolog_orthology_type | <biomart.Attribute name='hgfemale_homolog_orthology_type', display_name='Naked mole-rat female homology type', description=''>
hgfemale_homolog_perc_id | <biomart.Attribute name='hgfemale_homolog_perc_id', display_name='%id. target Naked mole-rat female gene identical to query gene', description=''>
hgfemale_homolog_perc_id_r1 | <biomart.Attribute name='hgfemale_homolog_perc_id_r1', display_name='%id. query gene identical to target Naked mole-rat female gene', description=''>
hgfemale_homolog_goc_score | <biomart.Attribute name='hgfemale_homolog_goc_score', display_name='Naked mole-rat female Gene-order conservation score', description=''>
hgfemale_homolog_wga_coverage | <biomart.Attribute name='hgfemale_homolog_wga_coverage', display_name='Naked mole-rat female Whole-genome alignment coverage', description=''>
hgfemale_homolog_orthology_confidence | <biomart.Attribute name='hgfemale_homolog_orthology_confidence', display_name='Naked mole-rat female orthology confidence [0 low, 1 high]', description=''>
hgmale_homolog_ensembl_gene | <biomart.Attribute name='hgmale_homolog_ensembl_gene', display_name='Naked mole-rat male gene stable ID', description=''>
hgmale_homolog_associated_gene_name | <biomart.Attribute name='hgmale_homolog_associated_gene_name', display_name='Naked mole-rat male gene name', description=''>
hgmale_homolog_ensembl_peptide | <biomart.Attribute name='hgmale_homolog_ensembl_peptide', display_name='Naked mole-rat male protein or transcript stable ID', description=''>
hgmale_homolog_chromosome | <biomart.Attribute name='hgmale_homolog_chromosome', display_name='Naked mole-rat male chromosome/scaffold name', description=''>
hgmale_homolog_chrom_start | <biomart.Attribute name='hgmale_homolog_chrom_start', display_name='Naked mole-rat male chromosome/scaffold start (bp)', description=''>
hgmale_homolog_chrom_end | <biomart.Attribute name='hgmale_homolog_chrom_end', display_name='Naked mole-rat male chromosome/scaffold end (bp)', description=''>
hgmale_homolog_canonical_transcript_protein | <biomart.Attribute name='hgmale_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
hgmale_homolog_subtype | <biomart.Attribute name='hgmale_homolog_subtype', display_name='Last common ancestor with Naked mole-rat male', description=''>
hgmale_homolog_orthology_type | <biomart.Attribute name='hgmale_homolog_orthology_type', display_name='Naked mole-rat male homology type', description=''>
hgmale_homolog_perc_id | <biomart.Attribute name='hgmale_homolog_perc_id', display_name='%id. target Naked mole-rat male gene identical to query gene', description=''>
hgmale_homolog_perc_id_r1 | <biomart.Attribute name='hgmale_homolog_perc_id_r1', display_name='%id. query gene identical to target Naked mole-rat male gene', description=''>
hgmale_homolog_goc_score | <biomart.Attribute name='hgmale_homolog_goc_score', display_name='Naked mole-rat male Gene-order conservation score', description=''>
hgmale_homolog_wga_coverage | <biomart.Attribute name='hgmale_homolog_wga_coverage', display_name='Naked mole-rat male Whole-genome alignment coverage', description=''>
hgmale_homolog_orthology_confidence | <biomart.Attribute name='hgmale_homolog_orthology_confidence', display_name='Naked mole-rat male orthology confidence [0 low, 1 high]', description=''>
oniloticus_homolog_ensembl_gene | <biomart.Attribute name='oniloticus_homolog_ensembl_gene', display_name='Nile tilapia gene stable ID', description=''>
oniloticus_homolog_associated_gene_name | <biomart.Attribute name='oniloticus_homolog_associated_gene_name', display_name='Nile tilapia gene name', description=''>
oniloticus_homolog_ensembl_peptide | <biomart.Attribute name='oniloticus_homolog_ensembl_peptide', display_name='Nile tilapia protein or transcript stable ID', description=''>
oniloticus_homolog_chromosome | <biomart.Attribute name='oniloticus_homolog_chromosome', display_name='Nile tilapia chromosome/scaffold name', description=''>
oniloticus_homolog_chrom_start | <biomart.Attribute name='oniloticus_homolog_chrom_start', display_name='Nile tilapia chromosome/scaffold start (bp)', description=''>
oniloticus_homolog_chrom_end | <biomart.Attribute name='oniloticus_homolog_chrom_end', display_name='Nile tilapia chromosome/scaffold end (bp)', description=''>
oniloticus_homolog_canonical_transcript_protein | <biomart.Attribute name='oniloticus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
oniloticus_homolog_subtype | <biomart.Attribute name='oniloticus_homolog_subtype', display_name='Last common ancestor with Nile tilapia', description=''>
oniloticus_homolog_orthology_type | <biomart.Attribute name='oniloticus_homolog_orthology_type', display_name='Nile tilapia homology type', description=''>
oniloticus_homolog_perc_id | <biomart.Attribute name='oniloticus_homolog_perc_id', display_name='%id. target Nile tilapia gene identical to query gene', description=''>
oniloticus_homolog_perc_id_r1 | <biomart.Attribute name='oniloticus_homolog_perc_id_r1', display_name='%id. query gene identical to target Nile tilapia gene', description=''>
oniloticus_homolog_goc_score | <biomart.Attribute name='oniloticus_homolog_goc_score', display_name='Nile tilapia Gene-order conservation score', description=''>
oniloticus_homolog_wga_coverage | <biomart.Attribute name='oniloticus_homolog_wga_coverage', display_name='Nile tilapia Whole-genome alignment coverage', description=''>
oniloticus_homolog_orthology_confidence | <biomart.Attribute name='oniloticus_homolog_orthology_confidence', display_name='Nile tilapia orthology confidence [0 low, 1 high]', description=''>
pmbairdii_homolog_ensembl_gene | <biomart.Attribute name='pmbairdii_homolog_ensembl_gene', display_name='Northern American deer mouse gene stable ID', description=''>
pmbairdii_homolog_associated_gene_name | <biomart.Attribute name='pmbairdii_homolog_associated_gene_name', display_name='Northern American deer mouse gene name', description=''>
pmbairdii_homolog_ensembl_peptide | <biomart.Attribute name='pmbairdii_homolog_ensembl_peptide', display_name='Northern American deer mouse protein or transcript stable ID', description=''>
pmbairdii_homolog_chromosome | <biomart.Attribute name='pmbairdii_homolog_chromosome', display_name='Northern American deer mouse chromosome/scaffold name', description=''>
pmbairdii_homolog_chrom_start | <biomart.Attribute name='pmbairdii_homolog_chrom_start', display_name='Northern American deer mouse chromosome/scaffold start (bp)', description=''>
pmbairdii_homolog_chrom_end | <biomart.Attribute name='pmbairdii_homolog_chrom_end', display_name='Northern American deer mouse chromosome/scaffold end (bp)', description=''>
pmbairdii_homolog_canonical_transcript_protein | <biomart.Attribute name='pmbairdii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pmbairdii_homolog_subtype | <biomart.Attribute name='pmbairdii_homolog_subtype', display_name='Last common ancestor with Northern American deer mouse', description=''>
pmbairdii_homolog_orthology_type | <biomart.Attribute name='pmbairdii_homolog_orthology_type', display_name='Northern American deer mouse homology type', description=''>
pmbairdii_homolog_perc_id | <biomart.Attribute name='pmbairdii_homolog_perc_id', display_name='%id. target Northern American deer mouse gene identical to query gene', description=''>
pmbairdii_homolog_perc_id_r1 | <biomart.Attribute name='pmbairdii_homolog_perc_id_r1', display_name='%id. query gene identical to target Northern American deer mouse gene', description=''>
pmbairdii_homolog_goc_score | <biomart.Attribute name='pmbairdii_homolog_goc_score', display_name='Northern American deer mouse Gene-order conservation score', description=''>
pmbairdii_homolog_wga_coverage | <biomart.Attribute name='pmbairdii_homolog_wga_coverage', display_name='Northern American deer mouse Whole-genome alignment coverage', description=''>
pmbairdii_homolog_orthology_confidence | <biomart.Attribute name='pmbairdii_homolog_orthology_confidence', display_name='Northern American deer mouse orthology confidence [0 low, 1 high]', description=''>
elucius_homolog_ensembl_gene | <biomart.Attribute name='elucius_homolog_ensembl_gene', display_name='Northern pike gene stable ID', description=''>
elucius_homolog_associated_gene_name | <biomart.Attribute name='elucius_homolog_associated_gene_name', display_name='Northern pike gene name', description=''>
elucius_homolog_ensembl_peptide | <biomart.Attribute name='elucius_homolog_ensembl_peptide', display_name='Northern pike protein or transcript stable ID', description=''>
elucius_homolog_chromosome | <biomart.Attribute name='elucius_homolog_chromosome', display_name='Northern pike chromosome/scaffold name', description=''>
elucius_homolog_chrom_start | <biomart.Attribute name='elucius_homolog_chrom_start', display_name='Northern pike chromosome/scaffold start (bp)', description=''>
elucius_homolog_chrom_end | <biomart.Attribute name='elucius_homolog_chrom_end', display_name='Northern pike chromosome/scaffold end (bp)', description=''>
elucius_homolog_canonical_transcript_protein | <biomart.Attribute name='elucius_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
elucius_homolog_subtype | <biomart.Attribute name='elucius_homolog_subtype', display_name='Last common ancestor with Northern pike', description=''>
elucius_homolog_orthology_type | <biomart.Attribute name='elucius_homolog_orthology_type', display_name='Northern pike homology type', description=''>
elucius_homolog_perc_id | <biomart.Attribute name='elucius_homolog_perc_id', display_name='%id. target Northern pike gene identical to query gene', description=''>
elucius_homolog_perc_id_r1 | <biomart.Attribute name='elucius_homolog_perc_id_r1', display_name='%id. query gene identical to target Northern pike gene', description=''>
elucius_homolog_goc_score | <biomart.Attribute name='elucius_homolog_goc_score', display_name='Northern pike Gene-order conservation score', description=''>
elucius_homolog_wga_coverage | <biomart.Attribute name='elucius_homolog_wga_coverage', display_name='Northern pike Whole-genome alignment coverage', description=''>
elucius_homolog_orthology_confidence | <biomart.Attribute name='elucius_homolog_orthology_confidence', display_name='Northern pike orthology confidence [0 low, 1 high]', description=''>
panubis_homolog_ensembl_gene | <biomart.Attribute name='panubis_homolog_ensembl_gene', display_name='Olive baboon gene stable ID', description=''>
panubis_homolog_associated_gene_name | <biomart.Attribute name='panubis_homolog_associated_gene_name', display_name='Olive baboon gene name', description=''>
panubis_homolog_ensembl_peptide | <biomart.Attribute name='panubis_homolog_ensembl_peptide', display_name='Olive baboon protein or transcript stable ID', description=''>
panubis_homolog_chromosome | <biomart.Attribute name='panubis_homolog_chromosome', display_name='Olive baboon chromosome/scaffold name', description=''>
panubis_homolog_chrom_start | <biomart.Attribute name='panubis_homolog_chrom_start', display_name='Olive baboon chromosome/scaffold start (bp)', description=''>
panubis_homolog_chrom_end | <biomart.Attribute name='panubis_homolog_chrom_end', display_name='Olive baboon chromosome/scaffold end (bp)', description=''>
panubis_homolog_canonical_transcript_protein | <biomart.Attribute name='panubis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
panubis_homolog_subtype | <biomart.Attribute name='panubis_homolog_subtype', display_name='Last common ancestor with Olive baboon', description=''>
panubis_homolog_orthology_type | <biomart.Attribute name='panubis_homolog_orthology_type', display_name='Olive baboon homology type', description=''>
panubis_homolog_perc_id | <biomart.Attribute name='panubis_homolog_perc_id', display_name='%id. target Olive baboon gene identical to query gene', description=''>
panubis_homolog_perc_id_r1 | <biomart.Attribute name='panubis_homolog_perc_id_r1', display_name='%id. query gene identical to target Olive baboon gene', description=''>
panubis_homolog_goc_score | <biomart.Attribute name='panubis_homolog_goc_score', display_name='Olive baboon Gene-order conservation score', description=''>
panubis_homolog_wga_coverage | <biomart.Attribute name='panubis_homolog_wga_coverage', display_name='Olive baboon Whole-genome alignment coverage', description=''>
panubis_homolog_orthology_confidence | <biomart.Attribute name='panubis_homolog_orthology_confidence', display_name='Olive baboon orthology confidence [0 low, 1 high]', description=''>
mdomestica_homolog_ensembl_gene | <biomart.Attribute name='mdomestica_homolog_ensembl_gene', display_name='Opossum gene stable ID', description=''>
mdomestica_homolog_associated_gene_name | <biomart.Attribute name='mdomestica_homolog_associated_gene_name', display_name='Opossum gene name', description=''>
mdomestica_homolog_ensembl_peptide | <biomart.Attribute name='mdomestica_homolog_ensembl_peptide', display_name='Opossum protein or transcript stable ID', description=''>
mdomestica_homolog_chromosome | <biomart.Attribute name='mdomestica_homolog_chromosome', display_name='Opossum chromosome/scaffold name', description=''>
mdomestica_homolog_chrom_start | <biomart.Attribute name='mdomestica_homolog_chrom_start', display_name='Opossum chromosome/scaffold start (bp)', description=''>
mdomestica_homolog_chrom_end | <biomart.Attribute name='mdomestica_homolog_chrom_end', display_name='Opossum chromosome/scaffold end (bp)', description=''>
mdomestica_homolog_canonical_transcript_protein | <biomart.Attribute name='mdomestica_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mdomestica_homolog_subtype | <biomart.Attribute name='mdomestica_homolog_subtype', display_name='Last common ancestor with Opossum', description=''>
mdomestica_homolog_orthology_type | <biomart.Attribute name='mdomestica_homolog_orthology_type', display_name='Opossum homology type', description=''>
mdomestica_homolog_perc_id | <biomart.Attribute name='mdomestica_homolog_perc_id', display_name='%id. target Opossum gene identical to query gene', description=''>
mdomestica_homolog_perc_id_r1 | <biomart.Attribute name='mdomestica_homolog_perc_id_r1', display_name='%id. query gene identical to target Opossum gene', description=''>
mdomestica_homolog_goc_score | <biomart.Attribute name='mdomestica_homolog_goc_score', display_name='Opossum Gene-order conservation score', description=''>
mdomestica_homolog_wga_coverage | <biomart.Attribute name='mdomestica_homolog_wga_coverage', display_name='Opossum Whole-genome alignment coverage', description=''>
mdomestica_homolog_orthology_confidence | <biomart.Attribute name='mdomestica_homolog_orthology_confidence', display_name='Opossum orthology confidence [0 low, 1 high]', description=''>
pabelii_homolog_ensembl_gene | <biomart.Attribute name='pabelii_homolog_ensembl_gene', display_name='Orangutan gene stable ID', description=''>
pabelii_homolog_associated_gene_name | <biomart.Attribute name='pabelii_homolog_associated_gene_name', display_name='Orangutan gene name', description=''>
pabelii_homolog_ensembl_peptide | <biomart.Attribute name='pabelii_homolog_ensembl_peptide', display_name='Orangutan protein or transcript stable ID', description=''>
pabelii_homolog_chromosome | <biomart.Attribute name='pabelii_homolog_chromosome', display_name='Orangutan chromosome/scaffold name', description=''>
pabelii_homolog_chrom_start | <biomart.Attribute name='pabelii_homolog_chrom_start', display_name='Orangutan chromosome/scaffold start (bp)', description=''>
pabelii_homolog_chrom_end | <biomart.Attribute name='pabelii_homolog_chrom_end', display_name='Orangutan chromosome/scaffold end (bp)', description=''>
pabelii_homolog_canonical_transcript_protein | <biomart.Attribute name='pabelii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pabelii_homolog_subtype | <biomart.Attribute name='pabelii_homolog_subtype', display_name='Last common ancestor with Orangutan', description=''>
pabelii_homolog_orthology_type | <biomart.Attribute name='pabelii_homolog_orthology_type', display_name='Orangutan homology type', description=''>
pabelii_homolog_perc_id | <biomart.Attribute name='pabelii_homolog_perc_id', display_name='%id. target Orangutan gene identical to query gene', description=''>
pabelii_homolog_perc_id_r1 | <biomart.Attribute name='pabelii_homolog_perc_id_r1', display_name='%id. query gene identical to target Orangutan gene', description=''>
pabelii_homolog_goc_score | <biomart.Attribute name='pabelii_homolog_goc_score', display_name='Orangutan Gene-order conservation score', description=''>
pabelii_homolog_wga_coverage | <biomart.Attribute name='pabelii_homolog_wga_coverage', display_name='Orangutan Whole-genome alignment coverage', description=''>
pabelii_homolog_orthology_confidence | <biomart.Attribute name='pabelii_homolog_orthology_confidence', display_name='Orangutan orthology confidence [0 low, 1 high]', description=''>
sorbicularis_homolog_ensembl_gene | <biomart.Attribute name='sorbicularis_homolog_ensembl_gene', display_name='Orbiculate cardinalfish gene stable ID', description=''>
sorbicularis_homolog_associated_gene_name | <biomart.Attribute name='sorbicularis_homolog_associated_gene_name', display_name='Orbiculate cardinalfish gene name', description=''>
sorbicularis_homolog_ensembl_peptide | <biomart.Attribute name='sorbicularis_homolog_ensembl_peptide', display_name='Orbiculate cardinalfish protein or transcript stable ID', description=''>
sorbicularis_homolog_chromosome | <biomart.Attribute name='sorbicularis_homolog_chromosome', display_name='Orbiculate cardinalfish chromosome/scaffold name', description=''>
sorbicularis_homolog_chrom_start | <biomart.Attribute name='sorbicularis_homolog_chrom_start', display_name='Orbiculate cardinalfish chromosome/scaffold start (bp)', description=''>
sorbicularis_homolog_chrom_end | <biomart.Attribute name='sorbicularis_homolog_chrom_end', display_name='Orbiculate cardinalfish chromosome/scaffold end (bp)', description=''>
sorbicularis_homolog_canonical_transcript_protein | <biomart.Attribute name='sorbicularis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sorbicularis_homolog_subtype | <biomart.Attribute name='sorbicularis_homolog_subtype', display_name='Last common ancestor with Orbiculate cardinalfish', description=''>
sorbicularis_homolog_orthology_type | <biomart.Attribute name='sorbicularis_homolog_orthology_type', display_name='Orbiculate cardinalfish homology type', description=''>
sorbicularis_homolog_perc_id | <biomart.Attribute name='sorbicularis_homolog_perc_id', display_name='%id. target Orbiculate cardinalfish gene identical to query gene', description=''>
sorbicularis_homolog_perc_id_r1 | <biomart.Attribute name='sorbicularis_homolog_perc_id_r1', display_name='%id. query gene identical to target Orbiculate cardinalfish gene', description=''>
sorbicularis_homolog_goc_score | <biomart.Attribute name='sorbicularis_homolog_goc_score', display_name='Orbiculate cardinalfish Gene-order conservation score', description=''>
sorbicularis_homolog_wga_coverage | <biomart.Attribute name='sorbicularis_homolog_wga_coverage', display_name='Orbiculate cardinalfish Whole-genome alignment coverage', description=''>
sorbicularis_homolog_orthology_confidence | <biomart.Attribute name='sorbicularis_homolog_orthology_confidence', display_name='Orbiculate cardinalfish orthology confidence [0 low, 1 high]', description=''>
ampachon_homolog_ensembl_gene | <biomart.Attribute name='ampachon_homolog_ensembl_gene', display_name='Pachon cavefish gene stable ID', description=''>
ampachon_homolog_associated_gene_name | <biomart.Attribute name='ampachon_homolog_associated_gene_name', display_name='Pachon cavefish gene name', description=''>
ampachon_homolog_ensembl_peptide | <biomart.Attribute name='ampachon_homolog_ensembl_peptide', display_name='Pachon cavefish protein or transcript stable ID', description=''>
ampachon_homolog_chromosome | <biomart.Attribute name='ampachon_homolog_chromosome', display_name='Pachon cavefish chromosome/scaffold name', description=''>
ampachon_homolog_chrom_start | <biomart.Attribute name='ampachon_homolog_chrom_start', display_name='Pachon cavefish chromosome/scaffold start (bp)', description=''>
ampachon_homolog_chrom_end | <biomart.Attribute name='ampachon_homolog_chrom_end', display_name='Pachon cavefish chromosome/scaffold end (bp)', description=''>
ampachon_homolog_canonical_transcript_protein | <biomart.Attribute name='ampachon_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ampachon_homolog_subtype | <biomart.Attribute name='ampachon_homolog_subtype', display_name='Last common ancestor with Pachon cavefish', description=''>
ampachon_homolog_orthology_type | <biomart.Attribute name='ampachon_homolog_orthology_type', display_name='Pachon cavefish homology type', description=''>
ampachon_homolog_perc_id | <biomart.Attribute name='ampachon_homolog_perc_id', display_name='%id. target Pachon cavefish gene identical to query gene', description=''>
ampachon_homolog_perc_id_r1 | <biomart.Attribute name='ampachon_homolog_perc_id_r1', display_name='%id. query gene identical to target Pachon cavefish gene', description=''>
ampachon_homolog_goc_score | <biomart.Attribute name='ampachon_homolog_goc_score', display_name='Pachon cavefish Gene-order conservation score', description=''>
ampachon_homolog_wga_coverage | <biomart.Attribute name='ampachon_homolog_wga_coverage', display_name='Pachon cavefish Whole-genome alignment coverage', description=''>
ampachon_homolog_orthology_confidence | <biomart.Attribute name='ampachon_homolog_orthology_confidence', display_name='Pachon cavefish orthology confidence [0 low, 1 high]', description=''>
cpbellii_homolog_ensembl_gene | <biomart.Attribute name='cpbellii_homolog_ensembl_gene', display_name='Painted turtle gene stable ID', description=''>
cpbellii_homolog_associated_gene_name | <biomart.Attribute name='cpbellii_homolog_associated_gene_name', display_name='Painted turtle gene name', description=''>
cpbellii_homolog_ensembl_peptide | <biomart.Attribute name='cpbellii_homolog_ensembl_peptide', display_name='Painted turtle protein or transcript stable ID', description=''>
cpbellii_homolog_chromosome | <biomart.Attribute name='cpbellii_homolog_chromosome', display_name='Painted turtle chromosome/scaffold name', description=''>
cpbellii_homolog_chrom_start | <biomart.Attribute name='cpbellii_homolog_chrom_start', display_name='Painted turtle chromosome/scaffold start (bp)', description=''>
cpbellii_homolog_chrom_end | <biomart.Attribute name='cpbellii_homolog_chrom_end', display_name='Painted turtle chromosome/scaffold end (bp)', description=''>
cpbellii_homolog_canonical_transcript_protein | <biomart.Attribute name='cpbellii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cpbellii_homolog_subtype | <biomart.Attribute name='cpbellii_homolog_subtype', display_name='Last common ancestor with Painted turtle', description=''>
cpbellii_homolog_orthology_type | <biomart.Attribute name='cpbellii_homolog_orthology_type', display_name='Painted turtle homology type', description=''>
cpbellii_homolog_perc_id | <biomart.Attribute name='cpbellii_homolog_perc_id', display_name='%id. target Painted turtle gene identical to query gene', description=''>
cpbellii_homolog_perc_id_r1 | <biomart.Attribute name='cpbellii_homolog_perc_id_r1', display_name='%id. query gene identical to target Painted turtle gene', description=''>
cpbellii_homolog_goc_score | <biomart.Attribute name='cpbellii_homolog_goc_score', display_name='Painted turtle Gene-order conservation score', description=''>
cpbellii_homolog_wga_coverage | <biomart.Attribute name='cpbellii_homolog_wga_coverage', display_name='Painted turtle Whole-genome alignment coverage', description=''>
cpbellii_homolog_orthology_confidence | <biomart.Attribute name='cpbellii_homolog_orthology_confidence', display_name='Painted turtle orthology confidence [0 low, 1 high]', description=''>
amelanoleuca_homolog_ensembl_gene | <biomart.Attribute name='amelanoleuca_homolog_ensembl_gene', display_name='Panda gene stable ID', description=''>
amelanoleuca_homolog_associated_gene_name | <biomart.Attribute name='amelanoleuca_homolog_associated_gene_name', display_name='Panda gene name', description=''>
amelanoleuca_homolog_ensembl_peptide | <biomart.Attribute name='amelanoleuca_homolog_ensembl_peptide', display_name='Panda protein or transcript stable ID', description=''>
amelanoleuca_homolog_chromosome | <biomart.Attribute name='amelanoleuca_homolog_chromosome', display_name='Panda chromosome/scaffold name', description=''>
amelanoleuca_homolog_chrom_start | <biomart.Attribute name='amelanoleuca_homolog_chrom_start', display_name='Panda chromosome/scaffold start (bp)', description=''>
amelanoleuca_homolog_chrom_end | <biomart.Attribute name='amelanoleuca_homolog_chrom_end', display_name='Panda chromosome/scaffold end (bp)', description=''>
amelanoleuca_homolog_canonical_transcript_protein | <biomart.Attribute name='amelanoleuca_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
amelanoleuca_homolog_subtype | <biomart.Attribute name='amelanoleuca_homolog_subtype', display_name='Last common ancestor with Panda', description=''>
amelanoleuca_homolog_orthology_type | <biomart.Attribute name='amelanoleuca_homolog_orthology_type', display_name='Panda homology type', description=''>
amelanoleuca_homolog_perc_id | <biomart.Attribute name='amelanoleuca_homolog_perc_id', display_name='%id. target Panda gene identical to query gene', description=''>
amelanoleuca_homolog_perc_id_r1 | <biomart.Attribute name='amelanoleuca_homolog_perc_id_r1', display_name='%id. query gene identical to target Panda gene', description=''>
amelanoleuca_homolog_goc_score | <biomart.Attribute name='amelanoleuca_homolog_goc_score', display_name='Panda Gene-order conservation score', description=''>
amelanoleuca_homolog_wga_coverage | <biomart.Attribute name='amelanoleuca_homolog_wga_coverage', display_name='Panda Whole-genome alignment coverage', description=''>
amelanoleuca_homolog_orthology_confidence | <biomart.Attribute name='amelanoleuca_homolog_orthology_confidence', display_name='Panda orthology confidence [0 low, 1 high]', description=''>
pkingsleyae_homolog_ensembl_gene | <biomart.Attribute name='pkingsleyae_homolog_ensembl_gene', display_name='Paramormyrops kingsleyae gene stable ID', description=''>
pkingsleyae_homolog_associated_gene_name | <biomart.Attribute name='pkingsleyae_homolog_associated_gene_name', display_name='Paramormyrops kingsleyae gene name', description=''>
pkingsleyae_homolog_ensembl_peptide | <biomart.Attribute name='pkingsleyae_homolog_ensembl_peptide', display_name='Paramormyrops kingsleyae protein or transcript stable ID', description=''>
pkingsleyae_homolog_chromosome | <biomart.Attribute name='pkingsleyae_homolog_chromosome', display_name='Paramormyrops kingsleyae chromosome/scaffold name', description=''>
pkingsleyae_homolog_chrom_start | <biomart.Attribute name='pkingsleyae_homolog_chrom_start', display_name='Paramormyrops kingsleyae chromosome/scaffold start (bp)', description=''>
pkingsleyae_homolog_chrom_end | <biomart.Attribute name='pkingsleyae_homolog_chrom_end', display_name='Paramormyrops kingsleyae chromosome/scaffold end (bp)', description=''>
pkingsleyae_homolog_canonical_transcript_protein | <biomart.Attribute name='pkingsleyae_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pkingsleyae_homolog_subtype | <biomart.Attribute name='pkingsleyae_homolog_subtype', display_name='Last common ancestor with Paramormyrops kingsleyae', description=''>
pkingsleyae_homolog_orthology_type | <biomart.Attribute name='pkingsleyae_homolog_orthology_type', display_name='Paramormyrops kingsleyae homology type', description=''>
pkingsleyae_homolog_perc_id | <biomart.Attribute name='pkingsleyae_homolog_perc_id', display_name='%id. target Paramormyrops kingsleyae gene identical to query gene', description=''>
pkingsleyae_homolog_perc_id_r1 | <biomart.Attribute name='pkingsleyae_homolog_perc_id_r1', display_name='%id. query gene identical to target Paramormyrops kingsleyae gene', description=''>
pkingsleyae_homolog_goc_score | <biomart.Attribute name='pkingsleyae_homolog_goc_score', display_name='Paramormyrops kingsleyae Gene-order conservation score', description=''>
pkingsleyae_homolog_wga_coverage | <biomart.Attribute name='pkingsleyae_homolog_wga_coverage', display_name='Paramormyrops kingsleyae Whole-genome alignment coverage', description=''>
pkingsleyae_homolog_orthology_confidence | <biomart.Attribute name='pkingsleyae_homolog_orthology_confidence', display_name='Paramormyrops kingsleyae orthology confidence [0 low, 1 high]', description=''>
sscrofa_homolog_ensembl_gene | <biomart.Attribute name='sscrofa_homolog_ensembl_gene', display_name='Pig gene stable ID', description=''>
sscrofa_homolog_associated_gene_name | <biomart.Attribute name='sscrofa_homolog_associated_gene_name', display_name='Pig gene name', description=''>
sscrofa_homolog_ensembl_peptide | <biomart.Attribute name='sscrofa_homolog_ensembl_peptide', display_name='Pig protein or transcript stable ID', description=''>
sscrofa_homolog_chromosome | <biomart.Attribute name='sscrofa_homolog_chromosome', display_name='Pig chromosome/scaffold name', description=''>
sscrofa_homolog_chrom_start | <biomart.Attribute name='sscrofa_homolog_chrom_start', display_name='Pig chromosome/scaffold start (bp)', description=''>
sscrofa_homolog_chrom_end | <biomart.Attribute name='sscrofa_homolog_chrom_end', display_name='Pig chromosome/scaffold end (bp)', description=''>
sscrofa_homolog_canonical_transcript_protein | <biomart.Attribute name='sscrofa_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sscrofa_homolog_subtype | <biomart.Attribute name='sscrofa_homolog_subtype', display_name='Last common ancestor with Pig', description=''>
sscrofa_homolog_orthology_type | <biomart.Attribute name='sscrofa_homolog_orthology_type', display_name='Pig homology type', description=''>
sscrofa_homolog_perc_id | <biomart.Attribute name='sscrofa_homolog_perc_id', display_name='%id. target Pig gene identical to query gene', description=''>
sscrofa_homolog_perc_id_r1 | <biomart.Attribute name='sscrofa_homolog_perc_id_r1', display_name='%id. query gene identical to target Pig gene', description=''>
sscrofa_homolog_goc_score | <biomart.Attribute name='sscrofa_homolog_goc_score', display_name='Pig Gene-order conservation score', description=''>
sscrofa_homolog_wga_coverage | <biomart.Attribute name='sscrofa_homolog_wga_coverage', display_name='Pig Whole-genome alignment coverage', description=''>
sscrofa_homolog_orthology_confidence | <biomart.Attribute name='sscrofa_homolog_orthology_confidence', display_name='Pig orthology confidence [0 low, 1 high]', description=''>
mnemestrina_homolog_ensembl_gene | <biomart.Attribute name='mnemestrina_homolog_ensembl_gene', display_name='Pig-tailed macaque gene stable ID', description=''>
mnemestrina_homolog_associated_gene_name | <biomart.Attribute name='mnemestrina_homolog_associated_gene_name', display_name='Pig-tailed macaque gene name', description=''>
mnemestrina_homolog_ensembl_peptide | <biomart.Attribute name='mnemestrina_homolog_ensembl_peptide', display_name='Pig-tailed macaque protein or transcript stable ID', description=''>
mnemestrina_homolog_chromosome | <biomart.Attribute name='mnemestrina_homolog_chromosome', display_name='Pig-tailed macaque chromosome/scaffold name', description=''>
mnemestrina_homolog_chrom_start | <biomart.Attribute name='mnemestrina_homolog_chrom_start', display_name='Pig-tailed macaque chromosome/scaffold start (bp)', description=''>
mnemestrina_homolog_chrom_end | <biomart.Attribute name='mnemestrina_homolog_chrom_end', display_name='Pig-tailed macaque chromosome/scaffold end (bp)', description=''>
mnemestrina_homolog_canonical_transcript_protein | <biomart.Attribute name='mnemestrina_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mnemestrina_homolog_subtype | <biomart.Attribute name='mnemestrina_homolog_subtype', display_name='Last common ancestor with Pig-tailed macaque', description=''>
mnemestrina_homolog_orthology_type | <biomart.Attribute name='mnemestrina_homolog_orthology_type', display_name='Pig-tailed macaque homology type', description=''>
mnemestrina_homolog_perc_id | <biomart.Attribute name='mnemestrina_homolog_perc_id', display_name='%id. target Pig-tailed macaque gene identical to query gene', description=''>
mnemestrina_homolog_perc_id_r1 | <biomart.Attribute name='mnemestrina_homolog_perc_id_r1', display_name='%id. query gene identical to target Pig-tailed macaque gene', description=''>
mnemestrina_homolog_goc_score | <biomart.Attribute name='mnemestrina_homolog_goc_score', display_name='Pig-tailed macaque Gene-order conservation score', description=''>
mnemestrina_homolog_wga_coverage | <biomart.Attribute name='mnemestrina_homolog_wga_coverage', display_name='Pig-tailed macaque Whole-genome alignment coverage', description=''>
mnemestrina_homolog_orthology_confidence | <biomart.Attribute name='mnemestrina_homolog_orthology_confidence', display_name='Pig-tailed macaque orthology confidence [0 low, 1 high]', description=''>
oprinceps_homolog_ensembl_gene | <biomart.Attribute name='oprinceps_homolog_ensembl_gene', display_name='Pika gene stable ID', description=''>
oprinceps_homolog_associated_gene_name | <biomart.Attribute name='oprinceps_homolog_associated_gene_name', display_name='Pika gene name', description=''>
oprinceps_homolog_ensembl_peptide | <biomart.Attribute name='oprinceps_homolog_ensembl_peptide', display_name='Pika protein or transcript stable ID', description=''>
oprinceps_homolog_chromosome | <biomart.Attribute name='oprinceps_homolog_chromosome', display_name='Pika chromosome/scaffold name', description=''>
oprinceps_homolog_chrom_start | <biomart.Attribute name='oprinceps_homolog_chrom_start', display_name='Pika chromosome/scaffold start (bp)', description=''>
oprinceps_homolog_chrom_end | <biomart.Attribute name='oprinceps_homolog_chrom_end', display_name='Pika chromosome/scaffold end (bp)', description=''>
oprinceps_homolog_canonical_transcript_protein | <biomart.Attribute name='oprinceps_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
oprinceps_homolog_subtype | <biomart.Attribute name='oprinceps_homolog_subtype', display_name='Last common ancestor with Pika', description=''>
oprinceps_homolog_orthology_type | <biomart.Attribute name='oprinceps_homolog_orthology_type', display_name='Pika homology type', description=''>
oprinceps_homolog_perc_id | <biomart.Attribute name='oprinceps_homolog_perc_id', display_name='%id. target Pika gene identical to query gene', description=''>
oprinceps_homolog_perc_id_r1 | <biomart.Attribute name='oprinceps_homolog_perc_id_r1', display_name='%id. query gene identical to target Pika gene', description=''>
oprinceps_homolog_goc_score | <biomart.Attribute name='oprinceps_homolog_goc_score', display_name='Pika Gene-order conservation score', description=''>
oprinceps_homolog_wga_coverage | <biomart.Attribute name='oprinceps_homolog_wga_coverage', display_name='Pika Whole-genome alignment coverage', description=''>
oprinceps_homolog_orthology_confidence | <biomart.Attribute name='oprinceps_homolog_orthology_confidence', display_name='Pika orthology confidence [0 low, 1 high]', description=''>
mmurdjan_homolog_ensembl_gene | <biomart.Attribute name='mmurdjan_homolog_ensembl_gene', display_name='Pinecone soldierfish gene stable ID', description=''>
mmurdjan_homolog_associated_gene_name | <biomart.Attribute name='mmurdjan_homolog_associated_gene_name', display_name='Pinecone soldierfish gene name', description=''>
mmurdjan_homolog_ensembl_peptide | <biomart.Attribute name='mmurdjan_homolog_ensembl_peptide', display_name='Pinecone soldierfish protein or transcript stable ID', description=''>
mmurdjan_homolog_chromosome | <biomart.Attribute name='mmurdjan_homolog_chromosome', display_name='Pinecone soldierfish chromosome/scaffold name', description=''>
mmurdjan_homolog_chrom_start | <biomart.Attribute name='mmurdjan_homolog_chrom_start', display_name='Pinecone soldierfish chromosome/scaffold start (bp)', description=''>
mmurdjan_homolog_chrom_end | <biomart.Attribute name='mmurdjan_homolog_chrom_end', display_name='Pinecone soldierfish chromosome/scaffold end (bp)', description=''>
mmurdjan_homolog_canonical_transcript_protein | <biomart.Attribute name='mmurdjan_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mmurdjan_homolog_subtype | <biomart.Attribute name='mmurdjan_homolog_subtype', display_name='Last common ancestor with Pinecone soldierfish', description=''>
mmurdjan_homolog_orthology_type | <biomart.Attribute name='mmurdjan_homolog_orthology_type', display_name='Pinecone soldierfish homology type', description=''>
mmurdjan_homolog_perc_id | <biomart.Attribute name='mmurdjan_homolog_perc_id', display_name='%id. target Pinecone soldierfish gene identical to query gene', description=''>
mmurdjan_homolog_perc_id_r1 | <biomart.Attribute name='mmurdjan_homolog_perc_id_r1', display_name='%id. query gene identical to target Pinecone soldierfish gene', description=''>
mmurdjan_homolog_goc_score | <biomart.Attribute name='mmurdjan_homolog_goc_score', display_name='Pinecone soldierfish Gene-order conservation score', description=''>
mmurdjan_homolog_wga_coverage | <biomart.Attribute name='mmurdjan_homolog_wga_coverage', display_name='Pinecone soldierfish Whole-genome alignment coverage', description=''>
mmurdjan_homolog_orthology_confidence | <biomart.Attribute name='mmurdjan_homolog_orthology_confidence', display_name='Pinecone soldierfish orthology confidence [0 low, 1 high]', description=''>
xmaculatus_homolog_ensembl_gene | <biomart.Attribute name='xmaculatus_homolog_ensembl_gene', display_name='Platyfish gene stable ID', description=''>
xmaculatus_homolog_associated_gene_name | <biomart.Attribute name='xmaculatus_homolog_associated_gene_name', display_name='Platyfish gene name', description=''>
xmaculatus_homolog_ensembl_peptide | <biomart.Attribute name='xmaculatus_homolog_ensembl_peptide', display_name='Platyfish protein or transcript stable ID', description=''>
xmaculatus_homolog_chromosome | <biomart.Attribute name='xmaculatus_homolog_chromosome', display_name='Platyfish chromosome/scaffold name', description=''>
xmaculatus_homolog_chrom_start | <biomart.Attribute name='xmaculatus_homolog_chrom_start', display_name='Platyfish chromosome/scaffold start (bp)', description=''>
xmaculatus_homolog_chrom_end | <biomart.Attribute name='xmaculatus_homolog_chrom_end', display_name='Platyfish chromosome/scaffold end (bp)', description=''>
xmaculatus_homolog_canonical_transcript_protein | <biomart.Attribute name='xmaculatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
xmaculatus_homolog_subtype | <biomart.Attribute name='xmaculatus_homolog_subtype', display_name='Last common ancestor with Platyfish', description=''>
xmaculatus_homolog_orthology_type | <biomart.Attribute name='xmaculatus_homolog_orthology_type', display_name='Platyfish homology type', description=''>
xmaculatus_homolog_perc_id | <biomart.Attribute name='xmaculatus_homolog_perc_id', display_name='%id. target Platyfish gene identical to query gene', description=''>
xmaculatus_homolog_perc_id_r1 | <biomart.Attribute name='xmaculatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Platyfish gene', description=''>
xmaculatus_homolog_goc_score | <biomart.Attribute name='xmaculatus_homolog_goc_score', display_name='Platyfish Gene-order conservation score', description=''>
xmaculatus_homolog_wga_coverage | <biomart.Attribute name='xmaculatus_homolog_wga_coverage', display_name='Platyfish Whole-genome alignment coverage', description=''>
xmaculatus_homolog_orthology_confidence | <biomart.Attribute name='xmaculatus_homolog_orthology_confidence', display_name='Platyfish orthology confidence [0 low, 1 high]', description=''>
oanatinus_homolog_ensembl_gene | <biomart.Attribute name='oanatinus_homolog_ensembl_gene', display_name='Platypus gene stable ID', description=''>
oanatinus_homolog_associated_gene_name | <biomart.Attribute name='oanatinus_homolog_associated_gene_name', display_name='Platypus gene name', description=''>
oanatinus_homolog_ensembl_peptide | <biomart.Attribute name='oanatinus_homolog_ensembl_peptide', display_name='Platypus protein or transcript stable ID', description=''>
oanatinus_homolog_chromosome | <biomart.Attribute name='oanatinus_homolog_chromosome', display_name='Platypus chromosome/scaffold name', description=''>
oanatinus_homolog_chrom_start | <biomart.Attribute name='oanatinus_homolog_chrom_start', display_name='Platypus chromosome/scaffold start (bp)', description=''>
oanatinus_homolog_chrom_end | <biomart.Attribute name='oanatinus_homolog_chrom_end', display_name='Platypus chromosome/scaffold end (bp)', description=''>
oanatinus_homolog_canonical_transcript_protein | <biomart.Attribute name='oanatinus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
oanatinus_homolog_subtype | <biomart.Attribute name='oanatinus_homolog_subtype', display_name='Last common ancestor with Platypus', description=''>
oanatinus_homolog_orthology_type | <biomart.Attribute name='oanatinus_homolog_orthology_type', display_name='Platypus homology type', description=''>
oanatinus_homolog_perc_id | <biomart.Attribute name='oanatinus_homolog_perc_id', display_name='%id. target Platypus gene identical to query gene', description=''>
oanatinus_homolog_perc_id_r1 | <biomart.Attribute name='oanatinus_homolog_perc_id_r1', display_name='%id. query gene identical to target Platypus gene', description=''>
oanatinus_homolog_goc_score | <biomart.Attribute name='oanatinus_homolog_goc_score', display_name='Platypus Gene-order conservation score', description=''>
oanatinus_homolog_wga_coverage | <biomart.Attribute name='oanatinus_homolog_wga_coverage', display_name='Platypus Whole-genome alignment coverage', description=''>
oanatinus_homolog_orthology_confidence | <biomart.Attribute name='oanatinus_homolog_orthology_confidence', display_name='Platypus orthology confidence [0 low, 1 high]', description=''>
umaritimus_homolog_ensembl_gene | <biomart.Attribute name='umaritimus_homolog_ensembl_gene', display_name='Polar bear gene stable ID', description=''>
umaritimus_homolog_associated_gene_name | <biomart.Attribute name='umaritimus_homolog_associated_gene_name', display_name='Polar bear gene name', description=''>
umaritimus_homolog_ensembl_peptide | <biomart.Attribute name='umaritimus_homolog_ensembl_peptide', display_name='Polar bear protein or transcript stable ID', description=''>
umaritimus_homolog_chromosome | <biomart.Attribute name='umaritimus_homolog_chromosome', display_name='Polar bear chromosome/scaffold name', description=''>
umaritimus_homolog_chrom_start | <biomart.Attribute name='umaritimus_homolog_chrom_start', display_name='Polar bear chromosome/scaffold start (bp)', description=''>
umaritimus_homolog_chrom_end | <biomart.Attribute name='umaritimus_homolog_chrom_end', display_name='Polar bear chromosome/scaffold end (bp)', description=''>
umaritimus_homolog_canonical_transcript_protein | <biomart.Attribute name='umaritimus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
umaritimus_homolog_subtype | <biomart.Attribute name='umaritimus_homolog_subtype', display_name='Last common ancestor with Polar bear', description=''>
umaritimus_homolog_orthology_type | <biomart.Attribute name='umaritimus_homolog_orthology_type', display_name='Polar bear homology type', description=''>
umaritimus_homolog_perc_id | <biomart.Attribute name='umaritimus_homolog_perc_id', display_name='%id. target Polar bear gene identical to query gene', description=''>
umaritimus_homolog_perc_id_r1 | <biomart.Attribute name='umaritimus_homolog_perc_id_r1', display_name='%id. query gene identical to target Polar bear gene', description=''>
umaritimus_homolog_goc_score | <biomart.Attribute name='umaritimus_homolog_goc_score', display_name='Polar bear Gene-order conservation score', description=''>
umaritimus_homolog_wga_coverage | <biomart.Attribute name='umaritimus_homolog_wga_coverage', display_name='Polar bear Whole-genome alignment coverage', description=''>
umaritimus_homolog_orthology_confidence | <biomart.Attribute name='umaritimus_homolog_orthology_confidence', display_name='Polar bear orthology confidence [0 low, 1 high]', description=''>
mochrogaster_homolog_ensembl_gene | <biomart.Attribute name='mochrogaster_homolog_ensembl_gene', display_name='Prairie vole gene stable ID', description=''>
mochrogaster_homolog_associated_gene_name | <biomart.Attribute name='mochrogaster_homolog_associated_gene_name', display_name='Prairie vole gene name', description=''>
mochrogaster_homolog_ensembl_peptide | <biomart.Attribute name='mochrogaster_homolog_ensembl_peptide', display_name='Prairie vole protein or transcript stable ID', description=''>
mochrogaster_homolog_chromosome | <biomart.Attribute name='mochrogaster_homolog_chromosome', display_name='Prairie vole chromosome/scaffold name', description=''>
mochrogaster_homolog_chrom_start | <biomart.Attribute name='mochrogaster_homolog_chrom_start', display_name='Prairie vole chromosome/scaffold start (bp)', description=''>
mochrogaster_homolog_chrom_end | <biomart.Attribute name='mochrogaster_homolog_chrom_end', display_name='Prairie vole chromosome/scaffold end (bp)', description=''>
mochrogaster_homolog_canonical_transcript_protein | <biomart.Attribute name='mochrogaster_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mochrogaster_homolog_subtype | <biomart.Attribute name='mochrogaster_homolog_subtype', display_name='Last common ancestor with Prairie vole', description=''>
mochrogaster_homolog_orthology_type | <biomart.Attribute name='mochrogaster_homolog_orthology_type', display_name='Prairie vole homology type', description=''>
mochrogaster_homolog_perc_id | <biomart.Attribute name='mochrogaster_homolog_perc_id', display_name='%id. target Prairie vole gene identical to query gene', description=''>
mochrogaster_homolog_perc_id_r1 | <biomart.Attribute name='mochrogaster_homolog_perc_id_r1', display_name='%id. query gene identical to target Prairie vole gene', description=''>
mochrogaster_homolog_goc_score | <biomart.Attribute name='mochrogaster_homolog_goc_score', display_name='Prairie vole Gene-order conservation score', description=''>
mochrogaster_homolog_wga_coverage | <biomart.Attribute name='mochrogaster_homolog_wga_coverage', display_name='Prairie vole Whole-genome alignment coverage', description=''>
mochrogaster_homolog_orthology_confidence | <biomart.Attribute name='mochrogaster_homolog_orthology_confidence', display_name='Prairie vole orthology confidence [0 low, 1 high]', description=''>
ocuniculus_homolog_ensembl_gene | <biomart.Attribute name='ocuniculus_homolog_ensembl_gene', display_name='Rabbit gene stable ID', description=''>
ocuniculus_homolog_associated_gene_name | <biomart.Attribute name='ocuniculus_homolog_associated_gene_name', display_name='Rabbit gene name', description=''>
ocuniculus_homolog_ensembl_peptide | <biomart.Attribute name='ocuniculus_homolog_ensembl_peptide', display_name='Rabbit protein or transcript stable ID', description=''>
ocuniculus_homolog_chromosome | <biomart.Attribute name='ocuniculus_homolog_chromosome', display_name='Rabbit chromosome/scaffold name', description=''>
ocuniculus_homolog_chrom_start | <biomart.Attribute name='ocuniculus_homolog_chrom_start', display_name='Rabbit chromosome/scaffold start (bp)', description=''>
ocuniculus_homolog_chrom_end | <biomart.Attribute name='ocuniculus_homolog_chrom_end', display_name='Rabbit chromosome/scaffold end (bp)', description=''>
ocuniculus_homolog_canonical_transcript_protein | <biomart.Attribute name='ocuniculus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ocuniculus_homolog_subtype | <biomart.Attribute name='ocuniculus_homolog_subtype', display_name='Last common ancestor with Rabbit', description=''>
ocuniculus_homolog_orthology_type | <biomart.Attribute name='ocuniculus_homolog_orthology_type', display_name='Rabbit homology type', description=''>
ocuniculus_homolog_perc_id | <biomart.Attribute name='ocuniculus_homolog_perc_id', display_name='%id. target Rabbit gene identical to query gene', description=''>
ocuniculus_homolog_perc_id_r1 | <biomart.Attribute name='ocuniculus_homolog_perc_id_r1', display_name='%id. query gene identical to target Rabbit gene', description=''>
ocuniculus_homolog_goc_score | <biomart.Attribute name='ocuniculus_homolog_goc_score', display_name='Rabbit Gene-order conservation score', description=''>
ocuniculus_homolog_wga_coverage | <biomart.Attribute name='ocuniculus_homolog_wga_coverage', display_name='Rabbit Whole-genome alignment coverage', description=''>
ocuniculus_homolog_orthology_confidence | <biomart.Attribute name='ocuniculus_homolog_orthology_confidence', display_name='Rabbit orthology confidence [0 low, 1 high]', description=''>
rnorvegicus_homolog_ensembl_gene | <biomart.Attribute name='rnorvegicus_homolog_ensembl_gene', display_name='Rat gene stable ID', description=''>
rnorvegicus_homolog_associated_gene_name | <biomart.Attribute name='rnorvegicus_homolog_associated_gene_name', display_name='Rat gene name', description=''>
rnorvegicus_homolog_ensembl_peptide | <biomart.Attribute name='rnorvegicus_homolog_ensembl_peptide', display_name='Rat protein or transcript stable ID', description=''>
rnorvegicus_homolog_chromosome | <biomart.Attribute name='rnorvegicus_homolog_chromosome', display_name='Rat chromosome/scaffold name', description=''>
rnorvegicus_homolog_chrom_start | <biomart.Attribute name='rnorvegicus_homolog_chrom_start', display_name='Rat chromosome/scaffold start (bp)', description=''>
rnorvegicus_homolog_chrom_end | <biomart.Attribute name='rnorvegicus_homolog_chrom_end', display_name='Rat chromosome/scaffold end (bp)', description=''>
rnorvegicus_homolog_canonical_transcript_protein | <biomart.Attribute name='rnorvegicus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
rnorvegicus_homolog_subtype | <biomart.Attribute name='rnorvegicus_homolog_subtype', display_name='Last common ancestor with Rat', description=''>
rnorvegicus_homolog_orthology_type | <biomart.Attribute name='rnorvegicus_homolog_orthology_type', display_name='Rat homology type', description=''>
rnorvegicus_homolog_perc_id | <biomart.Attribute name='rnorvegicus_homolog_perc_id', display_name='%id. target Rat gene identical to query gene', description=''>
rnorvegicus_homolog_perc_id_r1 | <biomart.Attribute name='rnorvegicus_homolog_perc_id_r1', display_name='%id. query gene identical to target Rat gene', description=''>
rnorvegicus_homolog_goc_score | <biomart.Attribute name='rnorvegicus_homolog_goc_score', display_name='Rat Gene-order conservation score', description=''>
rnorvegicus_homolog_wga_coverage | <biomart.Attribute name='rnorvegicus_homolog_wga_coverage', display_name='Rat Whole-genome alignment coverage', description=''>
rnorvegicus_homolog_orthology_confidence | <biomart.Attribute name='rnorvegicus_homolog_orthology_confidence', display_name='Rat orthology confidence [0 low, 1 high]', description=''>
vvulpes_homolog_ensembl_gene | <biomart.Attribute name='vvulpes_homolog_ensembl_gene', display_name='Red fox gene stable ID', description=''>
vvulpes_homolog_associated_gene_name | <biomart.Attribute name='vvulpes_homolog_associated_gene_name', display_name='Red fox gene name', description=''>
vvulpes_homolog_ensembl_peptide | <biomart.Attribute name='vvulpes_homolog_ensembl_peptide', display_name='Red fox protein or transcript stable ID', description=''>
vvulpes_homolog_chromosome | <biomart.Attribute name='vvulpes_homolog_chromosome', display_name='Red fox chromosome/scaffold name', description=''>
vvulpes_homolog_chrom_start | <biomart.Attribute name='vvulpes_homolog_chrom_start', display_name='Red fox chromosome/scaffold start (bp)', description=''>
vvulpes_homolog_chrom_end | <biomart.Attribute name='vvulpes_homolog_chrom_end', display_name='Red fox chromosome/scaffold end (bp)', description=''>
vvulpes_homolog_canonical_transcript_protein | <biomart.Attribute name='vvulpes_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
vvulpes_homolog_subtype | <biomart.Attribute name='vvulpes_homolog_subtype', display_name='Last common ancestor with Red fox', description=''>
vvulpes_homolog_orthology_type | <biomart.Attribute name='vvulpes_homolog_orthology_type', display_name='Red fox homology type', description=''>
vvulpes_homolog_perc_id | <biomart.Attribute name='vvulpes_homolog_perc_id', display_name='%id. target Red fox gene identical to query gene', description=''>
vvulpes_homolog_perc_id_r1 | <biomart.Attribute name='vvulpes_homolog_perc_id_r1', display_name='%id. query gene identical to target Red fox gene', description=''>
vvulpes_homolog_goc_score | <biomart.Attribute name='vvulpes_homolog_goc_score', display_name='Red fox Gene-order conservation score', description=''>
vvulpes_homolog_wga_coverage | <biomart.Attribute name='vvulpes_homolog_wga_coverage', display_name='Red fox Whole-genome alignment coverage', description=''>
vvulpes_homolog_orthology_confidence | <biomart.Attribute name='vvulpes_homolog_orthology_confidence', display_name='Red fox orthology confidence [0 low, 1 high]', description=''>
pnattereri_homolog_ensembl_gene | <biomart.Attribute name='pnattereri_homolog_ensembl_gene', display_name='Red-bellied piranha gene stable ID', description=''>
pnattereri_homolog_associated_gene_name | <biomart.Attribute name='pnattereri_homolog_associated_gene_name', display_name='Red-bellied piranha gene name', description=''>
pnattereri_homolog_ensembl_peptide | <biomart.Attribute name='pnattereri_homolog_ensembl_peptide', display_name='Red-bellied piranha protein or transcript stable ID', description=''>
pnattereri_homolog_chromosome | <biomart.Attribute name='pnattereri_homolog_chromosome', display_name='Red-bellied piranha chromosome/scaffold name', description=''>
pnattereri_homolog_chrom_start | <biomart.Attribute name='pnattereri_homolog_chrom_start', display_name='Red-bellied piranha chromosome/scaffold start (bp)', description=''>
pnattereri_homolog_chrom_end | <biomart.Attribute name='pnattereri_homolog_chrom_end', display_name='Red-bellied piranha chromosome/scaffold end (bp)', description=''>
pnattereri_homolog_canonical_transcript_protein | <biomart.Attribute name='pnattereri_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pnattereri_homolog_subtype | <biomart.Attribute name='pnattereri_homolog_subtype', display_name='Last common ancestor with Red-bellied piranha', description=''>
pnattereri_homolog_orthology_type | <biomart.Attribute name='pnattereri_homolog_orthology_type', display_name='Red-bellied piranha homology type', description=''>
pnattereri_homolog_perc_id | <biomart.Attribute name='pnattereri_homolog_perc_id', display_name='%id. target Red-bellied piranha gene identical to query gene', description=''>
pnattereri_homolog_perc_id_r1 | <biomart.Attribute name='pnattereri_homolog_perc_id_r1', display_name='%id. query gene identical to target Red-bellied piranha gene', description=''>
pnattereri_homolog_goc_score | <biomart.Attribute name='pnattereri_homolog_goc_score', display_name='Red-bellied piranha Gene-order conservation score', description=''>
pnattereri_homolog_wga_coverage | <biomart.Attribute name='pnattereri_homolog_wga_coverage', display_name='Red-bellied piranha Whole-genome alignment coverage', description=''>
pnattereri_homolog_orthology_confidence | <biomart.Attribute name='pnattereri_homolog_orthology_confidence', display_name='Red-bellied piranha orthology confidence [0 low, 1 high]', description=''>
ecalabaricus_homolog_ensembl_gene | <biomart.Attribute name='ecalabaricus_homolog_ensembl_gene', display_name='Reedfish gene stable ID', description=''>
ecalabaricus_homolog_associated_gene_name | <biomart.Attribute name='ecalabaricus_homolog_associated_gene_name', display_name='Reedfish gene name', description=''>
ecalabaricus_homolog_ensembl_peptide | <biomart.Attribute name='ecalabaricus_homolog_ensembl_peptide', display_name='Reedfish protein or transcript stable ID', description=''>
ecalabaricus_homolog_chromosome | <biomart.Attribute name='ecalabaricus_homolog_chromosome', display_name='Reedfish chromosome/scaffold name', description=''>
ecalabaricus_homolog_chrom_start | <biomart.Attribute name='ecalabaricus_homolog_chrom_start', display_name='Reedfish chromosome/scaffold start (bp)', description=''>
ecalabaricus_homolog_chrom_end | <biomart.Attribute name='ecalabaricus_homolog_chrom_end', display_name='Reedfish chromosome/scaffold end (bp)', description=''>
ecalabaricus_homolog_canonical_transcript_protein | <biomart.Attribute name='ecalabaricus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ecalabaricus_homolog_subtype | <biomart.Attribute name='ecalabaricus_homolog_subtype', display_name='Last common ancestor with Reedfish', description=''>
ecalabaricus_homolog_orthology_type | <biomart.Attribute name='ecalabaricus_homolog_orthology_type', display_name='Reedfish homology type', description=''>
ecalabaricus_homolog_perc_id | <biomart.Attribute name='ecalabaricus_homolog_perc_id', display_name='%id. target Reedfish gene identical to query gene', description=''>
ecalabaricus_homolog_perc_id_r1 | <biomart.Attribute name='ecalabaricus_homolog_perc_id_r1', display_name='%id. query gene identical to target Reedfish gene', description=''>
ecalabaricus_homolog_goc_score | <biomart.Attribute name='ecalabaricus_homolog_goc_score', display_name='Reedfish Gene-order conservation score', description=''>
ecalabaricus_homolog_wga_coverage | <biomart.Attribute name='ecalabaricus_homolog_wga_coverage', display_name='Reedfish Whole-genome alignment coverage', description=''>
ecalabaricus_homolog_orthology_confidence | <biomart.Attribute name='ecalabaricus_homolog_orthology_confidence', display_name='Reedfish orthology confidence [0 low, 1 high]', description=''>
mcaroli_homolog_ensembl_gene | <biomart.Attribute name='mcaroli_homolog_ensembl_gene', display_name='Ryukyu mouse gene stable ID', description=''>
mcaroli_homolog_associated_gene_name | <biomart.Attribute name='mcaroli_homolog_associated_gene_name', display_name='Ryukyu mouse gene name', description=''>
mcaroli_homolog_ensembl_peptide | <biomart.Attribute name='mcaroli_homolog_ensembl_peptide', display_name='Ryukyu mouse protein or transcript stable ID', description=''>
mcaroli_homolog_chromosome | <biomart.Attribute name='mcaroli_homolog_chromosome', display_name='Ryukyu mouse chromosome/scaffold name', description=''>
mcaroli_homolog_chrom_start | <biomart.Attribute name='mcaroli_homolog_chrom_start', display_name='Ryukyu mouse chromosome/scaffold start (bp)', description=''>
mcaroli_homolog_chrom_end | <biomart.Attribute name='mcaroli_homolog_chrom_end', display_name='Ryukyu mouse chromosome/scaffold end (bp)', description=''>
mcaroli_homolog_canonical_transcript_protein | <biomart.Attribute name='mcaroli_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mcaroli_homolog_subtype | <biomart.Attribute name='mcaroli_homolog_subtype', display_name='Last common ancestor with Ryukyu mouse', description=''>
mcaroli_homolog_orthology_type | <biomart.Attribute name='mcaroli_homolog_orthology_type', display_name='Ryukyu mouse homology type', description=''>
mcaroli_homolog_perc_id | <biomart.Attribute name='mcaroli_homolog_perc_id', display_name='%id. target Ryukyu mouse gene identical to query gene', description=''>
mcaroli_homolog_perc_id_r1 | <biomart.Attribute name='mcaroli_homolog_perc_id_r1', display_name='%id. query gene identical to target Ryukyu mouse gene', description=''>
mcaroli_homolog_goc_score | <biomart.Attribute name='mcaroli_homolog_goc_score', display_name='Ryukyu mouse Gene-order conservation score', description=''>
mcaroli_homolog_wga_coverage | <biomart.Attribute name='mcaroli_homolog_wga_coverage', display_name='Ryukyu mouse Whole-genome alignment coverage', description=''>
mcaroli_homolog_orthology_confidence | <biomart.Attribute name='mcaroli_homolog_orthology_confidence', display_name='Ryukyu mouse orthology confidence [0 low, 1 high]', description=''>
scerevisiae_homolog_ensembl_gene | <biomart.Attribute name='scerevisiae_homolog_ensembl_gene', display_name='Saccharomyces cerevisiae gene stable ID', description=''>
scerevisiae_homolog_associated_gene_name | <biomart.Attribute name='scerevisiae_homolog_associated_gene_name', display_name='Saccharomyces cerevisiae gene name', description=''>
scerevisiae_homolog_ensembl_peptide | <biomart.Attribute name='scerevisiae_homolog_ensembl_peptide', display_name='Saccharomyces cerevisiae protein or transcript stable ID', description=''>
scerevisiae_homolog_chromosome | <biomart.Attribute name='scerevisiae_homolog_chromosome', display_name='Saccharomyces cerevisiae chromosome/scaffold name', description=''>
scerevisiae_homolog_chrom_start | <biomart.Attribute name='scerevisiae_homolog_chrom_start', display_name='Saccharomyces cerevisiae chromosome/scaffold start (bp)', description=''>
scerevisiae_homolog_chrom_end | <biomart.Attribute name='scerevisiae_homolog_chrom_end', display_name='Saccharomyces cerevisiae chromosome/scaffold end (bp)', description=''>
scerevisiae_homolog_canonical_transcript_protein | <biomart.Attribute name='scerevisiae_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
scerevisiae_homolog_subtype | <biomart.Attribute name='scerevisiae_homolog_subtype', display_name='Last common ancestor with Saccharomyces cerevisiae', description=''>
scerevisiae_homolog_orthology_type | <biomart.Attribute name='scerevisiae_homolog_orthology_type', display_name='Saccharomyces cerevisiae homology type', description=''>
scerevisiae_homolog_perc_id | <biomart.Attribute name='scerevisiae_homolog_perc_id', display_name='%id. target Saccharomyces cerevisiae gene identical to query gene', description=''>
scerevisiae_homolog_perc_id_r1 | <biomart.Attribute name='scerevisiae_homolog_perc_id_r1', display_name='%id. query gene identical to target Saccharomyces cerevisiae gene', description=''>
scerevisiae_homolog_orthology_confidence | <biomart.Attribute name='scerevisiae_homolog_orthology_confidence', display_name='Saccharomyces cerevisiae orthology confidence [0 low, 1 high]', description=''>
platipinna_homolog_ensembl_gene | <biomart.Attribute name='platipinna_homolog_ensembl_gene', display_name='Sailfin molly gene stable ID', description=''>
platipinna_homolog_associated_gene_name | <biomart.Attribute name='platipinna_homolog_associated_gene_name', display_name='Sailfin molly gene name', description=''>
platipinna_homolog_ensembl_peptide | <biomart.Attribute name='platipinna_homolog_ensembl_peptide', display_name='Sailfin molly protein or transcript stable ID', description=''>
platipinna_homolog_chromosome | <biomart.Attribute name='platipinna_homolog_chromosome', display_name='Sailfin molly chromosome/scaffold name', description=''>
platipinna_homolog_chrom_start | <biomart.Attribute name='platipinna_homolog_chrom_start', display_name='Sailfin molly chromosome/scaffold start (bp)', description=''>
platipinna_homolog_chrom_end | <biomart.Attribute name='platipinna_homolog_chrom_end', display_name='Sailfin molly chromosome/scaffold end (bp)', description=''>
platipinna_homolog_canonical_transcript_protein | <biomart.Attribute name='platipinna_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
platipinna_homolog_subtype | <biomart.Attribute name='platipinna_homolog_subtype', display_name='Last common ancestor with Sailfin molly', description=''>
platipinna_homolog_orthology_type | <biomart.Attribute name='platipinna_homolog_orthology_type', display_name='Sailfin molly homology type', description=''>
platipinna_homolog_perc_id | <biomart.Attribute name='platipinna_homolog_perc_id', display_name='%id. target Sailfin molly gene identical to query gene', description=''>
platipinna_homolog_perc_id_r1 | <biomart.Attribute name='platipinna_homolog_perc_id_r1', display_name='%id. query gene identical to target Sailfin molly gene', description=''>
platipinna_homolog_goc_score | <biomart.Attribute name='platipinna_homolog_goc_score', display_name='Sailfin molly Gene-order conservation score', description=''>
platipinna_homolog_wga_coverage | <biomart.Attribute name='platipinna_homolog_wga_coverage', display_name='Sailfin molly Whole-genome alignment coverage', description=''>
platipinna_homolog_orthology_confidence | <biomart.Attribute name='platipinna_homolog_orthology_confidence', display_name='Sailfin molly orthology confidence [0 low, 1 high]', description=''>
oaries_homolog_ensembl_gene | <biomart.Attribute name='oaries_homolog_ensembl_gene', display_name='Sheep gene stable ID', description=''>
oaries_homolog_associated_gene_name | <biomart.Attribute name='oaries_homolog_associated_gene_name', display_name='Sheep gene name', description=''>
oaries_homolog_ensembl_peptide | <biomart.Attribute name='oaries_homolog_ensembl_peptide', display_name='Sheep protein or transcript stable ID', description=''>
oaries_homolog_chromosome | <biomart.Attribute name='oaries_homolog_chromosome', display_name='Sheep chromosome/scaffold name', description=''>
oaries_homolog_chrom_start | <biomart.Attribute name='oaries_homolog_chrom_start', display_name='Sheep chromosome/scaffold start (bp)', description=''>
oaries_homolog_chrom_end | <biomart.Attribute name='oaries_homolog_chrom_end', display_name='Sheep chromosome/scaffold end (bp)', description=''>
oaries_homolog_canonical_transcript_protein | <biomart.Attribute name='oaries_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
oaries_homolog_subtype | <biomart.Attribute name='oaries_homolog_subtype', display_name='Last common ancestor with Sheep', description=''>
oaries_homolog_orthology_type | <biomart.Attribute name='oaries_homolog_orthology_type', display_name='Sheep homology type', description=''>
oaries_homolog_perc_id | <biomart.Attribute name='oaries_homolog_perc_id', display_name='%id. target Sheep gene identical to query gene', description=''>
oaries_homolog_perc_id_r1 | <biomart.Attribute name='oaries_homolog_perc_id_r1', display_name='%id. query gene identical to target Sheep gene', description=''>
oaries_homolog_goc_score | <biomart.Attribute name='oaries_homolog_goc_score', display_name='Sheep Gene-order conservation score', description=''>
oaries_homolog_wga_coverage | <biomart.Attribute name='oaries_homolog_wga_coverage', display_name='Sheep Whole-genome alignment coverage', description=''>
oaries_homolog_orthology_confidence | <biomart.Attribute name='oaries_homolog_orthology_confidence', display_name='Sheep orthology confidence [0 low, 1 high]', description=''>
cvariegatus_homolog_ensembl_gene | <biomart.Attribute name='cvariegatus_homolog_ensembl_gene', display_name='Sheepshead minnow gene stable ID', description=''>
cvariegatus_homolog_associated_gene_name | <biomart.Attribute name='cvariegatus_homolog_associated_gene_name', display_name='Sheepshead minnow gene name', description=''>
cvariegatus_homolog_ensembl_peptide | <biomart.Attribute name='cvariegatus_homolog_ensembl_peptide', display_name='Sheepshead minnow protein or transcript stable ID', description=''>
cvariegatus_homolog_chromosome | <biomart.Attribute name='cvariegatus_homolog_chromosome', display_name='Sheepshead minnow chromosome/scaffold name', description=''>
cvariegatus_homolog_chrom_start | <biomart.Attribute name='cvariegatus_homolog_chrom_start', display_name='Sheepshead minnow chromosome/scaffold start (bp)', description=''>
cvariegatus_homolog_chrom_end | <biomart.Attribute name='cvariegatus_homolog_chrom_end', display_name='Sheepshead minnow chromosome/scaffold end (bp)', description=''>
cvariegatus_homolog_canonical_transcript_protein | <biomart.Attribute name='cvariegatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
cvariegatus_homolog_subtype | <biomart.Attribute name='cvariegatus_homolog_subtype', display_name='Last common ancestor with Sheepshead minnow', description=''>
cvariegatus_homolog_orthology_type | <biomart.Attribute name='cvariegatus_homolog_orthology_type', display_name='Sheepshead minnow homology type', description=''>
cvariegatus_homolog_perc_id | <biomart.Attribute name='cvariegatus_homolog_perc_id', display_name='%id. target Sheepshead minnow gene identical to query gene', description=''>
cvariegatus_homolog_perc_id_r1 | <biomart.Attribute name='cvariegatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Sheepshead minnow gene', description=''>
cvariegatus_homolog_goc_score | <biomart.Attribute name='cvariegatus_homolog_goc_score', display_name='Sheepshead minnow Gene-order conservation score', description=''>
cvariegatus_homolog_wga_coverage | <biomart.Attribute name='cvariegatus_homolog_wga_coverage', display_name='Sheepshead minnow Whole-genome alignment coverage', description=''>
cvariegatus_homolog_orthology_confidence | <biomart.Attribute name='cvariegatus_homolog_orthology_confidence', display_name='Sheepshead minnow orthology confidence [0 low, 1 high]', description=''>
pmexicana_homolog_ensembl_gene | <biomart.Attribute name='pmexicana_homolog_ensembl_gene', display_name='Shortfin molly gene stable ID', description=''>
pmexicana_homolog_associated_gene_name | <biomart.Attribute name='pmexicana_homolog_associated_gene_name', display_name='Shortfin molly gene name', description=''>
pmexicana_homolog_ensembl_peptide | <biomart.Attribute name='pmexicana_homolog_ensembl_peptide', display_name='Shortfin molly protein or transcript stable ID', description=''>
pmexicana_homolog_chromosome | <biomart.Attribute name='pmexicana_homolog_chromosome', display_name='Shortfin molly chromosome/scaffold name', description=''>
pmexicana_homolog_chrom_start | <biomart.Attribute name='pmexicana_homolog_chrom_start', display_name='Shortfin molly chromosome/scaffold start (bp)', description=''>
pmexicana_homolog_chrom_end | <biomart.Attribute name='pmexicana_homolog_chrom_end', display_name='Shortfin molly chromosome/scaffold end (bp)', description=''>
pmexicana_homolog_canonical_transcript_protein | <biomart.Attribute name='pmexicana_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
pmexicana_homolog_subtype | <biomart.Attribute name='pmexicana_homolog_subtype', display_name='Last common ancestor with Shortfin molly', description=''>
pmexicana_homolog_orthology_type | <biomart.Attribute name='pmexicana_homolog_orthology_type', display_name='Shortfin molly homology type', description=''>
pmexicana_homolog_perc_id | <biomart.Attribute name='pmexicana_homolog_perc_id', display_name='%id. target Shortfin molly gene identical to query gene', description=''>
pmexicana_homolog_perc_id_r1 | <biomart.Attribute name='pmexicana_homolog_perc_id_r1', display_name='%id. query gene identical to target Shortfin molly gene', description=''>
pmexicana_homolog_goc_score | <biomart.Attribute name='pmexicana_homolog_goc_score', display_name='Shortfin molly Gene-order conservation score', description=''>
pmexicana_homolog_wga_coverage | <biomart.Attribute name='pmexicana_homolog_wga_coverage', display_name='Shortfin molly Whole-genome alignment coverage', description=''>
pmexicana_homolog_orthology_confidence | <biomart.Attribute name='pmexicana_homolog_orthology_confidence', display_name='Shortfin molly orthology confidence [0 low, 1 high]', description=''>
saraneus_homolog_ensembl_gene | <biomart.Attribute name='saraneus_homolog_ensembl_gene', display_name='Shrew gene stable ID', description=''>
saraneus_homolog_associated_gene_name | <biomart.Attribute name='saraneus_homolog_associated_gene_name', display_name='Shrew gene name', description=''>
saraneus_homolog_ensembl_peptide | <biomart.Attribute name='saraneus_homolog_ensembl_peptide', display_name='Shrew protein or transcript stable ID', description=''>
saraneus_homolog_chromosome | <biomart.Attribute name='saraneus_homolog_chromosome', display_name='Shrew chromosome/scaffold name', description=''>
saraneus_homolog_chrom_start | <biomart.Attribute name='saraneus_homolog_chrom_start', display_name='Shrew chromosome/scaffold start (bp)', description=''>
saraneus_homolog_chrom_end | <biomart.Attribute name='saraneus_homolog_chrom_end', display_name='Shrew chromosome/scaffold end (bp)', description=''>
saraneus_homolog_canonical_transcript_protein | <biomart.Attribute name='saraneus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
saraneus_homolog_subtype | <biomart.Attribute name='saraneus_homolog_subtype', display_name='Last common ancestor with Shrew', description=''>
saraneus_homolog_orthology_type | <biomart.Attribute name='saraneus_homolog_orthology_type', display_name='Shrew homology type', description=''>
saraneus_homolog_perc_id | <biomart.Attribute name='saraneus_homolog_perc_id', display_name='%id. target Shrew gene identical to query gene', description=''>
saraneus_homolog_perc_id_r1 | <biomart.Attribute name='saraneus_homolog_perc_id_r1', display_name='%id. query gene identical to target Shrew gene', description=''>
saraneus_homolog_goc_score | <biomart.Attribute name='saraneus_homolog_goc_score', display_name='Shrew Gene-order conservation score', description=''>
saraneus_homolog_wga_coverage | <biomart.Attribute name='saraneus_homolog_wga_coverage', display_name='Shrew Whole-genome alignment coverage', description=''>
saraneus_homolog_orthology_confidence | <biomart.Attribute name='saraneus_homolog_orthology_confidence', display_name='Shrew orthology confidence [0 low, 1 high]', description=''>
mpahari_homolog_ensembl_gene | <biomart.Attribute name='mpahari_homolog_ensembl_gene', display_name='Shrew mouse gene stable ID', description=''>
mpahari_homolog_associated_gene_name | <biomart.Attribute name='mpahari_homolog_associated_gene_name', display_name='Shrew mouse gene name', description=''>
mpahari_homolog_ensembl_peptide | <biomart.Attribute name='mpahari_homolog_ensembl_peptide', display_name='Shrew mouse protein or transcript stable ID', description=''>
mpahari_homolog_chromosome | <biomart.Attribute name='mpahari_homolog_chromosome', display_name='Shrew mouse chromosome/scaffold name', description=''>
mpahari_homolog_chrom_start | <biomart.Attribute name='mpahari_homolog_chrom_start', display_name='Shrew mouse chromosome/scaffold start (bp)', description=''>
mpahari_homolog_chrom_end | <biomart.Attribute name='mpahari_homolog_chrom_end', display_name='Shrew mouse chromosome/scaffold end (bp)', description=''>
mpahari_homolog_canonical_transcript_protein | <biomart.Attribute name='mpahari_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mpahari_homolog_subtype | <biomart.Attribute name='mpahari_homolog_subtype', display_name='Last common ancestor with Shrew mouse', description=''>
mpahari_homolog_orthology_type | <biomart.Attribute name='mpahari_homolog_orthology_type', display_name='Shrew mouse homology type', description=''>
mpahari_homolog_perc_id | <biomart.Attribute name='mpahari_homolog_perc_id', display_name='%id. target Shrew mouse gene identical to query gene', description=''>
mpahari_homolog_perc_id_r1 | <biomart.Attribute name='mpahari_homolog_perc_id_r1', display_name='%id. query gene identical to target Shrew mouse gene', description=''>
mpahari_homolog_goc_score | <biomart.Attribute name='mpahari_homolog_goc_score', display_name='Shrew mouse Gene-order conservation score', description=''>
mpahari_homolog_wga_coverage | <biomart.Attribute name='mpahari_homolog_wga_coverage', display_name='Shrew mouse Whole-genome alignment coverage', description=''>
mpahari_homolog_orthology_confidence | <biomart.Attribute name='mpahari_homolog_orthology_confidence', display_name='Shrew mouse orthology confidence [0 low, 1 high]', description=''>
bsplendens_homolog_ensembl_gene | <biomart.Attribute name='bsplendens_homolog_ensembl_gene', display_name='Siamese fighting fish gene stable ID', description=''>
bsplendens_homolog_associated_gene_name | <biomart.Attribute name='bsplendens_homolog_associated_gene_name', display_name='Siamese fighting fish gene name', description=''>
bsplendens_homolog_ensembl_peptide | <biomart.Attribute name='bsplendens_homolog_ensembl_peptide', display_name='Siamese fighting fish protein or transcript stable ID', description=''>
bsplendens_homolog_chromosome | <biomart.Attribute name='bsplendens_homolog_chromosome', display_name='Siamese fighting fish chromosome/scaffold name', description=''>
bsplendens_homolog_chrom_start | <biomart.Attribute name='bsplendens_homolog_chrom_start', display_name='Siamese fighting fish chromosome/scaffold start (bp)', description=''>
bsplendens_homolog_chrom_end | <biomart.Attribute name='bsplendens_homolog_chrom_end', display_name='Siamese fighting fish chromosome/scaffold end (bp)', description=''>
bsplendens_homolog_canonical_transcript_protein | <biomart.Attribute name='bsplendens_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
bsplendens_homolog_subtype | <biomart.Attribute name='bsplendens_homolog_subtype', display_name='Last common ancestor with Siamese fighting fish', description=''>
bsplendens_homolog_orthology_type | <biomart.Attribute name='bsplendens_homolog_orthology_type', display_name='Siamese fighting fish homology type', description=''>
bsplendens_homolog_perc_id | <biomart.Attribute name='bsplendens_homolog_perc_id', display_name='%id. target Siamese fighting fish gene identical to query gene', description=''>
bsplendens_homolog_perc_id_r1 | <biomart.Attribute name='bsplendens_homolog_perc_id_r1', display_name='%id. query gene identical to target Siamese fighting fish gene', description=''>
bsplendens_homolog_goc_score | <biomart.Attribute name='bsplendens_homolog_goc_score', display_name='Siamese fighting fish Gene-order conservation score', description=''>
bsplendens_homolog_wga_coverage | <biomart.Attribute name='bsplendens_homolog_wga_coverage', display_name='Siamese fighting fish Whole-genome alignment coverage', description=''>
bsplendens_homolog_orthology_confidence | <biomart.Attribute name='bsplendens_homolog_orthology_confidence', display_name='Siamese fighting fish orthology confidence [0 low, 1 high]', description=''>
choffmanni_homolog_ensembl_gene | <biomart.Attribute name='choffmanni_homolog_ensembl_gene', display_name='Sloth gene stable ID', description=''>
choffmanni_homolog_associated_gene_name | <biomart.Attribute name='choffmanni_homolog_associated_gene_name', display_name='Sloth gene name', description=''>
choffmanni_homolog_ensembl_peptide | <biomart.Attribute name='choffmanni_homolog_ensembl_peptide', display_name='Sloth protein or transcript stable ID', description=''>
choffmanni_homolog_chromosome | <biomart.Attribute name='choffmanni_homolog_chromosome', display_name='Sloth chromosome/scaffold name', description=''>
choffmanni_homolog_chrom_start | <biomart.Attribute name='choffmanni_homolog_chrom_start', display_name='Sloth chromosome/scaffold start (bp)', description=''>
choffmanni_homolog_chrom_end | <biomart.Attribute name='choffmanni_homolog_chrom_end', display_name='Sloth chromosome/scaffold end (bp)', description=''>
choffmanni_homolog_canonical_transcript_protein | <biomart.Attribute name='choffmanni_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
choffmanni_homolog_subtype | <biomart.Attribute name='choffmanni_homolog_subtype', display_name='Last common ancestor with Sloth', description=''>
choffmanni_homolog_orthology_type | <biomart.Attribute name='choffmanni_homolog_orthology_type', display_name='Sloth homology type', description=''>
choffmanni_homolog_perc_id | <biomart.Attribute name='choffmanni_homolog_perc_id', display_name='%id. target Sloth gene identical to query gene', description=''>
choffmanni_homolog_perc_id_r1 | <biomart.Attribute name='choffmanni_homolog_perc_id_r1', display_name='%id. query gene identical to target Sloth gene', description=''>
choffmanni_homolog_goc_score | <biomart.Attribute name='choffmanni_homolog_goc_score', display_name='Sloth Gene-order conservation score', description=''>
choffmanni_homolog_wga_coverage | <biomart.Attribute name='choffmanni_homolog_wga_coverage', display_name='Sloth Whole-genome alignment coverage', description=''>
choffmanni_homolog_orthology_confidence | <biomart.Attribute name='choffmanni_homolog_orthology_confidence', display_name='Sloth orthology confidence [0 low, 1 high]', description=''>
catys_homolog_ensembl_gene | <biomart.Attribute name='catys_homolog_ensembl_gene', display_name='Sooty mangabey gene stable ID', description=''>
catys_homolog_associated_gene_name | <biomart.Attribute name='catys_homolog_associated_gene_name', display_name='Sooty mangabey gene name', description=''>
catys_homolog_ensembl_peptide | <biomart.Attribute name='catys_homolog_ensembl_peptide', display_name='Sooty mangabey protein or transcript stable ID', description=''>
catys_homolog_chromosome | <biomart.Attribute name='catys_homolog_chromosome', display_name='Sooty mangabey chromosome/scaffold name', description=''>
catys_homolog_chrom_start | <biomart.Attribute name='catys_homolog_chrom_start', display_name='Sooty mangabey chromosome/scaffold start (bp)', description=''>
catys_homolog_chrom_end | <biomart.Attribute name='catys_homolog_chrom_end', display_name='Sooty mangabey chromosome/scaffold end (bp)', description=''>
catys_homolog_canonical_transcript_protein | <biomart.Attribute name='catys_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
catys_homolog_subtype | <biomart.Attribute name='catys_homolog_subtype', display_name='Last common ancestor with Sooty mangabey', description=''>
catys_homolog_orthology_type | <biomart.Attribute name='catys_homolog_orthology_type', display_name='Sooty mangabey homology type', description=''>
catys_homolog_perc_id | <biomart.Attribute name='catys_homolog_perc_id', display_name='%id. target Sooty mangabey gene identical to query gene', description=''>
catys_homolog_perc_id_r1 | <biomart.Attribute name='catys_homolog_perc_id_r1', display_name='%id. query gene identical to target Sooty mangabey gene', description=''>
catys_homolog_goc_score | <biomart.Attribute name='catys_homolog_goc_score', display_name='Sooty mangabey Gene-order conservation score', description=''>
catys_homolog_wga_coverage | <biomart.Attribute name='catys_homolog_wga_coverage', display_name='Sooty mangabey Whole-genome alignment coverage', description=''>
catys_homolog_orthology_confidence | <biomart.Attribute name='catys_homolog_orthology_confidence', display_name='Sooty mangabey orthology confidence [0 low, 1 high]', description=''>
loculatus_homolog_ensembl_gene | <biomart.Attribute name='loculatus_homolog_ensembl_gene', display_name='Spotted gar gene stable ID', description=''>
loculatus_homolog_associated_gene_name | <biomart.Attribute name='loculatus_homolog_associated_gene_name', display_name='Spotted gar gene name', description=''>
loculatus_homolog_ensembl_peptide | <biomart.Attribute name='loculatus_homolog_ensembl_peptide', display_name='Spotted gar protein or transcript stable ID', description=''>
loculatus_homolog_chromosome | <biomart.Attribute name='loculatus_homolog_chromosome', display_name='Spotted gar chromosome/scaffold name', description=''>
loculatus_homolog_chrom_start | <biomart.Attribute name='loculatus_homolog_chrom_start', display_name='Spotted gar chromosome/scaffold start (bp)', description=''>
loculatus_homolog_chrom_end | <biomart.Attribute name='loculatus_homolog_chrom_end', display_name='Spotted gar chromosome/scaffold end (bp)', description=''>
loculatus_homolog_canonical_transcript_protein | <biomart.Attribute name='loculatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
loculatus_homolog_subtype | <biomart.Attribute name='loculatus_homolog_subtype', display_name='Last common ancestor with Spotted gar', description=''>
loculatus_homolog_orthology_type | <biomart.Attribute name='loculatus_homolog_orthology_type', display_name='Spotted gar homology type', description=''>
loculatus_homolog_perc_id | <biomart.Attribute name='loculatus_homolog_perc_id', display_name='%id. target Spotted gar gene identical to query gene', description=''>
loculatus_homolog_perc_id_r1 | <biomart.Attribute name='loculatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Spotted gar gene', description=''>
loculatus_homolog_goc_score | <biomart.Attribute name='loculatus_homolog_goc_score', display_name='Spotted gar Gene-order conservation score', description=''>
loculatus_homolog_wga_coverage | <biomart.Attribute name='loculatus_homolog_wga_coverage', display_name='Spotted gar Whole-genome alignment coverage', description=''>
loculatus_homolog_orthology_confidence | <biomart.Attribute name='loculatus_homolog_orthology_confidence', display_name='Spotted gar orthology confidence [0 low, 1 high]', description=''>
itridecemlineatus_homolog_ensembl_gene | <biomart.Attribute name='itridecemlineatus_homolog_ensembl_gene', display_name='Squirrel gene stable ID', description=''>
itridecemlineatus_homolog_associated_gene_name | <biomart.Attribute name='itridecemlineatus_homolog_associated_gene_name', display_name='Squirrel gene name', description=''>
itridecemlineatus_homolog_ensembl_peptide | <biomart.Attribute name='itridecemlineatus_homolog_ensembl_peptide', display_name='Squirrel protein or transcript stable ID', description=''>
itridecemlineatus_homolog_chromosome | <biomart.Attribute name='itridecemlineatus_homolog_chromosome', display_name='Squirrel chromosome/scaffold name', description=''>
itridecemlineatus_homolog_chrom_start | <biomart.Attribute name='itridecemlineatus_homolog_chrom_start', display_name='Squirrel chromosome/scaffold start (bp)', description=''>
itridecemlineatus_homolog_chrom_end | <biomart.Attribute name='itridecemlineatus_homolog_chrom_end', display_name='Squirrel chromosome/scaffold end (bp)', description=''>
itridecemlineatus_homolog_canonical_transcript_protein | <biomart.Attribute name='itridecemlineatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
itridecemlineatus_homolog_subtype | <biomart.Attribute name='itridecemlineatus_homolog_subtype', display_name='Last common ancestor with Squirrel', description=''>
itridecemlineatus_homolog_orthology_type | <biomart.Attribute name='itridecemlineatus_homolog_orthology_type', display_name='Squirrel homology type', description=''>
itridecemlineatus_homolog_perc_id | <biomart.Attribute name='itridecemlineatus_homolog_perc_id', display_name='%id. target Squirrel gene identical to query gene', description=''>
itridecemlineatus_homolog_perc_id_r1 | <biomart.Attribute name='itridecemlineatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Squirrel gene', description=''>
itridecemlineatus_homolog_goc_score | <biomart.Attribute name='itridecemlineatus_homolog_goc_score', display_name='Squirrel Gene-order conservation score', description=''>
itridecemlineatus_homolog_wga_coverage | <biomart.Attribute name='itridecemlineatus_homolog_wga_coverage', display_name='Squirrel Whole-genome alignment coverage', description=''>
itridecemlineatus_homolog_orthology_confidence | <biomart.Attribute name='itridecemlineatus_homolog_orthology_confidence', display_name='Squirrel orthology confidence [0 low, 1 high]', description=''>
mspicilegus_homolog_ensembl_gene | <biomart.Attribute name='mspicilegus_homolog_ensembl_gene', display_name='Steppe mouse gene stable ID', description=''>
mspicilegus_homolog_associated_gene_name | <biomart.Attribute name='mspicilegus_homolog_associated_gene_name', display_name='Steppe mouse gene name', description=''>
mspicilegus_homolog_ensembl_peptide | <biomart.Attribute name='mspicilegus_homolog_ensembl_peptide', display_name='Steppe mouse protein or transcript stable ID', description=''>
mspicilegus_homolog_chromosome | <biomart.Attribute name='mspicilegus_homolog_chromosome', display_name='Steppe mouse chromosome/scaffold name', description=''>
mspicilegus_homolog_chrom_start | <biomart.Attribute name='mspicilegus_homolog_chrom_start', display_name='Steppe mouse chromosome/scaffold start (bp)', description=''>
mspicilegus_homolog_chrom_end | <biomart.Attribute name='mspicilegus_homolog_chrom_end', display_name='Steppe mouse chromosome/scaffold end (bp)', description=''>
mspicilegus_homolog_canonical_transcript_protein | <biomart.Attribute name='mspicilegus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mspicilegus_homolog_subtype | <biomart.Attribute name='mspicilegus_homolog_subtype', display_name='Last common ancestor with Steppe mouse', description=''>
mspicilegus_homolog_orthology_type | <biomart.Attribute name='mspicilegus_homolog_orthology_type', display_name='Steppe mouse homology type', description=''>
mspicilegus_homolog_perc_id | <biomart.Attribute name='mspicilegus_homolog_perc_id', display_name='%id. target Steppe mouse gene identical to query gene', description=''>
mspicilegus_homolog_perc_id_r1 | <biomart.Attribute name='mspicilegus_homolog_perc_id_r1', display_name='%id. query gene identical to target Steppe mouse gene', description=''>
mspicilegus_homolog_goc_score | <biomart.Attribute name='mspicilegus_homolog_goc_score', display_name='Steppe mouse Gene-order conservation score', description=''>
mspicilegus_homolog_wga_coverage | <biomart.Attribute name='mspicilegus_homolog_wga_coverage', display_name='Steppe mouse Whole-genome alignment coverage', description=''>
mspicilegus_homolog_orthology_confidence | <biomart.Attribute name='mspicilegus_homolog_orthology_confidence', display_name='Steppe mouse orthology confidence [0 low, 1 high]', description=''>
gaculeatus_homolog_ensembl_gene | <biomart.Attribute name='gaculeatus_homolog_ensembl_gene', display_name='Stickleback gene stable ID', description=''>
gaculeatus_homolog_associated_gene_name | <biomart.Attribute name='gaculeatus_homolog_associated_gene_name', display_name='Stickleback gene name', description=''>
gaculeatus_homolog_ensembl_peptide | <biomart.Attribute name='gaculeatus_homolog_ensembl_peptide', display_name='Stickleback protein or transcript stable ID', description=''>
gaculeatus_homolog_chromosome | <biomart.Attribute name='gaculeatus_homolog_chromosome', display_name='Stickleback chromosome/scaffold name', description=''>
gaculeatus_homolog_chrom_start | <biomart.Attribute name='gaculeatus_homolog_chrom_start', display_name='Stickleback chromosome/scaffold start (bp)', description=''>
gaculeatus_homolog_chrom_end | <biomart.Attribute name='gaculeatus_homolog_chrom_end', display_name='Stickleback chromosome/scaffold end (bp)', description=''>
gaculeatus_homolog_canonical_transcript_protein | <biomart.Attribute name='gaculeatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
gaculeatus_homolog_subtype | <biomart.Attribute name='gaculeatus_homolog_subtype', display_name='Last common ancestor with Stickleback', description=''>
gaculeatus_homolog_orthology_type | <biomart.Attribute name='gaculeatus_homolog_orthology_type', display_name='Stickleback homology type', description=''>
gaculeatus_homolog_perc_id | <biomart.Attribute name='gaculeatus_homolog_perc_id', display_name='%id. target Stickleback gene identical to query gene', description=''>
gaculeatus_homolog_perc_id_r1 | <biomart.Attribute name='gaculeatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Stickleback gene', description=''>
gaculeatus_homolog_goc_score | <biomart.Attribute name='gaculeatus_homolog_goc_score', display_name='Stickleback Gene-order conservation score', description=''>
gaculeatus_homolog_wga_coverage | <biomart.Attribute name='gaculeatus_homolog_wga_coverage', display_name='Stickleback Whole-genome alignment coverage', description=''>
gaculeatus_homolog_orthology_confidence | <biomart.Attribute name='gaculeatus_homolog_orthology_confidence', display_name='Stickleback orthology confidence [0 low, 1 high]', description=''>
malbus_homolog_ensembl_gene | <biomart.Attribute name='malbus_homolog_ensembl_gene', display_name='Swamp eel gene stable ID', description=''>
malbus_homolog_associated_gene_name | <biomart.Attribute name='malbus_homolog_associated_gene_name', display_name='Swamp eel gene name', description=''>
malbus_homolog_ensembl_peptide | <biomart.Attribute name='malbus_homolog_ensembl_peptide', display_name='Swamp eel protein or transcript stable ID', description=''>
malbus_homolog_chromosome | <biomart.Attribute name='malbus_homolog_chromosome', display_name='Swamp eel chromosome/scaffold name', description=''>
malbus_homolog_chrom_start | <biomart.Attribute name='malbus_homolog_chrom_start', display_name='Swamp eel chromosome/scaffold start (bp)', description=''>
malbus_homolog_chrom_end | <biomart.Attribute name='malbus_homolog_chrom_end', display_name='Swamp eel chromosome/scaffold end (bp)', description=''>
malbus_homolog_canonical_transcript_protein | <biomart.Attribute name='malbus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
malbus_homolog_subtype | <biomart.Attribute name='malbus_homolog_subtype', display_name='Last common ancestor with Swamp eel', description=''>
malbus_homolog_orthology_type | <biomart.Attribute name='malbus_homolog_orthology_type', display_name='Swamp eel homology type', description=''>
malbus_homolog_perc_id | <biomart.Attribute name='malbus_homolog_perc_id', display_name='%id. target Swamp eel gene identical to query gene', description=''>
malbus_homolog_perc_id_r1 | <biomart.Attribute name='malbus_homolog_perc_id_r1', display_name='%id. query gene identical to target Swamp eel gene', description=''>
malbus_homolog_goc_score | <biomart.Attribute name='malbus_homolog_goc_score', display_name='Swamp eel Gene-order conservation score', description=''>
malbus_homolog_wga_coverage | <biomart.Attribute name='malbus_homolog_wga_coverage', display_name='Swamp eel Whole-genome alignment coverage', description=''>
malbus_homolog_orthology_confidence | <biomart.Attribute name='malbus_homolog_orthology_confidence', display_name='Swamp eel orthology confidence [0 low, 1 high]', description=''>
csyrichta_homolog_ensembl_gene | <biomart.Attribute name='csyrichta_homolog_ensembl_gene', display_name='Tarsier gene stable ID', description=''>
csyrichta_homolog_associated_gene_name | <biomart.Attribute name='csyrichta_homolog_associated_gene_name', display_name='Tarsier gene name', description=''>
csyrichta_homolog_ensembl_peptide | <biomart.Attribute name='csyrichta_homolog_ensembl_peptide', display_name='Tarsier protein or transcript stable ID', description=''>
csyrichta_homolog_chromosome | <biomart.Attribute name='csyrichta_homolog_chromosome', display_name='Tarsier chromosome/scaffold name', description=''>
csyrichta_homolog_chrom_start | <biomart.Attribute name='csyrichta_homolog_chrom_start', display_name='Tarsier chromosome/scaffold start (bp)', description=''>
csyrichta_homolog_chrom_end | <biomart.Attribute name='csyrichta_homolog_chrom_end', display_name='Tarsier chromosome/scaffold end (bp)', description=''>
csyrichta_homolog_canonical_transcript_protein | <biomart.Attribute name='csyrichta_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
csyrichta_homolog_subtype | <biomart.Attribute name='csyrichta_homolog_subtype', display_name='Last common ancestor with Tarsier', description=''>
csyrichta_homolog_orthology_type | <biomart.Attribute name='csyrichta_homolog_orthology_type', display_name='Tarsier homology type', description=''>
csyrichta_homolog_perc_id | <biomart.Attribute name='csyrichta_homolog_perc_id', display_name='%id. target Tarsier gene identical to query gene', description=''>
csyrichta_homolog_perc_id_r1 | <biomart.Attribute name='csyrichta_homolog_perc_id_r1', display_name='%id. query gene identical to target Tarsier gene', description=''>
csyrichta_homolog_goc_score | <biomart.Attribute name='csyrichta_homolog_goc_score', display_name='Tarsier Gene-order conservation score', description=''>
csyrichta_homolog_wga_coverage | <biomart.Attribute name='csyrichta_homolog_wga_coverage', display_name='Tarsier Whole-genome alignment coverage', description=''>
csyrichta_homolog_orthology_confidence | <biomart.Attribute name='csyrichta_homolog_orthology_confidence', display_name='Tarsier orthology confidence [0 low, 1 high]', description=''>
sharrisii_homolog_ensembl_gene | <biomart.Attribute name='sharrisii_homolog_ensembl_gene', display_name='Tasmanian devil gene stable ID', description=''>
sharrisii_homolog_associated_gene_name | <biomart.Attribute name='sharrisii_homolog_associated_gene_name', display_name='Tasmanian devil gene name', description=''>
sharrisii_homolog_ensembl_peptide | <biomart.Attribute name='sharrisii_homolog_ensembl_peptide', display_name='Tasmanian devil protein or transcript stable ID', description=''>
sharrisii_homolog_chromosome | <biomart.Attribute name='sharrisii_homolog_chromosome', display_name='Tasmanian devil chromosome/scaffold name', description=''>
sharrisii_homolog_chrom_start | <biomart.Attribute name='sharrisii_homolog_chrom_start', display_name='Tasmanian devil chromosome/scaffold start (bp)', description=''>
sharrisii_homolog_chrom_end | <biomart.Attribute name='sharrisii_homolog_chrom_end', display_name='Tasmanian devil chromosome/scaffold end (bp)', description=''>
sharrisii_homolog_canonical_transcript_protein | <biomart.Attribute name='sharrisii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sharrisii_homolog_subtype | <biomart.Attribute name='sharrisii_homolog_subtype', display_name='Last common ancestor with Tasmanian devil', description=''>
sharrisii_homolog_orthology_type | <biomart.Attribute name='sharrisii_homolog_orthology_type', display_name='Tasmanian devil homology type', description=''>
sharrisii_homolog_perc_id | <biomart.Attribute name='sharrisii_homolog_perc_id', display_name='%id. target Tasmanian devil gene identical to query gene', description=''>
sharrisii_homolog_perc_id_r1 | <biomart.Attribute name='sharrisii_homolog_perc_id_r1', display_name='%id. query gene identical to target Tasmanian devil gene', description=''>
sharrisii_homolog_goc_score | <biomart.Attribute name='sharrisii_homolog_goc_score', display_name='Tasmanian devil Gene-order conservation score', description=''>
sharrisii_homolog_wga_coverage | <biomart.Attribute name='sharrisii_homolog_wga_coverage', display_name='Tasmanian devil Whole-genome alignment coverage', description=''>
sharrisii_homolog_orthology_confidence | <biomart.Attribute name='sharrisii_homolog_orthology_confidence', display_name='Tasmanian devil orthology confidence [0 low, 1 high]', description=''>
tnigroviridis_homolog_ensembl_gene | <biomart.Attribute name='tnigroviridis_homolog_ensembl_gene', display_name='Tetraodon gene stable ID', description=''>
tnigroviridis_homolog_associated_gene_name | <biomart.Attribute name='tnigroviridis_homolog_associated_gene_name', display_name='Tetraodon gene name', description=''>
tnigroviridis_homolog_ensembl_peptide | <biomart.Attribute name='tnigroviridis_homolog_ensembl_peptide', display_name='Tetraodon protein or transcript stable ID', description=''>
tnigroviridis_homolog_chromosome | <biomart.Attribute name='tnigroviridis_homolog_chromosome', display_name='Tetraodon chromosome/scaffold name', description=''>
tnigroviridis_homolog_chrom_start | <biomart.Attribute name='tnigroviridis_homolog_chrom_start', display_name='Tetraodon chromosome/scaffold start (bp)', description=''>
tnigroviridis_homolog_chrom_end | <biomart.Attribute name='tnigroviridis_homolog_chrom_end', display_name='Tetraodon chromosome/scaffold end (bp)', description=''>
tnigroviridis_homolog_canonical_transcript_protein | <biomart.Attribute name='tnigroviridis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
tnigroviridis_homolog_subtype | <biomart.Attribute name='tnigroviridis_homolog_subtype', display_name='Last common ancestor with Tetraodon', description=''>
tnigroviridis_homolog_orthology_type | <biomart.Attribute name='tnigroviridis_homolog_orthology_type', display_name='Tetraodon homology type', description=''>
tnigroviridis_homolog_perc_id | <biomart.Attribute name='tnigroviridis_homolog_perc_id', display_name='%id. target Tetraodon gene identical to query gene', description=''>
tnigroviridis_homolog_perc_id_r1 | <biomart.Attribute name='tnigroviridis_homolog_perc_id_r1', display_name='%id. query gene identical to target Tetraodon gene', description=''>
tnigroviridis_homolog_goc_score | <biomart.Attribute name='tnigroviridis_homolog_goc_score', display_name='Tetraodon Gene-order conservation score', description=''>
tnigroviridis_homolog_wga_coverage | <biomart.Attribute name='tnigroviridis_homolog_wga_coverage', display_name='Tetraodon Whole-genome alignment coverage', description=''>
tnigroviridis_homolog_orthology_confidence | <biomart.Attribute name='tnigroviridis_homolog_orthology_confidence', display_name='Tetraodon orthology confidence [0 low, 1 high]', description=''>
ptaltaica_homolog_ensembl_gene | <biomart.Attribute name='ptaltaica_homolog_ensembl_gene', display_name='Tiger gene stable ID', description=''>
ptaltaica_homolog_associated_gene_name | <biomart.Attribute name='ptaltaica_homolog_associated_gene_name', display_name='Tiger gene name', description=''>
ptaltaica_homolog_ensembl_peptide | <biomart.Attribute name='ptaltaica_homolog_ensembl_peptide', display_name='Tiger protein or transcript stable ID', description=''>
ptaltaica_homolog_chromosome | <biomart.Attribute name='ptaltaica_homolog_chromosome', display_name='Tiger chromosome/scaffold name', description=''>
ptaltaica_homolog_chrom_start | <biomart.Attribute name='ptaltaica_homolog_chrom_start', display_name='Tiger chromosome/scaffold start (bp)', description=''>
ptaltaica_homolog_chrom_end | <biomart.Attribute name='ptaltaica_homolog_chrom_end', display_name='Tiger chromosome/scaffold end (bp)', description=''>
ptaltaica_homolog_canonical_transcript_protein | <biomart.Attribute name='ptaltaica_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ptaltaica_homolog_subtype | <biomart.Attribute name='ptaltaica_homolog_subtype', display_name='Last common ancestor with Tiger', description=''>
ptaltaica_homolog_orthology_type | <biomart.Attribute name='ptaltaica_homolog_orthology_type', display_name='Tiger homology type', description=''>
ptaltaica_homolog_perc_id | <biomart.Attribute name='ptaltaica_homolog_perc_id', display_name='%id. target Tiger gene identical to query gene', description=''>
ptaltaica_homolog_perc_id_r1 | <biomart.Attribute name='ptaltaica_homolog_perc_id_r1', display_name='%id. query gene identical to target Tiger gene', description=''>
ptaltaica_homolog_goc_score | <biomart.Attribute name='ptaltaica_homolog_goc_score', display_name='Tiger Gene-order conservation score', description=''>
ptaltaica_homolog_wga_coverage | <biomart.Attribute name='ptaltaica_homolog_wga_coverage', display_name='Tiger Whole-genome alignment coverage', description=''>
ptaltaica_homolog_orthology_confidence | <biomart.Attribute name='ptaltaica_homolog_orthology_confidence', display_name='Tiger orthology confidence [0 low, 1 high]', description=''>
hcomes_homolog_ensembl_gene | <biomart.Attribute name='hcomes_homolog_ensembl_gene', display_name='Tiger tail seahorse gene stable ID', description=''>
hcomes_homolog_associated_gene_name | <biomart.Attribute name='hcomes_homolog_associated_gene_name', display_name='Tiger tail seahorse gene name', description=''>
hcomes_homolog_ensembl_peptide | <biomart.Attribute name='hcomes_homolog_ensembl_peptide', display_name='Tiger tail seahorse protein or transcript stable ID', description=''>
hcomes_homolog_chromosome | <biomart.Attribute name='hcomes_homolog_chromosome', display_name='Tiger tail seahorse chromosome/scaffold name', description=''>
hcomes_homolog_chrom_start | <biomart.Attribute name='hcomes_homolog_chrom_start', display_name='Tiger tail seahorse chromosome/scaffold start (bp)', description=''>
hcomes_homolog_chrom_end | <biomart.Attribute name='hcomes_homolog_chrom_end', display_name='Tiger tail seahorse chromosome/scaffold end (bp)', description=''>
hcomes_homolog_canonical_transcript_protein | <biomart.Attribute name='hcomes_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
hcomes_homolog_subtype | <biomart.Attribute name='hcomes_homolog_subtype', display_name='Last common ancestor with Tiger tail seahorse', description=''>
hcomes_homolog_orthology_type | <biomart.Attribute name='hcomes_homolog_orthology_type', display_name='Tiger tail seahorse homology type', description=''>
hcomes_homolog_perc_id | <biomart.Attribute name='hcomes_homolog_perc_id', display_name='%id. target Tiger tail seahorse gene identical to query gene', description=''>
hcomes_homolog_perc_id_r1 | <biomart.Attribute name='hcomes_homolog_perc_id_r1', display_name='%id. query gene identical to target Tiger tail seahorse gene', description=''>
hcomes_homolog_goc_score | <biomart.Attribute name='hcomes_homolog_goc_score', display_name='Tiger tail seahorse Gene-order conservation score', description=''>
hcomes_homolog_wga_coverage | <biomart.Attribute name='hcomes_homolog_wga_coverage', display_name='Tiger tail seahorse Whole-genome alignment coverage', description=''>
hcomes_homolog_orthology_confidence | <biomart.Attribute name='hcomes_homolog_orthology_confidence', display_name='Tiger tail seahorse orthology confidence [0 low, 1 high]', description=''>
csemilaevis_homolog_ensembl_gene | <biomart.Attribute name='csemilaevis_homolog_ensembl_gene', display_name='Tongue sole gene stable ID', description=''>
csemilaevis_homolog_associated_gene_name | <biomart.Attribute name='csemilaevis_homolog_associated_gene_name', display_name='Tongue sole gene name', description=''>
csemilaevis_homolog_ensembl_peptide | <biomart.Attribute name='csemilaevis_homolog_ensembl_peptide', display_name='Tongue sole protein or transcript stable ID', description=''>
csemilaevis_homolog_chromosome | <biomart.Attribute name='csemilaevis_homolog_chromosome', display_name='Tongue sole chromosome/scaffold name', description=''>
csemilaevis_homolog_chrom_start | <biomart.Attribute name='csemilaevis_homolog_chrom_start', display_name='Tongue sole chromosome/scaffold start (bp)', description=''>
csemilaevis_homolog_chrom_end | <biomart.Attribute name='csemilaevis_homolog_chrom_end', display_name='Tongue sole chromosome/scaffold end (bp)', description=''>
csemilaevis_homolog_canonical_transcript_protein | <biomart.Attribute name='csemilaevis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
csemilaevis_homolog_subtype | <biomart.Attribute name='csemilaevis_homolog_subtype', display_name='Last common ancestor with Tongue sole', description=''>
csemilaevis_homolog_orthology_type | <biomart.Attribute name='csemilaevis_homolog_orthology_type', display_name='Tongue sole homology type', description=''>
csemilaevis_homolog_perc_id | <biomart.Attribute name='csemilaevis_homolog_perc_id', display_name='%id. target Tongue sole gene identical to query gene', description=''>
csemilaevis_homolog_perc_id_r1 | <biomart.Attribute name='csemilaevis_homolog_perc_id_r1', display_name='%id. query gene identical to target Tongue sole gene', description=''>
csemilaevis_homolog_goc_score | <biomart.Attribute name='csemilaevis_homolog_goc_score', display_name='Tongue sole Gene-order conservation score', description=''>
csemilaevis_homolog_wga_coverage | <biomart.Attribute name='csemilaevis_homolog_wga_coverage', display_name='Tongue sole Whole-genome alignment coverage', description=''>
csemilaevis_homolog_orthology_confidence | <biomart.Attribute name='csemilaevis_homolog_orthology_confidence', display_name='Tongue sole orthology confidence [0 low, 1 high]', description=''>
tbelangeri_homolog_ensembl_gene | <biomart.Attribute name='tbelangeri_homolog_ensembl_gene', display_name='Tree Shrew gene stable ID', description=''>
tbelangeri_homolog_associated_gene_name | <biomart.Attribute name='tbelangeri_homolog_associated_gene_name', display_name='Tree Shrew gene name', description=''>
tbelangeri_homolog_ensembl_peptide | <biomart.Attribute name='tbelangeri_homolog_ensembl_peptide', display_name='Tree Shrew protein or transcript stable ID', description=''>
tbelangeri_homolog_chromosome | <biomart.Attribute name='tbelangeri_homolog_chromosome', display_name='Tree Shrew chromosome/scaffold name', description=''>
tbelangeri_homolog_chrom_start | <biomart.Attribute name='tbelangeri_homolog_chrom_start', display_name='Tree Shrew chromosome/scaffold start (bp)', description=''>
tbelangeri_homolog_chrom_end | <biomart.Attribute name='tbelangeri_homolog_chrom_end', display_name='Tree Shrew chromosome/scaffold end (bp)', description=''>
tbelangeri_homolog_canonical_transcript_protein | <biomart.Attribute name='tbelangeri_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
tbelangeri_homolog_subtype | <biomart.Attribute name='tbelangeri_homolog_subtype', display_name='Last common ancestor with Tree Shrew', description=''>
tbelangeri_homolog_orthology_type | <biomart.Attribute name='tbelangeri_homolog_orthology_type', display_name='Tree Shrew homology type', description=''>
tbelangeri_homolog_perc_id | <biomart.Attribute name='tbelangeri_homolog_perc_id', display_name='%id. target Tree Shrew gene identical to query gene', description=''>
tbelangeri_homolog_perc_id_r1 | <biomart.Attribute name='tbelangeri_homolog_perc_id_r1', display_name='%id. query gene identical to target Tree Shrew gene', description=''>
tbelangeri_homolog_goc_score | <biomart.Attribute name='tbelangeri_homolog_goc_score', display_name='Tree Shrew Gene-order conservation score', description=''>
tbelangeri_homolog_wga_coverage | <biomart.Attribute name='tbelangeri_homolog_wga_coverage', display_name='Tree Shrew Whole-genome alignment coverage', description=''>
tbelangeri_homolog_orthology_confidence | <biomart.Attribute name='tbelangeri_homolog_orthology_confidence', display_name='Tree Shrew orthology confidence [0 low, 1 high]', description=''>
xtropicalis_homolog_ensembl_gene | <biomart.Attribute name='xtropicalis_homolog_ensembl_gene', display_name='Tropical clawed frog gene stable ID', description=''>
xtropicalis_homolog_associated_gene_name | <biomart.Attribute name='xtropicalis_homolog_associated_gene_name', display_name='Tropical clawed frog gene name', description=''>
xtropicalis_homolog_ensembl_peptide | <biomart.Attribute name='xtropicalis_homolog_ensembl_peptide', display_name='Tropical clawed frog protein or transcript stable ID', description=''>
xtropicalis_homolog_chromosome | <biomart.Attribute name='xtropicalis_homolog_chromosome', display_name='Tropical clawed frog chromosome/scaffold name', description=''>
xtropicalis_homolog_chrom_start | <biomart.Attribute name='xtropicalis_homolog_chrom_start', display_name='Tropical clawed frog chromosome/scaffold start (bp)', description=''>
xtropicalis_homolog_chrom_end | <biomart.Attribute name='xtropicalis_homolog_chrom_end', display_name='Tropical clawed frog chromosome/scaffold end (bp)', description=''>
xtropicalis_homolog_canonical_transcript_protein | <biomart.Attribute name='xtropicalis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
xtropicalis_homolog_subtype | <biomart.Attribute name='xtropicalis_homolog_subtype', display_name='Last common ancestor with Tropical clawed frog', description=''>
xtropicalis_homolog_orthology_type | <biomart.Attribute name='xtropicalis_homolog_orthology_type', display_name='Tropical clawed frog homology type', description=''>
xtropicalis_homolog_perc_id | <biomart.Attribute name='xtropicalis_homolog_perc_id', display_name='%id. target Tropical clawed frog gene identical to query gene', description=''>
xtropicalis_homolog_perc_id_r1 | <biomart.Attribute name='xtropicalis_homolog_perc_id_r1', display_name='%id. query gene identical to target Tropical clawed frog gene', description=''>
xtropicalis_homolog_goc_score | <biomart.Attribute name='xtropicalis_homolog_goc_score', display_name='Tropical clawed frog Gene-order conservation score', description=''>
xtropicalis_homolog_wga_coverage | <biomart.Attribute name='xtropicalis_homolog_wga_coverage', display_name='Tropical clawed frog Whole-genome alignment coverage', description=''>
xtropicalis_homolog_orthology_confidence | <biomart.Attribute name='xtropicalis_homolog_orthology_confidence', display_name='Tropical clawed frog orthology confidence [0 low, 1 high]', description=''>
spunctatus_homolog_ensembl_gene | <biomart.Attribute name='spunctatus_homolog_ensembl_gene', display_name='Tuatara gene stable ID', description=''>
spunctatus_homolog_associated_gene_name | <biomart.Attribute name='spunctatus_homolog_associated_gene_name', display_name='Tuatara gene name', description=''>
spunctatus_homolog_ensembl_peptide | <biomart.Attribute name='spunctatus_homolog_ensembl_peptide', display_name='Tuatara protein or transcript stable ID', description=''>
spunctatus_homolog_chromosome | <biomart.Attribute name='spunctatus_homolog_chromosome', display_name='Tuatara chromosome/scaffold name', description=''>
spunctatus_homolog_chrom_start | <biomart.Attribute name='spunctatus_homolog_chrom_start', display_name='Tuatara chromosome/scaffold start (bp)', description=''>
spunctatus_homolog_chrom_end | <biomart.Attribute name='spunctatus_homolog_chrom_end', display_name='Tuatara chromosome/scaffold end (bp)', description=''>
spunctatus_homolog_canonical_transcript_protein | <biomart.Attribute name='spunctatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
spunctatus_homolog_subtype | <biomart.Attribute name='spunctatus_homolog_subtype', display_name='Last common ancestor with Tuatara', description=''>
spunctatus_homolog_orthology_type | <biomart.Attribute name='spunctatus_homolog_orthology_type', display_name='Tuatara homology type', description=''>
spunctatus_homolog_perc_id | <biomart.Attribute name='spunctatus_homolog_perc_id', display_name='%id. target Tuatara gene identical to query gene', description=''>
spunctatus_homolog_perc_id_r1 | <biomart.Attribute name='spunctatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Tuatara gene', description=''>
spunctatus_homolog_goc_score | <biomart.Attribute name='spunctatus_homolog_goc_score', display_name='Tuatara Gene-order conservation score', description=''>
spunctatus_homolog_wga_coverage | <biomart.Attribute name='spunctatus_homolog_wga_coverage', display_name='Tuatara Whole-genome alignment coverage', description=''>
spunctatus_homolog_orthology_confidence | <biomart.Attribute name='spunctatus_homolog_orthology_confidence', display_name='Tuatara orthology confidence [0 low, 1 high]', description=''>
smaximus_homolog_ensembl_gene | <biomart.Attribute name='smaximus_homolog_ensembl_gene', display_name='Turbot gene stable ID', description=''>
smaximus_homolog_associated_gene_name | <biomart.Attribute name='smaximus_homolog_associated_gene_name', display_name='Turbot gene name', description=''>
smaximus_homolog_ensembl_peptide | <biomart.Attribute name='smaximus_homolog_ensembl_peptide', display_name='Turbot protein or transcript stable ID', description=''>
smaximus_homolog_chromosome | <biomart.Attribute name='smaximus_homolog_chromosome', display_name='Turbot chromosome/scaffold name', description=''>
smaximus_homolog_chrom_start | <biomart.Attribute name='smaximus_homolog_chrom_start', display_name='Turbot chromosome/scaffold start (bp)', description=''>
smaximus_homolog_chrom_end | <biomart.Attribute name='smaximus_homolog_chrom_end', display_name='Turbot chromosome/scaffold end (bp)', description=''>
smaximus_homolog_canonical_transcript_protein | <biomart.Attribute name='smaximus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
smaximus_homolog_subtype | <biomart.Attribute name='smaximus_homolog_subtype', display_name='Last common ancestor with Turbot', description=''>
smaximus_homolog_orthology_type | <biomart.Attribute name='smaximus_homolog_orthology_type', display_name='Turbot homology type', description=''>
smaximus_homolog_perc_id | <biomart.Attribute name='smaximus_homolog_perc_id', display_name='%id. target Turbot gene identical to query gene', description=''>
smaximus_homolog_perc_id_r1 | <biomart.Attribute name='smaximus_homolog_perc_id_r1', display_name='%id. query gene identical to target Turbot gene', description=''>
smaximus_homolog_goc_score | <biomart.Attribute name='smaximus_homolog_goc_score', display_name='Turbot Gene-order conservation score', description=''>
smaximus_homolog_wga_coverage | <biomart.Attribute name='smaximus_homolog_wga_coverage', display_name='Turbot Whole-genome alignment coverage', description=''>
smaximus_homolog_orthology_confidence | <biomart.Attribute name='smaximus_homolog_orthology_confidence', display_name='Turbot orthology confidence [0 low, 1 high]', description=''>
mgallopavo_homolog_ensembl_gene | <biomart.Attribute name='mgallopavo_homolog_ensembl_gene', display_name='Turkey gene stable ID', description=''>
mgallopavo_homolog_associated_gene_name | <biomart.Attribute name='mgallopavo_homolog_associated_gene_name', display_name='Turkey gene name', description=''>
mgallopavo_homolog_ensembl_peptide | <biomart.Attribute name='mgallopavo_homolog_ensembl_peptide', display_name='Turkey protein or transcript stable ID', description=''>
mgallopavo_homolog_chromosome | <biomart.Attribute name='mgallopavo_homolog_chromosome', display_name='Turkey chromosome/scaffold name', description=''>
mgallopavo_homolog_chrom_start | <biomart.Attribute name='mgallopavo_homolog_chrom_start', display_name='Turkey chromosome/scaffold start (bp)', description=''>
mgallopavo_homolog_chrom_end | <biomart.Attribute name='mgallopavo_homolog_chrom_end', display_name='Turkey chromosome/scaffold end (bp)', description=''>
mgallopavo_homolog_canonical_transcript_protein | <biomart.Attribute name='mgallopavo_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mgallopavo_homolog_subtype | <biomart.Attribute name='mgallopavo_homolog_subtype', display_name='Last common ancestor with Turkey', description=''>
mgallopavo_homolog_orthology_type | <biomart.Attribute name='mgallopavo_homolog_orthology_type', display_name='Turkey homology type', description=''>
mgallopavo_homolog_perc_id | <biomart.Attribute name='mgallopavo_homolog_perc_id', display_name='%id. target Turkey gene identical to query gene', description=''>
mgallopavo_homolog_perc_id_r1 | <biomart.Attribute name='mgallopavo_homolog_perc_id_r1', display_name='%id. query gene identical to target Turkey gene', description=''>
mgallopavo_homolog_goc_score | <biomart.Attribute name='mgallopavo_homolog_goc_score', display_name='Turkey Gene-order conservation score', description=''>
mgallopavo_homolog_wga_coverage | <biomart.Attribute name='mgallopavo_homolog_wga_coverage', display_name='Turkey Whole-genome alignment coverage', description=''>
mgallopavo_homolog_orthology_confidence | <biomart.Attribute name='mgallopavo_homolog_orthology_confidence', display_name='Turkey orthology confidence [0 low, 1 high]', description=''>
ptephrosceles_homolog_ensembl_gene | <biomart.Attribute name='ptephrosceles_homolog_ensembl_gene', display_name='Ugandan red Colobus gene stable ID', description=''>
ptephrosceles_homolog_associated_gene_name | <biomart.Attribute name='ptephrosceles_homolog_associated_gene_name', display_name='Ugandan red Colobus gene name', description=''>
ptephrosceles_homolog_ensembl_peptide | <biomart.Attribute name='ptephrosceles_homolog_ensembl_peptide', display_name='Ugandan red Colobus protein or transcript stable ID', description=''>
ptephrosceles_homolog_chromosome | <biomart.Attribute name='ptephrosceles_homolog_chromosome', display_name='Ugandan red Colobus chromosome/scaffold name', description=''>
ptephrosceles_homolog_chrom_start | <biomart.Attribute name='ptephrosceles_homolog_chrom_start', display_name='Ugandan red Colobus chromosome/scaffold start (bp)', description=''>
ptephrosceles_homolog_chrom_end | <biomart.Attribute name='ptephrosceles_homolog_chrom_end', display_name='Ugandan red Colobus chromosome/scaffold end (bp)', description=''>
ptephrosceles_homolog_canonical_transcript_protein | <biomart.Attribute name='ptephrosceles_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ptephrosceles_homolog_subtype | <biomart.Attribute name='ptephrosceles_homolog_subtype', display_name='Last common ancestor with Ugandan red Colobus', description=''>
ptephrosceles_homolog_orthology_type | <biomart.Attribute name='ptephrosceles_homolog_orthology_type', display_name='Ugandan red Colobus homology type', description=''>
ptephrosceles_homolog_perc_id | <biomart.Attribute name='ptephrosceles_homolog_perc_id', display_name='%id. target Ugandan red Colobus gene identical to query gene', description=''>
ptephrosceles_homolog_perc_id_r1 | <biomart.Attribute name='ptephrosceles_homolog_perc_id_r1', display_name='%id. query gene identical to target Ugandan red Colobus gene', description=''>
ptephrosceles_homolog_goc_score | <biomart.Attribute name='ptephrosceles_homolog_goc_score', display_name='Ugandan red Colobus Gene-order conservation score', description=''>
ptephrosceles_homolog_wga_coverage | <biomart.Attribute name='ptephrosceles_homolog_wga_coverage', display_name='Ugandan red Colobus Whole-genome alignment coverage', description=''>
ptephrosceles_homolog_orthology_confidence | <biomart.Attribute name='ptephrosceles_homolog_orthology_confidence', display_name='Ugandan red Colobus orthology confidence [0 low, 1 high]', description=''>
ngalili_homolog_ensembl_gene | <biomart.Attribute name='ngalili_homolog_ensembl_gene', display_name='Upper Galilee mountains blind mole rat gene stable ID', description=''>
ngalili_homolog_associated_gene_name | <biomart.Attribute name='ngalili_homolog_associated_gene_name', display_name='Upper Galilee mountains blind mole rat gene name', description=''>
ngalili_homolog_ensembl_peptide | <biomart.Attribute name='ngalili_homolog_ensembl_peptide', display_name='Upper Galilee mountains blind mole rat protein or transcript stable ID', description=''>
ngalili_homolog_chromosome | <biomart.Attribute name='ngalili_homolog_chromosome', display_name='Upper Galilee mountains blind mole rat chromosome/scaffold name', description=''>
ngalili_homolog_chrom_start | <biomart.Attribute name='ngalili_homolog_chrom_start', display_name='Upper Galilee mountains blind mole rat chromosome/scaffold start (bp)', description=''>
ngalili_homolog_chrom_end | <biomart.Attribute name='ngalili_homolog_chrom_end', display_name='Upper Galilee mountains blind mole rat chromosome/scaffold end (bp)', description=''>
ngalili_homolog_canonical_transcript_protein | <biomart.Attribute name='ngalili_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
ngalili_homolog_subtype | <biomart.Attribute name='ngalili_homolog_subtype', display_name='Last common ancestor with Upper Galilee mountains blind mole rat', description=''>
ngalili_homolog_orthology_type | <biomart.Attribute name='ngalili_homolog_orthology_type', display_name='Upper Galilee mountains blind mole rat homology type', description=''>
ngalili_homolog_perc_id | <biomart.Attribute name='ngalili_homolog_perc_id', display_name='%id. target Upper Galilee mountains blind mole rat gene identical to query gene', description=''>
ngalili_homolog_perc_id_r1 | <biomart.Attribute name='ngalili_homolog_perc_id_r1', display_name='%id. query gene identical to target Upper Galilee mountains blind mole rat gene', description=''>
ngalili_homolog_goc_score | <biomart.Attribute name='ngalili_homolog_goc_score', display_name='Upper Galilee mountains blind mole rat Gene-order conservation score', description=''>
ngalili_homolog_wga_coverage | <biomart.Attribute name='ngalili_homolog_wga_coverage', display_name='Upper Galilee mountains blind mole rat Whole-genome alignment coverage', description=''>
ngalili_homolog_orthology_confidence | <biomart.Attribute name='ngalili_homolog_orthology_confidence', display_name='Upper Galilee mountains blind mole rat orthology confidence [0 low, 1 high]', description=''>
csabaeus_homolog_ensembl_gene | <biomart.Attribute name='csabaeus_homolog_ensembl_gene', display_name='Vervet-AGM gene stable ID', description=''>
csabaeus_homolog_associated_gene_name | <biomart.Attribute name='csabaeus_homolog_associated_gene_name', display_name='Vervet-AGM gene name', description=''>
csabaeus_homolog_ensembl_peptide | <biomart.Attribute name='csabaeus_homolog_ensembl_peptide', display_name='Vervet-AGM protein or transcript stable ID', description=''>
csabaeus_homolog_chromosome | <biomart.Attribute name='csabaeus_homolog_chromosome', display_name='Vervet-AGM chromosome/scaffold name', description=''>
csabaeus_homolog_chrom_start | <biomart.Attribute name='csabaeus_homolog_chrom_start', display_name='Vervet-AGM chromosome/scaffold start (bp)', description=''>
csabaeus_homolog_chrom_end | <biomart.Attribute name='csabaeus_homolog_chrom_end', display_name='Vervet-AGM chromosome/scaffold end (bp)', description=''>
csabaeus_homolog_canonical_transcript_protein | <biomart.Attribute name='csabaeus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
csabaeus_homolog_subtype | <biomart.Attribute name='csabaeus_homolog_subtype', display_name='Last common ancestor with Vervet-AGM', description=''>
csabaeus_homolog_orthology_type | <biomart.Attribute name='csabaeus_homolog_orthology_type', display_name='Vervet-AGM homology type', description=''>
csabaeus_homolog_perc_id | <biomart.Attribute name='csabaeus_homolog_perc_id', display_name='%id. target Vervet-AGM gene identical to query gene', description=''>
csabaeus_homolog_perc_id_r1 | <biomart.Attribute name='csabaeus_homolog_perc_id_r1', display_name='%id. query gene identical to target Vervet-AGM gene', description=''>
csabaeus_homolog_goc_score | <biomart.Attribute name='csabaeus_homolog_goc_score', display_name='Vervet-AGM Gene-order conservation score', description=''>
csabaeus_homolog_wga_coverage | <biomart.Attribute name='csabaeus_homolog_wga_coverage', display_name='Vervet-AGM Whole-genome alignment coverage', description=''>
csabaeus_homolog_orthology_confidence | <biomart.Attribute name='csabaeus_homolog_orthology_confidence', display_name='Vervet-AGM orthology confidence [0 low, 1 high]', description=''>
neugenii_homolog_ensembl_gene | <biomart.Attribute name='neugenii_homolog_ensembl_gene', display_name='Wallaby gene stable ID', description=''>
neugenii_homolog_associated_gene_name | <biomart.Attribute name='neugenii_homolog_associated_gene_name', display_name='Wallaby gene name', description=''>
neugenii_homolog_ensembl_peptide | <biomart.Attribute name='neugenii_homolog_ensembl_peptide', display_name='Wallaby protein or transcript stable ID', description=''>
neugenii_homolog_chromosome | <biomart.Attribute name='neugenii_homolog_chromosome', display_name='Wallaby chromosome/scaffold name', description=''>
neugenii_homolog_chrom_start | <biomart.Attribute name='neugenii_homolog_chrom_start', display_name='Wallaby chromosome/scaffold start (bp)', description=''>
neugenii_homolog_chrom_end | <biomart.Attribute name='neugenii_homolog_chrom_end', display_name='Wallaby chromosome/scaffold end (bp)', description=''>
neugenii_homolog_canonical_transcript_protein | <biomart.Attribute name='neugenii_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
neugenii_homolog_subtype | <biomart.Attribute name='neugenii_homolog_subtype', display_name='Last common ancestor with Wallaby', description=''>
neugenii_homolog_orthology_type | <biomart.Attribute name='neugenii_homolog_orthology_type', display_name='Wallaby homology type', description=''>
neugenii_homolog_perc_id | <biomart.Attribute name='neugenii_homolog_perc_id', display_name='%id. target Wallaby gene identical to query gene', description=''>
neugenii_homolog_perc_id_r1 | <biomart.Attribute name='neugenii_homolog_perc_id_r1', display_name='%id. query gene identical to target Wallaby gene', description=''>
neugenii_homolog_goc_score | <biomart.Attribute name='neugenii_homolog_goc_score', display_name='Wallaby Gene-order conservation score', description=''>
neugenii_homolog_wga_coverage | <biomart.Attribute name='neugenii_homolog_wga_coverage', display_name='Wallaby Whole-genome alignment coverage', description=''>
neugenii_homolog_orthology_confidence | <biomart.Attribute name='neugenii_homolog_orthology_confidence', display_name='Wallaby orthology confidence [0 low, 1 high]', description=''>
bmutus_homolog_ensembl_gene | <biomart.Attribute name='bmutus_homolog_ensembl_gene', display_name='Wild yak gene stable ID', description=''>
bmutus_homolog_associated_gene_name | <biomart.Attribute name='bmutus_homolog_associated_gene_name', display_name='Wild yak gene name', description=''>
bmutus_homolog_ensembl_peptide | <biomart.Attribute name='bmutus_homolog_ensembl_peptide', display_name='Wild yak protein or transcript stable ID', description=''>
bmutus_homolog_chromosome | <biomart.Attribute name='bmutus_homolog_chromosome', display_name='Wild yak chromosome/scaffold name', description=''>
bmutus_homolog_chrom_start | <biomart.Attribute name='bmutus_homolog_chrom_start', display_name='Wild yak chromosome/scaffold start (bp)', description=''>
bmutus_homolog_chrom_end | <biomart.Attribute name='bmutus_homolog_chrom_end', display_name='Wild yak chromosome/scaffold end (bp)', description=''>
bmutus_homolog_canonical_transcript_protein | <biomart.Attribute name='bmutus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
bmutus_homolog_subtype | <biomart.Attribute name='bmutus_homolog_subtype', display_name='Last common ancestor with Wild yak', description=''>
bmutus_homolog_orthology_type | <biomart.Attribute name='bmutus_homolog_orthology_type', display_name='Wild yak homology type', description=''>
bmutus_homolog_perc_id | <biomart.Attribute name='bmutus_homolog_perc_id', display_name='%id. target Wild yak gene identical to query gene', description=''>
bmutus_homolog_perc_id_r1 | <biomart.Attribute name='bmutus_homolog_perc_id_r1', display_name='%id. query gene identical to target Wild yak gene', description=''>
bmutus_homolog_goc_score | <biomart.Attribute name='bmutus_homolog_goc_score', display_name='Wild yak Gene-order conservation score', description=''>
bmutus_homolog_wga_coverage | <biomart.Attribute name='bmutus_homolog_wga_coverage', display_name='Wild yak Whole-genome alignment coverage', description=''>
bmutus_homolog_orthology_confidence | <biomart.Attribute name='bmutus_homolog_orthology_confidence', display_name='Wild yak orthology confidence [0 low, 1 high]', description=''>
sldorsalis_homolog_ensembl_gene | <biomart.Attribute name='sldorsalis_homolog_ensembl_gene', display_name='Yellowtail amberjack gene stable ID', description=''>
sldorsalis_homolog_associated_gene_name | <biomart.Attribute name='sldorsalis_homolog_associated_gene_name', display_name='Yellowtail amberjack gene name', description=''>
sldorsalis_homolog_ensembl_peptide | <biomart.Attribute name='sldorsalis_homolog_ensembl_peptide', display_name='Yellowtail amberjack protein or transcript stable ID', description=''>
sldorsalis_homolog_chromosome | <biomart.Attribute name='sldorsalis_homolog_chromosome', display_name='Yellowtail amberjack chromosome/scaffold name', description=''>
sldorsalis_homolog_chrom_start | <biomart.Attribute name='sldorsalis_homolog_chrom_start', display_name='Yellowtail amberjack chromosome/scaffold start (bp)', description=''>
sldorsalis_homolog_chrom_end | <biomart.Attribute name='sldorsalis_homolog_chrom_end', display_name='Yellowtail amberjack chromosome/scaffold end (bp)', description=''>
sldorsalis_homolog_canonical_transcript_protein | <biomart.Attribute name='sldorsalis_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
sldorsalis_homolog_subtype | <biomart.Attribute name='sldorsalis_homolog_subtype', display_name='Last common ancestor with Yellowtail amberjack', description=''>
sldorsalis_homolog_orthology_type | <biomart.Attribute name='sldorsalis_homolog_orthology_type', display_name='Yellowtail amberjack homology type', description=''>
sldorsalis_homolog_perc_id | <biomart.Attribute name='sldorsalis_homolog_perc_id', display_name='%id. target Yellowtail amberjack gene identical to query gene', description=''>
sldorsalis_homolog_perc_id_r1 | <biomart.Attribute name='sldorsalis_homolog_perc_id_r1', display_name='%id. query gene identical to target Yellowtail amberjack gene', description=''>
sldorsalis_homolog_goc_score | <biomart.Attribute name='sldorsalis_homolog_goc_score', display_name='Yellowtail amberjack Gene-order conservation score', description=''>
sldorsalis_homolog_wga_coverage | <biomart.Attribute name='sldorsalis_homolog_wga_coverage', display_name='Yellowtail amberjack Whole-genome alignment coverage', description=''>
sldorsalis_homolog_orthology_confidence | <biomart.Attribute name='sldorsalis_homolog_orthology_confidence', display_name='Yellowtail amberjack orthology confidence [0 low, 1 high]', description=''>
tguttata_homolog_ensembl_gene | <biomart.Attribute name='tguttata_homolog_ensembl_gene', display_name='Zebra finch gene stable ID', description=''>
tguttata_homolog_associated_gene_name | <biomart.Attribute name='tguttata_homolog_associated_gene_name', display_name='Zebra finch gene name', description=''>
tguttata_homolog_ensembl_peptide | <biomart.Attribute name='tguttata_homolog_ensembl_peptide', display_name='Zebra finch protein or transcript stable ID', description=''>
tguttata_homolog_chromosome | <biomart.Attribute name='tguttata_homolog_chromosome', display_name='Zebra finch chromosome/scaffold name', description=''>
tguttata_homolog_chrom_start | <biomart.Attribute name='tguttata_homolog_chrom_start', display_name='Zebra finch chromosome/scaffold start (bp)', description=''>
tguttata_homolog_chrom_end | <biomart.Attribute name='tguttata_homolog_chrom_end', display_name='Zebra finch chromosome/scaffold end (bp)', description=''>
tguttata_homolog_canonical_transcript_protein | <biomart.Attribute name='tguttata_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
tguttata_homolog_subtype | <biomart.Attribute name='tguttata_homolog_subtype', display_name='Last common ancestor with Zebra finch', description=''>
tguttata_homolog_orthology_type | <biomart.Attribute name='tguttata_homolog_orthology_type', display_name='Zebra finch homology type', description=''>
tguttata_homolog_perc_id | <biomart.Attribute name='tguttata_homolog_perc_id', display_name='%id. target Zebra finch gene identical to query gene', description=''>
tguttata_homolog_perc_id_r1 | <biomart.Attribute name='tguttata_homolog_perc_id_r1', display_name='%id. query gene identical to target Zebra finch gene', description=''>
tguttata_homolog_goc_score | <biomart.Attribute name='tguttata_homolog_goc_score', display_name='Zebra finch Gene-order conservation score', description=''>
tguttata_homolog_wga_coverage | <biomart.Attribute name='tguttata_homolog_wga_coverage', display_name='Zebra finch Whole-genome alignment coverage', description=''>
tguttata_homolog_orthology_confidence | <biomart.Attribute name='tguttata_homolog_orthology_confidence', display_name='Zebra finch orthology confidence [0 low, 1 high]', description=''>
mzebra_homolog_ensembl_gene | <biomart.Attribute name='mzebra_homolog_ensembl_gene', display_name='Zebra mbuna gene stable ID', description=''>
mzebra_homolog_associated_gene_name | <biomart.Attribute name='mzebra_homolog_associated_gene_name', display_name='Zebra mbuna gene name', description=''>
mzebra_homolog_ensembl_peptide | <biomart.Attribute name='mzebra_homolog_ensembl_peptide', display_name='Zebra mbuna protein or transcript stable ID', description=''>
mzebra_homolog_chromosome | <biomart.Attribute name='mzebra_homolog_chromosome', display_name='Zebra mbuna chromosome/scaffold name', description=''>
mzebra_homolog_chrom_start | <biomart.Attribute name='mzebra_homolog_chrom_start', display_name='Zebra mbuna chromosome/scaffold start (bp)', description=''>
mzebra_homolog_chrom_end | <biomart.Attribute name='mzebra_homolog_chrom_end', display_name='Zebra mbuna chromosome/scaffold end (bp)', description=''>
mzebra_homolog_canonical_transcript_protein | <biomart.Attribute name='mzebra_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
mzebra_homolog_subtype | <biomart.Attribute name='mzebra_homolog_subtype', display_name='Last common ancestor with Zebra mbuna', description=''>
mzebra_homolog_orthology_type | <biomart.Attribute name='mzebra_homolog_orthology_type', display_name='Zebra mbuna homology type', description=''>
mzebra_homolog_perc_id | <biomart.Attribute name='mzebra_homolog_perc_id', display_name='%id. target Zebra mbuna gene identical to query gene', description=''>
mzebra_homolog_perc_id_r1 | <biomart.Attribute name='mzebra_homolog_perc_id_r1', display_name='%id. query gene identical to target Zebra mbuna gene', description=''>
mzebra_homolog_goc_score | <biomart.Attribute name='mzebra_homolog_goc_score', display_name='Zebra mbuna Gene-order conservation score', description=''>
mzebra_homolog_wga_coverage | <biomart.Attribute name='mzebra_homolog_wga_coverage', display_name='Zebra mbuna Whole-genome alignment coverage', description=''>
mzebra_homolog_orthology_confidence | <biomart.Attribute name='mzebra_homolog_orthology_confidence', display_name='Zebra mbuna orthology confidence [0 low, 1 high]', description=''>
drerio_homolog_ensembl_gene | <biomart.Attribute name='drerio_homolog_ensembl_gene', display_name='Zebrafish gene stable ID', description=''>
drerio_homolog_associated_gene_name | <biomart.Attribute name='drerio_homolog_associated_gene_name', display_name='Zebrafish gene name', description=''>
drerio_homolog_ensembl_peptide | <biomart.Attribute name='drerio_homolog_ensembl_peptide', display_name='Zebrafish protein or transcript stable ID', description=''>
drerio_homolog_chromosome | <biomart.Attribute name='drerio_homolog_chromosome', display_name='Zebrafish chromosome/scaffold name', description=''>
drerio_homolog_chrom_start | <biomart.Attribute name='drerio_homolog_chrom_start', display_name='Zebrafish chromosome/scaffold start (bp)', description=''>
drerio_homolog_chrom_end | <biomart.Attribute name='drerio_homolog_chrom_end', display_name='Zebrafish chromosome/scaffold end (bp)', description=''>
drerio_homolog_canonical_transcript_protein | <biomart.Attribute name='drerio_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
drerio_homolog_subtype | <biomart.Attribute name='drerio_homolog_subtype', display_name='Last common ancestor with Zebrafish', description=''>
drerio_homolog_orthology_type | <biomart.Attribute name='drerio_homolog_orthology_type', display_name='Zebrafish homology type', description=''>
drerio_homolog_perc_id | <biomart.Attribute name='drerio_homolog_perc_id', display_name='%id. target Zebrafish gene identical to query gene', description=''>
drerio_homolog_perc_id_r1 | <biomart.Attribute name='drerio_homolog_perc_id_r1', display_name='%id. query gene identical to target Zebrafish gene', description=''>
drerio_homolog_goc_score | <biomart.Attribute name='drerio_homolog_goc_score', display_name='Zebrafish Gene-order conservation score', description=''>
drerio_homolog_wga_coverage | <biomart.Attribute name='drerio_homolog_wga_coverage', display_name='Zebrafish Whole-genome alignment coverage', description=''>
drerio_homolog_orthology_confidence | <biomart.Attribute name='drerio_homolog_orthology_confidence', display_name='Zebrafish orthology confidence [0 low, 1 high]', description=''>
marmatus_homolog_ensembl_gene | <biomart.Attribute name='marmatus_homolog_ensembl_gene', display_name='Zig-zag eel gene stable ID', description=''>
marmatus_homolog_associated_gene_name | <biomart.Attribute name='marmatus_homolog_associated_gene_name', display_name='Zig-zag eel gene name', description=''>
marmatus_homolog_ensembl_peptide | <biomart.Attribute name='marmatus_homolog_ensembl_peptide', display_name='Zig-zag eel protein or transcript stable ID', description=''>
marmatus_homolog_chromosome | <biomart.Attribute name='marmatus_homolog_chromosome', display_name='Zig-zag eel chromosome/scaffold name', description=''>
marmatus_homolog_chrom_start | <biomart.Attribute name='marmatus_homolog_chrom_start', display_name='Zig-zag eel chromosome/scaffold start (bp)', description=''>
marmatus_homolog_chrom_end | <biomart.Attribute name='marmatus_homolog_chrom_end', display_name='Zig-zag eel chromosome/scaffold end (bp)', description=''>
marmatus_homolog_canonical_transcript_protein | <biomart.Attribute name='marmatus_homolog_canonical_transcript_protein', display_name='Query protein or transcript ID', description=''>
marmatus_homolog_subtype | <biomart.Attribute name='marmatus_homolog_subtype', display_name='Last common ancestor with Zig-zag eel', description=''>
marmatus_homolog_orthology_type | <biomart.Attribute name='marmatus_homolog_orthology_type', display_name='Zig-zag eel homology type', description=''>
marmatus_homolog_perc_id | <biomart.Attribute name='marmatus_homolog_perc_id', display_name='%id. target Zig-zag eel gene identical to query gene', description=''>
marmatus_homolog_perc_id_r1 | <biomart.Attribute name='marmatus_homolog_perc_id_r1', display_name='%id. query gene identical to target Zig-zag eel gene', description=''>
marmatus_homolog_goc_score | <biomart.Attribute name='marmatus_homolog_goc_score', display_name='Zig-zag eel Gene-order conservation score', description=''>
marmatus_homolog_wga_coverage | <biomart.Attribute name='marmatus_homolog_wga_coverage', display_name='Zig-zag eel Whole-genome alignment coverage', description=''>
marmatus_homolog_orthology_confidence | <biomart.Attribute name='marmatus_homolog_orthology_confidence', display_name='Zig-zag eel orthology confidence [0 low, 1 high]', description=''>
hsapiens_paralog_ensembl_gene | <biomart.Attribute name='hsapiens_paralog_ensembl_gene', display_name='Human paralogue gene stable ID', description=''>
hsapiens_paralog_associated_gene_name | <biomart.Attribute name='hsapiens_paralog_associated_gene_name', display_name='Human paralogue associated gene name', description=''>
hsapiens_paralog_ensembl_peptide | <biomart.Attribute name='hsapiens_paralog_ensembl_peptide', display_name='Human paralogue protein or transcript ID', description=''>
hsapiens_paralog_chromosome | <biomart.Attribute name='hsapiens_paralog_chromosome', display_name='Human paralogue chromosome/scaffold name', description=''>
hsapiens_paralog_chrom_start | <biomart.Attribute name='hsapiens_paralog_chrom_start', display_name='Human paralogue chromosome/scaffold start (bp)', description=''>
hsapiens_paralog_chrom_end | <biomart.Attribute name='hsapiens_paralog_chrom_end', display_name='Human paralogue chromosome/scaffold end (bp)', description=''>
hsapiens_paralog_canonical_transcript_protein | <biomart.Attribute name='hsapiens_paralog_canonical_transcript_protein', display_name='Paralogue query protein or transcript ID', description=''>
hsapiens_paralog_subtype | <biomart.Attribute name='hsapiens_paralog_subtype', display_name='Paralogue last common ancestor with Human', description=''>
hsapiens_paralog_orthology_type | <biomart.Attribute name='hsapiens_paralog_orthology_type', display_name='Human paralogue homology type', description=''>
hsapiens_paralog_perc_id | <biomart.Attribute name='hsapiens_paralog_perc_id', display_name='Paralogue %id. target Human gene identical to query gene', description=''>
hsapiens_paralog_perc_id_r1 | <biomart.Attribute name='hsapiens_paralog_perc_id_r1', display_name='Paralogue %id. query gene identical to target Human gene', description=''>
snp_ensembl_gene_id | <biomart.Attribute name='snp_ensembl_gene_id', display_name='', description=''>
snp_gene_stable_id_version | <biomart.Attribute name='snp_gene_stable_id_version', display_name='', description=''>
snp_gene_version | <biomart.Attribute name='snp_gene_version', display_name='', description=''>
snp_ensembl_transcript_id | <biomart.Attribute name='snp_ensembl_transcript_id', display_name='', description=''>
snp_transcript_stable_id_version | <biomart.Attribute name='snp_transcript_stable_id_version', display_name='', description=''>
snp_transcript_version | <biomart.Attribute name='snp_transcript_version', display_name='', description=''>
snp_ensembl_peptide_id | <biomart.Attribute name='snp_ensembl_peptide_id', display_name='', description=''>
snp_translation_stable_id_version | <biomart.Attribute name='snp_translation_stable_id_version', display_name='', description=''>
snp_peptide_version | <biomart.Attribute name='snp_peptide_version', display_name='', description=''>
sequence_canonical_transcript_id | <biomart.Attribute name='sequence_canonical_transcript_id', display_name='', description=''>
snp_chromosome_name | <biomart.Attribute name='snp_chromosome_name', display_name='', description=''>
snp_start_position | <biomart.Attribute name='snp_start_position', display_name='', description=''>
snp_end_position | <biomart.Attribute name='snp_end_position', display_name='', description=''>
snp_strand | <biomart.Attribute name='snp_strand', display_name='', description=''>
snp_band | <biomart.Attribute name='snp_band', display_name='', description=''>
snp_external_gene_name | <biomart.Attribute name='snp_external_gene_name', display_name='', description=''>
snp_external_gene_source | <biomart.Attribute name='snp_external_gene_source', display_name='', description=''>
snp_ensembl_CDS_length | <biomart.Attribute name='snp_ensembl_CDS_length', display_name='', description=''>
snp_ensembl_cDNA_length | <biomart.Attribute name='snp_ensembl_cDNA_length', display_name='', description=''>
snp_ensembl_peptide_length | <biomart.Attribute name='snp_ensembl_peptide_length', display_name='', description=''>
snp_transcript_count | <biomart.Attribute name='snp_transcript_count', display_name='', description=''>
snp_percentage_gc_content | <biomart.Attribute name='snp_percentage_gc_content', display_name='', description=''>
snp_description | <biomart.Attribute name='snp_description', display_name='', description=''>
variation_name | <biomart.Attribute name='variation_name', display_name='Variant name', description=''>
germ_line_variation_source | <biomart.Attribute name='germ_line_variation_source', display_name='Variant source', description=''>
source_description | <biomart.Attribute name='source_description', display_name='Variant source description', description=''>
allele | <biomart.Attribute name='allele', display_name='Variant alleles', description='Snp Allele'>
validated | <biomart.Attribute name='validated', display_name='Variant supporting evidence', description=''>
mapweight | <biomart.Attribute name='mapweight', display_name='Mapweight', description=''>
minor_allele | <biomart.Attribute name='minor_allele', display_name='Minor allele', description=''>
minor_allele_freq | <biomart.Attribute name='minor_allele_freq', display_name='Minor allele frequency', description=''>
minor_allele_count | <biomart.Attribute name='minor_allele_count', display_name='Minor allele count', description=''>
clinical_significance | <biomart.Attribute name='clinical_significance', display_name='Clinical significance', description=''>
transcript_location | <biomart.Attribute name='transcript_location', display_name='Transcript location (bp)', description='Start of Variant on Transcript in Chromosomal Coordinates'>
snp_chromosome_strand | <biomart.Attribute name='snp_chromosome_strand', display_name='Variant chromosome Strand', description='Orientation of the SNP on the Transcript'>
peptide_location | <biomart.Attribute name='peptide_location', display_name='Protein location (aa)', description='Location of the SNP in the peptide'>
chromosome_start | <biomart.Attribute name='chromosome_start', display_name='chromosome/scaffold position start (bp)', description='Start of the Variant on the chromosome/scaffold in Chromosomal Coordinates'>
chromosome_end | <biomart.Attribute name='chromosome_end', display_name='Chromosome/scaffold position end (bp)', description='End of the Variant on the chromosome/scaffold in Chromosomal Coordinates'>
polyphen_prediction_2076 | <biomart.Attribute name='polyphen_prediction_2076', display_name='PolyPhen prediction', description=''>
polyphen_score_2076 | <biomart.Attribute name='polyphen_score_2076', display_name='PolyPhen score', description=''>
sift_prediction_2076 | <biomart.Attribute name='sift_prediction_2076', display_name='SIFT prediction', description=''>
sift_score_2076 | <biomart.Attribute name='sift_score_2076', display_name='SIFT score', description=''>
distance_to_transcript_2076 | <biomart.Attribute name='distance_to_transcript_2076', display_name='Distance to transcript', description=''>
cds_start_2076 | <biomart.Attribute name='cds_start_2076', display_name='CDS start', description=''>
cds_end_2076 | <biomart.Attribute name='cds_end_2076', display_name='CDS end', description=''>
peptide_shift | <biomart.Attribute name='peptide_shift', display_name='Protein allele', description=''>
synonymous_status | <biomart.Attribute name='synonymous_status', display_name='Variant consequence', description=''>
allele_string_2076 | <biomart.Attribute name='allele_string_2076', display_name='Consequence specific allele', description=''>
snp_som_ensembl_gene_id | <biomart.Attribute name='snp_som_ensembl_gene_id', display_name='', description=''>
snp_som_gene_stable_id_version | <biomart.Attribute name='snp_som_gene_stable_id_version', display_name='', description=''>
snp_som_gene_version | <biomart.Attribute name='snp_som_gene_version', display_name='', description=''>
snp_som_ensembl_transcript_id | <biomart.Attribute name='snp_som_ensembl_transcript_id', display_name='', description=''>
snp_som_transcript_stable_id_version | <biomart.Attribute name='snp_som_transcript_stable_id_version', display_name='', description=''>
snp_som_transcript_version | <biomart.Attribute name='snp_som_transcript_version', display_name='', description=''>
snp_som_ensembl_peptide_id | <biomart.Attribute name='snp_som_ensembl_peptide_id', display_name='', description=''>
snp_som_translation_stable_id_version | <biomart.Attribute name='snp_som_translation_stable_id_version', display_name='', description=''>
snp_som_peptide_version | <biomart.Attribute name='snp_som_peptide_version', display_name='', description=''>
sequence_som_canonical_transcript_id | <biomart.Attribute name='sequence_som_canonical_transcript_id', display_name='', description=''>
snp_som_chromosome_name | <biomart.Attribute name='snp_som_chromosome_name', display_name='', description=''>
snp_som_start_position | <biomart.Attribute name='snp_som_start_position', display_name='', description=''>
snp_som_end_position | <biomart.Attribute name='snp_som_end_position', display_name='', description=''>
snp_som_strand | <biomart.Attribute name='snp_som_strand', display_name='', description=''>
snp_som_band | <biomart.Attribute name='snp_som_band', display_name='', description=''>
snp_som_external_gene_name | <biomart.Attribute name='snp_som_external_gene_name', display_name='', description=''>
snp_som_external_gene_source | <biomart.Attribute name='snp_som_external_gene_source', display_name='', description=''>
snp_som_ensembl_CDS_length | <biomart.Attribute name='snp_som_ensembl_CDS_length', display_name='', description=''>
snp_som_ensembl_cDNA_length | <biomart.Attribute name='snp_som_ensembl_cDNA_length', display_name='', description=''>
snp_som_ensembl_peptide_length | <biomart.Attribute name='snp_som_ensembl_peptide_length', display_name='', description=''>
snp_som_transcript_count | <biomart.Attribute name='snp_som_transcript_count', display_name='', description=''>
snp_som_percentage_gc_content | <biomart.Attribute name='snp_som_percentage_gc_content', display_name='', description=''>
snp_som_description | <biomart.Attribute name='snp_som_description', display_name='', description=''>
somatic_variation_name | <biomart.Attribute name='somatic_variation_name', display_name='Variant name', description=''>
somatic_source_name | <biomart.Attribute name='somatic_source_name', display_name='Variant source', description=''>
somatic_source_description | <biomart.Attribute name='somatic_source_description', display_name='Variant source description', description=''>
somatic_allele | <biomart.Attribute name='somatic_allele', display_name='Variant alleles', description='Snp Allele'>
somatic_validated | <biomart.Attribute name='somatic_validated', display_name='Variant supporting evidence', description=''>
somatic_mapweight | <biomart.Attribute name='somatic_mapweight', display_name='Mapweight', description=''>
somatic_transcript_location | <biomart.Attribute name='somatic_transcript_location', display_name='Transcript location (bp)', description='Start of Variant on Transcript in Chromosomal Coordinates'>
somatic_snp_chromosome_strand | <biomart.Attribute name='somatic_snp_chromosome_strand', display_name='Variant chromosome/scaffold strand', description='Orientation of the Variant on the Transcript'>
somatic_peptide_location | <biomart.Attribute name='somatic_peptide_location', display_name='Protein location (aa)', description='Location of the Variant in the peptide'>
somatic_chromosome_start | <biomart.Attribute name='somatic_chromosome_start', display_name='Chromosome/scaffold position start (bp)', description='Start of the Variant on the chromosome/scaffold in Chromosomal Coordinates'>
somatic_chromosome_end | <biomart.Attribute name='somatic_chromosome_end', display_name='Chromosome/scaffold position end (bp)', description='End of the Variant on the chromosome/scaffold in Chromosomal Coordinates'>
mart_transcript_variation_som__dm_distance_to_transcript_2076 | <biomart.Attribute name='mart_transcript_variation_som__dm_distance_to_transcript_2076', display_name='Distance to transcript', description=''>
somatic_cds_start_2076 | <biomart.Attribute name='somatic_cds_start_2076', display_name='CDS start', description=''>
somatic_cds_end_2076 | <biomart.Attribute name='somatic_cds_end_2076', display_name='CDS end', description=''>
somatic_synonymous_status | <biomart.Attribute name='somatic_synonymous_status', display_name='Variant consequence', description=''>
mart_transcript_variation_som__dm_allele_string_2076 | <biomart.Attribute name='mart_transcript_variation_som__dm_allele_string_2076', display_name='Consequence specific allele', description=''>
seq_edits | <biomart.Attribute name='seq_edits', display_name='seq_edits', description=''>
rna_seq_edits | <biomart.Attribute name='rna_seq_edits', display_name='rna_seq_edits', description=''>
codon_table_id | <biomart.Attribute name='codon_table_id', display_name='codon_table_id', description=''>
end_exon_id | <biomart.Attribute name='end_exon_id', display_name='End exon id', description=''>
start_exon_id | <biomart.Attribute name='start_exon_id', display_name='Start exon id', description=''>
exon_id | <biomart.Attribute name='exon_id', display_name='exon_id', description=''>
coding_end_offset | <biomart.Attribute name='coding_end_offset', display_name='Coding end offset', description=''>
coding_start_offset | <biomart.Attribute name='coding_start_offset', display_name='Coding start offset', description=''>
go | <biomart.Attribute name='go', display_name='Dbprimary acc', description=''>
transcript_id_key | <biomart.Attribute name='transcript_id_key', display_name='', description=''>
structure_gene_id | <biomart.Attribute name='structure_gene_id', display_name='', description=''>
structure_transcript_id | <biomart.Attribute name='structure_transcript_id', display_name='', description=''>
structure_rank | <biomart.Attribute name='structure_rank', display_name='', description=''>
transcript_exon_intron | <biomart.Attribute name='transcript_exon_intron', display_name='', description=''>
gene_exon_intron | <biomart.Attribute name='gene_exon_intron', display_name='', description=''>
transcript_flank | <biomart.Attribute name='transcript_flank', display_name='', description=''>
gene_flank | <biomart.Attribute name='gene_flank', display_name='', description=''>
coding_transcript_flank | <biomart.Attribute name='coding_transcript_flank', display_name='', description=''>
coding_gene_flank | <biomart.Attribute name='coding_gene_flank', display_name='', description=''>
5utr | <biomart.Attribute name='5utr', display_name='', description=''>
3utr | <biomart.Attribute name='3utr', display_name='3 UTR only', description=''>
gene_exon | <biomart.Attribute name='gene_exon', display_name='', description=''>
cdna | <biomart.Attribute name='cdna', display_name='', description=''>
coding | <biomart.Attribute name='coding', display_name='', description=''>
peptide | <biomart.Attribute name='peptide', display_name='', description=''>
upstream_flank | <biomart.Attribute name='upstream_flank', display_name='', description=''>
downstream_flank | <biomart.Attribute name='downstream_flank', display_name='', description=''>
sequence_gene_stable_id | <biomart.Attribute name='sequence_gene_stable_id', display_name='', description=''>
sequence_gene_stable_id_version | <biomart.Attribute name='sequence_gene_stable_id_version', display_name='', description=''>
sequence_description | <biomart.Attribute name='sequence_description', display_name='', description=''>
sequence_external_gene_name | <biomart.Attribute name='sequence_external_gene_name', display_name='', description=''>
sequence_external_source_name | <biomart.Attribute name='sequence_external_source_name', display_name='', description=''>
sequence_str_chrom_name | <biomart.Attribute name='sequence_str_chrom_name', display_name='', description=''>
sequence_gene_chrom_start | <biomart.Attribute name='sequence_gene_chrom_start', display_name='', description=''>
sequence_gene_chrom_end | <biomart.Attribute name='sequence_gene_chrom_end', display_name='', description=''>
sequence_gene_biotype | <biomart.Attribute name='sequence_gene_biotype', display_name='', description=''>
sequence_gene_version | <biomart.Attribute name='sequence_gene_version', display_name='', description=''>
sequence_family | <biomart.Attribute name='sequence_family', display_name='', description=''>
sequence_upi | <biomart.Attribute name='sequence_upi', display_name='', description=''>
sequence_uniprot_swissprot_accession | <biomart.Attribute name='sequence_uniprot_swissprot_accession', display_name='', description=''>
sequence_uniprot_sptrembl | <biomart.Attribute name='sequence_uniprot_sptrembl', display_name='', description=''>
sequence_uniprot_sptrembl_predicted | <biomart.Attribute name='sequence_uniprot_sptrembl_predicted', display_name='', description=''>
sequence_transcript_stable_id | <biomart.Attribute name='sequence_transcript_stable_id', display_name='', description=''>
sequence_transcript_stable_id_version | <biomart.Attribute name='sequence_transcript_stable_id_version', display_name='', description=''>
sequence_translation_stable_id | <biomart.Attribute name='sequence_translation_stable_id', display_name='', description=''>
sequence_translation_stable_id_version | <biomart.Attribute name='sequence_translation_stable_id_version', display_name='', description=''>
sequence_biotype | <biomart.Attribute name='sequence_biotype', display_name='', description=''>
sequence_transcript_biotype | <biomart.Attribute name='sequence_transcript_biotype', display_name='', description=''>
sequence_transcript_version | <biomart.Attribute name='sequence_transcript_version', display_name='', description=''>
sequence_peptide_version | <biomart.Attribute name='sequence_peptide_version', display_name='', description=''>
sequence_transcript_chrom_strand | <biomart.Attribute name='sequence_transcript_chrom_strand', display_name='', description=''>
sequence_transcript_chrom_start | <biomart.Attribute name='sequence_transcript_chrom_start', display_name='', description=''>
sequence_transcript_chrom_end | <biomart.Attribute name='sequence_transcript_chrom_end', display_name='', description=''>
sequence_transcription_start_site | <biomart.Attribute name='sequence_transcription_start_site', display_name='', description=''>
sequence_peptide_length | <biomart.Attribute name='sequence_peptide_length', display_name='', description=''>
sequence_transcript_length | <biomart.Attribute name='sequence_transcript_length', display_name='', description=''>
sequence_cdna_length | <biomart.Attribute name='sequence_cdna_length', display_name='', description=''>
5_utr_start | <biomart.Attribute name='5_utr_start', display_name="5' UTR start", description=''>
5_utr_end | <biomart.Attribute name='5_utr_end', display_name="5' UTR end", description=''>
3_utr_start | <biomart.Attribute name='3_utr_start', display_name="3' UTR start", description=''>
3_utr_end | <biomart.Attribute name='3_utr_end', display_name="3' UTR end", description=''>
sequence_exon_stable_id | <biomart.Attribute name='sequence_exon_stable_id', display_name='', description=''>
sequence_type | <biomart.Attribute name='sequence_type', display_name='', description=''>
sequence_exon_chrom_start | <biomart.Attribute name='sequence_exon_chrom_start', display_name='', description=''>
sequence_exon_chrom_end | <biomart.Attribute name='sequence_exon_chrom_end', display_name='', description=''>
sequence_exon_chrom_strand | <biomart.Attribute name='sequence_exon_chrom_strand', display_name='', description=''>
sequence_rank | <biomart.Attribute name='sequence_rank', display_name='', description=''>
sequence_phase | <biomart.Attribute name='sequence_phase', display_name='start phase', description=''>
sequence_end_phase | <biomart.Attribute name='sequence_end_phase', display_name='end phase', description=''>
sequence_cdna_coding_start | <biomart.Attribute name='sequence_cdna_coding_start', display_name='cDNA coding start', description=''>
sequence_cdna_coding_end | <biomart.Attribute name='sequence_cdna_coding_end', display_name='cDNA coding end', description=''>
sequence_genomic_coding_start | <biomart.Attribute name='sequence_genomic_coding_start', display_name='Genomic coding start', description=''>
sequence_genomic_coding_end | <biomart.Attribute name='sequence_genomic_coding_end', display_name='Genomic coding end', description=''>
sequence_constitutive | <biomart.Attribute name='sequence_constitutive', display_name='', description=''>
cds_length | <biomart.Attribute name='cds_length', display_name='CDS Length', description=''>
cds_start | <biomart.Attribute name='cds_start', display_name='CDS start', description=''>
cds_end | <biomart.Attribute name='cds_end', display_name='CDS end', description=''>
 
 
 
 
 - filters
 
 options | name
 --- | --- 
 'biotype'|<biomart.Filter name='biotype', type='list'>,
 'chromosomal_region'|<biomart.Filter name='chromosomal_region', type='text'>,
 'chromosome_name'|<biomart.Filter name='chromosome_name', type='text'>,
 'end'|<biomart.Filter name='end', type='text'>,
 'gene_id'|<biomart.Filter name='gene_id', type='text'>,
 'germ_line_variation_source'|<biomart.Filter name='germ_line_variation_source', type='list'>,
 'go_clos'|<biomart.Filter name='go_clos', type=''>,
 'go_evidence_code'|<biomart.Filter name='go_evidence_code', type='list'>,
 'go_name'|<biomart.Filter name='go_name', type=''>,
 'homolog_filters'|<biomart.Filter name='homolog_filters', type='boolean_list'>,
 'id_list_limit_microarray_filters'|<biomart.Filter name='id_list_limit_microarray_filters', type='id_list'>,
 'id_list_limit_protein_domain_filters'|<biomart.Filter name='id_list_limit_protein_domain_filters', type='id_list'>,
 'id_list_limit_xrefs_filters'|<biomart.Filter name='id_list_limit_xrefs_filters', type='id_list'>,
 'id_list_microarray_filters'|<biomart.Filter name='id_list_microarray_filters', type='boolean_list'>,
 'id_list_protein_domain_and_feature_filters'|<biomart.Filter name='id_list_protein_domain_and_feature_filters', type='boolean_list'>,
 'id_list_xrefs_filters'|<biomart.Filter name='id_list_xrefs_filters', type='boolean_list'>,
 'link_ensembl_gene_id'|<biomart.Filter name='link_ensembl_gene_id', type='text'>,
 'link_ensembl_transcript_stable_id'|<biomart.Filter name='link_ensembl_transcript_stable_id', type='text'>,
 'link_go_closure'|<biomart.Filter name='link_go_closure', type='text'>,
 'link_so_mini_closure'|<biomart.Filter name='link_so_mini_closure', type='list'>,
 'mane_select'|<biomart.Filter name='mane_select', type='boolean'>,
 'phenotype_description'|<biomart.Filter name='phenotype_description', type='list'>,
 'phenotype_source'|<biomart.Filter name='phenotype_source', type='list'>,
 'so_consequence_name'|<biomart.Filter name='so_consequence_name', type='list'>,
 'somatic_variation_source'|<biomart.Filter name='somatic_variation_source', type='list'>,
 'source'|<biomart.Filter name='source', type='list'>,
 'start'|<biomart.Filter name='start', type='text'>,
 'strand'|<biomart.Filter name='strand', type='text'>,
 'transcript_appris'|<biomart.Filter name='transcript_appris', type='boolean'>,
 'transcript_biotype'|<biomart.Filter name='transcript_biotype', type='list'>,
 'transcript_count_greater_than'|<biomart.Filter name='transcript_count_greater_than', type='text'>,
 'transcript_count_less_than'|<biomart.Filter name='transcript_count_less_than', type='text'>,
 'transcript_gencode_basic'|<biomart.Filter name='transcript_gencode_basic', type='boolean'>,
 'transcript_id'|<biomart.Filter name='transcript_id', type='text'>,
 'transcript_source'|<biomart.Filter name='transcript_source', type='list'>,
 'transcript_tsl'|<biomart.Filter name='transcript_tsl', type='boolean'>,
 'with_validated_snp'|<biomart.Filter name='with_validated_snp', type='boolean'>
