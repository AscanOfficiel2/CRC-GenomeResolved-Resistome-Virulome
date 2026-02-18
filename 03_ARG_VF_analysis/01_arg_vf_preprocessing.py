# -*- coding: utf-8 -*-

import pandas as pd
import re

# --------------------------------------------------
# 1. Load your TSV file
# --------------------------------------------------
df = pd.read_csv("arg_with_species.tsv", sep="\t")

# --------------------------------------------------
# 2. Function: Relabel species
# --------------------------------------------------
def relabel_species(tax):
    if pd.isna(tax):
        return "Unknown"

    genus_match = re.search(r"g__([^;]+)", tax)
    species_match = re.search(r"s__([^;]*)$", tax)

    if not genus_match:
        return "Unknown"

    genus = genus_match.group(1).strip()
    species = species_match.group(1).strip() if species_match else ""

    if species == "":
        return f"{genus} spp."
    else:
        return species

# --------------------------------------------------
# 3. Function: Extract gene name
# --------------------------------------------------
def extract_gene(gene_string):
    if pd.isna(gene_string):
        return "Unknown"

    # Look for gene=something
    match = re.search(r"gene=([^|]+)", gene_string)
    if match:
        return match.group(1).strip()
    return "Unknown"

# --------------------------------------------------
# 4. Apply functions
# --------------------------------------------------
df["Relabeled_Species"] = df["Species"].apply(relabel_species)
df["Gene_clean"] = df["Gene"].apply(extract_gene)

# --------------------------------------------------
# 5. Save processed file
# --------------------------------------------------
output = "arg_with_species_cleaned.tsv"
df.to_csv(output, sep="\t", index=False)

print("DONE! File saved as:", output)

import pandas as pd
import re

# ==========================================
# 1. LOAD FILE
# ==========================================
df = pd.read_csv("vfdb_with_species.tsv", sep="\t")


# ==========================================
# 2. SPECIES CLEANER
# ==========================================
def clean_species(taxon):
    if pd.isna(taxon):
        return "Unclassified spp."

    # extract genus
    g = re.search(r"g__([^;]+)", taxon)
    genus = g.group(1).strip() if g else "Unclassified"

    # extract species
    s = re.search(r"s__([^;]*)$", taxon)
    species = s.group(1).strip() if s else ""

    # if species missing -> genus spp.
    if species == "":
        return f"{genus} spp."
    else:
        return species


# ==========================================
# 3. GENE CLEANER
# ==========================================
def extract_gene(gene_str):
    if pd.isna(gene_str):
        return "Unknown"
    m = re.search(r"gene=([^|]+)", gene_str)
    return m.group(1).strip() if m else "Unknown"


# ==========================================
# 4. CATEGORY CLEANER
# ==========================================
def extract_category(gene_str):
    if pd.isna(gene_str):
        return "Unknown"

    m = re.search(r"category=([^|]+)", gene_str)
    if not m:
        return "Unknown"

    cat = m.group(1)

    # remove VFC codes: __VFC0272_ or _VFC0272_
    cat = re.sub(r"__?VFC\d+_?", "", cat)

    # underscores → spaces
    cat = cat.replace("_", " ")

    return cat.strip()


# ==========================================
# 5. APPLY CLEANERS
# ==========================================
df["Gene_clean"] = df["Gene"].apply(extract_gene)
df["Category_clean"] = df["Gene"].apply(extract_category)
df["Species_clean"] = df["Species"].apply(clean_species)


# ==========================================
# 6. SAVE CLEAN OUTPUT
# ==========================================
out = "vfdb_with_species_cleaned.tsv"
df.to_csv(out, sep="\t", index=False)

print("DONE! Saved:", out)
print(df.head(10).to_string(index=False))

"""Generate the ARG Mechanisms file

"""

# ==========================================================
#  LOAD DATA
# ==========================================================
import pandas as pd
import re

arg = pd.read_csv("arg_with_species.tsv", sep="\t")
meta = pd.read_csv("/content/metadata_MAGs.csv")
aro = pd.read_csv("aro_index.tsv", sep="\t")

# ==========================================================
#  CLEANING FUNCTIONS
# ==========================================================
def clean_taxonomy(tax):
    if pd.isna(tax):
        return None
    tax = re.sub(r'__+', '_', tax)
    tax = tax.replace(';_', ';').strip(';_ ')
    return tax

def extract_species_name(tax):
    if pd.isna(tax):
        return "Unclassified"
    species_match = re.search(r's_([\w\s\.\-\(\)]+)$', tax)
    genus_match   = re.search(r'g_([\w\s\.\-\(\)]+)', tax)
    family_match  = re.search(r'f_([\w\s\.\-\(\)]+)', tax)
    if species_match and species_match.group(1).strip():
        return species_match.group(1).strip()
    elif genus_match and genus_match.group(1).strip():
        return f"{genus_match.group(1).strip()} spp."
    elif family_match and family_match.group(1).strip():
        return family_match.group(1).strip()
    else:
        return "Unclassified"

# ==========================================================
# ARG CLEANING
# ==========================================================
arg['Gene_name'] = arg['Gene'].str.extract(r'gene=([^|]+)')
arg['Gene_name'] = arg['Gene_name'].str.strip()

arg['Genome'] = arg['Gene'].str.extract(r'sample=([\w\-\.\_]+)')
arg['Genome'] = arg['Genome'].str.strip()

arg['Species'] = arg['Species'].apply(clean_taxonomy)
arg['Species_name'] = arg['Species'].apply(extract_species_name)

arg.rename(columns={'Read Count': 'Read_Count'}, inplace=True, errors='ignore')

cols = ['Sample_ID', 'Genome', 'Gene_name', 'TPM', 'Read_Count', 'Species', 'Species_name']
arg = arg[cols].copy()

arg.dropna(subset=['Genome', 'Gene_name'], inplace=True)
arg.drop_duplicates(inplace=True)

# ==========================================================
# ARG MECHANISM ANNOTATION (with fuzzy matching)
# ==========================================================
aro.columns = aro.columns.str.strip().str.replace(" ", "_")

aro_subset = aro[['Model_Name', 'AMR_Gene_Family', 'Drug_Class',
                  'Resistance_Mechanism', 'CARD_Short_Name']].copy()

# ---- Exact match
arg = arg.merge(aro_subset, how="left", left_on="Gene_name", right_on="Model_Name")

# ---- Fuzzy match fallback
missing = arg[arg['AMR_Gene_Family'].isna()].copy()
mapped_rows = []

for gene in missing['Gene_name'].unique():
    matches = aro_subset[aro_subset['Model_Name'].str.contains(gene, case=False, regex=False)]
    if not matches.empty:
        m = matches.iloc[0]
        mapped_rows.append({
            'Gene_name': gene,
            'AMR_Gene_Family': m['AMR_Gene_Family'],
            'Drug_Class': m['Drug_Class'],
            'Resistance_Mechanism': m['Resistance_Mechanism'],
            'CARD_Short_Name': m['CARD_Short_Name']
        })

mapped_df = pd.DataFrame(mapped_rows)

if not mapped_df.empty:
    arg = arg.merge(mapped_df, on='Gene_name', how='left', suffixes=('', '_fuzzy'))
    for col in ['AMR_Gene_Family','Drug_Class','Resistance_Mechanism','CARD_Short_Name']:
        fcol = f"{col}_fuzzy"
        if fcol in arg.columns:
            arg[col] = arg[col].fillna(arg[fcol])
            arg.drop(columns=[fcol], inplace=True, errors='ignore')

# Rename columns
arg.rename(columns={
    "AMR_Gene_Family": "ARG_Gene_Family",
    "Drug_Class": "ARG_Drug_Class",
    "Resistance_Mechanism": "ARG_Mechanism",
    "CARD_Short_Name": "ARG_Short_Name"
}, inplace=True)

arg.fillna("Unknown", inplace=True)
arg['Data_source'] = 'ARG'

# ==========================================================
# MERGE WITH METADATA
# ==========================================================
meta.columns = meta.columns.str.strip()
meta = meta.rename(columns=lambda x: x.replace(" ", "_"))

arg_meta = arg.merge(meta, how='left', on='Sample_ID')

# ==========================================================
# SAVE RESULTS
# ==========================================================
arg.to_csv("ARG_clean.tsv", sep="\t", index=False)
arg_meta.to_csv("ARG_mechanism.tsv", sep="\t", index=False)
arg_meta.to_csv("Unified_ARG_Metadata_Mechanisms.tsv", sep="\t", index=False)

# ==========================================================
#  SUMMARY
# ==========================================================
print("\n ARG-ONLY Mechanism Integration Complete!")
print(f"ARG_clean rows: {len(arg)}")
print(f"ARG_mechanism rows: {len(arg_meta)}")
print(f"Unique ARG mechanisms: {arg_meta['ARG_Mechanism'].nunique()}")
print("\nSaved:")
print(" ARG_clean.tsv")
print(" ARG_mechanism.tsv")
print(" Unified_ARG_Metadata_Mechanisms.tsv")

"""Merge the mechanism to the Gene"""

import pandas as pd

# ============================================================
#  LOAD BOTH FILES
# ============================================================

df_mech = pd.read_csv("ARG_mechanism.tsv", sep="\t")
df_clean = pd.read_csv("arg_cleaned_final.csv")

print("Mechanism file:", df_mech.shape)
print("Clean ARG file:", df_clean.shape)

# Standardize column names
df_mech["Gene_name"] = df_mech["Gene_name"].astype(str).str.strip()
df_clean["Gene"]     = df_clean["Gene"].astype(str).str.strip()

# ============================================================
# SUBSET MECHANISM ANNOTATIONS
# ============================================================

annot = df_mech[[
    "Gene_name",
    "ARG_Gene_Family",
    "ARG_Drug_Class",
    "ARG_Mechanism"
]].drop_duplicates()

print("Annotation table:", annot.shape)

# ============================================================
# MERGE ON GENE NAME
# ============================================================

df_merged = df_clean.merge(
    annot,
    how="left",
    left_on="Gene",
    right_on="Gene_name"
)

df_merged.drop(columns=["Gene_name"], inplace=True)

print("Merged file:", df_merged.shape)

# ============================================================
# SAVE OUTPUT
# ============================================================

df_merged.to_csv("arg_cleaned_final_with_mechanisms.csv", index=False)

print("\n SUCCESS: Saved arg_cleaned_final_with_mechanisms.csv")

import pandas as pd

# ============================================================
#  LOAD BOTH FILES
# ============================================================

df_mech = pd.read_csv("ARG_mechanism.tsv", sep="\t")
df_clean = pd.read_csv("arg_cleaned_final.csv")

print("Mechanism file:", df_mech.shape)
print("Clean ARG file:", df_clean.shape)

# Ensure clean matching keys
df_mech["Gene_name"] = df_mech["Gene_name"].astype(str).str.strip()
df_clean["Gene"]     = df_clean["Gene"].astype(str).str.strip()

# ============================================================
# SELECT ONLY THE MECHANISM COLUMNS YOU NEED
# ============================================================

annot = df_mech[[
    "Gene_name",
    "ARG_Gene_Family",
    "ARG_Drug_Class",
    "ARG_Mechanism"
]].drop_duplicates()

print("Annotation table:", annot.shape)

# ============================================================
# EXACT MERGE (NO fuzzy match, NO transforms)
# ============================================================

df_merged = df_clean.merge(
    annot,
    how="left",
    left_on="Gene",
    right_on="Gene_name"
)

df_merged.drop(columns=["Gene_name"], inplace=True)

print("\nMerged file:", df_merged.shape)

# ============================================================
# REMOVE ROWS WHERE NO MECHANISM WAS FOUND
# ============================================================

df_merged = df_merged.dropna(
    subset=["ARG_Gene_Family", "ARG_Drug_Class", "ARG_Mechanism"],
    how="any"
)

print("After dropping NaN rows:", df_merged.shape)


# ============================================================
# SAVE OUTPUT
# ============================================================

df_merged.to_csv("arg_cleaned_final_with_mechanisms.csv", index=False)

print("\n SUCCESS: Saved arg_cleaned_final_with_mechanisms.csv")

###Clean the arg and vf file FROM ABOVE before continuing below:
#ARG should have columns: Sample_ID, Gene, RPKM, TPM, Read Count Species...
#VF should have columns: Sample_ID, Gene, RPKM, TPM, Read Count, Mechanism, Species

import pandas as pd

# ============================================================
# LOAD THE CSV FILE
# ============================================================

df = pd.read_csv("/content/vfdb_with_species_cleaned.tsv", sep="\t")   # <-- replace with your actual filename
print("Loaded file:", df.shape)

# Ensure TPM is numeric
df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce").fillna(0)

# Check required columns
required = {"Sample_ID", "Gene", "TPM"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

print("Columns OK.")


# ============================================================
# BUILD TPM MATRIX: Sample_ID × Gene
# ============================================================

TPM_matrix = df.pivot_table(
    index="Sample_ID",
    columns="Gene",
    values="TPM",
    aggfunc="sum",    # sum TPM if gene appears multiple times per sample
    fill_value=0
)

# Save matrix
TPM_matrix.to_csv("VF_TPM_Matrix.csv")

print("\nSaved VF_TPM_Matrix.csv →", TPM_matrix.shape)
print("Preview:")
print(TPM_matrix.head())

import pandas as pd

# ============================================================
# LOAD THE CSV FILE
# ============================================================

df = pd.read_csv("ARG_TPM_File.csv")   # <-- replace with your actual filename
print("Loaded file:", df.shape)

# Ensure TPM is numeric
df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce").fillna(0)

# Check required columns
required = {"Sample_ID", "Gene", "TPM"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

print("Columns OK.")


# ============================================================
# BUILD TPM MATRIX: Sample_ID × Gene
# ============================================================

TPM_matrix = df.pivot_table(
    index="Sample_ID",
    columns="Gene",
    values="TPM",
    aggfunc="sum",    # sum TPM if gene appears multiple times per sample
    fill_value=0
)

# Save matrix
TPM_matrix.to_csv("ARG_TPM_Matrix.csv")

print("\nSaved arg_TPM_Matrix.csv →", TPM_matrix.shape)
print("Preview:")
print(TPM_matrix.head())

##########################################################################

import pandas as pd

# ============================================================
# LOAD ARG + VF TPM MATRICES
# ============================================================

df_arg = pd.read_csv("ARG_TPM_Matrix.csv", index_col=0)
df_vf  = pd.read_csv("VF_TPM_Matrix.csv", index_col=0)

print("ARG matrix:", df_arg.shape)
print("VF matrix:", df_vf.shape)

# ============================================================
# RENAME COLUMNS WITH PREFIXES
# ============================================================

df_arg_prefixed = df_arg.copy()
df_arg_prefixed.columns = ["ARG_" + str(col) for col in df_arg.columns]

df_vf_prefixed = df_vf.copy()
df_vf_prefixed.columns = ["VF_" + str(col) for col in df_vf.columns]

# ============================================================
# ALIGN SAMPLE_ID INDEX BETWEEN MATRICES
# ============================================================

combined = df_arg_prefixed.join(df_vf_prefixed, how="outer")
combined = combined.fillna(0)

print("Combined matrix:", combined.shape)

# ============================================================
# SAVE OUTPUT
# ============================================================

combined.to_csv("Combined_ARG_VF_TPM_Matrix.csv")

print("\n SUCCESS: Saved Combined_ARG_VF_TPM_Matrix.csv")
