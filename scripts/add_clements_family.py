import os
import shutil

import pandas as pd


EPP_PATH = "epp_data_June2026.csv"
CLEMENTS_SRC = r"C:\Users\aneta\Downloads\Clements_v2025-October-2025.csv"
CLEMENTS_DST = r"data_raw\Clements_v2025-October-2025.csv"


os.makedirs("data_raw", exist_ok=True)
shutil.copy2(CLEMENTS_SRC, CLEMENTS_DST)

epp = pd.read_csv(EPP_PATH)
clements = pd.read_csv(CLEMENTS_SRC, low_memory=False)

species = clements[clements["category"].eq("species")].copy()
species["Clements_name_key"] = species["scientific name"].astype(str).str.replace(
    " ", "_", regex=False
)

family_full = species["family"].astype(str)
species["family"] = family_full.str.extract(r"^([^()]+)", expand=False).str.strip()
species["family_common"] = family_full.str.extract(r"\((.*)\)", expand=False)

taxonomy = species[
    ["Clements_name_key", "order", "family", "family_common"]
].drop_duplicates("Clements_name_key")

epp = epp.loc[
    :,
    [col for col in epp.columns if col not in ["order", "family", "family_common"]],
]

insert_at = list(epp.columns).index("Clements_name") + 1
merged = epp.merge(
    taxonomy,
    left_on="Clements_name",
    right_on="Clements_name_key",
    how="left",
    validate="many_to_one",
).drop(columns=["Clements_name_key"])

cols = list(merged.columns)
for col in ["family_common", "family", "order"]:
    cols.remove(col)
cols[insert_at:insert_at] = ["order", "family", "family_common"]
merged = merged[cols]

merged.to_csv(EPP_PATH, index=False)

print(f"rows={len(merged)}")
print(f"cols={len(merged.columns)}")
print(f"family_nonmissing={int(merged['family'].notna().sum())}")
print(f"families={merged['family'].nunique()}")
print(f"orders={merged['order'].nunique()}")
