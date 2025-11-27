import gzip
from collections import Counter
import matplotlib.pyplot as plt
import os

vcf_path = "results/variants/NA12878.vcf.gz"
os.makedirs("results/plots", exist_ok=True)

positions = []
types = Counter()

with gzip.open(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom, pos, vid, ref, alt = fields[:5]

        # chr20 dışındakileri at (genelde sadece chr20 var ama emin olalım)
        if chrom != "chr20":
            continue

        positions.append(int(pos))

        # multi-allelic (A,G,T vs) ile uğraşma, sunum için gerek yok
        if "," in alt:
            continue

        if len(ref) == 1 and len(alt) == 1:
            types["SNP"] += 1
        else:
            types["INDEL"] += 1

print(f"Number of variants used: {len(positions)}")
print(f"Counts by type: {types}")

# 1) Pozisyon dağılımı histogramı
plt.figure()
plt.hist(positions, bins=50)
plt.xlabel("Position on chr20")
plt.ylabel("Number of variants")
plt.title("Variant distribution along chr20 (NA12878)")
plt.tight_layout()
plt.savefig("results/plots/chr20_variant_distribution.png", dpi=200)

# 2) SNP vs INDEL bar grafiği
plt.figure()
labels = list(types.keys())
values = [types[k] for k in labels]
plt.bar(labels, values)
plt.ylabel("Number of variants")
plt.title("Variant types on chr20 (NA12878)")
plt.tight_layout()
plt.savefig("results/plots/variant_types_snp_vs_indel.png", dpi=200)

