import msprime
import gzip
from subprocess import run

ts = msprime.sim_ancestry(
    samples=500,
    ploidy=2,
    sequence_length=100_000_000,
    random_seed=2,
    recombination_rate=1e-8,
    population_size=10_000,
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=2)
with open("test.vcf", "w") as f:
    ts.write_vcf(f)

run("bcftools view test.vcf -Ob -o test.bcf", shell=True, check=True)
