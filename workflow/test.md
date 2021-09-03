from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "../config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
ids= [1,2]
genes= ["CFTR","CYP2C9","CYP4F2","CYP3A5","DPYD","NAT2","NUDT15","SLCO1B1","TPMT","G6PD"]

##### target rules #####
