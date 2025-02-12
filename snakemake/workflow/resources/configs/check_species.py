import os
import yaml

fname='species.yaml'
with open(fname,'r') as file:
	species_config = yaml.safe_load(file)

config = {"placeholder": 0}

keys = list(species_config.keys())

for SPECIES in keys:
	print(SPECIES)

	config.update({
        'GENOME': species_config[SPECIES]['GENOME'],
        'GTF': species_config[SPECIES]['GTF'],
        'ANNO_TAB': species_config[SPECIES]['ANNO_TAB'],
        'GSEA_DB_PATH': species_config[SPECIES]['GSEA_DB_PATH']
    })

	if not os.path.exists(config['GENOME']):
		print(f"{config['GENOME']} does not exist")
	if not os.path.exists(config['GTF']):
		print(f"{config['GTF']} does not exist")


print(config)
