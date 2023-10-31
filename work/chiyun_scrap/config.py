import collections

vcf_file_format = "/lustre/scratch124/gsu/legacy/pipelines/builds/pf_70_build/pf_70_internal_release/vcf/%s.pf7.vcf.gz" # Location of VCFs, allowing chrom to vary
samples_fn = '/lustre/scratch126/gsu/team112/pf7/ftp_20221021/Pf7_samples.txt'
ref_genome_fn = '/lustre/scratch124/gsu/legacy/pfalciparum/resources/Pfalciparum.genome.fasta'
gff_fn = '/lustre/scratch124/gsu/legacy/pfalciparum/resources/snpEff/data/Pfalciparum_GeneDB_Feb2020/Pfalciparum_replace_Pf3D7_MIT_v3_with_Pf_M76611.gff'

population_colours = collections.OrderedDict()

population_colours['AF-W']= "#e31a1c"
population_colours['AF-C'] = "#fd8d3c" 
population_colours['AF-NE'] = "#bb8129" 
population_colours['AF-E'] = "#fecc5c"

# locations with at least 25 samples (see 20210218_drm_frequency_over_time_pops_adminlevel1.ipynb) (borderline exceptions are marked on this notebook)
# I think we should be able to generate this with code, but for now this is fine
locations = [
    ('AF-W', 'Gabon', 'Wouleu-Ntem'),
    ('AF-W', 'Cameroon', 'Sud-Ouest'),
    ('AF-W', 'Nigeria', 'Lagos'),
    ('AF-W', 'Benin', 'Littoral'),
    ('AF-W', 'Benin', 'Atlantique'),
    ('AF-W', 'Ghana', 'Volta'),
    ('AF-W', 'Ghana', 'Greater Accra'),
    ('AF-W', 'Ghana', 'Upper East'),
    ('AF-W', 'Ghana', 'Central'),
    ('AF-W', 'Ghana', 'Ashanti'),
    ('AF-W', 'Ghana', 'Brong Ahafo'),
    ('AF-W', 'Burkina Faso', 'Haut-Bassins'),
    ('AF-W', 'Mali', 'Segou'),
    ('AF-W', 'Mali', 'Sikasso'),
    ('AF-W', 'Mali', 'Koulikoro'),
    ('AF-W', 'Mali', 'Bamako'),
    ('AF-W', 'Mali', 'Kayes'),
    
    
    ('AF-W', "Côte d'Ivoire", 'Abidjan'),
    ('AF-W', 'Mauritania', 'Hodh el Gharbi'),
    ('AF-W', 'Mauritania', 'Hodh ech Chargui'), 
    ('AF-W', 'Guinea', 'Nzerekore'),
    ('AF-W', 'Guinea', 'Faranah'),
    ('AF-W', 'Senegal', 'Sedhiou'),
    ('AF-W', 'Senegal', 'Dakar'),
    ('AF-W', 'Gambia', 'Upper River'),
    ('AF-W', 'Gambia', 'North Bank'),
    ('AF-W', 'Gambia', 'Western'),
    
    ('AF-C', 'DRC', 'Kinshasa'),

    ('AF-NE', 'Kenya', 'Kisumu'),
    ('AF-NE', 'Sudan', 'Khartoum'),

    ('AF-E', 'Kenya', 'Kilifi'),
    ('AF-E', 'Mozambique', 'Gaza'), 
    ('AF-E', 'Tanzania', 'Lindi'),
    ('AF-E', 'Tanzania', 'Tanga'),
    ('AF-E', 'Tanzania', 'Morogoro'),
    ('AF-E', 'Tanzania', 'Kagera'),
    ('AF-E', 'Tanzania', 'Kigoma'),
    ('AF-E', 'Malawi', 'Zomba'),   
    ('AF-E', 'Malawi', 'Chikwawa'),
    
#     ('AS-S-E', 'India', 'West Bengal'),
#     ('AS-S-E', 'India', 'Odisha'),

#     ('AS-S-FE', 'Bangladesh', 'Chittagong'),
#     ('AS-S-FE', 'India', 'Tripura'),

#     ('AS-SE-W', 'Thailand', 'Tak'),
#     ('AS-SE-W', 'Myanmar', 'Tanintharyi'),
#     ('AS-SE-W', 'Myanmar', 'Shan'),
#     ('AS-SE-W', 'Myanmar', 'Kayin'),
#     ('AS-SE-W', 'Myanmar', 'Kachin'),
#     ('AS-SE-W', 'Myanmar', 'Bago'),
#     ('AS-SE-W', 'Myanmar', 'Mandalay'),
#     ('AS-SE-W', 'Myanmar', 'Sagaing'),

#     ('AS-SE-E', 'Vietnam', 'Khanh Hoa'),
#     ('AS-SE-E', 'Vietnam', 'Ninh Thuan'),
#     ('AS-SE-E', 'Vietnam', 'Gia Lai'),
#     ('AS-SE-E', 'Vietnam', 'Dak Lak'),
#     ('AS-SE-E', 'Vietnam', 'Quang Nam'),
#     ('AS-SE-E', 'Vietnam', 'Dak Nong'),
#     ('AS-SE-E', 'Vietnam', 'Quang Tri'),
#     ('AS-SE-E', 'Vietnam', 'Binh Phuoc'),
#     ('AS-SE-E', 'Cambodia', 'Ratanakiri'),
#     ('AS-SE-E', 'Cambodia', 'Stueng Traeng'),
#     ('AS-SE-E', 'Cambodia', 'Preah Vihear'),
#     ('AS-SE-E', 'Cambodia', 'Pursat'),
#     ('AS-SE-E', 'Cambodia', 'Battambang'),
#     ('AS-SE-E', 'Cambodia', 'Pailin'),
    
#     ('AS-SE-E', 'Laos', 'Sekong'),
#     ('AS-SE-E', 'Laos', 'Attapeu'),
#     ('AS-SE-E', 'Laos', 'Salavan'),
#     ('AS-SE-E', 'Laos', 'Champasak'),
#     ('AS-SE-E', 'Laos', 'Savannakhet'),
#     ('AS-SE-E', 'Thailand', 'Sisakhet'),

#     ('OC-NG', 'Papua New Guinea', 'Milne Bay'),
#     ('OC-NG', 'Papua New Guinea', 'Madang'),
#     ('OC-NG', 'Papua New Guinea', 'East Sepik'),
#     ('OC-NG', 'Indonesia', 'Papua')
]       
