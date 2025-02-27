import pandas as pd

color_dict = {"fwd":"coral","rv":"mediumaquamarine","no_seq":"w"}

#species and genome names
species_dict={
    'hg38':["Human","Homo_sapiens"],
    'panTro5':["Chimp","Pan_troglodytes"],
    'panPan2':["Bonobo","Pan_paniscus"],
    'gorGor5':["Gorilla","Gorilla_gorilla_gorilla"],
    'ponAbe2':["Orangutan","Pongo_pygmaeus_abelii"],
    'nomLeu3':["Gibbon","Nomascus_leucogenys"],
    'rheMac8':["Rhesus","Macaca_mulatta"],
    'macFas5':["Crab-eating macaque","Macaca_fascicularis"],
    'macNem1':["Pig-tailed macaque","Macaca_nemestrina"],
    'cerAty1':["Sooty_mangabey", "Cercocebus_atys"],
    'papAnu3':["Baboon", "Papio_anubis"],
    'chlSab2':["Green_monkey","Chlorocebus_sabaeus"],
    'manLeu1':["Drill","Mandrillus_leucophaeus"],
    'nasLar1':["Proboscis_monkey","Nasalis_larvatus"], 
    'colAng1':["Angolan_colobus","Colobus_angolensis_palliatus"],
    'rhiRox1':["Golden_snub-nosed monkey","Rhinopithecus_roxellana"],
    'rhiBie1':["Black_snub-nosed monkey","Rhinopithecus_bieti"],
    "calJac3":["Marmoset","Callithrix_jacchus"],
    "saiBol1":["Squirrel_monkey","Saimiri_boliviensis"],
    "cebCap1":["White-faced_sapajou","Cebus_capucinus_imitator"],
    "aotNan1":["Ma's_night monkey","Aotus_nancymaae"],
    "tarSyr2":["Tarsier","Tarsius_syrichta"],
    "micMur3":["Mouse_lemur","Microcebus_murinus"],
    "proCoq1":["Coquerel's_sifaka","Propithecus_coquereli"],
    "eulMac1":["Black_lemur","Eulemur_macaco"],
    "eulFla1":["Sclater's_lemur","Eulemur_flavifrons"],
    "otoGar3":["Bushbaby","Otolemur_garnettii"],
    "mm10":["Mouse","Mus musculus"],
    "canFam3":["Dog","Canis_lupus_familiaris"],
    "dasNov3":["Armadillo","Dasypus novemcinctus"]
}

species_names_trivial=[v[0] for v in species_dict.values()]
species_names_scientific=[v[1] for v in species_dict.values()]

species_info_full = pd.DataFrame({"species":[],"species_trivial":[], "strand":[]})

for ss,st in zip(species_names_scientific,species_names_trivial):
    tmp_df = pd.DataFrame({"species":[ss,ss],"species_trivial":[st,st], "strand":["fwd","rv"]})
    species_info_full=pd.concat([species_info_full,tmp_df])


