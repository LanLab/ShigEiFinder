#!/usr/bin/env python
# Developer: Thanh Nguyen
import argparse
import os
import sys
import subprocess
import json
import re
import collections
from pkg_resources import resource_filename
import uuid
import shutil
from shigeifinder import __version__

def file_type(f, m):
    if m == 'a' and os.path.splitext(f)[-1] == ".fasta":
      return True
    elif m == 'a' and os.path.splitext(f)[-1] == ".fa":
      return True
    elif m == 'a' and os.path.splitext(f)[-1] == ".fna":
      return True
    elif m == 'r':
        if f.endswith('.fastq.gz') or f.endswith('.fastq'):
            return True
    return False

def get_json_data():
    try:
        json_file = resource_filename("shigeifinder","genes.json")
    except:
        curr_dir = os.path.dirname(os.path.realpath(__file__))
        json_file = os.path.join(curr_dir, 'genes.json')
    with open(json_file) as f:
        data = json.load(f)
    return data

def get_db_data():
    try:
        db_file = resource_filename("shigeifinder","resources/genes.fasta")
    except:
        curr_dir = os.path.dirname(os.path.realpath(__file__))
        db_file = os.path.join(curr_dir, 'resources/genes.fasta')
    # sys.stderr.write("db file location: " + db_file)
    return db_file

def ipaH_detect(genes):
    if 'ipaH' in genes:
        return True
    return False

def SB13_detect(genes):
    if 'CSB13_gene' in genes:
        return 'SB13'
    elif 'CSB13A_gene' in genes:
        return 'SB13-atypical'
    return 'Unknown'

def plasmid_genes(genes):
    plasmid = ["acp", "icsA (virG)","icsA","icsB","ipaA","ipaB","ipaC","ipaD","ipaJ","ipgA","ipgB1","ipgC","ipgD","ipgE","ipgF","mxiA","mxiC","mxiD","mxiE","mxiG","mxiH","mxiI","mxiJ","mxiK","mxiL","mxiM","mxiN","spa13","spa15","spa24","spa29","spa32","spa33","spa40","spa47","spa9","virA","virB","virF"]
    inter = list(set(plasmid) & set(genes))
    return len(inter)

def make_tmp_folder(args):
    uid = str(uuid.uuid1())
    fulltmp = args.tmpdir + "/shigeifinder_"+uid
    if os.path.exists(fulltmp):
        pass
        # shutil.rmtree(reads2run)
        # os.mkdir(reads2run)

    else:
        os.mkdir(fulltmp)
    return fulltmp

def type_cluster3(genes, data):
    gtr = ["gtrI", "gtrIC", "gtrII", "gtrIV", "gtrV", "gtrX"]
    opt = ["optIII", "optII"]
    oac = ["oac", "oac1b", "oacB", "oacC", "oacD"]
    g_list = ["SF1-5_wzx"]
    type_list = {}

    if "SBP_wzx" in genes.keys() and "SBP_wzy" in genes.keys():
        type_list["SB22"] = ["SBP_wzx", "SBP_wzy"]
    elif all(item in genes.keys() for item in g_list):
        if not any(item in genes.keys() for item in gtr) and not any(item in genes.keys() for item in opt) and not any(item in genes.keys() for item in oac):
            type_list["SFY"] = g_list
        else:
            if not any(item in genes.keys() for item in gtr):
                if all(item in genes.keys() for item in data["C3"]["SF3b"]["v1"]):
                    type_list["SF3b"] = data["C3"]["SF3b"]["v1"]
                elif all(item in genes.keys() for item in data["C3"]["SFY"]["v1"]):
                    type_list["SFY"] = data["C3"]["SFY"]["v1"]
                elif all(item in genes.keys() for item in data["C3"]["SFY"]["v2"]):
                    type_list["SFY"] = data["C3"]["SFY"]["v2"]
                elif all(item in genes.keys() for item in data["C3"]["SFY"]["v3"]):
                    type_list["SFY"] = data["C3"]["SFY"]["v3"]
                elif all(item in genes.keys() for item in data["C3"]["SFYv"]["v1"]):
                    type_list["SFYv"] = data["C3"]["SFYv"]["v1"]
                elif all(item in genes.keys() for item in data["C3"]["SFYv"]["v2"]):
                    type_list["SFYv"] = data["C3"]["SFYv"]["v2"]
            else:
                for s in data["C3"]:
                    if s != "cluster-genes" and s != "SF3b" and s != "SFY" and s != "SFYv":
                        for t in data["C3"][s]:
                            if all(item in genes.keys() for item in data["C3"][s][t]):
                                type_list[s] = data["C3"][s][t]

    if len(type_list) > 1:
        if "SFX" in type_list.keys():
            del type_list["SFX"]
        if "SF1a" in type_list.keys() and "SF1b" in type_list.keys():
            del type_list["SF1a"]
        if "SF1a" in type_list.keys() and "SF1c (7a)" in type_list.keys():
            del type_list["SF1a"]
        if "SF1a" in type_list.keys() and "SF1d" in type_list.keys():
            del type_list["SF1a"]
        if "SF2a" in type_list.keys() and "SF2b" in type_list.keys():
            del type_list["SF2a"]
        if "SF3a" in type_list.keys() and "SF2b" in type_list.keys():
            del type_list["SF3a"]
        if "SF3a" in type_list.keys() and "SF5b" in type_list.keys():
            del type_list["SF3a"]
        if "SF3b" in type_list.keys() and "SF3a" in type_list.keys():
            del type_list["SF3b"]
        if "SF4a" in type_list.keys() and "SF4b" in type_list.keys():
            del type_list["SF4a"]
        elif "SF4a" in type_list.keys() and "SF4av" in type_list.keys():
            del type_list["SF4a"]
        elif "SF4b" in type_list.keys() and "SF4bv" in type_list.keys():
            del type_list["SF4b"]
        elif "SF4av" in type_list.keys() and "SF4bv" in type_list.keys():
            del type_list["SF4av"]
        if "SF5a" in type_list.keys() and "SF5b" in type_list.keys():
            del type_list["SF5a"]
        if "SFX" in type_list.keys() and "SF3a" in type_list.keys():
            del type_list["SFX"]
        if "SFX" in type_list.keys() and "SFXv (4c)" in type_list.keys():
            del type_list["SFX"]
        if "SFX" in type_list.keys() and "SF1d" in type_list.keys():
            del type_list["SFX"]
        if "SFX" in type_list.keys() and "SF2b" in type_list.keys():
            del type_list["SFX"]
        if "SFX" in type_list.keys() and "SF5b" in type_list.keys():
            del type_list["SFX"]
        if "SFY" in type_list.keys() and "SFYv" in type_list.keys():
            del type_list["SFY"]
        if "SFY" in type_list.keys() and "SF5a" in type_list.keys():
            del type_list["SFY"]
        if "SFY" in type_list.keys() and "SF3a" in type_list.keys():
            del type_list["SFY"]
        if "SFY" in type_list.keys() and "SF2b" in type_list.keys():
            del type_list["SFY"]
        if "SFY" in type_list.keys() and "SF2a" in type_list.keys():
            del type_list["SFY"]
        if "SFY" in type_list.keys() and "SF1a" in type_list.keys():
            del type_list["SFY"]

    elif len(type_list) == 0:
        type_list["SF Unknown"] = list(set(genes).intersection(g_list))
        type_list["SF Unknown"].extend(set(genes).intersection(gtr))
        type_list["SF Unknown"].extend(set(genes).intersection(opt))
        type_list["SF Unknown"].extend(set(genes).intersection(oac))

    return type_list  

def cluster3_notes(info):
    genes = ["SF1-5_wzx","oac", "oac1b", "oacB", "oacC", "oacD", "optIII", "optII", "gtrI", "gtrIC", "gtrII", "gtrIV", "gtrV", "gtrX"]
    if "SF Unknown" in info.keys():
        absent = list(set(genes)-set(info["SF Unknown"]))
        return "Genes found: " + ",".join(info["SF Unknown"]) + ";Not found: " + ",".join(absent) 
    else:
        found = []
        for types, seq in info.items():
            found.extend(seq)
    return "Genes found: " + ",".join(sorted(list(set(found))))

def sb610_present(genes):
    count10 = 0
    count6 = 0
    for g in genes:
        if "SB10_wzy" == g or "SB10_wzx" == g:
            count10 += 1
        elif "SB6_wzy" == g or "SB6_wzx" == g:
            count6 += 1
    if count10 > 0 and count6 > 0:
        return True
    return False

def determine_serotype(genes_, cluster):
    shigella_clusters = ["C1", "C2", "C3", "CSS", "CSB12", "CSB13", "CSB13-atypical", "CSD1", "CSD8", "CSD10"]
    # Check for the all cluster-specific genes found in genes
    if ',' in cluster:
        clusters = cluster.split(',')
        serotypes = []
        for c in clusters:
            stype, note = determine_serotype(genes_, c)
            serotypes.append(stype)
        return ','.join(serotypes), ""

    genes = o_filter(cluster, genes_)
    type_list = {}
    data = get_json_data()
    
    if cluster == 'C3':
        type_list = type_cluster3(genes, data)
    elif cluster == 'C1' and sb610_present(genes):
        return "SB10", "Found a SB10 SNP"
    else:
        for s in data[cluster]:
            if s != "cluster-genes":
                if all(item in genes.keys() for item in data[cluster][s]):
                    type_list[s] = data[cluster][s]
    if len(type_list) == 1:
        if list(type_list.keys())[0] == "SF Unknown":
            return list(type_list.keys())[0], cluster3_notes(type_list)
        return list(type_list.keys())[0], ""
    elif len(type_list) > 1:
        if "EIEC H7" in type_list.keys():
            del type_list["EIEC H7"]
        if "EIEC O164:H7" in type_list.keys() and "EIEC O164" in type_list.keys():
            del type_list["EIEC O164"]
        if "EIEC O152:H30" in type_list.keys() and "EIEC O152" in type_list.keys():
            del type_list["EIEC O152"]
        if "EIEC O28ac" in type_list.keys() and "EIEC O28ac:H7" in type_list.keys():
            del type_list["EIEC O28ac"]
        if "SB1" in type_list.keys() and "SB20" in type_list.keys():
            del type_list["SB1"]
        if cluster == "C3":
            return '/'.join(type_list.keys()), cluster3_notes(type_list)
        return '/'.join(type_list.keys()), ""
    elif len(type_list) == 0:
        # Check for antigen matches
        note = ""
        antigen = antigen_search(genes)
        if cluster in shigella_clusters:
            shigella_antigens = get_oantigens(genes, cluster, False)
            if len(shigella_antigens) > 0:
                return '/'.join(shigella_antigens), note
            note = "Atypical Antigen in Shigella Cluster."
        
        if "EIEC " + antigen in data[cluster].keys():
            return "EIEC " + antigen, note
        
        return antigen, note
    return "Untypeable", ""

def o_filter(cluster, genes):
    # Removing the antigens that don't belong 
    o = get_oantigens(genes, cluster, False)
    shigella_clusters = ["C1", "C2", "C3", "CSS", "CSB12", "CSB13", "CSB13-atypical", "CSD1", "CSD8", "CSD10"]
    # o={}
    new_remove = []
    for g in genes:
        if 'wzx' in g:
            if 'O' not in g and cluster not in shigella_clusters:
                new_remove.append(g)
    # Removing wzx genes that are least likely for oantigen
    if len(o) > 1:
        #using wfep:
        if "O124" in o and "O164" in o:
            if "O124_wfep" in genes.keys():
                if genes["O124_wfep"] != 0:
                    new_remove.append("O164_wzx")
                    if "O164_wzy" in genes.keys():
                        new_remove.append("O164_wzy")
            else:
                new_remove.append("O124_wzx")
                if "O124_wzy" in genes.keys():
                    new_remove.append("O124_wzy")
        if "SS" in o and cluster != "CSS":
            new_remove.append("SS_wzx")
            
    return delete_genes(new_remove, genes)

def delete_genes(remove, genes):
    for g in remove:
        del genes[g]
    return genes

def antigen_search(genes):
    o_antigens = ["O105_wzy","O105_wzx","O112ab_wzx","O112ab_wzy","O112ac_wzy","O112ac_wzx","O121_wzx","O121_wzy",\
        "O124_wzx","O124_wzy","O129_wzx","O129_wzy","O13_wzx","O13_wzy","O130_wzx","O130_wzy","O135_wzx","O135_wzy",\
        "O143_wzx","O143_wzy","O147_wzx","O147_wzy","O149_wzx","O149_wzy","O152_wzx","O152_wzy","O16_wzx","O16_wzy",\
        "O167_wzx","O167_wzy","O180_wzx","O180_wzy","O183_wzy","O183_wzx","O25_wzx","O25_wzy","O32_wzx","O32_wzy",\
        "O40_wzy","O40_wzx","O53_wzx","O53_wzy","O7_wzx","O7_wzy","O79_wzy","O79_wzx","O8_wzy","O8_wzx","O8_wzm","O8_wzt"]
    h = {}
    ha = '0'
    o = {}
    oa = '0'
    for g in genes:
        if re.search(r'(.*)\_(.*)', g) is not None:
            temp = re.search(r'(.*)\_(.*)', g).group(1)
            if "fliC" in g:
                if re.search(r'(.*)\_(.*)\_.*', g) is not None:
                    temp = re.search(r'(.*)\_fliC\_.*', g).group(1)
                h[temp] = genes[g]
            elif g in o_antigens:
                o[temp] = genes[g]
    
    if len(h) > 0: 
        ha = '/'.join(h.keys())
    if len(o) > 0:
        oa ='/'.join(o.keys())
    return oa + ":" + ha

def gene_rename(gene):
    if gene.startswith('C'):
        return gene
    elif re.search(r'SF\-(.*)\_(.*)', gene) is not None:
        temp = re.search(r'SF\-(.*)\_(.*)', gene)
        return temp.group(1)
    return gene

def get_oantigens(genes, cluster, search):
    o = {}
    o_list = []
    shigella_clusters = ["C1", "C2", "C3", "CSS", "CSB12", "CSB13", "CSB13-atypical", "CSD1", "CSD8", "CSD10"]
    if "," in cluster:
        clusters = cluster.split(',')
        for c in clusters:
            o_list.extend(get_oantigens(genes, c, search))
    elif "CSP" in cluster and cluster != "CSP53":
        data = get_json_data()
        serotype = data["sporadic"][cluster]["serotype"]
        if ":" in serotype:
            antigens = serotype.split(':')
            return [antigens[0]]
        elif ":" not in serotype and "O" in serotype:
            return [serotype]
    elif cluster == "Shigella/EIEC Unclustered":
        for g in genes:
            if '_wz' in g:
                name = re.search(r'(.*)\_(wz.*)', g).group(1)
                if name in o:
                    o[name].append(g)
                    continue
                o[name] = [ g ]
    else:
        cluster_antigens = oantigen_cluster_specific(cluster)
        for g in genes:
            if '_wz' in g:
                name = re.search(r'(.*)\_(wz.*)', g).group(1)
                if cluster in shigella_clusters and "O" in name and not search:
                    continue
                elif cluster not in shigella_clusters and "O" not in name:
                    continue
                elif name not in cluster_antigens and not search:
                    continue

                if name in o:
                    o[name].append(g)
                    continue
                elif cluster != "C3" and g == "SF1-5_wzx":
                    continue
                # elif cluster != "C3" and g == "SF1-5_wzy":
                #     continue
                elif cluster == "C3" and "SB" in g and "SBP" not in g:
                    continue
                o[name] = [ g ]

    for antigen in o:
        if len(o[antigen]) > 0:
            o_list.append(antigen)

    return o_list

def oantigen_cluster_specific(cluster):
    if cluster == 'Shigella/EIEC Unclustered':
        return []

    antigens = set()
    single_eclusters = ["C7", "C8", "C9", "C10"]
    single_sclusters = ["CSD1", "CSD8", "CSD10", "CSB12", "CSB13", "CSS"]
    # genes conbinations from json file
    data = get_json_data()

    if "," in cluster:
        for c in cluster.split(','):
            antigens.update(oantigen_cluster_specific(c))
    elif "CSP" in cluster or cluster == "Unknown" or cluster == "Shigella/EIEC Unclustered":
        return antigens
    elif cluster == "C3":
        return [ "SF1-5", "SBP" ]
    elif cluster in single_eclusters:
        for key in data[cluster]:
            if key != "cluster-genes":
                oantigen = re.search(r'EIEC (O.*):.*', key).group(1)
                antigens.add(oantigen)
    elif cluster in single_sclusters:
        oantigen = re.search(r'C(S.*)', cluster).group(1)
        return [oantigen]
    elif cluster == "CSB13-atypical":
        return ["SB13"]
    else:
        for genes in data[cluster]:
            if genes != "cluster-genes":
                for g in data[cluster][genes]:
                    if "_wzx" in g or "_wzy" in g:
                        oantigen = re.search(r'(.*)\_wz.*', g).group(1)
                        antigens.add(oantigen)
    return list(antigens)

def get_hantigens(genes, cluster, search):
    h_list = set()
    shigella_clusters = ["C1", "C2", "C3", "CSS", "CSB12", "CSB13", "CSB13-atypical", "CSD1", "CSD8", "CSD10"]
    if cluster in shigella_clusters and not search:
        return []
    elif "," in cluster:
        clusters = cluster.split(',')
        for c in clusters:
            h_list.update(get_hantigens(genes, c, search))
    elif "CSP" in cluster and cluster != "CSP53":
        data = get_json_data()
        serotype = data["sporadic"][cluster]["serotype"]
        if ":" in serotype:
            antigens = serotype.split(':')
            return [antigens[1]]
        elif ":" not in serotype and "H" in serotype:
            return [serotype]

    for g in genes:
        if '_fliC' in g:
            name = re.search(r'(.*)\_fliC.*', g).group(1)
            antigen_list = hantigens_cluster_specific(cluster)
            if name not in antigen_list and not search:
                continue
            h_list.add(name)
    return list(h_list)

def hantigens_cluster_specific(cluster):
    antigens = set()
    shigella_clusters = ["C1", "C2", "C3", "CSS", "CSB12", "CSB13", "CSB13-atypical", "CSD1", "CSD8", "CSD10"]
    
    data = get_json_data()
    
    if cluster in shigella_clusters or "CSP" in cluster or cluster == "Unknown" or cluster == "Shigella/EIEC Unclustered":
        return antigens

    if ',' in cluster:
        clusters = cluster.split(',')
        for c in clusters:
            antigens.update(hantigens_cluster_specific(c))
    else:
        for types in data[cluster]:
            if types != "cluster-genes" and re.search(r'EIEC O.*:(H.*)', types):
                hantigen = re.search(r'EIEC O.*:(H.*)', types).group(1)
                antigens.add(hantigen)
    return list(antigens)

def blastn_cleanup(blast):
    blast.remove('')
    genes_set = {}
    for line in blast:
        info = line.rstrip().split('\t')
        gene = gene_rename(info[0])
        start = int(info[3])
        end = int(info[4])
        perc_identity = float(info[5])
       
        if start > end:
            save = start
            start = end
            end = save

        if gene in genes_set.keys():
            genes_set[gene]['positions'].extend(list(range(start, end+1)))
            genes_set[gene]['match'] = len(set(genes_set[gene]['positions']))
            if "SB10_" in gene or "SB6_" in gene:
                # print("yessss", gene, perc_identity, start, end)
                if "wzx" in gene and 904 in list(range(start, end+1)):
                    genes_set[gene]['pident'] = perc_identity
                elif "wzy" in gene and 141 in list(range(start, end+1)):
                    genes_set[gene]['pident'] = perc_identity
            elif "O124_wfep" == gene and 429 in list(range(start, end+1)) and 430 in list(range(start, end+1)):
                genes_set[gene]['pident'] = perc_identity
            elif perc_identity > genes_set[gene]['pident']:
                genes_set[gene]['pident'] = perc_identity
            continue
        genes_set[gene] = {}
        genes_set[gene]['slength'] = float(info[1])
        genes_set[gene]['match'] = float(info[2])
        genes_set[gene]['positions'] = list(range(start, end+1))
        genes_set[gene]['pident'] = perc_identity

    return lcoverage_filter(genes_set)

def lcoverage_filter(genes):
    genes_set = {}
    for gene in genes:
        len_coverage = 100*genes[gene]['match']/genes[gene]['slength']
        if gene == 'ipaH' and len_coverage > 10:
            genes_set[gene] = len_coverage
        elif 'C1_gene_4' == gene and genes[gene]['match'] > 253:
            genes_set[gene] = len_coverage
        elif 'C1_gene_2' == gene and genes[gene]['match'] > 250.99:
            genes_set[gene] = len_coverage
        elif 'C2_gene_1' == gene and genes[gene]['match'] > 392:
            genes_set[gene] = len_coverage
        elif gene.startswith('C') and len_coverage > 49.99: # THIS IS FOR CLUSTER SPECFIC
            genes_set[gene] = len_coverage
        elif gene == 'O124_wfep':
            if 429 not in genes[gene]['positions'] and 430 not in genes[gene]['positions']:
                genes_set[gene] = 0
            elif genes[gene]['pident'] != 100:
                continue
            elif len_coverage > 10:
                genes_set[gene] = len_coverage
        elif not gene.startswith('C') and len_coverage >= 50:
            genes_set[gene] = len_coverage
        
    sb610_blast_filter(genes, genes_set)

    return h_filter(genes_set)

def sb610_blast_filter(genes, genes_set):
    if "SB6_wzx" in genes_set.keys() and "SB10_wzx" in genes_set.keys():
        if genes_set["SB6_wzx"] == genes_set["SB10_wzx"]:
            if 904 not in genes["SB10_wzx"]["positions"]:
                del genes_set["SB10_wzx"]
            elif genes["SB6_wzx"]["pident"] != 100 and genes["SB10_wzx"]["pident"] == 100:
                del genes_set["SB6_wzx"]
            elif genes["SB6_wzx"]["pident"] == 100 and genes["SB10_wzx"]["pident"] != 100:
                del genes_set["SB10_wzx"]
            elif genes["SB6_wzx"]["pident"] > genes["SB10_wzx"]["pident"]:
                del genes_set["SB10_wzx"]
            elif genes["SB6_wzx"]["pident"] < genes["SB10_wzx"]["pident"]:
                del genes_set["SB6_wzx"]
    
    if "SB6_wzy" in genes_set.keys() and "SB10_wzy" in genes_set.keys():
        if genes_set["SB6_wzy"] == genes_set["SB10_wzy"]:
            if 141 not in genes["SB10_wzy"]["positions"]:
                del genes_set["SB10_wzy"]
            elif genes["SB6_wzy"]["pident"] != 100 and genes["SB10_wzy"]["pident"] == 100:
                del genes_set["SB6_wzy"]
            elif genes["SB6_wzy"]["pident"] == 100 and genes["SB10_wzy"]["pident"] != 100:
                del genes_set["SB10_wzy"]
            elif genes["SB6_wzy"]["pident"] > genes["SB10_wzy"]["pident"]:
                del genes_set["SB10_wzy"]
            elif genes["SB6_wzy"]["pident"] < genes["SB10_wzy"]["pident"]:
                del genes_set["SB6_wzy"]
    return genes

def h_filter(genes):
    remove = []
    hant = {}
    for g, l in genes.items():
        if "_fliC" in g:
            hant[g] = l

    if len(hant) > 0:
        max_len = hant[min(hant, key=(lambda k: abs(hant[k]-100)))]
        for a, l in hant.items():
            if l != max_len:
                remove.append(a)
    delete_genes(remove, genes)
    return h_duplicate_remove(genes)

def h_duplicate_remove(genes):
    remove = []
    hant = {}
    for g, l in genes.items():
        if "_fliC" in g:
            name = re.search(r'(.*\_fliC).*', g).group(1)
            hant[name] = l
            if name in genes.keys():
                continue
            remove.append(g)

    genes.update(hant)
    return delete_genes(remove, genes)

def sb610_snps(mpileup):
    result = {'wzx':'T', 'wzy':'A'}
    wzx = collections.Counter()
    wzy = collections.Counter()
    for d in mpileup:
        info = d.split('\t')
        bases = info[4].lower()
        if "wzx" in info[0]:
            wzx = collections.Counter(bases)
        elif "wzy" in info[0]:
            wzy = collections.Counter(bases)

    common_wzx = [x[0] for x in wzx.most_common(1)]
    common_wzy = [x[0] for x in wzy.most_common(1)]
    if "g" in common_wzx:
       result['wzx'] = 'G'
    
    if "g" in common_wzy:
       result['wzy'] = 'G'
    return result

def wfep_indel(mpileup):
    result = {'429':'A', '430':'A'}
    p429 = collections.Counter()
    p430 = collections.Counter()
    for d in mpileup:
        info = d.split('\t')
        bases = info[4].lower()
        if "O124_wfep" in info[0] and "429" in info[1]:
            p429 = collections.Counter(bases)
        elif "O124_wfep" in info[0] and "430" in info[1]:
            p430 = collections.Counter(bases)

    common_429 = [x[0] for x in p429.most_common(1)]
    common_430 = [x[0] for x in p430.most_common(1)]
    if "*" in common_429:
       result['429'] = '*'
    
    if "*" in common_430:
       result['430'] = '*'
    
    return result

def get_oantigen_geneids(data):
    ogenels = []
    for c in data:
        if c != "sporadic":
            for genesetname,geneset in data[c].items():
                if genesetname != "cluster-genes" and isinstance(geneset,list):
                    for gene in geneset:
                        if gene not in ogenels:
                            ogenels.append(gene)
                elif genesetname != "cluster-genes":
                    for subgenesetname,subgeneset in geneset.items():
                        for gene in subgeneset:
                            if gene not in ogenels:
                                ogenels.append(gene)
    return ogenels

def mapping_mode(bam, mpileup,non_ipah_cut,args):
    genes_set = {}
    depth_cut = mapping_depth_cutoff(bam)
    sb610_snp = sb610_snps(mpileup)
    wfep = wfep_indel(mpileup)
    data = get_json_data()
    oantigenids = get_oantigen_geneids(data)
    for line in bam:
        info = line.split('\t')
        if len(info) > 1 and '#rname' not in line:
            gene = gene_rename(info[0])
            lenperc = float(info[5]) # Meant to the percentage length cov calculated by samtools coverage
            meandepth = float(info[6])
            geneperccov = 100*meandepth/depth_cut
            # If the mapping ratio is < 10%
            if gene == 'ipaH' and geneperccov < args.ipaH_depth:
                continue
            elif gene in oantigenids and geneperccov < args.o_depth:
                continue
            elif gene != 'ipaH' and gene not in oantigenids and geneperccov < args.depth: #'group' in gene and
                continue

            if gene == 'ipaH' and lenperc > 10:
                genes_set[gene] = lenperc
            elif 'C1_gene_4' == gene and float(info[4]) > 253:
                genes_set[gene] = lenperc
            elif 'C1_gene_2' == gene and float(info[4]) > 250.99:
                genes_set[gene] = lenperc
            elif 'C2_gene_1' == gene and float(info[4]) > 392:
                genes_set[gene] = lenperc
            elif gene.startswith('C') and lenperc > 49.99: # THIS IS FOR CLUSTER SPECFIC
                genes_set[gene] = lenperc
            elif not gene.startswith('C') and lenperc >= 50:
                genes_set[gene] = lenperc
    if "SB6_wzx" in genes_set.keys() and "SB10_wzx" in genes_set.keys():
        if sb610_snp["wzx"] == "G":
            del genes_set["SB6_wzx"]
        else:
            del genes_set["SB10_wzx"]
    
    if "SB6_wzy" in genes_set.keys() and "SB10_wzy" in genes_set.keys():
        if sb610_snp["wzy"] == "G":
            del genes_set["SB6_wzy"]
        else:
            del genes_set["SB10_wzy"]

    if "O124_wfep" in genes_set.keys():
        if wfep['429'] == "*" and wfep['430'] == "*":
            del genes_set["O124_wfep"]

    return h_filter(genes_set)

def map_depth_ratios(bam):
    genes_set = []
    depth_cut = mapping_depth_cutoff(bam)
    genes_set.append('average 7 HS genes: ' + str(depth_cut))  
    for line in bam:
        info = line.split('\t')
        if len(info) > 1 and '#rname' not in line:
            gene = gene_rename(info[0])
            meandepth = float(info[6])
            ratio = 100*meandepth/depth_cut
            genes_set.append(gene + '\t' + str(ratio))
    return genes_set

def mapping_depth_cutoff(bam):
    mlst = ["NC_000913.3:recA", "NC_000913.3:purA", "NC_000913.3:mdh", "NC_000913.3:icd", "NC_000913.3:gyrB", "NC_000913.3:fumC", "NC_000913.3:adk"]
    depth = 0

    for line in bam:
        info = line.split('\t')
        if len(info) > 1 and '#rname' not in line:
            gene = info[0]
            meandepth = float(info[6])
            if gene in mlst:
                depth += meandepth
    return depth/7

def determine_cluster(genes):
    data = get_json_data()

    cluster_list = {}
    # Check for the all cluster-specific genes found in genes
    for s in data:
        if s != "sporadic":
            if all(item in genes.keys() for item in data[s]["cluster-genes"]):
                cluster_list[s] = data[s]["cluster-genes"]
    # Need to check if there is a sporadic(SP) cluster gene present if yes then it is the SP
    sporadic = sporadic_clusters(genes.keys())
    if sporadic != "Unknown":
        return "sporadic:" + sporadic
    if len(cluster_list) == 1:
        return list(cluster_list.keys())[0]
    elif 'C1' in cluster_list.keys() and 'CSB12' in cluster_list.keys():
        return 'CSB12'
    elif 'C5' in cluster_list.keys() and 'C3' in cluster_list.keys():
        return 'C3'
    elif 'C5' in cluster_list.keys() and 'C8' in cluster_list.keys():
        return 'C8'
    elif 'C1' in cluster_list.keys() and 'CSB12' in cluster_list.keys():
        return 'CSB12'
    elif 'C1' in cluster_list.keys() and 'C2' in cluster_list.keys():
        return 'C2'
    elif 'C1' in cluster_list.keys() and 'CSD1' in cluster_list.keys():
        return 'CSD1'
    elif 'C2' in cluster_list.keys() and 'CSS' in cluster_list.keys():
        return 'C2'   
    elif 'CSB12' in cluster_list.keys() and 'CSD1' in cluster_list.keys():
        return 'CSB12'
    elif len(cluster_list) > 1:
        return ','.join(cluster_list.keys())
    return "Unknown Cluster"

def sporadic_clusters(genes):
    data = get_json_data()

    # Check for the all cluster-specific genes found in genes
    for s in data["sporadic"]:
        if all(item in genes for item in data["sporadic"][s]["cluster-genes"]):
            return s
    return "Unknown"

def sporadic_serotype(cluster):
    data = get_json_data()
    return data["sporadic"][cluster]["serotype"]

def string_result(res, output):
    result = res['sample'] + "\t" + res['ipaH'] + "\t" + str(res['plasmid']) + "\t" + res['cluster'] + "\t" + \
        res['serotype'] + "\t" + ','.join(res['oantigens']) + "\t" + ','.join(res['hantigens']) + "\t" + res['notes']

    if output:
        result += "\n"

    return result

def get_gene_type(gene):
    gene_type = ""

    mlst = ["NC_000913.3:recA", "NC_000913.3:purA", "NC_000913.3:mdh", "NC_000913.3:icd", "NC_000913.3:gyrB", "NC_000913.3:fumC", "NC_000913.3:adk"]
    plasmid = ["acp", "icsA (virG)","icsA","icsB","ipaA","ipaB","ipaC","ipaD","ipaJ","ipgA","ipgB1","ipgC","ipgD","ipgE","ipgF","mxiA","mxiC","mxiD","mxiE","mxiG","mxiH","mxiI","mxiJ","mxiK","mxiL","mxiM","mxiN","spa13","spa15","spa24","spa29","spa32","spa33","spa40","spa47","spa9","virA","virB","virF"]
    if gene in mlst:
        gene_type = "House Keeping"
    elif "_fliC" in gene or "_wz" in gene:
        gene_type = "O/H-antigen Specific"
    elif gene in plasmid:
        gene_type = "Virulence Plasmid"
    elif gene.startswith('C'):
        gene_type = "Cluster-Specific"

    return gene_type

def run_mapping(args,dir,r1,r2,threads):
    # Run mapping to gather genes
    genesdb = get_db_data()
    coverage_mapped = []

    if args.single_end:
        name = r1.replace(".fastq","")
        name = name.replace(".gz", "")
        name = os.path.basename(name)
        bam_file = args.fulltmp + "/" + name + '.bam'
        qry1 = "bwa mem -t {thread} {gdb} {r1}  | samtools sort -@ {thread} -O bam -o {bamfile} - " \
               "&& samtools index {bamfile}".format(gdb=genesdb,r1=r1,r2=r2,thread=threads,bamfile=bam_file)
    else:
        name = re.search(r'(.*)\_.*\.fastq.*', r1).group(1)
        name = os.path.basename(name)
        bam_file = args.fulltmp + "/" + name + '.bam'
        qry1 = "bwa mem -t {thread} {gdb} {r1} {r2}  | samtools sort -@ {thread} -O bam -o {bamfile} - " \
               "&& samtools index {bamfile}".format(gdb=genesdb,r1=r1,r2=r2,thread=threads,bamfile=bam_file)
    try:
        mapping = subprocess.check_output(qry1, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        print(exc.output.decode('ascii'))
    
    qry2 = 'samtools coverage ' + bam_file
    try:
        coverage_mapped = subprocess.check_output(qry2, shell=True, stderr=subprocess.STDOUT).decode("ascii").split('\n')
    except subprocess.CalledProcessError as exc:
        print(exc.output)
    
    positions = ["SB6_wzx:904-904", "SB6_wzy:141-141", "SB10_wzx:904-904", "SB10_wzy:141-141", "O124_wfep:429-430"]
    mpileups = []
    for pos in positions:
        qry3 = 'samtools mpileup -r ' + pos + ' ' + bam_file
        try:
            mpileup = subprocess.check_output(qry3, shell=True, stderr=subprocess.STDOUT).decode("ascii").split('\n')
            for m in mpileup:
                if "[" not in m:
                    mpileups.append(m)
            mpileups.remove('')
        except subprocess.CalledProcessError as exc:
            print(exc.output)

    results = []
    for line in coverage_mapped:
        info = line.split('\t')
        if len(info) > 1 and info[5] != '0':
            results.append(line)

    os.remove(bam_file)
    os.remove(bam_file + '.bai')
    return results, mpileups
    
def run_blast(dir, fileA):
    # Get all genes for ipaH & cluster genes
    genesdb = get_db_data()

    qry = 'blastn -db ' + genesdb + ' -outfmt "6 sseqid slen length sstart send pident" -perc_identity 80 -query ' + \
        fileA
    blast_hits = subprocess.check_output(qry, shell=True, stderr=subprocess.STDOUT)
    blast_hits = blast_hits.decode("ascii").split('\n')
    return blast_hits

def three_properties(genes):
    ipah = False
    virplas = False
    clustgenes = False
    if ipaH_detect(genes.keys()):
        ipah = True
    
    if plasmid_genes(genes.keys()) >= 26:
        virplas = True
    
    if determine_cluster(genes) != "Unknown Cluster":
        clustgenes = True

    return ipah,virplas,clustgenes

def run_typing(dir, files, mode,args):
    threads = str(args.t)
    hits = args.hits
    ratios = args.dratio
    output = args.output
    args.fulltmp = make_tmp_folder(args)
    result = {}
    result['notes'] = ""
    if mode == "r":
        if args.single_end:
            hit_results, depths = run_mapping(args,dir,files,"", threads)
            genes = mapping_mode(hit_results, depths, 10, args)
            name = os.path.basename(files)
            name = name.replace(".fastq", "")
            name = name.replace(".gz", "")
            result['sample'] = name
        else:
            hit_results, depths = run_mapping(args,dir,files[0], files[1], threads)
            genes = mapping_mode(hit_results, depths, 10, args)
            name = os.path.basename(files[1])
            result['sample'] = re.search(r'(.*)\_.*\.fastq[\.gz]?', name).group(1)
    else:
        hit_results = run_blast(dir,files)
        genes = blastn_cleanup(hit_results)
        # print(genes)
        result['sample'] = os.path.basename(files).split('.')[0]
    shutil.rmtree(args.fulltmp)
    # Check presecnce of ipaH, cluster & plasmid
    # if not ipaH_detect(genes.keys()) and determine_cluster(genes) == "Unknown Cluster" and plasmid_genes(genes.keys()) < 26:
    ipah, virplas, clustgenes = three_properties(genes)
    if not ipah and SB13_detect(genes.keys()) == 'Unknown':
        result['ipaH'] = '-'
        result['plasmid'] = plasmid_genes(genes.keys())
        result['cluster'] = 'Not Shigella/EIEC'
        result['serotype'] = 'Not Shigella/EIEC'
        result['oantigens'] = []
        result['hantigens'] = []
        result['notes'] = 'ipaH negative, cluster-specific genes negative and less than 26 plasmid genes'
    else:
        # Check for the ipaH gene, if negative check for SB13 genes
        if not ipah and SB13_detect(genes.keys()) != 'Unknown':
            result['ipaH'] = '-'
            result['serotype'] = SB13_detect(genes.keys())
            result['cluster'] = "C"+SB13_detect(genes.keys())
            result['plasmid'] = plasmid_genes(genes.keys())
        else:
            result['ipaH'] = '+'
            result['plasmid'] = plasmid_genes(genes.keys())
            # Identify Cluster
            cluster = determine_cluster(genes)
            if "sporadic" in cluster:
                result['cluster'] = cluster.split(':')[1]
                result['serotype'] = sporadic_serotype(cluster.split(':')[1])
            elif cluster == "Unknown Cluster":
                result['cluster'] = "Shigella/EIEC Unclustered"
                result['serotype'] = antigen_search(genes)
                if mode == "r":
                    geneslowcov = mapping_mode(hit_results, depths, 1, args)
                    clusterlowcov = determine_cluster(geneslowcov)
                    serotypelowcov = antigen_search(geneslowcov)
                    if clusterlowcov != "Unknown Cluster":
                        result['notes'] = "Possible contamination by or low specific gene coverage of Shigella/EIEC strain in cluster {}, serotype {}. Some or all cluster specific gene to housekeeping genes ratios < 10%.".format(clusterlowcov,serotypelowcov)
            else:
                result['cluster'] = cluster
                # Determine Serotype
                result['serotype'], result['notes'] = determine_serotype(genes, cluster)
        search = False
        if "Atypical" in result['notes']:
            search = True
        result['oantigens'] = get_oantigens(genes, result['cluster'], search)
        result['hantigens'] = get_hantigens(genes, result['cluster'], search)
    
    # Print result
    if output:
        outp = open(output, "a+")
        outp.write(string_result(result,output))
        outp.close()
    else:
        print(string_result(result,output))

    if output:
        outp = open(output, "a+")
        if hits:
            outp.write("--------------- GENE SET ---------------\n")
            outp.write("#gene\tlength_coverage/depth\tgene_type\n")
            for key, item in genes.items():
                outp.write(key + '\t' + str(item)+ '\t' + get_gene_type(key) + "\n")
            if mode == "a":
                outp.write("---------- BLAST GENE HITS ----------\n")
                outp.write("#seqid\tslen\tlength\tsstart\tsend\tpident\n")
            else:
                outp.write("---------- BWA/SAMTOOLS GENE HITS ----------\n")
            outp.write('\n'.join(hit_results)+"\n")
            outp.write("----------------------------------------\n")

        if ratios:
            outp.write(
                "--------------- RATIOS OF DEPTH SPECFIC GENES TO AVERAGE DEPTH OF 7 HOUSE KEEPING GENES---------------\n")
            for g in map_depth_ratios(hit_results):
                outp.write(str(g)+"\n")
            outp.write(
                "--------------------------------------------------------------------------------------------------------------\n")
        outp.close()
    else:
        if hits:
            print("--------------- GENE SET ---------------")
            print("#gene\tlength_coverage/depth\tgene_type")
            for key, item in genes.items():
                print(key + '\t' + str(item) + '\t' + get_gene_type(key))
            if mode == "a":
                print("---------- BLAST GENE HITS ----------")
                print("#seqid\tslen\tlength\tsstart\tsend\tpident")
            else:
                print("---------- BWA/SAMTOOLS GENE HITS ----------")
            print('\n'.join(hit_results))
            print("----------------------------------------")

        if ratios:
            print("--------------- RATIOS OF DEPTH SPECIFIC GENES TO AVERAGE DEPTH OF 7 HOUSE KEEPING GENES---------------")
            for g in map_depth_ratios(hit_results):
                print(g)
            print("--------------------------------------------------------------------------------------------------------------")


def check_deps(checkonly,args):
    depslist = ["bwa","samtools","blastn"]
    f = 0
    for dep in depslist:
        rc = subprocess.call(['which', dep],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if rc == 0:
            sys.stderr.write(f'{dep:<10}:{"installed":<10}\n')
        else:
            sys.stderr.write(f'{dep:<10}:{"missing in path, Please install ":<10}{dep}\n')
            f+=1
    db_path = resource_filename("shigeifinder","resources/genes.fasta")
    if os.path.exists(db_path):
        sys.stderr.write(f'{"resources/genes.fasta":<10}:{"file present":<10}\n')
    else:
        sys.stderr.write(f'{db_path:<10}:{"not at expected path, possible broken installation?":<10}\n')
        f+=1
    json_path = resource_filename("shigeifinder","genes.json")
    if os.path.exists(json_path):
        sys.stderr.write(f'{"genes.json":<10}:{"file present":<10}\n')
    else:
        sys.stderr.write(f'{json_path:<10}:{"not at expected path, possible broken installation?":<10}\n')
        f+=1

    if f > 0:
        sys.stderr.write("One or more dependencies are missing or install paths broken.\n")
        sys.exit(1)
    else:
        sys.stderr.write("All dependencies present.\n")
        if checkonly:
            sys.exit(0)
        else:
            return

def write_header(args):
    outpath = args.output
    if not args.noheader:
        if outpath:
            outp = open(outpath, "w")
            outp.write("#SAMPLE\tipaH\tVIRULENCE_PLASMID\tCLUSTER\tSEROTYPE\tO_ANTIGEN\tH_ANTIGEN\tNOTES\n")
            outp.close()
        else:
            print("#SAMPLE\tipaH\tVIRULENCE_PLASMID\tCLUSTER\tSEROTYPE\tO_ANTIGEN\tH_ANTIGEN\tNOTES")

def main():
    parser = argparse.ArgumentParser(
        usage='\nAssembly fasta input/s:\n ShigeiFinder.py -i <input_data1> <input_data2> ... OR\n ShigeiFinder.py -i <directory/*> \nPaired end raw read fastq(.gz) input/s:\n ShigeiFinder.py -r -i <Read1> <Read2> OR \n ShigeiFinder.py -r -i <directory/*> \nSingle end raw read fastq(.gz) input/s:\n ShigeiFinder.py -r --single_end -i <Reads> OR \n ShigeiFinder.py -r --single_end -i <directory/*>\n')
    parser.add_argument("-i", nargs="+", help="<string>: path/to/input_data")
    parser.add_argument("-r", action='store_true', help="Add flag if file is raw reads.")
    parser.add_argument("-t", type=int, default='4', help="number of threads. Default 4.")
    parser.add_argument("--single_end", action='store_true', help="Add flag if raw reads are single end rather than paired.")
    parser.add_argument("--hits", action='store_true', help="To show the blast/alignment hits")
    parser.add_argument("--dratio", action='store_true', help="To show the depth ratios of cluster-specific genes to House Keeping genes")
    parser.add_argument("--update_db", action='store_true', help="Add flag if you added new sequences to genes database.")
    parser.add_argument("--output",
                        help="output file to write to (if not used writes to stdout)")
    parser.add_argument("--check", action='store_true', help="Check that dependencies are installed and required databases are available.")
    parser.add_argument("--o_depth", type=float, default=1.0, help="When using reads as input the minimum depth percentage relative to genome average for positive O antigen gene call (default 1.0).")
    parser.add_argument("--ipaH_depth", type=float,
                        help="When using reads as input the minimum depth percentage relative to genome average "
                             "for positive ipaH gene call (default 1.0).", default=1.0)
    parser.add_argument("--depth", type=float,
                        help="When using reads as input the minimum read depth for non ipaH/Oantigen gene to be called (default 10.0).", default=10.0)
    parser.add_argument("--tmpdir", help="temporary folder to use for intermediate files",
                        default=".")
    parser.add_argument("--noheader", help="do not print output header",
                        action='store_true')
    parser.add_argument("-v","--version", help="Print version information.", action='version', version=f"shigeifinder {__version__}")
    args = parser.parse_args()

    # Directory current script is in
    dir = os.path.dirname(os.path.realpath(__file__))

    #run and get the intermediate files for blast and bwa (makes updating the genes db easier)
    if args.update_db:
        dbpath = get_db_data()
        subprocess.run("bwa index " + dbpath + " >/dev/null 2>&1", shell=True)
        subprocess.run('makeblastdb -in ' + dbpath + ' -parse_seqids -blastdb_version 5 -title "Shigella/EIEC DB" -dbtype nucl >/dev/null 2>&1', shell=True)
        sys.exit()

    if args.check:
        check_deps(True,args)

    if args.dratio and not args.r:
        parser.error("-dratio requires -r. Only applies for raw reads.")
    if not args.i:
        parser.error("-i is required")

    # if not args.check:
    #     check_deps(False,args)


    if len(sys.argv) == 1:
        os.system("python3 " + dir + "/ShigeiFinder.py -h")
    else:
        mode = 'a'
        if args.r:
            mode = 'r'
        for files in args.i:
            if "*" in args.i[0]:
                dir1 = files.replace("*", "")
                if not os.path.isdir(dir1):
                    sys.exit('Invalid Directory Input! Directory: ' + dir1)
                break
            else:
                if not os.path.isfile(files):
                    sys.exit('Invalid Input File(s)! File Not Found:' + files)
                if not file_type(files, mode):
                    sys.exit('Incorrect File Type! File:' + files)
        if mode == 'r':
            # Run Raw Reads version
            # Check that there is 2 Reads inputed
            if (len(args.i) < 2 or len(args.i) % 2 != 0) and not args.single_end:
                if "*" in args.i[0] and len(args.i) == 1:
                    samples = set()
                    for files in os.listdir(dir1):
                        if files.endswith(".fastq.gz"):
                            path = dir1 + files
                            samples.add(path)
                    # print(samples)
                    reads = sorted(samples)
                    if len(reads) % 2 != 0:
                        sys.exit('Missing Input File(s)!! if you are using single end reads use the --single_end flag')
                    i = 0
                    write_header(args)
                    while i < len(reads):
                        f = [reads[i], reads[i+1]]
                        run_typing(dir, f, mode, args)
                        i += 2
                    sys.exit()
                else:
                    sys.exit('Missing Input File(s)!! if you are using single end reads use the --single_end flag')
            elif len(args.i) > 2 and not args.single_end:
                files = sorted(args.i)
                # print(files)
                i = 0
                write_header(args)
                while i < len(args.i):
                    f = [files[i], files[i+1]]
                    run_typing(dir, f, mode, args)
                    i += 2
                    # print(i, f)
                sys.exit()
            elif args.single_end:
                write_header(args)
                files=list(args.i)
                i=0
                while i < len(args.i):
                    f = files[i]
                    run_typing(dir, f, mode, args)
                    i += 1
                sys.exit()
            write_header(args)
            run_typing(dir, args.i, mode, args)
        else:
            # Run assembled genome version
            write_header(args)
            if "*" in args.i[0]:
                list_files = os.listdir(dir1)
                for f in list_files:
                    path = dir1 + "/" + f
                    if file_type(path, mode):
                        run_typing(dir, path, mode, args)
            elif len(args.i) > 1:
                for f in args.i:
                    run_typing(dir, f, mode, args)
            else:
                run_typing(dir, args.i[0], mode, args)
            
if __name__ == '__main__':
    main()
