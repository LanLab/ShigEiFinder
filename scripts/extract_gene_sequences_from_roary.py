#!/usr/bin/env python
# Developer: Michael Payne
import argparse
from Bio import SeqIO
import csv



def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,usage='used to identify '
                                                                                                  'clade specific gene '
                                                                                                  'sets from roary '
                                                                                                  'output files')
    parser.add_argument("-p","--presence", help="presence absence csv from roary output")
    parser.add_argument("-s", "--specificgenes", help="file containing cluster and specific gene set ids")
    parser.add_argument("-c", "--clusterassigns", help="assignments of genomes to clusters")
    parser.add_argument("-g", "--genesfolder", help="folder containing fasta files (.fasta suffix) of each genomes genes")
    parser.add_argument("-o", "--outprefix", help="prefix for output files")
    args = parser.parse_args()

    return args

def main():
    args = parseargs()

    csvfile = open(args.presence, "r")

    reader = csv.reader(csvfile)

    ser_group = args.specificgenes

    srr_ser = args.clusterassigns

    genomefolder = args.genesfolder

    outpath = open(args.outprefix + "_specific_gene_info.txt","w")

    outpath.write("Group\tSerovar\tGenome\tGeneID\n")

    ingroups = open(ser_group,"r").read().splitlines()

    groups = {}
    for line in ingroups[1:]:
        col = line.split("\t")
        if col[0] not in groups:
            groups[col[0]] = [col[1]]
        else:
            groups[col[0]].append(col[1])


    sero = {}
    insero = open(srr_ser,"r").read().splitlines()

    for i in insero[1:]:
        col = i.split("\t")
        if col[1] not in sero:
            sero[col[1]] = [col[0]]
        else:
            sero[col[1]].append(col[0])

    print("done_1")
    c=0
    for i in reader:
        if c==0:
            header1 = i
            header2 = [x.split("_")[-1] for x in header1]
            break

    # header = p_a[0]
    indexes = {}
    for serovar in sero:
        print(serovar)
        for acc in sero[serovar]:
            mat = acc
            ind = header2.index(mat)
            if serovar not in indexes:
                indexes[serovar] = [ind]
            else:
                indexes[serovar].append(ind)

    outd = {}

    print("done_2")
    c=0
    for line in reader:
        if c != 0:
            for serovar in groups:
                outd[serovar] = []
                if line[0] in groups[serovar]:
                    print(serovar)
                    inds = indexes[serovar]
                    m = 0
                    for pos in inds:
                        if line[pos] != "" and m == 0:
                            grp = line[0]
                            genome = header1[pos]
                            gene = line[pos]
                            m=1
                            outd[serovar].append((grp,genome,gene))
                            genome_path = genomefolder+"/"+genome+"/"+genome+".fasta"
                            ingenome = SeqIO.parse(genome_path,"fasta")


                            outpath.write(grp+'\t'+serovar+"\t"+genome+"\t"+gene+"\n")

                            if "___" in gene:
                                short_gene = gene.split("___")[0]
                            else:
                                short_gene = gene.split("(")[0]


                            gene_out = args.outprefix + "_" + serovar+"-"+grp+".fasta"

                            for seq in ingenome:
                                if seq.id == short_gene:
                                    SeqIO.write(seq,gene_out,"fasta")
        else:
            c+=1

    outpath.close()






