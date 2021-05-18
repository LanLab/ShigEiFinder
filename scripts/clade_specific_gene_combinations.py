#!/usr/bin/env python
# Developer: Michael Payne
import argparse
import pandas
from itertools import combinations
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import sys
from os import remove

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,usage='used to identify '
                                                                                                  'clade specific gene '
                                                                                                  'sets from roary '
                                                                                                  'output files')
    parser.add_argument("-p","--presence", help="presence absence csv from roary output")
    parser.add_argument("-a", "--assign", help="clade assignments for isolates corresponding to -p")
    parser.add_argument("--falseneg", help="percentage false negatives allowed (rounds up)", default=0, type=int)
    parser.add_argument("-c", "--clade", help="clade for which the gene/s are being selected")
    parser.add_argument("--maxfp", help="maximum percentage of false positive to allow", default=10)
    parser.add_argument("--combsize", help="number of loci to use in combination to get specific sets", default=2,type=int)
    parser.add_argument("--reportno", help="number of loci sets to report", default=10, type=int)
    parser.add_argument("--minalign", help="minimum alignment length to report BLAST result", default=50, type=int)
    parser.add_argument("--maxcombfp", help="maximum total number of false positives for a set", default=20, type=int)
    parser.add_argument("--cpus", help="number of threads to use (for BLAST)", default=2, type=int)
    parser.add_argument("--mingenesize", help="minimum length for a gene to be used (bp)", default=200, type=int)
    parser.add_argument("-o", "--output", help="output file prefix")
    parser.add_argument("-d", "--blastdb", help="database of genomes to blast")
    parser.add_argument("-g", "--genomefolder", help="folder containing fasta files of genes")
    args = parser.parse_args()

    # make output file
    args.output = args.output + args.clade + "_fp" + str(args.maxfp) + ".txt"


    return args


def run_blast(query_seq, blastdb, pident,args):
    """

    :param query_seq: query sequence - can be multiple fasta seqs
    :param locus_db: blastdb path
    :param wordsize: blast word size
    :return: returns list of blast results
    """
    tmp_out = args.clade + "_fp" + str(args.maxfp) + "_tmp_blast.txt"
    cline = NcbiblastnCommandline(
        query=query_seq,
        db=blastdb,
        perc_identity=pident,
        out=tmp_out,
        outfmt=6,
        max_target_seqs=100000,
        num_threads=2,
        task="blastn")
    try:
        stdout, stderr = cline()
    except Exception as e:
        print(e)
        sys.exit()

    blast_records = open(tmp_out,"r").read().splitlines()

    remove(tmp_out)

    """
    blast_records structure: 
    list of results (if multifasta input, one result per fasta seq) 

    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'
    """

    return blast_records


def countgenes(x):
    if x == "":
        return 0
    else:
        res = len(list(x.split("\t")))
        return res

def not0(x):
    if x > 0:
        return 1
    else:
        return 0

def loadPresenceAbsense(args):
    inp = pandas.read_csv(args.presence,delimiter=",", dtype=str,index_col=0)
    inp.fillna("", inplace=True)
    inp.drop(inp.columns[[0,1,2,3,4,5,6,7,8,9,10,11,12]], axis=1, inplace=True)
    return inp

def get_clade_strainls(args):
    inf = open(args.assign,"r").read().splitlines()
    outls = []
    all = []
    genometoclade = {}
    for i in inf[1:]:
        col = i.split("\t")
        all.append(col[0])
        genometoclade[col[0]] = col[1]
        if col[1] == args.clade:
            outls.append(col[0])
    return outls,all,genometoclade

def groupls_keep_in_clade(args, dataframe,strainls):

    strainmin = int(len(strainls) * ((100 - int(args.falseneg)) / 100))

    cladedf = dataframe[strainls].T

    cladedfcounts = cladedf[(cladedf > 0) & (cladedf < 3)].count()

    cladedfsub = cladedfcounts[cladedfcounts >= strainmin]

    # get false negative strains for each gene

    falseneg = cladedfsub[cladedfsub < len(strainls)]

    falseneginf = cladedf[falseneg.to_dict().keys()]

    negdict = {}
    for col in falseneginf:
        subdf = falseneginf[col]
        lista = subdf[subdf == 0].index.to_list()
        negdict[col] = lista

    pass_groups = cladedfsub.index.tolist()
    print(pass_groups)

    return pass_groups, cladedf, negdict

def get_groups_under_fp_lim(args, dataframe, strainls):

    strainmax = len(strainls) * int((int(args.maxfp) / 100))

    noncladedf = dataframe[strainls].T

    noncladedfcounts = noncladedf[noncladedf > 0].count()

    noncladedfsub = noncladedfcounts[noncladedfcounts <= strainmax]

    pass_groups = noncladedfsub.index.tolist()

    return pass_groups, noncladedf

def get_candidate_genes(args,presabsdf,cladeStrains,allstrains):
    """
    steps:
    1 get presence or absense of genes in each strain from "presence"
    2 group strains into clusters as defined by "assign"
    3 get genes that are present in "clade" with "falseneg" allowable missing
    4 if any are perfect then report as single gene targets
    """
    presabsdf_orig = pandas.DataFrame(presabsdf)

    presabsdf = pandas.DataFrame(presabsdf.applymap(countgenes))
    pass_groups_fn, cladedf, negdict = groupls_keep_in_clade(args,presabsdf, cladeStrains)

    nonclade = list(set(allstrains).difference(set(cladeStrains)))

    pass_groups_fp, noncladedf = get_groups_under_fp_lim(args,presabsdf,nonclade)

    pass_groups = list(set(pass_groups_fn).intersection(set(pass_groups_fp)))

    clade_rep = cladeStrains[0]
    genomepath = "{0}/{1}/{1}.ffn".format(args.genomefolder, clade_rep)
    genomeseq = SeqIO.parse(genomepath,"fasta")

    grouptogene = {}
    genetogroup = {}
    passgenes = []

    for groupid in pass_groups:
        geneid = presabsdf_orig.loc[groupid,clade_rep]
        geneid = geneid.split("\t")[0][:-3]
        grouptogene[groupid] = geneid
        genetogroup[geneid] = groupid
        passgenes.append(grouptogene[groupid])

    pass_genes = [x.id for x in genomeseq if (x.id in passgenes and len(x.seq) > args.mingenesize)]

    pass_groups = [genetogroup[x] for x in pass_genes]

    return pass_groups, cladedf, noncladedf,genomeseq,grouptogene,genetogroup, negdict


def pairHasPrevResultSubset(pair,results):

    for i in results:
        if all(x in pair for x in i):
            return True
    return False



def try_combinations(args,candidates,noncladedf,cladeStrains,allstrains):
    """
    steps:
    1 get all unique pairwise combinations of genes
    2 get false positive set for both genes
    3 if pair has no common false positives report as useful combination
    4 if no successful pairs progress to triples, quads etc
    """

    result = 0
    results = []
    for limit in range(1,args.combsize+1):
        if len(results) < args.reportno:
            # print(limit)
            comb = combinations(candidates, limit)
            comb = sorted(comb)
            for pair in comb:
                if len(results) < args.reportno:
                    pairls = list(pair)
                    if not pairHasPrevResultSubset(pairls,results):
                        pairdf = noncladedf[pairls]

                        pairdf = pairdf.applymap(not0)
                        pairdf["sumRows"] = pairdf.sum(axis=1)
                        strains_with_all = pairdf[pairdf["sumRows"]==len(pairls)].index.tolist()
                        no_strains_with_all = len(strains_with_all)

                        if no_strains_with_all == 0:
                            print(limit,pairls,no_strains_with_all,strains_with_all)

                            results.append(pairls)
    return results

def getseqs_blast_returnfps(args,candidate_comb,noncladels,grouptogene,genetogroup,cladels):

    clade_rep = cladels[0]
    genomepath = "{0}/{1}/{1}.ffn".format(args.genomefolder, clade_rep)
    genomeseq = SeqIO.parse(genomepath, "fasta")

    genels = []
    for comb in candidate_comb:
        for group in comb:
            if grouptogene[group] not in genels:
                genels.append(grouptogene[group])

    print(genels)

    # too_short = [x for x in genomeseq if len(x.seq) >= args.mingenesize]
    toblast = [x for x in genomeseq if x.id in genels]

    genes2blast = args.clade + "_fp" + str(args.maxfp) + "_genes2blast_tmp.fasta"

    SeqIO.write(toblast,genes2blast,"fasta")

    blastres = run_blast(genes2blast,args.blastdb,args.cpus,80,args)


    remove(genes2blast)

    fp = {x:[] for x in genels}
    for i in blastres:
        col = i.split("\t")
        gene = col[0]

        hit_strain = col[1].split("_")[0]
        align_len = int(col[3])
        if hit_strain in noncladels and align_len > args.minalign:
            if hit_strain not in fp[gene]:
                fp[gene].append(hit_strain)

    return fp



def check_combs(false_pos,potential_combinations,grouptogene,args,genometoclade, negdict, cladestrains):
    outf = open(args.output, "w")

    # remove combinations with too small genes

    falseneg_genes = [x for x in negdict if len(negdict[x]) != 0]

    nofn_combs = []
    fn_combs = []

    for comb in potential_combinations:
        if any([x in falseneg_genes for x in comb]):
            fn_combs.append(comb)
        else:
            nofn_combs.append(comb)

    potential_combinations = nofn_combs + fn_combs

    for comb in potential_combinations:
        combfp = ""
        total_fn = []
        fn_details = []
        for group in comb:
            gene = grouptogene[group]
            fpstrains = false_pos[gene]
            if combfp == "":
                combfp = fpstrains
            else:
                newbfp = []
                for strain in combfp:
                    if strain in fpstrains:
                        newbfp.append(strain)
                combfp = list(newbfp)
            if group in negdict:
                for isolate in negdict[group]:
                    if isolate not in total_fn:
                        total_fn.append(isolate)
                if len(negdict[group]) > 0:
                    isolates = negdict[group]
                    details_string = "{}:{}".format(group,",".join(isolates))
                    fn_details.append(details_string)




        if len(combfp) <= args.maxcombfp:
            fpclades = [genometoclade[x] for x in combfp]
            combfpclade = [x[0]+"_"+x[1] for x in zip(combfp,fpclades)]
            ## fn output: \t total_isolates_with_fn\tgroupid:comma_sep_list;groupid:comma_sep_list.....
            outf.write(",".join(comb) + "\t" + str(len(combfp)) + "\t" + ",".join(combfpclade) + "\t" + str(len(cladestrains))+ "\t" + str(len(total_fn)) + "\t" + ";".join(fn_details)+ "\n")

    outf.close()





def main():

    args = parseargs()

    presabsdf = loadPresenceAbsense(args)  # loads presence absence into dataframe and transforms cells into gene counts

    cladeStrains, allstrains,genometoclade = get_clade_strainls(args)  # gets lists of all strains and those in the clade of interest

    # gets list of genes with < args.falseneg percentage false negatives - by default this is 0
    # also creates presence absence dataframes of clade and non clade isolates
    candidates, cladedf, noncladedf,genomeseq,grouptogene,genetogroup, negdict = get_candidate_genes(args,presabsdf,cladeStrains,allstrains)

    print(len(candidates),candidates)

    # for the list of candidates each combination of genes is tested to see whether all are found in the same false
    # positive strain. If a set of genes are never all found together then that set of genes is reported as a result
    # the size of the gene combinations starts at 1 for the whole list and increases progressively
    # at each size, successful sets of genes are reported until the total number of reported sets equals args.reportno
    # additionally if a successful set of 2 genes(for example) is found within a subsequent set of 3 genes that three
    # gene set will be excluded because the additional gene provides no benefit
    potential_combinations = try_combinations(args,candidates,noncladedf,cladeStrains,allstrains)

    noncladeStrains = list(set(allstrains).difference(set(cladeStrains)))

    false_pos = getseqs_blast_returnfps(args,potential_combinations,noncladeStrains,grouptogene,genetogroup,cladeStrains)

    check_combs(false_pos,potential_combinations,grouptogene,args, genometoclade, negdict,cladeStrains)

    """
    get sequences for each of the candidates (from clade of interest - args.clade)
    blast gene against all isolates
    get hits for isolates not in clade
    check combinations for no overlap in hits
    
    """




if __name__ == '__main__':
    main()