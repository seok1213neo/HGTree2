import os
import sys
from multiprocessing import Pool, Queue
from queue import Queue
import threading
import pandas as pd
from Bio import SeqIO
import re
from ete3 import Tree
import glob

from datetime import datetime,timedelta
import time

#유저 인풋 변수로 받기
#archaea 인지 bacteria 인지 확인

USRinput = os.sys.argv[1]
processing_kingdom = os.sys.argv[2] # archaea or bacteria

id_value = os.sys.argv[3] #blast의 id 필터 값
e_value = os.sys.argv[4]
coverage_value = os.sys.argv[5] #blast의 coverage 필터 값

ranger_mode = os.sys.argv[6]

#example) python main.py assembly.fasta bacteria 80 0.00001 80 standard

USRQ_path= '/disk8/3.YS/01.HGT/07.USRQ' #유저 쿼리를 진행할 path
db_path='/disk8/3.YS/01.HGT/00.DATA' # 필요한 데이터들이 있는 path
program_path='/disk8/3.YS/00.programs' # 필요한 프로그램들이 있는 path

bacteria_list_txt = '{}/bacteria_list_species6.txt'.format(db_path) #사용된 bacteria genome의 정리파일
archaea_list_txt = '{}/archaea_list_species4.txt'.format(db_path) #사용된 archaea genome의 정리파일

prot_groups_with_hgt='{}/00.protDB/{}/{}/groups_with_hgt'.format(db_path,processing_kingdom,ranger_mode)
rrna_groups_with_hgt='{}/00.16srrna/{}/{}/groups_with_hgt'.format(db_path,processing_kingdom,ranger_mode)


hgt_db_path = '{0}/00.protDB/{1}/{2}/{2}_db.dmnd'.format(db_path,processing_kingdom,ranger_mode)
hgt_db_fasta_path = '{0}/00.protDB/{1}/{2}/{2}_db.fna'.format(db_path,processing_kingdom,ranger_mode)
rrna_path='{}/00.16srrna/{}'.format(db_path,processing_kingdom)
prokka_path='{}/00.prokka'.format(USRQ_path)
rgi_path='{}/09.rgi'.format(USRQ_path)
vf_path = '{}/10.vf'.format(USRQ_path)
pfam_path = '{}/11.pfam'.format(USRQ_path)
prodigal_path='{}/01.prodigal'.format(USRQ_path)
rnammer_path="{}/02.rnammer".format(USRQ_path)
reci_dia_path='{}/03.reci_diamond'.format(USRQ_path)
data_as_db="{}/03.reci_diamond/2.DATA_as_DB".format(USRQ_path)
query_as_db="{}/03.reci_diamond/1.USRQ_as_DB".format(USRQ_path)

new_prot_group_path='{}/04.groups_with_USRQ/2.protein'.format(USRQ_path)
new_rrna_group_path='{}/04.groups_with_USRQ/1.rrna'.format(USRQ_path)
USRQ_rrna_path='{}/02.rnammer'.format(USRQ_path)
prot_aligned_path='{}/05.aligned_groups/2.protein'.format(USRQ_path)
rrna_aligned_path='{}/05.aligned_groups/1.rrna'.format(USRQ_path)

#fasttree
fasttree_result='{}/06.fasttree'.format(USRQ_path)
species_tree_path='{}/1.rrna/2.outgrouped'.format(fasttree_result)
gene_tree_path='{}/2.protein'.format(fasttree_result)
merged_tree_path='{}/3.merged_nosupport'.format(fasttree_result)

#programs
rgi = '{}/CARD/rgi-master/rgi'.format(program_path)
vf_db = '{}/vfdb_prot.dmnd'.format(db_path)
pfam_files = '{}/PfamScan'.format(program_path)
pfam_scan = '{}/PfamScan/pfam_scan.pl'.format(program_path)
prokka_program_path = '{}/prokka-master/bin/prokka'.format(program_path)
ranger = '{}/5.Ranger/CorePrograms'.format(program_path)

#ranger
usrq_ranger_path = '{}/07.ranger'.format(USRQ_path)
opts_path='{}/07.ranger/1.opts'.format(USRQ_path)
opted_tree_path='{}/07.ranger/2.opt_trees'.format(USRQ_path)
ranger_result_path='{}/07.ranger/3.ranger_dtl'.format(USRQ_path)
ranger_result_path_T4='{}/07.ranger/3.ranger_dtl_t4'.format(USRQ_path)
aggregate_result_path='{}/07.ranger/4.aggregate_ranger'.format(USRQ_path)

bac_list={}
bacteria_species_list = {}

with open (bacteria_list_txt,'r')as baclist:
    for line in baclist:
        GCF=line.split("\t")[0]
        name=line.split("\t")[1].strip("\n")
        bac_list[GCF]=name

        species=line.split("\t")[2].strip()
        bacteria_species_list[GCF]=species

#내가 가진 아키아들의 GCF 아이디와 이름을 딕셔너리로

arc_list={}
archaea_species_list = {}

with open (archaea_list_txt,'r')as arclist:
    for line in arclist:
        GCF=line.split("\t")[0]
        name=line.split("\t")[1].strip("\n")
        arc_list[GCF]=name

        species=line.split("\t")[2].strip()
        archaea_species_list[GCF]=species

#내가 가진 아키아들의 GCF 아이디와 이름을 딕셔너리로

#유저의 파일 rename

def rename(usrq):
    start_time = time.monotonic()
    if (os.path.isfile('{}/USRQ.fasta'.format(USRQ_path))):
        os.system('rm {}/USRQ.fasta'.format(USRQ_path))

    contig_i=0

    if ( not os.path.isfile('{}/USRQ.fasta'.format(USRQ_path))):
         with open(usrq,'r') as USR_inpt:
            for line in USR_inpt:
                if line.startswith(">"):
                    contig_i+=1
                    with open ('{}/USRQ.fasta'.format(USRQ_path), "a")as USR_output:
                        USR_output.write(">USRQ_chr{}\n".format(contig_i))
                else:
                    with open ('{}/USRQ.fasta'.format(USRQ_path), "a")as USR_output:
                        USR_output.write(line)
    end_time = time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = end_time-start_time)))

def prokka_usrq():
    s=time.monotonic()
    os.chdir(USRQ_path)
    if not os.path.isdir(USRQ_path+'/00.prokka'):
        os.system('mkdir -p {}/00.prokka'.format(USRQ_path))
    
    os.system('prokka ./USRQ.fasta --outdir 00.prokka --prefix USRQ --force --quiet'.format(prokka_program_path))
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#위에서 rename된 USRQ.fasta를 prodigal을 사용하여 프로틴 코딩 부분을 확인.

def prodigal_usrq():
    s=time.monotonic()

    if not os.path.isdir('{}/01.prodigal'.format(USRQ_path)):
        os.system('mkdir -p {}/01.prodigal'.format(USRQ_path))

    os.system('prodigal -i {0}/USRQ.fasta -a {1}/USRQ_prodigal.faa -o {1}/USRQ_coordinate.txt -f sco -q'.format(USRQ_path,prodigal_path))
    os.system('prodigal -i {0}/USRQ.fasta -o {1}/USRQ.gff -f gff -q'.format(USRQ_path,prodigal_path))

    USR_protein=0
    USRQ_prodigal_dict={}

    with open ("{}/USRQ_prodigal.faa".format(prodigal_path), "r")as USRQ_prodigal:
        for line in USRQ_prodigal: # >USRQ_chr1_1
            if line.startswith(">"):
                USR_protein+=1
                USRQ_prot_id="_".join(line.split(" #")[0].split("_")[0:2])
                prot_no=line.split(" #")[0].split("_")[-1]

                USRQ_prot_id = USRQ_prot_id+"_"+prot_no

                USRQ_prodigal_dict.setdefault(USRQ_prot_id, "")

                with open ("{}/USRQ_prodigal2.faa".format(prodigal_path), "a")as USRQ_prodigal_new:
                    USRQ_prodigal_new.write(USRQ_prot_id+"\n")
                continue
            else:
                line = line.replace('*','')
                line = line.replace('J','X')
                line = line.replace('O','X')
                line = line.replace('U','X')

                with open ("{}/USRQ_prodigal2.faa".format(prodigal_path), "a")as USRQ_prodigal_new:
                    USRQ_prodigal_new.write(line)
                USRQ_prodigal_dict[USRQ_prot_id] += line.strip()
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#os.system ('{}/rnammer/rnammer -S bac -m ssu -f USRQ_rrna.fa < {}/USRQ.fasta'.format(USRQ_path))
#rnammer가 안돼서 barnnap으로 한 것.

def rrna_usrq():
    s=time.monotonic()
    os.chdir(USRQ_path)
    if not os.path.isdir("./02.rnammer"):
        os.system('mkdir 02.rnammer')

    os.system('barrnap -o {0}/USRQ_rrna.fa < {1}/USRQ.fasta --quiet'.format(rnammer_path,USRQ_path))

    with open ('{}/USRQ_rrna.fa'.format(rnammer_path),'r')as barnnap:
        lines = barnnap.read().split(">")[1].split("\n")[1]
        with open ('{}/USRQ_rrna2.fa'.format(rnammer_path),'a')as result:
            result.write(">USRQ\n"+lines)
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#alignment를 위한 threading

class Threading(threading.Thread):
        def __init__(self,queue):
                threading.Thread.__init__(self)
                self.queue = queue
        def run(self):
                while True:
                    systemstr = self.queue.get()
#                     print (systemstr)
                    os.system(systemstr)

                    self.queue.task_done()

#diamond를 통해서 reciprocal blast 해주기.
#처음은 우리 데이터를 db화해서

def DIAMOND_data_as_db(identity, e_value, coverage, ranger_mode): #80, 0.00001, 80
    s=time.monotonic()
    if not os.path.isdir(data_as_db):
        os.system('mkdir -p {}/03.reci_diamond/2.DATA_as_DB'.format(USRQ_path))
    
    data_as_db_result = '{}/2.DATA_as_DB/USRQ_{}_data_db.tsv'.format(reci_dia_path,ranger_mode)
    
    systemstr="{0}/diamond blastp --quiet --id {4} -e {5} --query-cover {6} --subject-cover {6} --threads 10 --db {1} --out {2} --outfmt tab --query {3}/01.prodigal/USRQ_prodigal2.faa".format(program_path, hgt_db_path, data_as_db_result, USRQ_path,identity, e_value, coverage)
    
    os.system(systemstr)
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))
    
#유저퀘리룰 디비화 해서
def DIAMOND_usrq_as_db(identity, e_value, coverage, ranger_mode): #80, 0.00001, 80
    s=time.monotonic()
    if not os.path.isdir(query_as_db):
        os.system('mkdir -p {}/03.reci_diamond/1.USRQ_as_DB'.format(USRQ_path))
    
    usrq_as_db_result = '{}/1.USRQ_as_DB/USRQ_{}_usrq_db.tsv'.format(reci_dia_path,ranger_mode)
   
    os.system('{}/diamond makedb --in {}/USRQ_prodigal2.faa --db {}/USRQ_db --quiet'.format(program_path, prodigal_path, reci_dia_path))

    systemstr="{0}/diamond blastp --quiet --id {5} -e {6} --query-cover {7} --subject-cover {7} --threads 10 --db {1}/USRQ_db.dmnd --out {4} --outfmt tab --query {3}".format(program_path, reci_dia_path, query_as_db, hgt_db_fasta_path, usrq_as_db_result, identity, e_value, coverage)
    
    os.system(systemstr)
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))
                                                
#결과 파일 둘을 읽어 USRQ 프로틴 마다 해당되는 그룹,프로틴을 모두 딕셔너리로 넣고, 그 중 id가 제일 높은 것을 확인
#id가 같다면 그 다음 순서였던 프로틴이 결과가 되도록.

def filter_diamond_result():
    
    usrq_db = open('{}/USRQ_{}_usrq_db.tsv'.format(query_as_db,ranger_mode), 'r')
    usrq_db = usrq_db.readlines()
    data_db = open('{}/USRQ_{}_data_db.tsv'.format(data_as_db,ranger_mode), 'r')
    data_db = data_db.readlines()

    best_hit_temp = {}

    for line in usrq_db:
        USRQ = line.split('\t')[1]
        prot = line.split('\t')[0]
        id = line.split('\t')[2]

        best_hit_temp.setdefault(USRQ, {})
        best_hit_temp[USRQ].setdefault(prot, id)

    for line in data_db:
        USRQ = line.split('\t')[0]
        prot = line.split('\t')[1]
        id = line.split('\t')[2]

        best_hit_temp.setdefault(USRQ, {})
        best_hit_temp[USRQ].setdefault(prot, id)

    best_hit = {}

    for usrq_protein in best_hit_temp.keys():
        for db_prot in best_hit_temp[usrq_protein]:
            if best_hit_temp[usrq_protein][db_prot] == max(best_hit_temp[usrq_protein].values()):

                best_hit.setdefault(usrq_protein, db_prot)
    
    return best_hit

#filter diamond result에서 나온 best hit을 토대로 새로운 그룹 만들어주기

def new_groups_with_USRQ():
    s=time.monotonic()
    if not os.path.isdir(new_prot_group_path):
        os.system('mkdir -p {}'.format(new_prot_group_path))

    if not os.path.isdir(new_rrna_group_path):
        os.system('mkdir -p {}'.format(new_rrna_group_path))

    best_hit = filter_diamond_result()

    USRQ_prots_seqs = {}

    USRQ_prodigal = open('{}/USRQ_prodigal2.faa'.format(prodigal_path), 'r')
    USRQ_prodigal = USRQ_prodigal.read().split('>')[1:]

    for prot in USRQ_prodigal:
        prot_number = prot.split('\n')[0].strip()
        prot_seq = ''.join(prot.split('\n')[1:]).strip()

        USRQ_prots_seqs[prot_number] = prot_seq

    USRQ_rrna_seq = open('{}/USRQ_rrna2.fa'.format(rnammer_path), 'r')
    USRQ_rrna_seq = USRQ_rrna_seq.read()

    for key in best_hit:
        group = best_hit[key].split('|')[0]

        with open ('{}/{}.fsa'.format(prot_groups_with_hgt,group), 'r') as before:
            before = ''.join(list(set(before.readlines())))
            with open ('{}/{}.fsa'.format(new_prot_group_path,group), 'a') as after:
                after.write(before +'>{}\n{}'.format(key,USRQ_prots_seqs[key]))

        with open ('{}/{}.fsa'.format(rrna_groups_with_hgt,group), 'r') as before2:
            before2 = ''.join(list(set(before2.readlines())))
            with open ('{}/{}.fsa'.format(new_rrna_group_path,group), 'a') as after2:
                after2.write(before2 +'{}'.format(USRQ_rrna_seq))
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#USRQ가 붙어서 새로 만들어진 그룹들을 얼라인 해줘야함

def align_new_groups():
    s=time.monotonic()
    prot_aligned_path='{}/05.aligned_groups/2.protein'.format(USRQ_path)
    rrna_aligned_path='{}/05.aligned_groups/1.rrna'.format(USRQ_path)

    if not os.path.isdir(prot_aligned_path) or not os.path.isdir(rrna_aligned_path):
        os.system('mkdir -p {}'.format(prot_aligned_path))
        os.system('mkdir -p {}'.format(rrna_aligned_path))

    flist=[]

    new_prots=os.listdir(new_prot_group_path)
    new_rrnas=os.listdir(new_rrna_group_path)

    for rrna in new_rrnas:
        group=rrna.split('.')[0]
        systemstr1 = ("clustalo --threads=1 -i {}/{} -o {}/{}.msa -v --outfmt=fa".format(new_rrna_group_path,rrna,rrna_aligned_path,group))
        systemstr2 = ("clustalo --threads=1 -i {}/{} -o {}/{}.msa -v --outfmt=fa".format(new_prot_group_path,rrna,prot_aligned_path,group))
        flist.append(systemstr1)
        flist.append(systemstr2)

    queue = Queue()
    for i in range(5):
        t = Threading(queue)
        t.setDaemon(True)
        t.start()
    for fname in flist:
        queue.put(fname)

    queue.join()
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#얼라인 된 애들은 다시 fasttree로 newick form을 만들어 준다.

def fasttree_groups():
    s=time.monotonic()
    if not os.path.isdir(fasttree_result):
        os.system('mkdir -p {}/1.rrna/1.unrooted'.format(fasttree_result))
        os.system('mkdir -p {}/1.rrna/2.outgrouped'.format(fasttree_result))
        os.system('mkdir -p {}/2.protein'.format(fasttree_result))
        os.system('mkdir -p {}/3.merged_nosupport'.format(fasttree_result))

    rrna_aligns=glob.glob('{}/*.msa'.format(rrna_aligned_path))
    
    for rrna in rrna_aligns:
        group=rrna.split('/')[-1].split(".")[0]
        prot = rrna.replace(rrna_aligned_path, prot_aligned_path)
        os.system('FastTree -quiet -gtr -nt -nopr {} > {}/1.rrna/1.unrooted/{}.txt'.format(rrna,fasttree_result,group))
        os.system('nw_reroot {0}/1.rrna/1.unrooted/{1}.txt SaCer > {0}/1.rrna/2.outgrouped/{1}.txt'.format(fasttree_result,group))
        os.system('FastTree -quiet -nopr {} > {}/2.protein/{}.txt'.format(prot,fasttree_result,group))
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))


#위에 얼라인 된 애들을 하나의 파일로 다시 만들어 준다.
# 하나로 통합 하면서, support value는 없애준다.

def merge_trees():
    s=time.monotonic()
    species_trees=glob.glob('{}/*.txt'.format(species_tree_path))
    
    for species_tree in species_trees:
        group=species_tree.split('/')[-1].split('.')[0]
        with open (species_tree,'r')as spec_tree:
            spec_tree=spec_tree.read().split(":")
            for dd in spec_tree:
                if ")" in dd:
                    dd_new = dd.replace(dd.split(")")[1],"")
                    spec_tree[spec_tree.index(dd)]=dd_new
            spec_tree_nosupport = ":".join(spec_tree)

            with open ('{}/{}.txt'.format(gene_tree_path,group),'r')as gene_tree:
                gene_tree = gene_tree.read().split(":")
                for dd in gene_tree:
                    if ")" in dd:
                        dd_new = dd.replace(dd.split(")")[1],"")
                        gene_tree[gene_tree.index(dd)]=dd_new
                gene_tree_nosupport = ":".join(gene_tree)

                with open ('{}/{}.txt'.format(merged_tree_path,group),'a')as merg_tree:
                    merg_tree.write(spec_tree_nosupport+";\n"+gene_tree_nosupport+";")
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#하나로 합쳐진 트리 파일을 레인저로 돌린다.
# 1 - optroot

def optroot():
    s=time.monotonic()
    if not os.path.isdir('{}/07.ranger'.format(USRQ_path)):
        os.system('mkdir -p {}'.format(opts_path))
        os.system('mkdir -p {}'.format(opted_tree_path))
        
        os.system('mkdir -p {}'.format(ranger_result_path))
        os.system('mkdir -p {}'.format(aggregate_result_path))

    merged_trees=os.listdir(merged_tree_path)

    for tree in merged_trees:
        group=tree.split(".")[0]
        os.system('{0}/OptRoot.linux -q -i {1}/{2} -o {3}/{2}'.format(ranger,merged_tree_path,tree,opts_path))
        #해당 디렉토리에는 gene tree들의 가능한 옵티멀하게 rooted 된 애들이 숙숙 나온다.
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#옾티멀한 gene tree들을 다시 맞는 group의 species tree랑 붙여줘야 한다.

def merge_opts():
    s=time.monotonic()
    opts_options=os.listdir(opts_path)

    for opts in opts_options:
        with open ('{}/{}'.format(opts_path, opts), 'r')as opt:
            opt_gene_trees=opt.readlines()[3:-4]

            for gene_tree in opt_gene_trees:
                i = opt_gene_trees.index(gene_tree)+1
                with open ('{}/{}'.format(merged_tree_path,opts),'r') as species_tree:
                    speceis_tree=species_tree.readlines()[0]
                    with open ('{}/{}_{}'.format(opted_tree_path,opts,i),'a') as result:
                        result.write(speceis_tree+gene_tree)

                    #같은 스피시스트리 같은 진트리지만 topology가 다른 애들이 그룹은 같게, i 는 다르게 한 디렉토리에 다 들어감
                    # 07.ranger/2.opt_tree 에 넣었다고 가정
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

# 2_1 - Ranger-DTL            

def ranger_dtl_T3():
    s=time.monotonic()
    opt_trees=os.listdir(opted_tree_path)

    for tree in opt_trees:
        group=tree.split('.')[0]
        i=tree.split('_')[-1]
        os.system('{}/Ranger-DTL.linux -q -i {}/{} -o {}/USRQ.{}_{}'.format(ranger,opted_tree_path,tree,ranger_result_path,group,i))
        #위에서 optimal하다고 한 gene tree의 숫자만큼 ranger result가 나옴
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

# ranger T4 돌리면서 결과들 이름 바꿔주기

# for group in groups:
def ranger_dtl_T4():
    s=time.monotonic()
    groups = os.listdir(opted_tree_path)
    
    for group in groups:
    
        GROUP = group.split(".")[0]
        Group = GROUP+'.'
        
        num = int(group.split("_")[-1])

        group_trees = [i for i in groups if Group in i]
        
        numbers =[]

        for a in group_trees:
            number = int(a.split("_")[-1])
            numbers.append(number)

        numbers.sort()
        try:
            i = max(numbers)

            NUMBER = num+i

            os.system ('{0}/Ranger-DTL.linux -q -T 4 -i {1}/{4} -o {3}/USRQ.{2}_{5}'.format(ranger,opted_tree_path,GROUP,ranger_result_path,group,NUMBER))

        except:
            print (group)
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

# 3 - Aggregate_ranger 그 여러가지 리설트들을 aggregate_ranger로 하나로 만들어줌

#위에서 USRQ.group_xx_i로 명명하여 USRQ.group_xx_로 묶어주었으니 묶어진 애들 

def aggregate_ranger():
    s=time.monotonic()
    merged_trees=os.listdir(merged_tree_path)

    for tree in merged_trees:
        group=tree.split(".")[0]
        os.system('{0}/AggregateRanger_recipient.linux {1}/USRQ.{2}_ > {3}/{2}.txt'.format(ranger,ranger_result_path,group,aggregate_result_path))
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

def parse_hgt_bac_1():
    s=time.monotonic()
    aggregate_results = os.listdir(aggregate_result_path)
    for file in aggregate_results:
        group = file.split('.')[0]
        file = open ('{}/{}'.format(aggregate_result_path,file), 'r')
        file_ = file.readlines()
        try:
            if float(file_[-3].split(' ')[-1]) >= 90 and float(file_[-4].split(' ')[-1]) >= 90:
                result = open ('{}/hgt_result_USRQ.txt'.format(usrq_ranger_path), 'a')
                hgt_event = []
                gene = open('{}/{}.txt_1'.format(opted_tree_path,group), 'r')
                gene_tree=gene.readlines()[1].split(",")
                for line in file_:
                    if 'recipient -->' in line:# and 'USRQ' in line:
                        donor = line.split("mapping --> ")[1].split(',')[0].strip()
                        recip = line.split("recipient --> ")[1].split(',')[0].strip()
                        if 0 < len(donor) <= 3  and 'USRQ' in recip: # n -> USRQ
                            tree = open ('{}/USRQ.{}_1'.format(ranger_result_path,group),'r')
                            tree_ = Tree(tree.readlines()[4].strip("\n"), format=1)
                            nodes = list(tree_.traverse())
                            x = list([i for i in nodes if '{}'.format(donor) in i][-1]) 
                            node_names=[]
                            for y in x:   # 해당 nxx 에서 노드별로 진트리와 비교해서 어떤 진들과 관련이 있는지 확인
                                for a in gene_tree:
                                    if y.name in a:
                                        a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                        GCFid = 'GCF_'+a.split("_")[0]
                                        protein = 'WP_'+a.split("_")[1]
                                        newname = bac_list[GCFid]
                                        a = GCFid+'|'+newname+'|'+protein
                                        node_names.append(a)
                            RECIP=''
                            for b in gene_tree:
                                    if recip in b:
                                        b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                        USRQ= b.split("_")[0]
                                        chr_protein = "_".join(b.split("_")[1:3])
                                        b = USRQ+'|'+chr_protein
                                        RECIP+=b
                            donor=",".join(node_names)
                            hgt_event.append('{}-->{}'.format(donor, RECIP))
                            tree.close()
                        elif 0 < len(recip) <= 3 and 'USRQ' in donor: # USRQ -> n
                            tree = open ('{}/USRQ.{}_1'.format(ranger_result_path,group),'r')
                            tree_ = Tree(tree.readlines()[4].strip("\n"), format=1)
                            nodes = list(tree_.traverse())
                            x = list([i for i in nodes if '{}'.format(recip) in i][-1])
                            node_names=[]
                            for y in x:
                                for a in gene_tree:
                                    if y.name in a:
                                        a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                        GCFid = 'GCF_'+a.split("_")[0]
                                        protein = 'WP_'+a.split("_")[1]
                                        newname = bac_list[GCFid]
                                        a = GCFid+'|'+newname+'|'+protein
                                        node_names.append(a)
                            DONOR = ''
                            for b in gene_tree:
                                    if donor in b:
                                        b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                        USRQ= b.split("_")[0]
                                        chr_protein = "_".join(b.split("_")[1:3])
                                        b = USRQ+'|'+chr_protein
                                        DONOR+=b
                            recip=",".join(node_names)
                            hgt_event.append('{}-->{}'.format(DONOR, recip))
                            tree.close()
                        elif 0 < len(donor)<=3 and 0 < len(recip)<=3: # n -> n
                            tree = open ('{}/USRQ.{}_1'.format(ranger_result_path,group),'r')
                            tree2 = Tree(tree.readlines()[4].strip("\n"), format=1)
                            nodes = list(tree2.traverse())
                            x = list([i for i in nodes if '{}'.format(donor) in i][-1])
                            x2 = list([i for i in nodes if '{}'.format(recip) in i][-1])
                            node_names=[]
                            for y in x:
                                for a in gene_tree:
                                    if y.name in a:
                                        if 'USRQ' in y.name:
                                            a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                            USRQ= a.split("_")[0]
                                            chr_protein = "_".join(a.split("_")[1:3])
                                            a = USRQ+'|'+chr_protein
                                            node_names.append(a)
                                        else:
                                            a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                            GCFid = 'GCF_'+a.split("_")[0]
                                            protein = 'WP_'+a.split("_")[1]
                                            newname = bac_list[GCFid]
                                            a = GCFid+'|'+newname+'|'+protein
                                            node_names.append(a)
                            donor=",".join(node_names)
                            node_names2=[]
                            for y2 in x2:
                                for b in gene_tree:
                                    if y2.name in b:
                                        if 'USRQ' in y2.name:
                                            b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                            USRQ= b.split("_")[0]
                                            chr_protein = "_".join(b.split("_")[1:3])
                                            b = USRQ+'|'+chr_protein
                                            node_names2.append(b)
                                        else:
                                            b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                            GCFid = 'GCF_'+b.split("_")[0]
                                            protein = 'WP_'+b.split("_")[1]
                                            newname = bac_list[GCFid]
                                            b = GCFid+'|'+newname+'|'+protein
                                            node_names2.append(b)
                            recip=",".join(node_names2)
                            hgt_event.append('{}-->{}'.format(donor, recip))
                            tree.close()
                        elif len(donor)==0 or len(recip)==0:
                            pass
#                 print (hgt_event)
                if len(hgt_event)>0:
                    result.write(group+"\t")
                    result.write("\t".join(hgt_event) + '\n')
#                     print ("\n".join(hgt_event))

                gene.close()
                result.close()
            else:
                file.close()
            file.close()
        except:
            print ('kkk')
            continue
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

def parse_hgt_arc_1():
    s=time.monotonic()
    aggregate_results = os.listdir(aggregate_result_path)
    for file in aggregate_results:
        group = file.split('.')[0]
        file = open ('{}/{}'.format(aggregate_result_path,file), 'r')
        file_ = file.readlines()
        try:
            if float(file_[-3].split(' ')[-1]) >= 90 and float(file_[-4].split(' ')[-1]) >= 90:
                result = open ('{}/hgt_result_USRQ.txt'.format(usrq_ranger_path), 'a')
                hgt_event = []
                gene = open('{}/{}.txt_1'.format(opted_tree_path,group), 'r')
                gene_tree=gene.readlines()[1].split(",")
                for line in file_:
                    if 'recipient -->' in line: #and 'USRQ' in line:
                        donor = line.split("mapping --> ")[1].split(',')[0].strip()
                        recip = line.split("recipient --> ")[1].split(',')[0].strip()
                        if 0 < len(donor) <= 3  and 'USRQ' in recip: # n -> USRQ
                            tree = open ('{}/USRQ.{}_1'.format(ranger_result_path,group),'r')
                            tree_ = Tree(tree.readlines()[4].strip("\n"), format=1)
                            nodes = list(tree_.traverse())
                            x = list([i for i in nodes if '{}'.format(donor) in i][-1])
                            node_names=[]
                            for y in x:   # 해당 nxx 에서 노드별로 진트리와 비교해서 어떤 진들과 관련이 있는지 확인
                                for a in gene_tree:
                                    if y.name in a:
                                        a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                        newname = a.split("_")[0]
                                        protein = 'WP_'+a.split("_")[1]
                                        GCFid = arc_list[newname]
                                        a = GCFid+'|'+newname+'|'+protein
                                        node_names.append(a)
#                             print(node_names)
                            RECIP=''
                            for b in gene_tree:
                                    if recip in b:
                                        b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                        USRQ= b.split("_")[0]
                                        chr_protein = "_".join(b.split("_")[1:3])
                                        b = USRQ+'|'+chr_protein
                                        RECIP+=b
                            donor=",".join(node_names)
                            hgt_event.append('{}-->{}'.format(donor, RECIP))
                            tree.close()
                        elif 0 < len(recip) <= 3 and 'USRQ' in donor: # USRQ -> n
                            tree = open ('{}/USRQ.{}_1'.format(ranger_result_path,group),'r')
                            tree_ = Tree(tree.readlines()[4].strip("\n"), format=1)
                            nodes = list(tree_.traverse())
                            x = list([i for i in nodes if '{}'.format(recip) in i][-1])
                            node_names=[]
                            for y in x:
                                for a in gene_tree:
                                    if y.name in a:
                                        a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                        newname = a.split("_")[0]
                                        protein = 'WP_'+a.split("_")[1]
                                        GCFid = arc_list[newname]
                                        a = GCFid+'|'+newname+'|'+protein
                                        node_names.append(a)
#                             print(node_names)
                            DONOR = ''
                            for b in gene_tree:
                                    if donor in b:
                                        b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                        USRQ= b.split("_")[0]
                                        chr_protein = "_".join(b.split("_")[1:3])
                                        b = USRQ+'|'+chr_protein
                                        DONOR+=b
                            recip=",".join(node_names)
                            hgt_event.append('{}-->{}'.format(DONOR, recip))
                            tree.close()
                        elif 0 < len(donor)<=3 and 0 < len(recip)<=3: # n -> n
                            tree = open ('{}/USRQ.{}_1'.format(ranger_result_path,group),'r')
                            tree2 = Tree(tree.readlines()[4].strip("\n"), format=1)
                            nodes = list(tree2.traverse())
                            x = list([i for i in nodes if '{}'.format(donor) in i][-1])
                            x2 = list([i for i in nodes if '{}'.format(recip) in i][-1])
                            node_names=[]
                            for y in x:
                                for a in gene_tree:
                                    if y.name in a:
                                        if 'USRQ' in y.name:
                                            a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                            USRQ= a.split("_")[0]
                                            chr_protein = "_".join(a.split("_")[1:3])
                                            a = USRQ+'|'+chr_protein
                                            node_names.append(a)
                                        else:
                                            a = a.replace('(', '').replace(')', '').replace(';', '').strip()
                                            newname = a.split("_")[0]
                                            protein = 'WP_'+a.split("_")[1]
                                            GCFid = arc_list[newname]
                                            a = GCFid+'|'+newname+'|'+protein
                                            node_names.append(a)
                            donor=",".join(node_names)
                            node_names2=[]
                            for y2 in x2:
                                for b in gene_tree:
                                    if y2.name in b:
                                        if 'USRQ' in y2.name:
                                            b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                            USRQ= b.split("_")[0]
                                            chr_protein = "_".join(b.split("_")[1:3])
                                            b = USRQ+'|'+chr_protein
                                            node_names2.append(b)
                                        else:
                                            b = b.replace('(', '').replace(')', '').replace(';', '').strip()
                                            newname = b.split("_")[0]
                                            protein = 'WP_'+b.split("_")[1]
                                            GCFid = arc_list[newname]
                                            b = GCFid+'|'+newname+'|'+protein
                                            node_names2.append(b)
                            recip=",".join(node_names2)
                            hgt_event.append('{}-->{}'.format(donor, recip))
                            tree.close()
                        elif len(donor)==0 or len(recip)==0:
                            pass
                if len(hgt_event)>0:
                    result.write(group+"\t")
                    result.write("\t".join(hgt_event) + '\n')
#                     print ("\n".join(hgt_event))
                gene.close()
                result.close()
            else:
                file.close()
            file.close()
        except:
            continue
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))


#gtdbtk를 통해 USRQ의 species를 predict, result file에서 같은 species는 제외해준다
def gtdbtk():
    s=time.monotonic()
    if not os.path.isdir('{}/08.gtdbtk'.format(USRQ_path)):
        os.system('mkdir -p {}/08.gtdbtk'.format(USRQ_path))

#     os.system('export GTDBTK_DATA_PATH=/disk8/3.YS/01.HGT/03.GTDBTK/GTDBTk-master/release202')
    os.system('cp {0}/USRQ.fasta {0}/08.gtdbtk/USRQ.fna'.format(USRQ_path))
    os.system('gtdbtk classify_wf --genome_dir {0}/08.gtdbtk --out_dir {0}/08.gtdbtk --cpus 4 --pplacer_cpus 4'.format(USRQ_path))
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

def USRQ_gtdbtk_bac():
    with open ('{}/08.gtdbtk/gtdbtk.bac120.summary.tsv'.format(USRQ_path), 'r')as gtdbtk_result:
        gtdbtk_result = gtdbtk_result.readlines()[1]
        species = gtdbtk_result.split("\t")[1].split(";")[-1].split("__")[-1]
        
        Species = species.split(" ")
        if '_' in Species[1]:
            species = Species[0]+'_'+Species[1].split("_")[0]
        elif '_' in Species[0]:
            species = Species[0].split("_")[0]+'_'+Species[1]
        elif '_' in Species[0] and '_' in Species[1]:
            species = Species[0].split("_")[0]+'_'+Species[1].split("_")[0]
        else:
            species = species.replace(" ", "_")
        
        bacteria_species_list["USRQ"]=species

        return species

def USRQ_gtdbtk_arc():
    with open ('{}/08.gtdbtk/gtdbtk.ar122.summary.tsv'.format(USRQ_path), 'r')as gtdbtk_result:
        gtdbtk_result = gtdbtk_result.readlines()[1]
        species = gtdbtk_result.split("\t")[1].split(";")[-1].split("__")[-1]

        Species = species.split(" ")
        if '_' in Species[1]:
            species = Species[0]+'_'+Species[1].split("_")[0]
        elif '_' in Species[0]:
            species = Species[0].split("_")[0]+'_'+Species[1]
        elif '_' in Species[0] and '_' in Species[1]:
            species = Species[0].split("_")[0]+'_'+Species[1].split("_")[0]
        else:
            species = species.replace(" ", "_")

        bacteria_species_list["USRQ"]=species

        return species

def parse_hgt_bac_2():
    s=time.monotonic()
    USRQ_result = open('{}/hgt_result_USRQ.txt'.format(usrq_ranger_path), 'r')
    
    USRQ_species = USRQ_gtdbtk_bac()
    
    for line in USRQ_result:
        new_line = []

        group = line.split("\t")[0]
        new_line.append(group)
        events = line.split("\t")[1:]

        for event in events :
            if 'USRQ' in event:
                donor = event.split("-->")[0]
                recip = event.split("-->")[1]

                if "," not in donor and "," not in recip: # 1 --> 1
                    donor_species = bacteria_species_list[donor.split("|")[0]]
                    recip_species = bacteria_species_list[recip.split("|")[0]]
                    if donor_species == recip_species:
                        continue
                    else:
                        new_line.append(event)

                elif ',' in donor and ',' not in recip: # n --> 1
                    new_donor = []
                    subjects = donor.split(",")
                    for subject in subjects:
                        donor_species = bacteria_species_list[subject.split("|")[0]]
                        if donor_species == bacteria_species_list['USRQ']: #USRQ_species
                            continue
                        else:
                            new_donor.append(subject.strip())
                    if len(new_donor)>0:
                        new_line.append('{0}-->{1}'.format(",".join(new_donor),recip.strip()))

                elif ',' in recip and ',' not in donor: # 1 --> n
                    new_recip = []
                    subjects = recip.split(",")
                    for subject in subjects:
                        recip_species = bacteria_species_list[subject.split("|")[0]]
                        if recip_species == bacteria_species_list['USRQ']: #USRQ_species
                            continue
                        else:
                            new_recip.append(subject.strip())
                    if len(new_recip)>0:
                        new_line.append('{1}-->{0}'.format(",".join(new_recip),donor.strip()))

                elif ',' in recip and ',' in donor: # n --> n
                    subjects1 = recip.split(",")
                    subjects2 = donor.split(",")
                    recip_spe=[]
                    donor_spe=[]

                    for subject1 in subjects1:
                        recip_spe.append(bacteria_species_list[subject1.split("|")[0]])
                    for subject2 in subjects2:
                        donor_spe.append(bacteria_species_list[subject2.split("|")[0]])

                    if set(recip_spe) == set(donor_spe):
                        continue 
                    else:
                        new_line.append(event.strip())

        if len(new_line)>1:
            with open ('{}/hgt_result_USRQ2.txt'.format(usrq_ranger_path), 'a')as result:
                result.write('\t'.join(new_line)+'\n')
    USRQ_result.close()
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

def parse_hgt_arc_2():

    USRQ_result = open('{}/hgt_result_USRQ.txt'.format(usrq_ranger_path), 'r')
    
    USRQ_species = USRQ_gtdbtk_arc()
    
    for line in USRQ_result:
        new_line = []

        group = line.split("\t")[0]
        new_line.append(group)
        events = line.split("\t")[1:]

        for event in events :
            if 'USRQ' in event:
                donor = event.split("-->")[0]
                recip = event.split("-->")[1]

                if "," not in donor and "," not in recip: # 1 --> 1
                    donor_species = archaea_species_list[donor.split("|")[0]]
                    recip_species = archaea_species_list[recip.split("|")[0]]
                    if donor_species == recip_species:
                        continue
                    else:
                        new_line.append(event)

                elif ',' in donor and ',' not in recip: # n --> 1
                    new_donor = []
                    subjects = donor.split(",")
                    for subject in subjects:
                        donor_species = archaea_species_list[subject.split("|")[0]]
                        if donor_species == archaea_species_list['USRQ']: #USRQ_species
                            continue
                        else:
                            new_donor.append(subject.strip())
                    if len(new_donor)>0:
                        new_line.append('{0}-->{1}'.format(",".join(new_donor),recip.strip()))

                elif ',' in recip and ',' not in donor: # 1 --> n
                    new_recip = []
                    subjects = recip.split(",")
                    for subject in subjects:
                        recip_species = archaea_species_list[subject.split("|")[0]]
                        if recip_species == archaea_species_list['USRQ']: #USRQ_species
                            continue
                        else:
                            new_recip.append(subject.strip())
                    if len(new_recip)>0:
                        new_line.append('{1}-->{0}'.format(",".join(new_recip),donor.strip()))

                elif ',' in recip and ',' in donor: # n --> n
                    subjects1 = recip.split(",")
                    subjects2 = donor.split(",")
                    recip_spe=[]
                    donor_spe=[]

                    for subject1 in subjects1:
                        recip_spe.append(archaea_species_list[subject1.split("|")[0]])
                    for subject2 in subjects2:
                        donor_spe.append(archaea_species_list[subject2.split("|")[0]])

                    if set(recip_spe) == set(donor_spe):
                        continue 
                    else:
                        new_line.append(event.strip())

        if len(new_line)>1:
            with open ('{}/hgt_result_USRQ2.txt'.format(usrq_ranger_path), 'a')as result:
                result.write('\t'.join(new_line)+'\n')
    USRQ_result.close()

def putative_hgt_genes():
    s=time.monotonic()
    put_hgts = []
    with open ('{}/hgt_result_USRQ2.txt'.format(usrq_ranger_path), 'r')as hgt_result:
        for line in hgt_result:
            donor = line.split('\t')[1].split("-->")[0].split(",")
            recip = line.split('\t')[1].split("-->")[1].split(",")
            
            for a in donor:
                if 'USRQ' in a:
                    put_hgts.append(a.replace("|", '_').strip())
            for b in recip:
                if 'USRQ' in b:
                    put_hgts.append(b.replace("|", '_').strip())
                    
    put_hgts = set(put_hgts)
    sequence_info = open('{}/USRQ_prodigal2.faa'.format(prodigal_path),'r')
    sequence_info = sequence_info.read().split(">")[1:]
    
    put_hgts_seq = []
    
    for put_hgt in put_hgts:
        for a in sequence_info:
            if put_hgt == a.split('\n')[0]:
                a = ''.join(a.split("\n")[1:])
                put_hgts_seq.append(">"+put_hgt+'\n'+a+'\n')
    
    put_hgts_seq = set(put_hgts_seq)
    
    with open ('{}/putative_hgt_genes.faa'.format(USRQ_path), 'a') as puthgt:
        puthgt.write(''.join(put_hgts_seq))
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#vfdb
def USRQ_vf():
    s=time.monotonic()
    if not os.path.isdir('{}/10.vf'.format(USRQ_path)):
        os.system('mkdir {}/10.vf'.format(USRQ_path))

    vf_db = '{}/vfdb_prot.dmnd'.format(db_path)
    os.system('{0}/diamond blastp --quiet --threads 2 --db {1} --out {2}/USRQ_VF.tsv --outfmt tab --query {3}/putative_hgt_genes.faa --id 80'.format(program_path,vf_db,vf_path,USRQ_path))

    USRQ_VF = {}

    with open ('{}/USRQ_VF.tsv'.format(vf_path), 'r') as blastresult:
        for line in blastresult:
            USRQ_id = line.split('\t')[0]
            USRQ_VFid = line.split('\t')[1].split('(')[0]
            USRQ_VF.setdefault(USRQ_id, [])
            USRQ_VF[USRQ_id].append(USRQ_VFid)

    for USRQid in USRQ_VF.keys():
        vfid = ",".join(set(USRQ_VF[USRQid]))
        with open ('{}/USRQ_VF.txt'.format(vf_path),'a') as vfresult:
            vfresult.write(USRQid+"\t"+vfid+"\n")
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#USRQ_vf 파일 파싱

def parse_vf(USRQ_id):
    
    USRQ_VF = {}
    with open ('{}/USRQ_VF.txt'.format(vf_path),'r') as vfresult:
        for line in vfresult:
            USRQ_prot_id = line.split("\t")[0]
            USRQ_prot_vf = line.split("\t")[1].strip()

            USRQ_VF[USRQ_prot_id] = USRQ_prot_vf
    
    return USRQ_VF[USRQ_id]
    
#CARD_rgi
def USRQ_rgi():
    s=time.monotonic()
    if not os.path.isdir("{}/09.rgi".format(USRQ_path)):
        os.system('mkdir {}/09.rgi'.format(USRQ_path))

    os.chdir('{}'.format(rgi_path))
    os.system('{} load --card_json {}/CARD/card.json --local'.format(rgi,program_path))
#     os.system('{} main -i {}/putative_hgt_genes.faa -o {}/USRQ_rgi -t protein -n 4 --local --clean'.format(rgi,USRQ_path,rgi_path))
    os.system('{} main -i {}/USRQ.ffn -o {}/USRQ_rgi -t contig -n 4 --local --clean'.format(rgi,prokka_path,rgi_path))

    with open ('{}/USRQ_rgi.txt'.format(rgi_path), 'r') as rgiresult:
        rgiresult = rgiresult.readlines()[1:]
        for line in rgiresult:
            USRQid = line.split('\t')[0]
            ARO = 'ARO:'+line.split('\t')[10]

            with open ('{}/USRQ_rgi_final.txt'.format(rgi_path),'a')as rgi_final:
                rgi_final.write(USRQid+'\t'+ARO+'\n')
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#USRQ_rgi 파일 파싱
def parse_rgi(USRQ_id):
    USRQ_RGI = {}
    with open ('{}/USRQ_rgi_final.txt'.format(rgi_path),'r') as rgiresult:
        for line in rgiresult:
            USRQ_prot_id = line.split("\t")[0]
            USRQ_prot_rgi = line.split("\t")[1].strip()

            USRQ_RGI[USRQ_prot_id] = USRQ_prot_rgi
    
    return USRQ_RGI[USRQ_id]

#pfam_scan
def USRQ_pfam():
    s=time.monotonic()
    if not os.path.isdir("{}/11.pfam".format(USRQ_path)):
        os.system('mkdir {}/11.pfam'.format(USRQ_path))

    os.system('{0} -fasta {1}/putative_hgt_genes.faa -dir {2} > {3}/USRQ_pfam.txt'.format(pfam_scan,USRQ_path,pfam_files,pfam_path))

    pfam_dict = {}

    with open ('{}/USRQ_pfam.txt'.format(pfam_path),'r') as pfam_result:
        pfam_result = pfam_result.readlines()[28:]
        for line in pfam_result:
            words = line.split(" ")
            for word in words:
                if 'USRQ_' in word:
                    prot_id = word
                elif 'PF' in word and '.' in word:
                    pf_id = word

            pfam_dict.setdefault(prot_id, [])
            pfam_dict[prot_id].append(pf_id)

    for USRQ_prot_id in pfam_dict.keys():
        pfam_id = ','.join(pfam_dict[USRQ_prot_id])
        with open ('{}/USRQ_pfam_final.txt'.format(pfam_path),'a') as pfam_final:
            pfam_final.write(USRQ_prot_id+'\t'+pfam_id+'\n')
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#usrq_pfam 파일 파싱

def parse_pfam(USRQ_id):
    USRQ_PFAM = {}
    with open ('{}/USRQ_pfam_final.txt'.format(pfam_path),'r') as vfresult:
        for line in vfresult:
            USRQ_prot_id = line.split("\t")[0]
            USRQ_prot_pfam = line.split("\t")[1].strip()

            USRQ_PFAM[USRQ_prot_id] = USRQ_prot_pfam
    
    return USRQ_PFAM[USRQ_id]

#statistics 만드는데 활용될 gbk parsing 함수

def annotating_gene(i):
    prokka_gbk = '{}/USRQ.gbk'.format(prokka_path)
    genome_record = SeqIO.parse(prokka_gbk, "genbank")

    for records in genome_record:
        for feature in records.features:
            if feature.type == 'CDS' and i in feature.location:
                
                dd = {}
                
                try:
                    dd['gene'] = feature.qualifiers['gene'][0]
                except:
                    dd['gene'] = 'NA'
                try:
                    dd['product'] = feature.qualifiers['product'][0]
                except:
                    dd['product'] = 'NA'
                return dd


#statistics file 만들기 파일 1개

def USRQ_HGT_statistics_bac():
    s=time.monotonic()
    prokka_gbk = '{}/00.prokka/USRQ.gbk'.format(prokka_path)
    hgt_resultfile = open('{}/hgt_result_USRQ2.txt'.format(usrq_ranger_path),'r')
    predicted_hgt = hgt_resultfile.read().count("-->")

    with open ('{}/01.prodigal/USRQ.gff'.format(USRQ_path), 'r') as USRQ_gff:
        lines = USRQ_gff.readlines()
        USRQ_prot_num = len(lines)

    recip_number = 0
    donor_number = 0
    
    result_donor = []
    result_recip = []
    
    with open ('{}/07.ranger/hgt_result_USRQ2.txt'.format(USRQ_path), 'r') as hgt_result:
        for line in hgt_result:
            group = line.split("\t")[0]
            events = line.split("\t")[1:]

            for event in events:
                recip = event.split("-->")[1]
                donor = event.split("-->")[0]

                if 'USRQ' in recip:
                    recip_number+=1
                    recip_prot = (a for a in recip.split(",") if 'USRQ' in a)
                    
                    donor = event.split("-->")[0]
                    donors = []
                    for a in donor.split(","):
                        spe = bacteria_species_list[a.split("|")[0]]
                        spe = spe.replace(" ", "_")
                        a = a.split("|")[0]+"|"+spe+'|'+a.split("|")[1]+"|"+a.split("|")[2]
                        donors.append(a.strip())

                    USRQ_id = list(recip_prot)[0].strip().replace('|','_')
                    USRQ_i = int(USRQ_id.split("_")[-1])-1

                    with open ("{}/01.prodigal/USRQ.gff".format(USRQ_path), 'r') as gff_file:
                        line = gff_file.readlines()[USRQ_i+3]
                        annotation = int(line.split("\t")[4])-1

                        annotated_gene = annotating_gene(annotation)
                        try:
                            gene_name = annotated_gene['gene']
                        except:
                            gene_name = 'NA'
                        try:
                            product_name = annotated_gene['gene']
                        except:
                            product_name = 'NA'
                        donors = ','.join(donors)
                        
                        try:
                            pfam_id = parse_pfam(USRQ_id)
                        except:
                            pfam_id = "NA"
                        try:
                            vf_id = parse_vf(USRQ_id)
                        except:
                            vf_id = 'NA'
                        try:
                            ar_id = parse_rgi(USRQ_id)
                        except:
                            ar_id = 'NA'
                        
                        result_recip.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(USRQ_id, gene_name, product_name,donors,pfam_id,vf_id,ar_id))
               
                elif 'USRQ' in donor:
                    donor_number+=1
                    donor_prot = (a for a in donor.split(",") if 'USRQ' in a)
                    
                    recip = event.split("-->")[1]
                    recips = []
                    for a in recip.split(","):
                        spe = bacteria_species_list[a.split("|")[0]]
                        spe = spe.replace(" ", "_")
                        a = a.split("|")[0]+"|"+spe+'|'+a.split("|")[1]+"|"+a.split("|")[2]
                        recips.append(a.strip())

                    USRQ_id = list(donor_prot)[0].strip().replace('|','_')
                    USRQ_i = int(USRQ_id.split("_")[-1])-1

                    with open ("{}/01.prodigal/USRQ.gff".format(USRQ_path), 'r') as gff_file:
                        line = gff_file.readlines()[USRQ_i+3]
                        annotation = int(line.split("\t")[4])-1

                        annotated_gene = annotating_gene(annotation)
                        try:
                            gene_name = annotated_gene['gene']
                        except:
                            gene_name = 'NA'
                        try:
                            product_name = annotated_gene['gene']
                        except:
                            product_name = 'NA'
                        recipients = ','.join(recips)
                        
                        try:
                            pfam_id = parse_pfam(USRQ_id)
                        except:
                            pfam_id = "NA"
                        try:
                            vf_id = parse_vf(USRQ_id)
                        except:
                            vf_id = 'NA'
                        try:
                            ar_id = parse_rgi(USRQ_id)
                        except:
                            ar_id = 'NA'
                        
                        result_donor.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(USRQ_id, gene_name, product_name,recipients,pfam_id,vf_id,ar_id))
                            
    result_complete = open ('{}/USRQ_HGT_info.txt'.format(USRQ_path), 'a')
    result_complete.write('###The total number or predicted HGT genes is {}.\n'.format(recip_number+donor_number))
    result_complete.write('###{} genes are predicted to be received, {} genes are predicted to be donated.\n'.format(recip_number, donor_number))
    result_complete.write('###{0:.2f}% of qeury genome is predicted to be consisted of HGT genes.\n'.format((recip_number+donor_number)/USRQ_prot_num*100))
    result_complete.write('###{0:.2f}% received and {1:.2f}% donated.\n'.format(recip_number/USRQ_prot_num*100, donor_number/USRQ_prot_num*100))
    
    result_complete.write('###USRQ_as_recip\tgene\tproduct\tdonors\tpfam\tVF\tAR\n')
    result_complete.write(''.join(result_recip))
    result_complete.write('###USRQ_as_donor\tgene\tproduct\trecipients\tpfam\tVF\tAR\n')
    result_complete.write(''.join(result_donor))
    
    result_complete.close()
    
    e=time.monotonic()
    print ('duration = {}'.format(timedelta(seconds = e-s)))

#statistics file 만들기 파일 1개 for archaea

def USRQ_HGT_statistics_arc():
    prokka_gbk = '{}/00.prokka/USRQ.gbk'.format(prokka_path)
    hgt_resultfile = open('{}/hgt_result_USRQ2.txt'.format(usrq_ranger_path),'r')
    predicted_hgt = hgt_resultfile.read().count("-->")

    with open ('{}/01.prodigal/USRQ.gff'.format(USRQ_path), 'r') as USRQ_gff:
        lines = USRQ_gff.readlines()
        USRQ_prot_num = len(lines)

    recip_number = 0
    donor_number = 0
    
    result_donor = []
    result_recip = []
    
    with open ('{}/07.ranger/hgt_result_USRQ2.txt'.format(USRQ_path), 'r') as hgt_result:
        for line in hgt_result:
            group = line.split("\t")[0]
            events = line.split("\t")[1:]

            for event in events:
                recip = event.split("-->")[1]
                donor = event.split("-->")[0]

                if 'USRQ' in recip:
                    recip_number+=1
                    recip_prot = (a for a in recip.split(",") if 'USRQ' in a)
                    
                    donor = event.split("-->")[0]
                    donors = []
                    for a in donor.split(","):
                        spe = archaea_species_list[a.split("|")[0]]
                        spe = spe.replace(" ", "_")
                        a = a.split("|")[0]+"|"+spe+'|'+a.split("|")[1]+"|"+a.split("|")[2]
                        donors.append(a.strip())

                    USRQ_id = list(recip_prot)[0].strip().replace('|','_')
                    USRQ_i = int(USRQ_id.split("_")[-1])-1

                    with open ("{}/01.prodigal/USRQ.gff".format(USRQ_path), 'r') as gff_file:
                        line = gff_file.readlines()[USRQ_i+3]
                        annotation = int(line.split("\t")[4])-1

                        annotated_gene = annotating_gene(annotation)
                        try:
                            gene_name = annotated_gene['gene']
                        except:
                            gene_name = 'NA'
                        try:
                            product_name = annotated_gene['gene']
                        except:
                            product_name = 'NA'
                        donors = ','.join(donors)
                        
                        try:
                            pfam_id = parse_pfam(USRQ_id)
                        except:
                            pfam_id = "NA"
                        
                        result_recip.append('{}\t{}\t{}\t{}\t{}\n'.format(USRQ_id, gene_name, product_name,donors,pfam_id))
               
                elif 'USRQ' in donor:
                    donor_number+=1
                    donor_prot = (a for a in donor.split(",") if 'USRQ' in a)
                    
                    recip = event.split("-->")[1]
                    recips = []
                    for a in recip.split(","):
                        spe = archaea_species_list[a.split("|")[0]]
                        spe = spe.replace(" ", "_")
                        a = a.split("|")[0]+"|"+spe+'|'+a.split("|")[1]+"|"+a.split("|")[2]
                        recips.append(a.strip())

                    USRQ_id = list(donor_prot)[0].strip().replace('|','_')
                    USRQ_i = int(USRQ_id.split("_")[-1])-1

                    with open ("{}/01.prodigal/USRQ.gff".format(USRQ_path), 'r') as gff_file:
                        line = gff_file.readlines()[USRQ_i+3]
                        annotation = int(line.split("\t")[4])-1

                        annotated_gene = annotating_gene(annotation)
                        try:
                            gene_name = annotated_gene['gene']
                        except:
                            gene_name = 'NA'
                        try:
                            product_name = annotated_gene['gene']
                        except:
                            product_name = 'NA'
                        recipients = ','.join(recips)
                        
                        try:
                            pfam_id = parse_pfam(USRQ_id)
                        except:
                            pfam_id = "NA"
                        
                        result_donor.append('{}\t{}\t{}\t{}\t{}\n'.format(USRQ_id, gene_name, product_name,recipients,pfam_id))
                            
    result_complete = open ('{}/HGT_USRQ_stats.txt'.format(USRQ_path), 'a')
    result_complete.write('###The total number or predicted HGT genes is {}.\n'.format(recip_number+donor_number))
    result_complete.write('###{} genes are predicted to be received, {} genes are predicted to be donated.\n'.format(recip_number, donor_number))
    result_complete.write('###{0:.2f}% of qeury genome is predicted to be consisted of HGT genes.\n'.format((recip_number+donor_number)/USRQ_prot_num*100))
    result_complete.write('###{0:.2f}% received and {1:.2f}% donated.\n'.format(recip_number/USRQ_prot_num*100, donor_number/USRQ_prot_num*100))
    
    result_complete.write('###USRQ_as_recip\tgene\tproduct\tdonors\tpfam\n')
    result_complete.write(''.join(result_recip))
    result_complete.write('###USRQ_as_donor\tgene\tproduct\trecipients\tpfam\n')
    result_complete.write(''.join(result_donor))
    
    result_complete.close()

def main(USRinput, processing_kingdom, identity, e_value, coverage, ranger_mode):
    if processing_kingdom == 'bacteria':
        start=time.monotonic()
        rename(USRinput)
        prokka_usrq()
        print('prokka done')
        prodigal_usrq()
        print('prodigal done')
        rrna_usrq()
        print('barnnap done')
        DIAMOND_data_as_db(identity, e_value, coverage, ranger_mode)
        print('diamond 1 done')
        DIAMOND_usrq_as_db(identity, e_value, coverage, ranger_mode)
        print('diamond 2 done')
#         filter_diamond_result()
        new_groups_with_USRQ()
        print('new groups done')
        align_new_groups()
        print('align groups done')
        fasttree_groups()
        print('fasttree done')
        merge_trees()
        optroot()
        print('tree optimizing done')
        merge_opts()
        ranger_dtl_T3()
        print('ranger T3 done')
        if ranger_mode == 'strict':
            ranger_dtl_T4()
            print('ranger T4 done')
        aggregate_ranger()
        print('aggregate done')
        
        parse_hgt_bac_1()
        gtdbtk()
        print('gtdbtk done')
        USRQ_gtdbtk_bac()
        parse_hgt_bac_2()
        putative_hgt_genes()
        USRQ_vf()
        print('vf done')
#         parse_vf(USRQ_id)
        USRQ_rgi()
        print('CARD_rgi done')
#         parse_rgi(USRQ_id)
        USRQ_pfam()
        print('pfam scan done')
#         parse_pfam(USRQ_id)
#         annotating_gene(i)
        USRQ_HGT_statistics_bac()
        print('HGT2 done')
        
        end=time.monotonic()
        print ()
        print ('total duration = {}'.format(timedelta(seconds = end-start)))
        
    elif processing_kingdom == 'archaea':
        rename(usrq)
        prokka_usrq()
        prodigal_usrq()
        rrna_usrq()
        DIAMOND_data_as_db(identity, e_value, coverage)
        DIAMOND_usrq_as_db(identity, e_value, coverage)
        filter_diamond_result()
        new_groups_with_USRQ()
        align_new_groups()
        fasttree_groups()
        merge_trees()
        optroot()
        merge_opts()
        ranger_dtl_T3()
        
        if ranger_mode == 'strict':
            ranger_dtl_T4()
            
        aggregate_ranger()

        parse_hgt_arc_1()
        gtdbtk()
        USRQ_gtdbtk_arc()
        parse_hgt_arc_2()
        putative_hgt_genes()
        USRQ_pfam()
        parse_pfam(USRQ_id)
        annotating_gene(i)
        USRQ_HGT_statistics_arc()
    
    else:
        print ('Kingdom not valid')

main(USRinput, processing_kingdom, id_value, e_value, coverage_value, ranger_mode)
