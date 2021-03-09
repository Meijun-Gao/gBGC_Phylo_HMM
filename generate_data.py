"""
 use msHOT algo
 generate positive control, base recombination rate=Len/1000; base mutation rate=1
 when hotspot_num=0, recombination and mutation rate are base rr and mr.
"""
import sys, getopt
import os
import linecache
from simulation_src.myFun import change_seq_format
import ete3
from simulation_src.myFun import newick2ms


# hot_rho_set = [0, 1, 4, 6, 8, 10]      # multiple the base rho; control recombination rate in msHOT
# hot_dScale_set = [0, 1, 4, 6, 8, 10]   # multiple the base dscale; control mutation rate in seq-gen
# NumTaxa = 4
# modelname = 'HKY'
# hotspot_num = 1
# hot_rho = hot_rho_set[5]  # recombination rate       N0=(1/4)10^6 , Rr=10^-8 neighbor base pair
# hot_dScale = hot_dScale_set[5]
# FactualPath = '../original-data/msHOT/dif_hky2/taxa' + str(NumTaxa) + '/Btest'


def main(argv):
    usage = '''
    python generate_hotspot.py --numtaxa=4 --modelname='HKY' --hotspot_num=1 
    --rr_multiple=10 --mr_multiple=10 --path=''
    Options:
    -h show the basic usage for simulation data generation
    --numtaxa 4 or 5;   the length of the sequences are 5000, 2000 for 4 and 5 taxa, respectively
    --modelname HKY or GTR; default HKY
    --hotspot_num 0 or 1 or 2
        when hotspot_num=0, --rr_multiple --mr_multiple do not need  
    --hot_rho recombination rate (rr) multiple based on the rr in background
    --hot_dScale mutation rate (mr) multiple based on the mr in background
    --path PATH, where to put output files
    
    '''
    # default
    hot_rho = 10
    hot_dScale = 10
    if len(argv) <= 2:
        print(usage)
        sys.exit()

    try:
        opts, args = getopt.getopt(argv[1:], "-h", ["modelname=", "numtaxa=", "hotspot_num=", "rr_multiple=",
                                                    "mr_multiple=", "path="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt == "--modelname":
            modelname = arg
        elif opt == "--numtaxa":
            NumTaxa = int(arg)
        elif opt == "--hotspot_num":
            hotspot_num = int(arg)
        elif opt == "--rr_multiple":
            hot_rho = float(arg)
        elif opt == "--mr_multiple":
            hot_dScale = float(arg)
        elif opt == "--path":
            FactualPath = arg

    if NumTaxa == 4:
        numSites = 5000
        Tree_num = 3
    else:
        numSites = 2000
        Tree_num = 15
    # parameters   4taxa: 5000;  5taxa: 2000
    NumTrials = 3   # default 20
    base_rho = numSites/1000     # 1  for 1000
    base_dScale = 1  # default 1 scale     0.5, 1, 2
    fileContent = linecache.getlines('simulation_src/tree_files/taxa' + str(NumTaxa) + '_modeltree.txt')    # file contains all species trees
    fileContent2 = linecache.getlines('simulation_src/tree_files/taxa' + str(NumTaxa) + '_dif_topo.txt')   # used for infer
    os.mkdir(FactualPath)  # make a new directory named this one
    ########
    if hotspot_num == 1:
        if NumTaxa == 4:   # 4 taxa hotspot location (2000, 4000)
            begin = 2000
            end = 4000
        else:              # 5 taxa hotspot location (500, 13000)
            begin = 500
            end = 1300
    else:
        if NumTaxa == 4:   # 4 taxa 2 hotspots (1000, 3000) (4000,4500)
            begin = 1000
            end = 3000
            begin2 = 4000
            end2 = 4500
        else:              # 5 taxa 2 hotspots (500, 1200) (1500,1800)
            begin = 500
            end = 1200
            begin2 = 1500
            end2 = 1800

    base_freq = '0.267 0.200 0.213 0.320'
    base_freq_hot = '0.162 0.347 0.341 0.150'   # ACGT; AC;AG;AT;CG;CT;GT
    gtr_mat = '0.852139 2.726123 0.670220 0.742474 2.932312 1.000000'
    gtr_mat_hot = '0.994225 2.382991 0.603187 0.358835 2.312774 1.000000'
    tstv = 3.709560/2.0   # 1.85478
    tstv_hot = 4.116433/2.0  # 2.0582165
    I_pop = ' '.join(['1'] * NumTaxa)
    # generate data
    for trialCount in range(NumTrials):
        if hotspot_num == 1:
            idx_list = [begin, end, numSites]
            seqlen_list = [begin, end - begin, numSites - end]
            BP_list = [[], [], []]
            seq_list = [[], [], []]
        else:  # 2 hotspot
            idx_list = [begin, end, begin2, end2, numSites]
            seqlen_list = [begin, end - begin, begin2 - end, end2 - begin2, numSites - end2]
            BP_list = [[], [], [], [], []]
            seq_list = [[], [], [], [], []]

        path1 = FactualPath + '/REPLICATE-R' + str(trialCount + 1)
        os.mkdir(path1)  # put final genetrees and sequence
        # get random tree
        temp_tree = fileContent[trialCount]  # i speciestree all info  fix the parental tree
        oldspeciestree = temp_tree[temp_tree.index('('):temp_tree.index(';') + 1]  # i speciestree tree string
        speciestree = oldspeciestree.replace("t", "")
        mssplit = newick2ms(speciestree)
        # generate local gene trees according to species tree
        if hotspot_num == 0:
            os.system('''simulation_src/ms %s 1 -I %s %s -r %s %s %s -T > \
                            %s/REPLICATE-R%s/local-genetrees''' % (
                NumTaxa, NumTaxa, I_pop, base_rho, numSites, mssplit, FactualPath, trialCount + 1))

            localpath = FactualPath + '/REPLICATE-R' + str(trialCount + 1)
            t = open(localpath + '/local-genetrees', 'r').readlines()[4:]
            b = ''.join(t)
            p_value = len(t)  # get the number of lines in local-genetrees

            with open(localpath + '/local-genetrees', "w+", encoding="utf-8") as f:
                f.write(b)
            if modelname == 'GTR':
                os.system('simulation_src/seq-gen -m GTR -f %s -r %s -of -s %s -l %s -p %s\
                               < %s/REPLICATE-R%s/local-genetrees > \
                               %s/REPLICATE-R%s/sequence_file' % (
                    base_freq, gtr_mat, base_dScale, numSites, p_value, FactualPath, trialCount + 1, FactualPath,
                    trialCount + 1))
            else:  # hky
                os.system('simulation_src/seq-gen -m HKY -f %s -t %s -of -s %s -l %s -p %s\
                                       < %s/REPLICATE-R%s/local-genetrees > \
                                       %s/REPLICATE-R%s/sequence_file' % (
                    base_freq, tstv, base_dScale, numSites, p_value, FactualPath, trialCount + 1, FactualPath,
                    trialCount + 1))

            SequenceFile = localpath + '/sequence_file'
            Current_SeqFile2 = localpath + '/sequence.fasta'
            change_seq_format(SequenceFile, Current_SeqFile2)
            os.system('''rm %s/sequence_file*''' % localpath)

        else:
            if hotspot_num == 1:
                os.system('''simulation_src/msHOT/msHOT %s 1 -r %s %s -v 1 %s %s %s -I %s %s %s -T > \
                %s/REPLICATE-R%s/local-genetrees''' % (NumTaxa, base_rho, numSites, begin, end, hot_rho,
                                                       NumTaxa, I_pop, mssplit, FactualPath, trialCount + 1))
            else:  # 2 hot spots
                os.system('''simulation_src/msHOT/msHOT %s 1 -r %s %s -v 2 %s %s %s %s %s %s -I %s %s %s -T > \
                            %s/REPLICATE-R%s/local-genetrees''' % (NumTaxa, base_rho, numSites, begin, end, hot_rho,
                                                                   begin2, end2, hot_rho,
                                                                   NumTaxa, I_pop, mssplit, FactualPath, trialCount + 1))

            localpath = FactualPath + '/REPLICATE-R' + str(trialCount + 1)
            t = open(localpath + '/local-genetrees', 'r').readlines()[4:]
            b = ''.join(t)
            p_value = len(t)  # get the number of lines in local-genetrees
            with open(localpath + '/local-genetrees', "w+", encoding="utf-8") as f:
                f.write(b)

            # generate different tree files
            genetreeLine = open(localpath + '/local-genetrees', 'r').readlines()
            treenum = 0  # use to record current index (the head)
            head_treenum = 0  # use to record current index (the head)
            t = 0
            for k in range(len(genetreeLine)):
                num_end_ind = genetreeLine[k].index(']')
                current_num = int(genetreeLine[k][1:num_end_ind])
                treenum += current_num

                if treenum <= idx_list[t]:
                    BP_list[t].append(genetreeLine[k])
                    if treenum == idx_list[t]:
                        t += 1
                    head_treenum = treenum
                elif idx_list[t] < treenum <= idx_list[t+1]:   #
                    temp_num1 = idx_list[t] - head_treenum
                    changeline1 = '[' + str(temp_num1) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t].append(changeline1)
                    temp_num2 = current_num - temp_num1
                    changeline2 = '[' + str(temp_num2) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t+1].append(changeline2)
                    if treenum == idx_list[t+1]:
                        t += 1
                    t += 1

                    head_treenum = treenum
                elif idx_list[t+1] < treenum <= idx_list[t+2]:
                    temp_num1 = idx_list[t] - head_treenum
                    changeline1 = '[' + str(temp_num1) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t].append(changeline1)
                    changeline2 = '[' + str(idx_list[t+1] - idx_list[t]) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t+1].append(changeline2)
                    temp_num2 = current_num - temp_num1 - (idx_list[t+1] - idx_list[t])
                    changeline3 = '[' + str(temp_num2) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t+2].append(changeline3)
                    if treenum == idx_list[t+2]:
                        t += 1
                    t += 2
                    head_treenum = treenum

                elif idx_list[t+2] < treenum <= idx_list[t+3]:
                    temp_num1 = idx_list[t] - head_treenum
                    changeline1 = '[' + str(temp_num1) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t].append(changeline1)
                    changeline2 = '[' + str(idx_list[t+1] - idx_list[t]) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t+1].append(changeline2)
                    changeline3 = '[' + str(idx_list[t + 2] - idx_list[t+1]) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t + 2].append(changeline3)
                    temp_num2 = current_num - temp_num1 - (idx_list[t+1] - idx_list[t]) - (idx_list[t + 2] - idx_list[t+1])
                    changeline4 = '[' + str(temp_num2) + ']' + genetreeLine[k][num_end_ind + 1:]
                    BP_list[t+3].append(changeline3)
                    if treenum == idx_list[t+3]:
                        t += 1
                    t += 3
                    head_treenum = treenum
                else:
                    print('error happen, maybe ')

            for x in range(len(BP_list)):
                p_value = len(BP_list[x])
                bk = ''.join(BP_list[x])
                with open(localpath + '/BP_trees' + str(x), "w") as f:
                    f.write(bk)
                if x % 2 == 0:
                    dScale = base_dScale
                    base_freq_t = base_freq
                    gtr_mat_t = gtr_mat
                    tstv_t = tstv
                else:   # hotspot
                    dScale = hot_dScale
                    base_freq_t = base_freq_hot
                    gtr_mat_t = gtr_mat_hot
                    tstv_t = tstv_hot
                if modelname == 'GTR':
                    os.system('simulation_src/seq-gen -m GTR -f %s -r %s -of -s %s -l %s -p %s\
                                               < %s/REPLICATE-R%s/BP_trees%s > \
                                               %s/REPLICATE-R%s/sequence_file%s' % (
                        base_freq_t, gtr_mat_t,
                        dScale, seqlen_list[x], p_value, FactualPath, trialCount + 1, x, FactualPath, trialCount + 1, x))
                elif modelname == 'HKY':
                    os.system('simulation_src/seq-gen -m HKY -f %s -t %s -of -s %s -l %s -p %s\
                                                               < %s/REPLICATE-R%s/BP_trees%s > \
                                                               %s/REPLICATE-R%s/sequence_file%s' % (
                        base_freq_t, tstv_t,
                        dScale, seqlen_list[x], p_value, FactualPath, trialCount + 1, x, FactualPath, trialCount + 1, x))

            # change the original seq file
            for x in range(len(BP_list)):
                change_seq_format(localpath + '/sequence_file' + str(x), localpath + '/tempseq_' + str(x))
            # assembly the sequences
            for x in range(len(BP_list)):
                seq_list[x] = open(localpath + '/tempseq_' + str(x)).readlines()
            seqfile = ''
            for j in range(NumTaxa):
                seqfile += '>' + str(j+1) + '\n'
                for x in range(len(seq_list)):
                    for s in range(1, len(seq_list[x]), 2):
                        if str(j + 1) in seq_list[x][s-1]:
                            seqfile += seq_list[x][s][:-1]

                seqfile += '\n'
            with open(localpath + '/sequence.fasta', "w", encoding="utf-8") as f:
                f.write(seqfile)

            os.system('''rm %s/BP_*''' % localpath)
            os.system('''rm %s/sequence_file*''' % localpath)
            os.system('''rm %s/tempseq_*''' % localpath)

    # generate all different unrooted gene trees
    for trialCount in range(1, NumTrials + 1):
        path = FactualPath + '/REPLICATE-R' + str(trialCount)
        GeneTreesLines = open(path + '/local-genetrees', 'r').readlines()
        new_local_trees = ''
        for line in GeneTreesLines:
            num_end_ind = line.index(']')
            num_GeneTree = int(line[1:num_end_ind])
            Current_GeneTree = line[num_end_ind + 1:]
            for i in range(num_GeneTree):
                new_local_trees += Current_GeneTree

        # generate gene trees in each site, for debug
        # with open(path + '/new-local-genetrees', "w+", encoding="utf-8") as f:
        #     f.write(new_local_trees)

        ttt = []
        for i in range(len(GeneTreesLines)):
            begin_i = GeneTreesLines[i].index(']')
            temp_str = GeneTreesLines[i][begin_i + 1:-1]
            ttt.append(temp_str)

        Treelist = []
        for i in range(len(ttt)):
            temp_t = ete3.Tree(ttt[i])
            if i == 0:
                temp_x = temp_t.write(format=5) + '\n'
                Treelist.append(temp_x)
            else:
                flag = 0
                for t in range(len(Treelist)):
                    temp_t2 = ete3.Tree(Treelist[t])
                    if temp_t.robinson_foulds(temp_t2, unrooted_trees=True)[0] != 0:
                        flag += 1

                if flag == len(Treelist):
                    temp_x = temp_t.write() + '\n'
                    Treelist.append(temp_x)

        if len(Treelist) < Tree_num:
            for x in range(len(fileContent2)):
                current_tree = fileContent2[x][:-1]
                current_tree = current_tree.replace("t", "")
                tree1 = ete3.Tree(current_tree)
                flag = 0
                for t in range(len(Treelist)):
                    tree2 = ete3.Tree(Treelist[t])
                    if tree1.robinson_foulds(tree2, unrooted_trees=True)[0] != 0:
                        flag += 1
                if flag == len(Treelist):
                    Treelist.append(current_tree + '\n')
        print('The number of trees:' + str(len(Treelist)))
        final_treestr = ''.join(Treelist)
        with open(path + '/all-genetrees', "w+") as f:
            f.write(final_treestr)

    print('data simulation finished')


if __name__ == "__main__":
    main(sys.argv)
