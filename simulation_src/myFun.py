import re
import ete3
import numpy as np


def change_seq_format(original_path, save_path):
    # delete >
    old_str1 = '>'
    old_str2 = '\n'
    new_str = ""
    file_data2 = ""  # for recHMM
    with open(original_path, "r", encoding="utf-8") as f:
        for line in f:
            if old_str1 in line:
                line2 = line.replace(old_str1, old_str2 + '>')
            else:
                line2 = line.replace(old_str2, new_str)

            file_data2 += line2
        file_data2 += '\n'
    file_data2 = file_data2[1:]
    with open(save_path, "w+", encoding="utf-8") as f:
        f.write(file_data2)


def newick2ms(speciestree):
    # 将 newick格式的树转变成 ms 中split的语句
    import linecache
    import re
    stack = []
    mssplit = ''  # str
    newstring = ''  # str
    current_i = 0  # record current index in speciestree
    for char in speciestree:
        if char == '(':
            stack.append(char)
            newstring = newstring + char
        elif char == ')':
            stack.pop()
            bracket_list = [m.start() for m in re.finditer('\(', newstring)]  # all ( index in newstring
            start_ind = bracket_list.pop()
            end_ind = len(newstring)
            tempstr = newstring[start_ind + 1:end_ind]  # extract content in the bracket, not include bracket
            numstr = list(map(str, re.split('[,]', tempstr)))  # list(str), 每个taxa为一个item
            # write ms command and rewrite trees in newstring
            tempsort = []
            CA_taxa = []
            for ss in numstr:
                numseq = list(map(float, re.split('[:]', ss)))  # each item is float, taxa:t1:t2
                tempsort.append(numseq)
            new_numstr = sorted(tempsort, key=lambda x: x[1])  # 按高度排序，总是低的合并到高的上
            temptime = new_numstr[0]  # 括号内合并的时间肯定是相同的，这里选第一个的时间
            if len(temptime) == 2:
                time = temptime[1]  # float
            else:
                time = temptime[1] + temptime[2]
            for tt in new_numstr:
                CA_taxa.append(tt[0])  # 升序 common ancestor taxa合并顺序 [3,2,1] 3->2, 2->1 ,合并到最后一个上
            for i in range(len(CA_taxa) - 1):  # 同源taxa number not fixed， 循环输出
                mssplit += ' -ej ' + str(time) + ' ' + str(int(CA_taxa[i])) + ' ' + str(int(CA_taxa[i + 1]))

            newstring = newstring[:start_ind] + str(int(CA_taxa[len(CA_taxa) - 1])) + ':' + str(time)

        else:  # char except ( , )
            newstring = newstring + char

    return mssplit







