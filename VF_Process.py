from VF_variables import *
#  Parameters you might want to change.
#	ESA switch threshold is "PercentRESA * MeanAllESA" or the (Minimum sequence ESA)+0.00001 so that there is at least one ESA switch.
PercentRESA = 0.75
#
inlist = open('ListFiles.txt', 'r')
Files = inlist.read()
FileList = Files.rsplit()
inlist.close()
for iFL in range(len(FileList)):
    infile = open(FileList[iFL]+'_In.txt', 'r')
    RespList = infile.read()
    DataList = RespList.rsplit()
    infile.close()
    #
    outfile = open('output_' + FileList[iFL]+'_Out.txt', 'w')
    #
    AnimalIndex = []
    import math
    import statistics
    LogFreq = []
    Syllabs = []
    Typical = []
    N_UNKN = 0
    #
    Subject = DataList[0]		# 120626NNAUR
    Sub_Age = DataList[1]		# 61
    Sub_Edu = DataList[2]		# 10
    WordStart = 3
    #
    for ia in range(WordStart, len(DataList)):
        Looking = True
        for i in range(len(VF_Animals)):
            if ((DataList[ia] == VF_Animals[i]) | (('UNKNOWN_'+DataList[ia]) == VF_Animals[i])):
                AnimalIndex = AnimalIndex + [i]
                Looking = False
                break
        if Looking:
            i = len(VF_Animals)
            AnimalIndex = AnimalIndex + [len(VF_Animals)]
            VF_Animals = VF_Animals + ['UNKNOWN_'+DataList[ia]]
            N_UNKN = N_UNKN + 1
            WordFreq = WordFreq + [1099]
            Syllables = Syllables + [0]
            Typicality = Typicality + [0.0022]
            TroyCat = TroyCat + [('none', 'none', 'none', 'none')]
        LogFreq = LogFreq + [math.log10(WordFreq[i])]
        Syllabs = Syllabs + [Syllables[i]]
        Typical = Typical + [Typicality[i]]

    #
    Dups = []
    N_Dups = 0
    for i in range(len(AnimalIndex)):
        DUP = 0
        if i > 0:
            for idup in range(i):
                if DataList[idup+WordStart] == DataList[i+WordStart]:
                    DUP = 1
        Dups = Dups + [DUP]
        N_Dups = N_Dups + DUP

    N_Totl = len(AnimalIndex)
    N_Uniq = N_Totl - N_Dups - N_UNKN
    outfile.write('Number of animals:  ')
    outfile.write(format(N_Totl, '>3d'))
    outfile.write('\nNumber correct:      ')
    outfile.write(format(N_Uniq, '>3d'))
    outfile.write('\nRepetitions:        ')
    outfile.write(format(N_Dups, '>3d'))
    outfile.write('\nUnknown Animals:    ')
    outfile.write(format(N_UNKN, '>3d'))
    outfile.write('\n\nMean Log Freq:      ')
    outfile.write(format(statistics.mean(LogFreq), '10.4f'))
    outfile.write('  (')
    outfile.write(format(statistics.stdev(LogFreq), '6.4f'))
    outfile.write(')')
    outfile.write('\nMean Syllables:     ')
    outfile.write(format(statistics.mean(Syllabs), '10.4f'))
    outfile.write('  (')
    outfile.write(format(statistics.stdev(Syllabs), '6.4f'))
    outfile.write(')')
    outfile.write('\nMean Typicality     ')
    outfile.write(format(statistics.mean(Typical), '10.4f'))
    outfile.write('  (')
    outfile.write(format(statistics.stdev(Typical), '6.4f'))
    outfile.write(')')
    #
    #
    ESA = []
    #
    for i in range(len(AnimalIndex)):
        if i == 0:
            ESA = []
            SumESA = 0.0
        elif ('UNKNOWN' in VF_Animals[AnimalIndex[i]]) | ('UNKNOWN' in VF_Animals[AnimalIndex[i-1]]):
            ESA = ESA + [0.0000000]
        else:
            OneESA = ESA_2s[AnimalIndex[i]][AnimalIndex[i-1]]
            ESA = ESA + [OneESA]
            SumESA = SumESA + OneESA

    #
    MeanESA = SumESA / len(ESA)
    #
    N_AllESA = 0
    SumAllESA = 0.0
    #
    tmp = 0
    for i in range(0, (len(AnimalIndex)-1)):
        for j in range(i+1, len(AnimalIndex)):
            if (Dups[i] == 1) | (Dups[j] == 1):
                tmp = tmp + 1
            elif ('UNKNOWN' in VF_Animals[AnimalIndex[i]]) | ('UNKNOWN' in VF_Animals[AnimalIndex[j]]):
                N_AllESA = N_AllESA + 1
            else:
                OneESA = ESA_2s[AnimalIndex[i]][AnimalIndex[j]]
                N_AllESA = N_AllESA + 1
                SumAllESA = SumAllESA + OneESA

    #
    if (N_AllESA == 0):
        MeanAllESA = 0.0
        SOI = 0.0
    else:
        MeanAllESA = SumAllESA / N_AllESA
        if MeanAllESA == 0.0:
            SOI = 0.0
        else:
            SOI = MeanESA / MeanAllESA
    #
    outfile.write('\n\nSequence Mean ESA:  ')
    outfile.write(format(MeanESA, '10.4f'))
    outfile.write('\nTotal Mean ESA:     ')
    outfile.write(format(MeanAllESA, '10.4f'))
    outfile.write('\nSem Org Index (SOI):')
    outfile.write(format(SOI, '10.4f'))
    #
    ESAsort = []		# Need independent ESA list to sort, for quartiles.
    for i in range(len(ESA)):
        ESAsort = ESAsort + [ESA[i]]

    #
    ESAsort.sort()
    ESA_Threshold = PercentRESA * MeanAllESA  # ESA switch if ESA < ESA_Threshold
    if ESA_Threshold < ESAsort[0]:
        ESA_Threshold = ESAsort[0] + 0.000001
    #
    ESAsw = []  # Initialize ESA switch array.
    ESAcs = []  # Initialize ESA cluster array.
    for i in range(len(AnimalIndex)-1):
        ESAcs = ESAcs + [' ']
        if ESA[i] < ESA_Threshold:
            ESAsw = ESAsw + [1]
        else:
            ESAsw = ESAsw + [0]

    # There is one ESAcs for each word, but one less ESAsw
    ESAcs = ESAcs + [' ']
    iiESAsw = len(ESAsw) - 1
    ii = 1
    #  SumCS and NumCS to compute mean multi-word cluster size.
    SumCS = 0.0
    NumCS = 0.0
    for i in range(iiESAsw, -1, -1):
        if (ESAsw[i] == 0):
            ii = ii + 1
        else:
            #		print (i, ii)
            ESAcs[i+1] = str(ii)
            if (ii > 1):
                SumCS = SumCS + ii
                NumCS = NumCS + 1
            ii = 1
    ESAcs[0] = str(ii)
    if (ii > 1):
        SumCS = SumCS + ii
        NumCS = NumCS + 1
    #
    if NumCS == 0:
        NumCS = 1
    sumESAsw = sum(ESAsw)
    if sumESAsw == 0:
        sumESAsw = 1
    outfile.write('\n\nESA switches:       ')
    outfile.write(format(sum(ESAsw), '>4d'))
    outfile.write('  (')
    outfile.write(format(N_Totl/sumESAsw, '5.2f'))
    outfile.write(' words/switch)')
    outfile.write('  (')
    outfile.write(format(SumCS/NumCS, '5.2f'))
    outfile.write(' words/cluster)')
    #
    TROYsw = []  # Initialize TROY switch array.
    TROYcs = []  # Initialize TROY cluster array.
    for i in range(1, len(AnimalIndex)):
        j = i - 1
        ai1 = AnimalIndex[j]
        ai2 = AnimalIndex[i]
        tsw = 1
        for it in range(0, 4):
            for jt in range(0, 4):
                if ((TroyCat[ai1][it] is not 'none') & (TroyCat[ai2][jt] is not 'none') & (TroyCat[ai1][it] is TroyCat[ai2][jt])):
                    #print ('not a switch')
                    tsw = 0
        TROYsw = TROYsw + [tsw]
        TROYcs = TROYcs + [' ']
        #print (i, j, tsw)

    # There is one TROYcs for each word, but one less TROYsw
    TROYcs = TROYcs + [' ']
    iiTROYsw = len(TROYsw) - 1
    ii = 1
    #  SumCS and NumCS to compute mean multi-word cluster size.
    SumCS = 0.0
    NumCS = 0.0
    for i in range(iiTROYsw, -1, -1):
        if (TROYsw[i] == 0):
            ii = ii + 1
        else:
            #print (i, ii)
            TROYcs[i+1] = str(ii)
            if (ii > 1):
                SumCS = SumCS + ii
                NumCS = NumCS + 1
            ii = 1

    TROYcs[0] = str(ii)
    if (ii > 1):
        SumCS = SumCS + ii
        NumCS = NumCS + 1

    #
    if NumCS == 0:
        NumCS = 1
    sumTROYsw = sum(TROYsw)
    if sumTROYsw == 0:
        sumTROYsw = 1
    outfile.write('\n\nTROY switches:      ')
    outfile.write(format(sum(TROYsw), '>4d'))
    outfile.write('  (')
    outfile.write(format(N_Totl/sumTROYsw, '5.2f'))
    outfile.write(' words/switch)')
    outfile.write('  (')
    outfile.write(format(SumCS/NumCS, '5.2f'))
    outfile.write(' words/cluster)')

    for i in range(len(AnimalIndex)):
        ai = AnimalIndex[i]
        n = i + 1
        DUP = 0
        if i == 0:
            outfile.write(
                '\n\n\nOrder\tWord\t   Rep\t Freq\t Syll\tTypic\tESA\tESA-sw\tESA-CS\tT-cat\t    T-sw\tT-CS\n')
        else:
            for idup in range(i):
                if DataList[idup+WordStart] == DataList[i+WordStart]:
                    DUP = 1

        LocCat = TroyCat[ai][0]
        if (TroyCat[ai][1] != 'none'):
            LocCat = LocCat + ',' + TroyCat[ai][1]
        if (TroyCat[ai][2] != 'none'):
            LocCat = LocCat + ',' + TroyCat[ai][2]
        if (TroyCat[ai][3] != 'none'):
            LocCat = LocCat + ',' + TroyCat[ai][3]
        outfile.write(format(n, '>2d'))
        outfile.write(format(VF_Animals[ai], '^17s'))
        outfile.write(format(DUP, '>2d'))
        outfile.write(format(WordFreq[ai], '8d'))
        outfile.write(format(Syllables[ai], '7d'))
        outfile.write(format(Typicality[ai], '10.4f'))
        if i > 0:
            outfile.write(format(ESA[i-1], '8.4f'))
        else:
            outfile.write('\t')
        if i > 0:
            outfile.write(format(ESAsw[i-1], '6d'))
        else:
            outfile.write('\t')
        outfile.write('\t')
        outfile.write(format(ESAcs[i], '>4s'))
        outfile.write('\t')
        outfile.write(format(LocCat, '12s'))
    #	outfile.write('\t')
        if i > 0:
            outfile.write(format(TROYsw[i-1], '>3d'))
            outfile.write('\t')
        else:
            outfile.write('\t')
        outfile.write('\t')
        outfile.write(format(TROYcs[i], '>3s'))
        outfile.write('\n\n')

    #   Write next-to-last line in by-subject analysis format.
    outfile.write('Subj')
    outfile.write('\t\t')
    outfile.write('Age')
    outfile.write('\t')
    outfile.write('Edu')
    outfile.write('\t')
    outfile.write('Words')
    outfile.write('\t')
    outfile.write('Corr')
    outfile.write('\t')
    outfile.write('Dups')
    outfile.write('\t')
    outfile.write('NonAnml')
    outfile.write('\t')
    outfile.write('MlogF')
    outfile.write('\t')
    outfile.write('SDlogF')
    outfile.write('\t')
    outfile.write('Msyls')
    outfile.write('\t')
    outfile.write('SDsyls')
    outfile.write('\t')
    outfile.write('Mtyp')
    outfile.write('\t')
    outfile.write('SDtyp')
    outfile.write('\t')
    outfile.write('OrdESA')
    outfile.write('\t')
    outfile.write('AllESA')
    outfile.write('\t')
    outfile.write('SOI')
    outfile.write('\t')
    outfile.write('ESAsw')
    outfile.write('\t')
    outfile.write('Wds/sw')
    outfile.write('\t')
    outfile.write('ESA-CS')
    outfile.write('\t')
    outfile.write('TROYsw')
    outfile.write('\t')
    outfile.write('Wds/sw')
    outfile.write('\t')
    outfile.write('TROY-CS')
    outfile.write('\n')

    #   Write last line in by-subject analysis format.
    outfile.write(Subject)
    outfile.write('\t')
    outfile.write(Sub_Age)
    outfile.write('\t')
    outfile.write(Sub_Edu)
    outfile.write('\t')
    outfile.write(format(N_Totl, 'd'))
    outfile.write('\t')
    outfile.write(format(N_Uniq, 'd'))
    outfile.write('\t')
    outfile.write(format(N_Dups, 'd'))
    outfile.write('\t')
    outfile.write(format(N_UNKN, 'd'))
    outfile.write('\t')
    outfile.write(format(statistics.mean(LogFreq), '6.4f'))
    outfile.write('\t')
    outfile.write(format(statistics.stdev(LogFreq), '6.4f'))
    outfile.write('\t')
    outfile.write(format(statistics.mean(Syllabs), '6.4f'))
    outfile.write('\t')
    outfile.write(format(statistics.stdev(Syllabs), '6.4f'))
    outfile.write('\t')
    outfile.write(format(statistics.mean(Typical), '6.4f'))
    outfile.write('\t')
    outfile.write(format(statistics.stdev(Typical), '6.4f'))
    outfile.write('\t')
    outfile.write(format(MeanESA, '6.4f'))
    outfile.write('\t')
    outfile.write(format(MeanAllESA, '6.4f'))
    outfile.write('\t')
    outfile.write(format(SOI, '6.4f'))
    outfile.write('\t')
    outfile.write(format(sum(ESAsw), 'd'))
    outfile.write('\t')
    outfile.write(format(N_Totl/sumESAsw, '5.2f'))
    outfile.write('\t')
    outfile.write(format(SumCS/NumCS, '5.2f'))
    outfile.write('\t')
    outfile.write(format(sum(TROYsw), 'd'))
    outfile.write('\t')
    outfile.write(format(N_Totl/sumTROYsw, '5.2f'))
    outfile.write('\t')
    outfile.write(format(SumCS/NumCS, '5.2f'))
    outfile.write('\n')
    outfile.close()
