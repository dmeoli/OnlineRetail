import copy

# Dataset structure.
# The dataset is a list of sequences.
# Each element in a sequence is a couple of:
# 	- items
#	- timestamp
# In the algorithm a sequence is a list of two dimensional tuples, 
# hence an element is a tuple.
# Each tuple is composed of a list, representing an itemset
# and a positive integer representing the timestamp.
# [ [([ ... ], t), ... , ([ ... ], t)],
#	... 
#	[([ ... ], t), ... , ([ ... ], t)] ]
#
# In this first implemenetation we mine all the frequent subsequence 
# with the maxspan, mingap and maxgap constraint.
# We say that a sequence v with time information supports a sequence s if:
#	- s is a subsequence of v
#	- if s is a k-sequence with k>1 then 
#	  the time constraints of s with respect to v are satisfied
# We say that a sequence s is a frequent-k-sequence if the support of s
# is higher then a threshold.
# Since the maxgap constrain violates the Apriori principale then we use the notion of
# contiguous subsequence.
# In particuare the new apriori principale says that "if a k-sequence is frequent then
# all its contiguous (k-1)-subsequence are frequent".
# We say that a sequence s is a candidate-k-sequence if every contiguous (k-1)-subsequence of 
# s is a frequent sequence.


# Returns True if subSequence is a subsequence of mainSequence
def isSubsequence(mainSequence, subSequence):
    subSequenceClone = list(subSequence)  # clone the sequence, because we will alter it
    return isSubsequenceRecursive(mainSequence, subSequenceClone)  # start recursion

def isSubsequenceRecursive(mainSequence, subSequenceClone, start=0):
    # Check if empty: End of recursion, all itemsets have been found
    if (not subSequenceClone):
        return True
    # retrieves element of the subsequence and removes is from subsequence
    firstElem = set(subSequenceClone.pop(0))
    # Search for the first itemset...
    # Note that mainSequence = [([ ... ], t), ... , ([ ... ], t)]
    # so mainSequence[i] = ([ ... ], t).
    # So in order to make the code work we must ignore the timestamp (add [0]).
    for i in range(start, len(mainSequence)):
        if (set(mainSequence[i][0]).issuperset(firstElem)):
            # and recurse
            return isSubsequenceRecursive(mainSequence, subSequenceClone, i + 1)
    return False


# Given a sequence and a dataset, count the frequence of the sequence in the dataset
# for seq in dataset
#	if if isSubsequence(seq, sequence)
#		return 1
# if sequence is a subsequence of seq then increase couting.
# Note that seq = [([ ... ], t), ... , ([ ... ], t)], we have to deal with the timestamp
def countFreq(sequence, dataset):
	return sum(1 
		for seq in dataset 
			if isSubsequence(seq, sequence))

# Merge two (k-1)-sequence into a k-sequence
def generateCandidatesForPair(cand1, cand2):
    cand1Clone = copy.deepcopy(cand1)
    cand2Clone = copy.deepcopy(cand2)
    
    # drop the leftmost item from cand1:
    if (len(cand1[0]) == 1):
    	# Se il 'leftmost element' contiene un solo item 
    	# (quindi l'unico da essere scartato).
        cand1Clone.pop(0)
    else:
    	# Se il 'leftmost element' contiene più di u item 
    	# (quindi l'ultimo item).
        cand1Clone[0] = cand1Clone[0][1:]
    
    # drop the rightmost item from cand2:
    if (len(cand2[-1]) == 1):
    	# Se il 'rightmost element' contiene un solo item 
    	# (quindi l'unico da essere scartato).
        cand2Clone.pop(-1)
    else:
    	# Se il 'rightmostmost element' contiene più di u item 
    	# (quindi l'ultimo item).
        cand2Clone[-1] = cand2Clone[-1][:-1]

    # if the result is not the same, then we dont need to join
    # and the function returns an empty list, indicating that the sequences
    # can't be merged. Otherwise returns the merged sequence.
    if not cand1Clone == cand2Clone:
        return []
    else:
        newCandidate = copy.deepcopy(cand1)
        if (len(cand2[-1]) == 1):
            newCandidate.append(cand2[-1])
        else:
            newCandidate[-1].extend(cand2[-1][-1])
        return newCandidate

# Generate k-sequences from a set of (k-1)-sequences.
def generateCandidates(lastLevelCandidates):
    k = sequenceSize(lastLevelCandidates[0]) + 1

    if (k == 2):
    	flatShortCandidates = [item 
    		for sublist2 in lastLevelCandidates 
    			for sublist1 in sublist2 
    				for item in sublist1]

        result = [[[a, b]] 
        	for a in flatShortCandidates 
        		for b in flatShortCandidates 
        			if b > a]

        result.extend([[[a], [b]] 
        	for a in flatShortCandidates 
        		for b in flatShortCandidates])

        return result
    else:
        candidates = []

        for i in range(0, len(lastLevelCandidates)):
            for j in range(0, len(lastLevelCandidates)):
                newCand = generateCandidatesForPair(lastLevelCandidates[i], lastLevelCandidates[j])
                if (not newCand == []):
                    candidates.append(newCand)
        candidates.sort()
        return candidates

def generateDirectSubsequences(sequence):
    result = []
    for i, itemset in enumerate(sequence):
        if (len(itemset) == 1):
            sequenceClone = copy.deepcopy(sequence)
            sequenceClone.pop(i)
            result.append(sequenceClone)
        else:
            for j in range(len(itemset)):
                sequenceClone = copy.deepcopy(sequence)
                sequenceClone[i].pop(j)
                result.append(sequenceClone)
    return result

# Return True if mseq supports cseq, False otherwire.
# Remember that:
#	- mseq = [([ ... ], t), ... , ([ ... ], t)]
#	- cseq = [[ ... ], ... , [ ... ]]
# mseq support cseq if cseq is a subsequence of mseq and and the time contraints hold. 
def supports(mseq, cseq, maxspan, mingap, maxgap):
    i = 0
    for j in range(len(mseq)):
        if set(mseq[j][0]).issuperset(set(cseq[i])):
            min_t = mseq[j][1]
            i += 1
            
            # Special case if cseq is a sequence of one element
            if i == len(cseq):
                return True
            
            prev_t = mseq[j][1]

            for itemset, t in mseq[j+1 :]:

            	# The mingap constraint is violated
            	# The time difference between the current and the previous itemset
            	# must be geater than
            	if not (t - prev_t > mingap):
            		i = 0
            		break

            	# The mingap constraint is violated
            	# The time difference between the current and the previous itemset
            	# must be geater than
            	if not (t - prev_t <= maxgap):
            		i = 0
            		break

            	# The maxspan constraint is violated
                if t - min_t > maxspan:
                    i = 0
                    break 

                if set(itemset).issuperset(set(cseq[i])):
                    i += 1

                # The whole sequence is found satisfying all the time constraints
                if i == len(cseq):
                    return True
    return False

def countSupport(dataset, cseq, maxspan, mingap. maxgap):
	support_count = 0
	for seq in dataset:
		if supports(seq, cseq, maxspan):
			support_count += 1
	return support_count

# Apriori pseudocode
#
# k = 1
# F1 = {generate all the frequent 1-sequences} (si usa solo una base frequentistica)
# while Fk is not empty
#	k = k + 1
#	Ck = {generate all the candidate-k-sequences from F(k-1)}
#	for each t in the dataset
#		for each c in Ck
#			if t supports c (cioè è una sottosequenza e rispetta il time contraint)
#				increase the support of c
#	Fk = {filter out all the sequences with a support lower than the threshold}
# return the union of all the Fk
def apriori_tc(dataset, maxspan, minsup, mingap, maxgap, verbose=False):
	# Overall is a list of lists.
	# An element of Overall is a list of tuples.
	# Each tuple is a couple of a list and a positive integer.
	# The list is a frequent-k-sequence and the integer is its support.
	# Overall = [
	#	[(freq-1, sup), ... , (freq-1, sup)],
	#	... ,
	#	[(freq-k, sup), ... , (freq-k, sup)] ]
	# Note that a k-sequence is also a list of lists.
	# freq-k = [[...], ... , [...]]
	Overall = []

	# Extract all the frequent-1-sequences
	itemsInDataset = sorted(set([item 
		for sublist1 in dataset # sublist1 = [([ ... ], t), ... , ([ ... ], t)]
			for sublist2 in sublist1 # sublist2 = ([ ... ], t)
				for item in sublist2[0] ])) # Must add the zero in order to pick the list.

	# Generate sequences of single items. 
	# Remeber, a sequence is a list of lists, hence a 1-sequence has the form: [ [item] ]
    singleItemSequences = [ [[item]] for item in itemsInDataset]

    # Count the frequence of the sequence i in the dataset and extract the frequent ones
    singleItemCounts = [(i, countFreq(i, dataset)) 
    	for i in singleItemSequences 
    		if countSupport(i, dataset) >= minsup]

    # Overall[0] = F1
    Overall.append(singleItemCounts)

    if verbose:
        print("Result, lvl 1: " + str(Overall[0]))

    k = 1
    while (True):

    	# if Fk is empty
        if not Overall[k - 1]:
            break

        # candidatesLastLevel = F(k-1)
    	candidatesLastLevel = [x[0] for x in Overall[k - 1]]

    	# Generate a set of k-sequences from the frequent (k-1)-sequences
    	candidatesGenerated = generateCandidates(candidatesLastLevel)

    	# Generate all the (k-1)-subsequence from a sequence candidate cand.
    	# If each direct subsequence is also a frequent-(k-1)-sequnce 
    	# (if is contained in candidatesLastLevel) then the k-sequence 
    	# is a candidate-k-sequence.
    	candidatesPruned = [cand 
    		for cand in candidatesGenerated 
    			if all(x in candidatesLastLevel 
    				for x in generateDirectSubsequences(cand))]

    	# Genetare a list of tuples (seq, sup) where sup is the support count
    	# of the k-sequence seq in the dataset.
    	# Remember that the support count must also consider the time constraint.
    	candidatesCounts = [(i, countSupport(dataset, i, maxspan, mingap, maxgap)) for i in candidatesPruned]

    	# Extract all the sequences with a support count grater then minsup
    	resultLvl = [(i, count) for (i, count) in candidatesCounts if (count >= minsup)]
        if verbose:
            print("Candidates generated, lvl " + str(k + 1) + ": " + str(candidatesGenerated))
            print("Candidates pruned, lvl " + str(k + 1) + ": " + str(candidatesPruned))
            print("Result, lvl " + str(k + 1) + ": " + str(resultLvl))
        Overall.append(resultLvl)
        k = k + 1

    # Remeber, the last element of Overall is empty so it must be removed
    Overall = Overall[:-1]
    
    # Rember the structure of Overall
    # [ [(freq-1, sup), ... , (freq-1, sup)]
    # 	...
    #	[(freq-k sup), ... , (freq-k, sup)] ]
 	# This last list comprehension transform Overll into
    # [ (freq-1, sup), ... , (freq-1, sup),
    #	(freq-2, sup), ... , (freq-2, sup),
    #	...
    #	(freq-k sup), ... , (freq-k, sup)]
	Overall = [item 
		for sublist in Overall 
			for item in sublist]
    return Overall