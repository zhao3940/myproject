#!/usr/bin/env python3
# Author: Dr. Robert Heckendorn
# Computer Science Department
# University of Idaho
# 2017
#
# ID3 decision tree algorithm
#
import sys
from math import *
from optparse import OptionParser

Epsilon = 1E-10


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#
# IO support for reading from stdin and writing to stdout
#

# read in a classifier problem
def readProblem():
    global Features, FeatureList, FeatureValues, Data

    FeatureList = []  # list of all the features including Ans
    FeatureValues = {}  # potential values of all the features, even ones not in examples
    Data = []  # example classification data

    # read number of features
    numFeatures = int(sys.stdin.readline())

    # read in features and answer which must be called: Ans
    for i in range(0, numFeatures + 1):
        line = sys.stdin.readline().strip()
        fields = line.split()
        FeatureList.append(fields[0])
        FeatureValues[fields[0]] = fields[2:]  # dictionary value is a list

    # read number of samples
    numSamples = int(sys.stdin.readline())

    # read in example classifications
    for line in sys.stdin.readlines():
        fields = line.split()
        sample = {}
        for i in range(0, len(FeatureList)):
            sample[FeatureList[i]] = fields[i]
        Data.append(sample)


# write out indented classifier tree
amountIndent = 3 * " "


def printDTree(tree):
    printDTreeAux("", tree)


def printDTreeAux(indent, tree):
    name = tree[0]
    d = tree[1]
    if type(d) is dict:
        for v in FeatureValues[name]:
            print(indent + name + "=" + v)
            printDTreeAux(indent + amountIndent, d[v])
    else:
        print(indent + d)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#
# Decision Tree:  list either [Feature, Value] or [Feature, dictionary of value:tree]
# Author: Robert B. Heckendorn, University of Idaho
#
# The OUTPUT:
#
# ["Ans", "Yes"]
# or
# ["Price", {"$" : tree1, "$$" : tree2, "$$$" : tree1}]
#
#
# DATA: a list of dictionary for each training case
#
# [{"Ans" : "Yes", "Alt" : "Yes", "Price" : "$$" ... },
#  ... ]
#
# FEATUREVALUES: dictionary of Feature:[list of values]
#
# {"Ans" : ["Yes", "No"], "Alt" : ["Yes", "No"], "Price" : ["$", "$$", "$$$"] ... }
#

#
# select(data, feature, value) -> subset of data as list of dictionaries that have feature=value
# count(data, feature, value) -> number of cases in data in which feature has value
# entropy(data, feature) -> [0, 1], feature is generally "Ans" for ID3


# list of the items in data that have feature equal to value
def select(data, feature, value):
    return [item for item in data if item[feature] == value]


# count how many items in the data have feature equal to value
def count(data, feature, value):
    num = 0
    for d in data:
        if d[feature] == value: num += 1
    return num


# what is the entropy of a question about feature?
# sum the entropy over the possible values of the feature.
def entropy(data, feature):
    lens=len(data)
    entropys = 0.0
    if lens==0:
        return entropys
    for v in FeatureValues[feature]:
        num=count(data,"Ans",v)
        pro=num/float(lens)
        if pro != 0:
            entropys-=pro*(log(pro,2))
    return entropys


# current entropy - expected entropy after getting info about feature
# entropy(data, "Ans") - sum_{v=featurevalues} p_v * entropy(select(data, feature, v), "Ans")
def gain(data, feature):
    lens=len(data)
    sum=0.0
    if lens == 0:
        sum=0.0
    else:
        for i in FeatureValues[feature]:
            p_v=float(count(data,feature,i))/lens
            sum+=p_v*entropy(select(data,feature,i),"Ans")
    return entropy(data,"Ans")-sum


# If there one and only one value for the given feature in given data
# If not return None
def isOneLabel(data, feature):
    old = None
    for i in data:
        new = i["Ans"]
        if new == old or old == None:
            old = new
        else:
            return None
    return old

# select the most popular Ans value left in the data for the constraints
# up to now.
def maxAns(data):
    maxans=0
    ans=None
    for i in FeatureValues["Ans"]:
        if count(data,"Ans",i)>maxans or maxans==0:
            maxans=count(data,"Ans",i)
            ans=i
    return ans
# this is the ID3 algorithm
def ID3BuildTree(data, availableFeatures):
# only one value for the Ans feature at this point?
    if isOneLabel(data, availableFeatures)!= None:
        return ["Ans",isOneLabel(data, availableFeatures)]

# ran out of discriminating features
    if len(availableFeatures)==0:
        print("***  out of features ***")
        return ["Ans",maxAns(data)]

# pick maximum information gain
    else:
        bestFeature = None
        bestGain = None
        for feature in availableFeatures:
            g = gain(data, feature)
            print("GAIN: ", feature, ":", round(g, 4));
            if bestGain == None or g > bestGain + Epsilon:
                bestGain = g
                bestFeature = feature
                bestList = [feature]
            elif abs(g - bestGain) < Epsilon:
                bestList.append(feature)
        print("BEST:", round(bestGain, 4), bestList);
        print()

# recursively construct tree on return
    treeLeaves = {}  # start with empty dictionary
    availableFeatures = availableFeatures[:]
    availableFeatures.remove(bestFeature)

#??? something IN this loop to not build a subtree if data is empty for any feature value
    if  len(data)==0:
        return ["Ans",maxAns(data)]
    for v in FeatureValues[bestFeature]:
        treeLeaves[v] = ID3BuildTree(select(data, bestFeature, v), availableFeatures)  # recurse

        return [bestFeature, treeLeaves]  # list of best feature and dictionary of trees

# read the tree data structure from the treeFile
def getTree(treeFile):
    try:
        inn = open(treeFile, "r")
        treeTxt = inn.readline()
        inn.close()
        return eval(treeTxt)
    except:
        print("ERROR: unable to open", treeFile)
        return None

# use the tree data structure to compute an Ans for each data line
def ID3DecodeTree(tree, data):
    for d in data:
        printDataLine(FeatureList, d)
        dtAns = getAns(tree, d)  # here is where we compute the Ans from tree
        print(dtAns, end=" ")
        if d["Ans"] == dtAns:
            print("Equal")
        else:
            print("NotEqual")

# print the data line only with no newline
def printDataLine(features, d):
    for f in features:
        print(f + ":" + d[f], end="   ")

def getAns(tree, d):
    name = tree[0]
    tem=tree[1]
    if type(tem) is dict:
        value=d[name]
        del d[name]
        return (getAns(tem[value],d))
    else:
        return tree[1]

def main():
    # parse the command line args
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="treeFile",
                      help="The decision tree")
    (options, args) = parser.parse_args()
    treeFile = options.treeFile

    # read the problem
    readProblem()

    # if -t option then use tree data structure otherwise compute tree data structure
    if options.treeFile:
        tree = getTree(options.treeFile)
        if tree:
            ID3DecodeTree(tree, Data)
    else:
        FeatureList.remove("Ans")
        tree = ID3BuildTree(Data, FeatureList)
        printDTree(tree)
        print(tree)

main()
