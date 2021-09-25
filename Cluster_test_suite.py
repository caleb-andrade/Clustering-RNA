"""
Created on Thu Nov 12 16:27:41 2015
Created on Wed Nov 18 21:06:13 2015

@author: Caleb Andrade
"""
from Cluster import Cluster, codonAbundance

def readFastq(filename):
    """
    Parse read and quality strings from a FASTQ file with sequencing reads.
    @author: Ben Langmead & Jacob Pritt.
    
    Input: file path
    Output: A list of reads, the list of qualities
    """
    sequences = []
    qualities = []
    
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() #read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
            
    return sequences, qualities
    

def kmerHashMap(reads, k):
    """
    Create a hash map between kmers and readings. 
    Note: Keys are kmers. Values are sets of those reads containing them.
    
    Input: list of reads, length of kmers.
    Output: hash map.
    """
    kmers_dict = {}
    # loop through all reads
    for i in range(len(reads)):
        # loop read's bases, except for the last k, to obtain its kmers
        for j in range(1+len(reads[i])-k):
            kmer = reads[i][j:k+j]
            if kmers_dict.has_key(kmer):
                kmers_dict[kmer].add(i)
            else:
                kmers_dict[kmer] = set([i])
    
    return kmers_dict

    
# Following are reads of a portion of the sequences of three species
reads1 = readFastq('ERR266411_1.for_asm.fastq')[0][:10]
reads2 = readFastq('ERR037900_1.first1000.fastq')[0][:10]

print "\n************************************************"
print "\nFollowing a sample of data to test cluster class: \n"    
print "     reads                  reads' abundance"
read_abundance1 = codonAbundance(reads1)
for i in range(len(read_abundance1)):
    print reads1[i], " ", read_abundance1[i]

reads = []
for i in range(10):
    reads.append((i, reads1[i]))

cluster1 = Cluster(reads)
cluster4 = Cluster([reads[0]])
cluster5 = Cluster([reads[1]])

print "\n*************** TESTING print ******************"
print "\nCluster1: ", cluster1
print "\nCluster4: ", cluster4
print "\nCluster5: ", cluster5

print "\n*************** TESTING getIDs *****************"
print cluster1.getIDs()

print "\n*************** TESTING getSize ****************"
print cluster1.getSize()

print "\n************ TESTING getCodonAbundance ***********"
for item in cluster1.getCodonAbundance():
    print item

print "\n************** TESTING getAvgAbundance *************"
print cluster1.getAvgAbundance()

print "\n************** TESTING getSumAbundance *************"
print cluster1.getSumAbundance()


print "\n************************************************"
print "\nAnother sample of data to test cluster mutability: \n"    
print "     reads                  reads' abundance"
read_abundance2 = codonAbundance(reads2)
for i in range(len(read_abundance2)):
    print reads2[i], " ", read_abundance2[i]

reads = []
for i in range(10):
    reads.append((i+10, reads2[i]))

cluster2 = Cluster(reads)

print "\n****************** TESTING copy ****************"
cluster3 = cluster1.copy()
print "\ncopy cluster1 into cluster3: "
print cluster3

print "\n*************** TESTING isEqualTo ***************"
print "\nIs cluster1 equal to cluster3?", cluster1.isEqualTo(cluster3)
print "\nIs cluster1 equal to cluster2?", cluster1.isEqualTo(cluster2)
# Uncommenting the following line raises an exception
#print "\nIs cluster1 equal to None?", cluster1.isEqualTo(None)

print "\n*************** TESTING update ******************"
print "\nFirst, lets reset average abundance vector: "
cluster1.avg_abundance_vectors = 64*[0]
print cluster1
cluster1.update()
print "\nAfter update: "
print cluster1

print "\n*************** TESTING distance ***************"
print "\nLet's test distance between cluster1 and cluster2: \n"
print "Cluster1's: ", cluster1.getAvgAbundance()
print "Cluster2's: ", cluster2.getAvgAbundance()
print "\nManhattan: ", cluster1.distance(cluster2)
print "Euclidean: ", cluster2.distance(cluster1, dim = 2)

print "\n*************** TESTING merge ******************"
print "\nLet's merge cluster1 with its copy, cluster3: "
cluster1.mergeClusters(cluster3)
print cluster1
print "\nLet's now merge cluster1 with cluster2: "
cluster1.mergeClusters(cluster2)
print cluster1

print "\n*************** TESTING error ******************"
print "\nCluster4: ", cluster4
print "\nCluster5: ", cluster5
print "\nMerge cluster4 & cluster5"
cluster4.mergeClusters(cluster5)
print cluster4
print cluster4.getAvgAbundance()

