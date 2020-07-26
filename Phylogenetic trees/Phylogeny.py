from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyTreeConstructor
from Bio import AlignIO, Phylo


def textToFasta(fileName):
    fileName = fileName[:-4]
    fastaFormat = open(fileName + '.Fasta', 'w+')
    with open(fileName + '.txt') as f:
        sequence = ""
        for fileLine in f:
            if fileLine[0] == "#":
                fastaFormat.write(sequence)
                sequence = ">" + fileLine[1:]
            else:
                sequence = sequence + fileLine
    fastaFormat.close()
    return fileName


def fastaToPhy(fileName):
    sequences = AlignIO.parse(fileName + '.fasta', "fasta")
    AlignIO.write(sequences, fileName + '.phy', 'phylip')
    return AlignIO.read(fileName + '.phy', 'phylip')


filePath = input("Enter the text file path:")
alignment = fastaToPhy(textToFasta(filePath))
distanceCalculator = DistanceCalculator('identity')
distanceMatrix = distanceCalculator.get_distance(alignment)

print('Distance Matrix:')
print(distanceMatrix)
print()

distanceConstructor = DistanceTreeConstructor()
##############################UPGMA##############################
UPGMATree = distanceConstructor.upgma(distanceMatrix)
Phylo.draw(UPGMATree)
print('UPGMA Phylogenetic Tree:')
Phylo.draw_ascii(UPGMATree)

##############################NJ##############################
NJTree = distanceConstructor.nj(distanceMatrix)
Phylo.draw(NJTree)
print('Neighbor joining Phylogenetic Tree:')
Phylo.draw_ascii(NJTree)

##############################PARSIMONY##############################
parsimonyScorer = Phylo.TreeConstruction.ParsimonyScorer()
searcher = Phylo.TreeConstruction.NNITreeSearcher(parsimonyScorer)
parsimonyConstructor = ParsimonyTreeConstructor(searcher, NJTree)
parsimonyTree = parsimonyConstructor.build_tree(alignment)
print('Parsimony Phylogenetic Tree:')
Phylo.draw_ascii(parsimonyTree)
Phylo.draw(parsimonyTree)
