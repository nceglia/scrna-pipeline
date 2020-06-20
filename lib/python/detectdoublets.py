from singlecellexperiment import SingleCellExperiment
import scrublet
import sys
import numpy

sce = sys.argv[1]
csv = sys.argv[2]

sce = SingleCellExperiment.fromRData(sce)

counts = numpy.transpose(sce.assays["counts"])

scrub = scrublet.Scrublet(counts)

doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.call_doublets(threshold=0.25)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

output = open(csv,"w")
output.write("Score,Prediction\n")
for score, prediction in zip(doublet_scores, predicted_doublets):
    output.write("{},{}\n".format(score, prediction))
output.close()