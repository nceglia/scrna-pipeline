from singlecellexperiment import SingleCellExperiment
import scrublet
import sys
import numpy

sce = sys.argv[1]
csv = sys.argv[2]
print(sce)
sce = SingleCellExperiment.fromRData(sce)
print("here")
counts = numpy.transpose(sce.assays["counts"])
scrub = scrublet.Scrublet(counts)

doublet_scores, predicted_doublets = scrub.scrub_doublets(n_prin_comps=10)

scrub.call_doublets(threshold=0.35)
doublet_scores, predicted_doublets = scrub.scrub_doublets(n_prin_comps=10)


output = open(csv,"w")
output.write("Score,Prediction\n")
for score in doublet_scores:
    if score < 0.35:
        prediction = "False"
    else:
        prediction = "True"
    output.write("{},{}\n".format(score, prediction))
output.close()
