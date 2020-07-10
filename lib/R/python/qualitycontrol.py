import container
import os

class QualityControl(object):

    def __init__(self, matrix, output, filtered_matrix):
        self.matrix = matrix
        self.output = output
        self.filtered_matrix = filtered_matrix
        self.current = os.path.split(os.path.realpath(__file__))[0]
    
    def container(self):
        script = os.path.join(self.current,"../R/qc.R")
        return container.Container(["Rscript", script, self.matrix, self.output, self.filtered_matrix])

