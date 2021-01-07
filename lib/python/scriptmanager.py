import os

class Script(object):
    def __init__(self, name):
        self.name = name
        self.args = dict()
        self.outs = dict()

    def set_path(self, path):
        self.path = path
        ext = os.path.splitext(path)[1]
        if ext == ".py":
            self.interpreter = "python3"
        elif ext == ".R":
            self.interpreter = "Rscript"
        else:
            self.interpreter = "sh"

    def set_args(self, position, value):
        self.args[position] = value

    def set_outs(self, position, value):
        self.outs[position] = value

    def set_interpreter(self, interpreter):
        self.interpreter = interpreter

class ScriptManager(object):

    def __init__(self):
        self.current = os.path.split(os.path.realpath(__file__))[0]

    def assigncelltypes(self):
        script = Script("assigncelltypes")
        script.set_args(1, "merged_batch_h5ad")
        script.set_args(2, "marker_csv")
        script.set_outs(3, "celltype_csv")
        script.set_path(os.path.join(self.current,"../R/cellassign.R"))
        return script