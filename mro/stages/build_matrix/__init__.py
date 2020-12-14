import container
import scriptmanager
import martian

__MRO__ = '''
stage BUILD_MATRIX(
    in  path inventory,
    in  path image,
    in  string runtime,
    out h5 project_matrix,
    out path matrix,
    out path features,
    out path barcodes,
    src py "stages/build_matrix",
) using (
    threads = 40,
    mem_gb  = 10,
)
'''

def main(args, outs):
    project_matrix = "matrix.h5"
    outs.project_matrix = martian.make_path(project_matrix)
    scripts = scriptmanager.ScriptManager()
    script = scripts.buildmatrix()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

