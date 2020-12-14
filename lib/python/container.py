import subprocess

class Container(object):

    def __init__(self):
        self.runtime = None
        self.container_command = None

    def set_runtime(self, runtime):
        if runtime == "singularity":
            self.container_command = "singularity exec --bind /juno/work/shah/:/juno/work/shah/ --bind /work/shah/:/work/shah".split()
            self.runtime = runtime
        elif runtime == "docker":
            self.container_command = "docker run --mount type=bind,source=/Users/ceglian/,target=/Users/ceglian/ -w $(pwd)"
            self.runtime = runtime
        else:
            self.container_command = ""
            self.runtime = "local"

    def set_image(self, image):
        self.image = image

    def run(self, script, args, outs):
        cmd = [script.interpreter, script.path]
        cmd = cmd + [None for _ in range(len(script.args.keys())+len(script.outs.keys()))]
        print(len(cmd), cmd)
        for position, name in script.args.items():
            print(position, name)
            cmd[position+1] = str(getattr(args, name))
        for position, name in script.outs.items():
            cmd[position+1] = str(getattr(outs, name))
        full_cmd = self.container_command + [self.image] + cmd
        subprocess.check_output(full_cmd)
