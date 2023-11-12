import NvTK

NvTK.load_module("../BenchmarksInManuscript/BaselineModel.py")

print(dir(NvTK.BaselineModel))

print(NvTK.BaselineModel.BaselineEpigCNN(output_size=10))
