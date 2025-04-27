import pyvista as pv

mesh1 = pv.read("/home/ilya/Downloads/unzipping7/Reflect_70001_fcps1G_n.obj")
mesh2 = pv.read("new_ChangedModels/Reflect_70001_fcps1G.obj")

plotter = pv.Plotter(shape=(1, 2))  # 1 row, 2 columns

# Left subplot (Model 1)
plotter.subplot(0, 0)
plotter.add_mesh(mesh1, color="red")
plotter.add_title("Model 1")

# Right subplot (Model 2)
plotter.subplot(0, 1)
plotter.add_mesh(mesh2, color="blue")
plotter.add_title("Model 2")

plotter.show()
