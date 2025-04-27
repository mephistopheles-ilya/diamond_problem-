import pyvista as pv

file1 = "/home/ilya/Downloads/unzipping/Reflect_70001_fcps1G_n.obj"
#file1 = "ChangedModels/Reflect_70001.obj"
file2 = "new_ChangedModels/Reflect_70001_fcps1G.obj"

mesh1 = pv.read(file1)
mesh2 = pv.read(file2)

plotter = pv.Plotter(shape=(1, 2))  # 1 row, 2 columns

# Left subplot (Model 1)
plotter.subplot(0, 0)
plotter.add_mesh(mesh1, color="green", show_edges=True)
plotter.add_title(file1, font_size=5)

# Right subplot (Model 2)
plotter.subplot(0, 1)
plotter.add_mesh(mesh2, color="red", show_edges=True)
plotter.add_title(file2, font_size=5)

plotter.show()
