import pyvista as pv

#file1 = "/home/ilya/Downloads/unzipping7/Reflect_70001_fcps1G_n.obj"
#file2 = "new_ChangedModels/Reflect_70001_fcps1G.obj"

#file1 = "/home/ilya/Downloads/unzipping1/Reflect_10001_fcps3d3G_n.obj"
#file2 = "new_ChangedModels/Reflect_10001_fcps3d3G.obj"

file1 = "/home/ilya/Downloads/unzipping3/Reflect_30002_fcps3a1d1G_n.obj"
file2 = "new_ChangedModels/Reflect_30002_fcps3a1d1G.obj"

#file1 = "/home/ilya/Downloads/unzipping4/Reflect_40001_fcs3a1d1G_n.obj"
#file2 = "new_ChangedModels/Reflect_40001_fcs3a1d1G.obj"

#file1 = "/home/ilya/Downloads/unzipping6/Reflect_60001_fps3a3d1G_n.obj"
#file2 = "new_ChangedModels/Reflect_60001_fps3a3d1G.obj"

#file1 = "/home/ilya/Downloads/unzipping7/Reflect_70001_fcs1a2G_n.obj"
#file2 = "new_ChangedModels/Reflect_70001_fcs1a2G.obj"



mesh1 = pv.read(file1)
mesh2 = pv.read(file2)

plotter = pv.Plotter(shape=(1, 2))  # 1 row, 2 columns

# Left subplot (Model 1)
plotter.subplot(0, 0)
plotter.add_mesh(mesh1, color="green", show_edges=True)
plotter.add_title(file1, font_size=5)
plotter.camera_position = 'xy'

# Right subplot (Model 2)
plotter.subplot(0, 1)
plotter.add_mesh(mesh2, color="red", show_edges=True)
plotter.add_title(file2, font_size=5)
plotter.camera_position = 'xy'

plotter.show()
